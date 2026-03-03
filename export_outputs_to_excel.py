"""
Export PhysiCell outputs to time-aligned CSV summaries.

This script scans PhysiCell `output*.xml` files, parses each timestep via `pyMCDS`,
and writes multiple CSV tables for downstream analysis and plotting.

Input modes (exactly one required):
- `--folder PATH`
    Read `output*.xml` files from a folder tree.
- `--zip PATH`
    Read from a zip archive. The zip is either extracted to a temporary directory
    (default) or to a persistent folder with the same zip base name.

Flags:
- `--zip PATH`
    Path to zip containing PhysiCell outputs.
- `--folder PATH`
    Path to extracted output folder containing `output*.xml` files.
- `--extract-zip-to-folder` / `--extract-zuip-to-folder`
    When used with `--zip`, extract to `<zip_parent>/<zip_stem>` and process there.
    The `--extract-zuip-to-folder` form is kept as a compatibility alias.
- `--outdir NAME`
    Output directory for CSVs. If `NAME` contains `TIMESTAMP`, it is replaced with
    current datetime (`YYYYmmdd_HHMMSS`).
- `--domain-x FLOAT`, `--domain-y FLOAT`, `--domain-z FLOAT`
    Optional domain dimensions used for density calculations.
- `--every-nth N`
    Process every Nth XML timestep (N >= 1). Use this to speed up long runs.
- `--max-files N`
    Optional cap on number of selected XML files after sampling.
- `--workers N`
    Number of parallel worker processes (N >= 1). Use >1 to parse files in parallel.
- `--fast-mode`
    Skip expensive microenvironment, spatial, and attribute stats for faster runtime.

Primary output CSVs:
- `summary_by_time.csv`
- `density_by_time.csv`
- `cell_type_counts_long.csv`
- `cell_phase_counts_long.csv`
- `cell_type_phase_counts.csv`
- `cell_type_counts_wide.csv`
- `cell_phase_counts_wide.csv`
- `spatial_stats_by_type.csv`
- `region_counts.csv`
- `cell_attribute_stats.csv`
- `microenvironment_stats.csv`
- `time_aligned_long.csv`
- `metadata.csv`

Examples:
1) Basic run from extracted folder
   python export_outputs_to_excel.py --folder simulation_results_2 --outdir Physicell_results_2

2) Run directly from zip and persist extraction folder
   python export_outputs_to_excel.py --zip simulation_results_2.zip --extract-zip-to-folder --outdir Physicell_results_2

3) Speed-focused run (sample timesteps + parallel + fast mode)
   python export_outputs_to_excel.py --folder simulation_results_2 --every-nth 5 --workers 6 --fast-mode --outdir Physicell_results_fast

4) Bounded smoke test run
   python export_outputs_to_excel.py --folder simulation_results_2 --max-files 10 --every-nth 2 --workers 2 --outdir Physicell_results_smoke

Notes:
- The script attempts `pcdl.pyMCDS` first, then falls back to `beta/pyMCDS.py`.
- On Windows with process-heavy workloads, start with `--workers 2` or `--workers 4`
  and increase gradually based on memory and CPU behavior.
"""

import argparse
import concurrent.futures
import importlib.util
import inspect
import os
import re
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd


def load_pyMCDS():
    try:
        from pcdl import pyMCDS as pcdl_pyMCDS
        return pcdl_pyMCDS
    except Exception as exc:
        script_dir = Path(__file__).resolve().parent
        fallback_path = script_dir / "beta" / "pyMCDS.py"
        if not fallback_path.exists():
            raise RuntimeError(
                f"Failed to import pyMCDS from pcdl ({exc}) and fallback file not found at: {fallback_path}"
            ) from exc

        spec = importlib.util.spec_from_file_location("local_pyMCDS", fallback_path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Unable to load fallback module from: {fallback_path}") from exc

        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if not hasattr(module, "pyMCDS"):
            raise RuntimeError(f"Fallback module does not define pyMCDS: {fallback_path}") from exc

        print(f"Warning: using fallback parser from {fallback_path} due to pcdl import error: {exc}")
        return module.pyMCDS




OUTPUT_RE = re.compile(r"output(\d+)\.xml$", re.IGNORECASE)
TYPE_COLUMNS = ["cell_type_name", "cell_type"]
PHASE_COLUMNS = ["cycle_phase", "current_phase", "phase", "state"]
POSITION_COLUMNS = ["position_x", "position_y", "position_z"]

_WORKER_PYMCDS = None
_WORKER_CONSTRUCTOR_MODE = None


def resolve_output_files(root: Path) -> list[Path]:
    candidates = list(root.rglob("output*.xml"))

    def sort_key(path: Path):
        match = OUTPUT_RE.search(path.name)
        if match:
            return (0, int(match.group(1)))
        return (1, path.name)

    return sorted(candidates, key=sort_key)


def safe_value(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def output_dir_name(raw_name: str) -> Path:
    name = raw_name
    if "TIMESTAMP" in name:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        name = name.replace("TIMESTAMP", timestamp)
    return Path(name)


def zip_extract_folder(zip_path: Path) -> Path:
    return zip_path.parent / zip_path.stem


def create_mcds(pyMCDS_class, xml_name: str, output_path: str):
    init_fn = getattr(pyMCDS_class, "__init__", None)
    if init_fn is not None:
        try:
            params = inspect.signature(init_fn).parameters
            if "xml_file_name" in params and "output_path" in params:
                return pyMCDS_class(xml_file_name=xml_name, output_path=output_path)
            if "xml_file" in params and "output_path" in params:
                return pyMCDS_class(xml_file=xml_name, output_path=output_path)
            if "xml_file" in params:
                return pyMCDS_class(xml_file=xml_name)
        except Exception:
            pass

    for constructor in (
        lambda: pyMCDS_class(xml_file_name=xml_name, output_path=output_path),
        lambda: pyMCDS_class(xml_file=xml_name, output_path=output_path),
        lambda: pyMCDS_class(xml_name, output_path),
        lambda: pyMCDS_class(xml_name),
    ):
        try:
            return constructor()
        except TypeError:
            continue

    raise TypeError("Unable to construct pyMCDS with known argument signatures.")


def resolve_mcds_constructor_mode(pyMCDS_class) -> str:
    init_fn = getattr(pyMCDS_class, "__init__", None)
    if init_fn is not None:
        try:
            params = inspect.signature(init_fn).parameters
            if "xml_file_name" in params and "output_path" in params:
                return "xml_file_name_output_path"
            if "xml_file" in params and "output_path" in params:
                return "xml_file_output_path"
            if "xml_file" in params:
                return "xml_file_only"
        except Exception:
            pass

    for mode, constructor in (
        ("xml_file_name_output_path", lambda: pyMCDS_class(xml_file_name="__dummy__.xml", output_path=".")),
        ("xml_file_output_path", lambda: pyMCDS_class(xml_file="__dummy__.xml", output_path=".")),
        ("positional_output_path", lambda: pyMCDS_class("__dummy__.xml", ".")),
        ("positional_only", lambda: pyMCDS_class("__dummy__.xml")),
    ):
        try:
            constructor()
            return mode
        except TypeError:
            continue
        except Exception:
            return mode

    return "unknown"


def create_mcds_with_mode(pyMCDS_class, constructor_mode: str, xml_name: str, output_path: str):
    if constructor_mode == "xml_file_name_output_path":
        return pyMCDS_class(xml_file_name=xml_name, output_path=output_path)
    if constructor_mode == "xml_file_output_path":
        return pyMCDS_class(xml_file=xml_name, output_path=output_path)
    if constructor_mode == "xml_file_only":
        return pyMCDS_class(xml_file=xml_name)
    if constructor_mode == "positional_output_path":
        return pyMCDS_class(xml_name, output_path)
    if constructor_mode == "positional_only":
        return pyMCDS_class(xml_name)

    return create_mcds(pyMCDS_class, xml_name, output_path)


def sample_xml_files(xml_files: list[Path], every_nth: int, max_files: int | None) -> list[Path]:
    selected = xml_files[::every_nth]
    if max_files is not None:
        selected = selected[:max_files]
    return selected


def _process_xml_file(
    pyMCDS_class,
    constructor_mode: str,
    xml_file: Path,
    domain_dims: tuple[float | None, float | None, float | None],
    include_microenv: bool,
    include_spatial: bool,
    include_attribute: bool,
):
    output_path = str(xml_file.parent)
    xml_name = xml_file.name

    mcds = create_mcds_with_mode(pyMCDS_class, constructor_mode, xml_name, output_path)
    time_min = safe_value(mcds.get_time())
    cell_df = mcds.get_cell_df()

    type_col = choose_column(cell_df, TYPE_COLUMNS)
    phase_col = choose_column(cell_df, PHASE_COLUMNS)

    summary_rows = []
    type_rows = []
    phase_rows = []
    type_phase_rows = []
    micro_rows = []
    spatial_rows = []
    region_rows = []
    attribute_rows = []
    density_rows = []

    summary_row = {
        "xml_file": xml_name,
        "time_min": time_min,
        "total_cells": int(len(cell_df)),
    }
    if "dead" in cell_df.columns:
        dead_count = int(cell_df["dead"].sum())
        summary_row["dead_cells"] = dead_count
        summary_row["live_cells"] = int(len(cell_df)) - dead_count
    summary_rows.append(summary_row)

    type_counts = count_by_column(cell_df, type_col)
    if not type_counts.empty:
        for _, row in type_counts.iterrows():
            type_rows.append(
                {
                    "xml_file": xml_name,
                    "time_min": time_min,
                    "cell_type": safe_value(row[type_col]),
                    "count": int(row["count"]),
                }
            )

    phase_counts = count_by_column(cell_df, phase_col)
    if not phase_counts.empty:
        for _, row in phase_counts.iterrows():
            phase_rows.append(
                {
                    "xml_file": xml_name,
                    "time_min": time_min,
                    "phase": safe_value(row[phase_col]),
                    "count": int(row["count"]),
                }
            )

    if type_col and phase_col:
        type_phase_counts = (
            cell_df.groupby([type_col, phase_col], dropna=False)
            .size()
            .reset_index(name="count")
        )
        for _, row in type_phase_counts.iterrows():
            type_phase_rows.append(
                {
                    "xml_file": xml_name,
                    "time_min": time_min,
                    "cell_type": safe_value(row[type_col]),
                    "phase": safe_value(row[phase_col]),
                    "count": int(row["count"]),
                }
            )

    if include_microenv:
        micro_stats = compute_microenv_stats(mcds)
        if not micro_stats.empty:
            micro_stats.insert(0, "xml_file", xml_name)
            micro_stats.insert(1, "time_min", time_min)
            micro_rows.extend(micro_stats.to_dict("records"))

    if type_col and include_spatial:
        spatial_stats, region_counts = compute_spatial_stats(cell_df, type_col)
        if not spatial_stats.empty:
            spatial_stats.insert(0, "xml_file", xml_name)
            spatial_stats.insert(1, "time_min", time_min)
            spatial_rows.extend(spatial_stats.to_dict("records"))

        if not region_counts.empty:
            region_counts.insert(0, "xml_file", xml_name)
            region_counts.insert(1, "time_min", time_min)
            region_counts = region_counts.rename(columns={type_col: "cell_type"})
            region_rows.extend(region_counts.to_dict("records"))

    if type_col and include_attribute:
        attribute_stats = compute_attribute_stats(cell_df, type_col)
        if not attribute_stats.empty:
            attribute_stats.insert(0, "xml_file", xml_name)
            attribute_stats.insert(1, "time_min", time_min)
            attribute_stats = attribute_stats.rename(columns={type_col: "cell_type"})
            attribute_rows.extend(attribute_stats.to_dict("records"))

    density_stats = compute_density(cell_df, domain_dims)
    if not density_stats.empty:
        density_stats.insert(0, "xml_file", xml_name)
        density_stats.insert(1, "time_min", time_min)
        density_rows.extend(density_stats.to_dict("records"))

    return {
        "type_col": type_col,
        "phase_col": phase_col,
        "summary_rows": summary_rows,
        "type_rows": type_rows,
        "phase_rows": phase_rows,
        "type_phase_rows": type_phase_rows,
        "micro_rows": micro_rows,
        "spatial_rows": spatial_rows,
        "region_rows": region_rows,
        "attribute_rows": attribute_rows,
        "density_rows": density_rows,
    }


def _worker_init(constructor_mode: str):
    global _WORKER_PYMCDS, _WORKER_CONSTRUCTOR_MODE
    _WORKER_PYMCDS = load_pyMCDS()
    _WORKER_CONSTRUCTOR_MODE = constructor_mode


def _worker_process(
    xml_file_str: str,
    domain_dims: tuple[float | None, float | None, float | None],
    include_microenv: bool,
    include_spatial: bool,
    include_attribute: bool,
):
    xml_file = Path(xml_file_str)
    return _process_xml_file(
        _WORKER_PYMCDS,
        _WORKER_CONSTRUCTOR_MODE,
        xml_file,
        domain_dims,
        include_microenv,
        include_spatial,
        include_attribute,
    )


def choose_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def count_by_column(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    if column_name is None or column_name not in df.columns:
        return pd.DataFrame()

    counts = df[column_name].value_counts(dropna=False).rename_axis(column_name).reset_index(name="count")
    return counts


def compute_microenv_stats(mcds) -> pd.DataFrame:
    stats_rows = []

    try:
        substrates = mcds.get_substrate_names()
    except Exception:
        substrates = []

    for substrate in substrates:
        values = None
        try:
            values = mcds.get_concentration(substrate)
        except Exception:
            try:
                values = mcds.get_substrate(substrate)
            except Exception:
                values = None

        if values is None:
            continue

        values = np.asarray(values).astype(float)
        stats_rows.append(
            {
                "substrate": substrate,
                "mean": safe_value(np.nanmean(values)),
                "min": safe_value(np.nanmin(values)),
                "max": safe_value(np.nanmax(values)),
                "std": safe_value(np.nanstd(values)),
                "p10": safe_value(np.nanpercentile(values, 10)),
                "p50": safe_value(np.nanpercentile(values, 50)),
                "p90": safe_value(np.nanpercentile(values, 90)),
            }
        )

    return pd.DataFrame(stats_rows)


def compute_spatial_stats(cell_df: pd.DataFrame, type_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not all(col in cell_df.columns for col in POSITION_COLUMNS[:2]):
        return pd.DataFrame(), pd.DataFrame()

    has_z = "position_z" in cell_df.columns
    positions = cell_df[["position_x", "position_y"] + (["position_z"] if has_z else [])].copy()
    center = positions.mean(axis=0)
    deltas = positions - center

    if has_z:
        r = np.sqrt(deltas["position_x"] ** 2 + deltas["position_y"] ** 2 + deltas["position_z"] ** 2)
    else:
        r = np.sqrt(deltas["position_x"] ** 2 + deltas["position_y"] ** 2)

    cell_df = cell_df.copy()
    cell_df["radius"] = r

    group = cell_df.groupby(type_col, dropna=False)
    stats = group[POSITION_COLUMNS[:2] + (["position_z"] if has_z else [])].agg(["mean", "std"])
    stats.columns = [f"{col}_{stat}" for col, stat in stats.columns]

    radius_stats = group["radius"].agg(["mean", "std", "min", "max"])
    radius_stats.columns = [f"radius_{stat}" for stat in radius_stats.columns]

    spatial_stats = pd.concat([stats, radius_stats], axis=1).reset_index()

    r_max = float(r.max()) if len(r) else 0.0
    if r_max <= 0:
        cell_df["region"] = "core"
    else:
        frac = r / r_max
        bins = [0.0, 0.25, 0.5, 0.75, 1.0]
        labels = ["core", "inner", "outer", "rim"]
        cell_df["region"] = pd.cut(frac, bins=bins, labels=labels, include_lowest=True)

    region_counts = (
        cell_df.groupby([type_col, "region"], dropna=False, observed=False)
        .size()
        .reset_index(name="count")
    )
    return spatial_stats, region_counts


def compute_attribute_stats(cell_df: pd.DataFrame, type_col: str) -> pd.DataFrame:
    numeric_cols = [
        col
        for col in cell_df.columns
        if pd.api.types.is_numeric_dtype(cell_df[col])
        and col not in {"position_x", "position_y", "position_z", "ID", "cell_ID", "cell_id"}
    ]
    if not numeric_cols:
        return pd.DataFrame()

    group = cell_df.groupby(type_col, dropna=False)
    stats = group[numeric_cols].agg(["mean", "std", "min", "max"])
    stats = stats.stack(level=0, future_stack=True).reset_index()
    stats = stats.rename(columns={"level_1": "attribute"})
    return stats


def compute_density(cell_df: pd.DataFrame, domain_dims: tuple[float | None, float | None, float | None]) -> pd.DataFrame:
    x_dim, y_dim, z_dim = domain_dims
    rows = []

    if x_dim and y_dim:
        if z_dim and z_dim > 0:
            volume = x_dim * y_dim * z_dim
            density = len(cell_df) / volume if volume > 0 else np.nan
            rows.append({"density_source": "domain_volume", "density": safe_value(density)})
        else:
            area = x_dim * y_dim
            density = len(cell_df) / area if area > 0 else np.nan
            rows.append({"density_source": "domain_area", "density": safe_value(density)})

    if all(col in cell_df.columns for col in POSITION_COLUMNS[:2]):
        x_extent = cell_df["position_x"].max() - cell_df["position_x"].min()
        y_extent = cell_df["position_y"].max() - cell_df["position_y"].min()
        if "position_z" in cell_df.columns:
            z_extent = cell_df["position_z"].max() - cell_df["position_z"].min()
            volume = x_extent * y_extent * z_extent
            density = len(cell_df) / volume if volume > 0 else np.nan
            rows.append({"density_source": "bbox_volume", "density": safe_value(density)})
        else:
            area = x_extent * y_extent
            density = len(cell_df) / area if area > 0 else np.nan
            rows.append({"density_source": "bbox_area", "density": safe_value(density)})

    return pd.DataFrame(rows)


def format_eta(seconds: float) -> str:
    if seconds < 0 or np.isnan(seconds):
        return "--:--"
    seconds = int(round(seconds))
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours > 0:
        return f"{hours:02d}:{minutes:02d}:{sec:02d}"
    return f"{minutes:02d}:{sec:02d}"


def render_progress(idx: int, total: int, start_time: float) -> None:
    elapsed = time.time() - start_time
    rate = elapsed / max(idx, 1)
    eta = rate * (total - idx)
    bar_width = 30
    filled = int(bar_width * idx / max(total, 1))
    bar = "#" * filled + "-" * (bar_width - filled)
    msg = f"[{bar}] {idx}/{total} | elapsed {format_eta(elapsed)} | eta {format_eta(eta)}"
    print(msg, end="\r", file=sys.stdout, flush=True)


def main():
    pyMCDS = load_pyMCDS()
    constructor_mode = resolve_mcds_constructor_mode(pyMCDS)

    parser = argparse.ArgumentParser(description="Export PhysiCell outputs to CSV summaries.")
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--zip", dest="zip_path", help="Path to zip with output*.xml files")
    input_group.add_argument("--folder", dest="folder_path", help="Path to folder with output*.xml files")
    parser.add_argument(
        "--extract-zip-to-folder",
        "--extract-zuip-to-folder",
        action="store_true",
        help="When used with --zip, extract zip contents into a folder named like the zip and parse from there.",
    )
    parser.add_argument(
        "--outdir",
        dest="out_dir",
        default="physicell_csvs_TIMESTAMP",
        help="Output directory name. Use TIMESTAMP to inject a timestamp.",
    )
    parser.add_argument("--domain-x", type=float, default=None, help="Domain size X (same units as positions)")
    parser.add_argument("--domain-y", type=float, default=None, help="Domain size Y (same units as positions)")
    parser.add_argument("--domain-z", type=float, default=None, help="Domain size Z (optional for 3D)")
    parser.add_argument("--every-nth", type=int, default=1, help="Process every Nth output file (default: 1)")
    parser.add_argument("--max-files", type=int, default=None, help="Optional cap on number of output files to process")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel worker processes for XML parsing (default: 1)",
    )
    parser.add_argument(
        "--fast-mode",
        action="store_true",
        help="Skip expensive microenvironment, spatial, and attribute stats for faster export.",
    )

    args = parser.parse_args()
    if args.every_nth < 1:
        raise ValueError("--every-nth must be >= 1")
    if args.max_files is not None and args.max_files < 1:
        raise ValueError("--max-files must be >= 1 when provided")
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")

    out_dir = output_dir_name(args.out_dir)

    if args.zip_path:
        input_path = Path(args.zip_path)
        if not input_path.exists():
            raise FileNotFoundError(f"Zip not found: {input_path}")
        use_zip = True
    else:
        input_path = Path(args.folder_path)
        if not input_path.exists():
            raise FileNotFoundError(f"Folder not found: {input_path}")
        use_zip = False

    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Starting export from: {input_path}")
    print(f"Writing CSVs to: {out_dir.resolve()}")

    summary_rows = []
    type_rows = []
    phase_rows = []
    type_phase_rows = []
    micro_rows = []
    spatial_rows = []
    region_rows = []
    attribute_rows = []
    density_rows = []

    type_col_used = None
    phase_col_used = None

    domain_dims = (args.domain_x, args.domain_y, args.domain_z)
    include_microenv = not args.fast_mode
    include_spatial = not args.fast_mode
    include_attribute = not args.fast_mode

    # Determine search path - extract zip if needed
    temp_context = tempfile.TemporaryDirectory() if use_zip and not args.extract_zip_to_folder else None
    try:
        if use_zip:
            if args.extract_zip_to_folder:
                extract_path = zip_extract_folder(input_path)
                extract_path.mkdir(parents=True, exist_ok=True)
                print(f"Extracting zip to: {extract_path.resolve()}")
            else:
                extract_path = Path(temp_context.name)
                print("Extracting zip...")

            with ZipFile(input_path, "r") as zf:
                zf.extractall(extract_path)

            search_path = extract_path
        else:
            search_path = input_path

        xml_files = resolve_output_files(search_path)
        if not xml_files:
            raise FileNotFoundError("No output*.xml files found.")

        original_count = len(xml_files)
        xml_files = sample_xml_files(xml_files, args.every_nth, args.max_files)
        if not xml_files:
            raise FileNotFoundError("No output*.xml files selected after sampling settings.")

        print(f"Found {original_count} output XML files.")
        print(f"Processing {len(xml_files)} file(s) with every_nth={args.every_nth}, workers={args.workers}.")
        if args.fast_mode:
            print("Fast mode enabled: skipping microenvironment/spatial/attribute stats.")

        start_time = time.time()
        total_files = len(xml_files)
        if args.workers == 1:
            for idx, xml_file in enumerate(xml_files, start=1):
                render_progress(idx, total_files, start_time)
                result = _process_xml_file(
                    pyMCDS,
                    constructor_mode,
                    xml_file,
                    domain_dims,
                    include_microenv,
                    include_spatial,
                    include_attribute,
                )

                if type_col_used is None and result["type_col"] is not None:
                    type_col_used = result["type_col"]
                if phase_col_used is None and result["phase_col"] is not None:
                    phase_col_used = result["phase_col"]

                summary_rows.extend(result["summary_rows"])
                type_rows.extend(result["type_rows"])
                phase_rows.extend(result["phase_rows"])
                type_phase_rows.extend(result["type_phase_rows"])
                micro_rows.extend(result["micro_rows"])
                spatial_rows.extend(result["spatial_rows"])
                region_rows.extend(result["region_rows"])
                attribute_rows.extend(result["attribute_rows"])
                density_rows.extend(result["density_rows"])
        else:
            max_workers = min(args.workers, total_files, (os.cpu_count() or 1))
            print(f"Parallel mode enabled with {max_workers} workers.")

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=_worker_init,
                initargs=(constructor_mode,),
            ) as executor:
                futures = {
                    executor.submit(
                        _worker_process,
                        str(xml_file),
                        domain_dims,
                        include_microenv,
                        include_spatial,
                        include_attribute,
                    ): xml_file
                    for xml_file in xml_files
                }

                completed = 0
                for future in concurrent.futures.as_completed(futures):
                    completed += 1
                    render_progress(completed, total_files, start_time)
                    result = future.result()

                    if type_col_used is None and result["type_col"] is not None:
                        type_col_used = result["type_col"]
                    if phase_col_used is None and result["phase_col"] is not None:
                        phase_col_used = result["phase_col"]

                    summary_rows.extend(result["summary_rows"])
                    type_rows.extend(result["type_rows"])
                    phase_rows.extend(result["phase_rows"])
                    type_phase_rows.extend(result["type_phase_rows"])
                    micro_rows.extend(result["micro_rows"])
                    spatial_rows.extend(result["spatial_rows"])
                    region_rows.extend(result["region_rows"])
                    attribute_rows.extend(result["attribute_rows"])
                    density_rows.extend(result["density_rows"])

        print("".ljust(80), end="\r")
        print("Processing complete.")
    finally:
        if temp_context:
            temp_context.cleanup()

    summary_df = pd.DataFrame(summary_rows)
    type_df = pd.DataFrame(type_rows)
    phase_df = pd.DataFrame(phase_rows)
    type_phase_df = pd.DataFrame(type_phase_rows)
    micro_df = pd.DataFrame(micro_rows)
    spatial_df = pd.DataFrame(spatial_rows)
    region_df = pd.DataFrame(region_rows)
    attribute_df = pd.DataFrame(attribute_rows)
    density_df = pd.DataFrame(density_rows)

    type_wide_df = pd.DataFrame()
    if not type_df.empty:
        type_wide_df = (
            type_df.pivot_table(index="time_min", columns="cell_type", values="count", aggfunc="sum")
            .fillna(0)
            .reset_index()
        )

    phase_wide_df = pd.DataFrame()
    if not phase_df.empty:
        phase_wide_df = (
            phase_df.pivot_table(index="time_min", columns="phase", values="count", aggfunc="sum")
            .fillna(0)
            .reset_index()
        )

    long_parts = []
    if not summary_df.empty:
        long_parts.append(
            summary_df[["time_min", "total_cells"]]
            .rename(columns={"total_cells": "value"})
            .assign(metric="total_cells", group="all")
            [["time_min", "metric", "group", "value"]]
        )
    if not type_df.empty:
        long_parts.append(
            type_df[["time_min", "cell_type", "count"]]
            .rename(columns={"cell_type": "group", "count": "value"})
            .assign(metric="cell_type")
            [["time_min", "metric", "group", "value"]]
        )
    if not phase_df.empty:
        long_parts.append(
            phase_df[["time_min", "phase", "count"]]
            .rename(columns={"phase": "group", "count": "value"})
            .assign(metric="cell_phase")
            [["time_min", "metric", "group", "value"]]
        )
    long_df = pd.concat(long_parts, ignore_index=True) if long_parts else pd.DataFrame(
        columns=["time_min", "metric", "group", "value"]
    )

    metadata_rows = [
        {"key": "type_column", "value": type_col_used or ""},
        {"key": "phase_column", "value": phase_col_used or ""},
        {"key": "domain_x", "value": args.domain_x or ""},
        {"key": "domain_y", "value": args.domain_y or ""},
        {"key": "domain_z", "value": args.domain_z or ""},
        {"key": "every_nth", "value": args.every_nth},
        {"key": "max_files", "value": args.max_files or ""},
        {"key": "workers", "value": args.workers},
        {"key": "fast_mode", "value": args.fast_mode},
        {"key": "constructor_mode", "value": constructor_mode},
    ]
    metadata_df = pd.DataFrame(metadata_rows)

    print("Writing CSV files...")
    summary_df.to_csv(out_dir / "summary_by_time.csv", index=False)
    density_df.to_csv(out_dir / "density_by_time.csv", index=False)
    type_df.to_csv(out_dir / "cell_type_counts_long.csv", index=False)
    phase_df.to_csv(out_dir / "cell_phase_counts_long.csv", index=False)
    type_phase_df.to_csv(out_dir / "cell_type_phase_counts.csv", index=False)
    type_wide_df.to_csv(out_dir / "cell_type_counts_wide.csv", index=False)
    phase_wide_df.to_csv(out_dir / "cell_phase_counts_wide.csv", index=False)
    spatial_df.to_csv(out_dir / "spatial_stats_by_type.csv", index=False)
    region_df.to_csv(out_dir / "region_counts.csv", index=False)
    attribute_df.to_csv(out_dir / "cell_attribute_stats.csv", index=False)
    micro_df.to_csv(out_dir / "microenvironment_stats.csv", index=False)
    long_df.to_csv(out_dir / "time_aligned_long.csv", index=False)
    metadata_df.to_csv(out_dir / "metadata.csv", index=False)

    print(f"Wrote CSV files to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()

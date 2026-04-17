"""
Add regimen-end-relative time columns to existing curve-fitter CSV outputs.

This script updates both:
- curve_fitter_table.csv
- curve_fitter_highlights.csv

inside an existing export folder, without rerunning export_outputs_to_excel.py.

Behavior:
- Adds/overwrites:
    regimen_end_time_min
    time_from_regimen_end_min = time_min - regimen_end_time_min
- If metadata.csv exists, updates/creates:
    regimen_end_minutes_input
    regimen_end_minutes_aligned

Alignment formula (matches export_outputs_to_excel.py):
    aligned = regimen_end_minutes_input - skip_initial_minutes + time_axis_applied_shift_min

If skip/shift keys are missing in metadata, defaults are used:
- skip_initial_minutes -> 0
- time_axis_applied_shift_min -> min(time_min) from curve_fitter_table.csv

Usage:
    python add_regimen_time_column.py --folder WP2S1_CSV --regimen-end-minutes 5760

If --regimen-end-minutes is omitted, the script tries metadata.csv keys in this order:
1) regimen_end_minutes_aligned (already aligned)
2) regimen_end_minutes_input (then aligns using skip/shift)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd


def _read_metadata_map(metadata_path: Path) -> dict[str, str]:
    if not metadata_path.exists():
        return {}

    meta_df = pd.read_csv(metadata_path)
    if "key" not in meta_df.columns or "value" not in meta_df.columns:
        return {}

    out: dict[str, str] = {}
    for _, row in meta_df.iterrows():
        key = str(row["key"]).strip()
        value = "" if pd.isna(row["value"]) else str(row["value"]).strip()
        if key:
            out[key] = value
    return out


def _to_float_or_none(value: str | None) -> Optional[float]:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _atomic_write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp_path, index=False)
    tmp_path.replace(path)


def _update_metadata(metadata_path: Path, updates: dict[str, float]) -> None:
    if metadata_path.exists():
        meta_df = pd.read_csv(metadata_path)
        if "key" not in meta_df.columns or "value" not in meta_df.columns:
            meta_df = pd.DataFrame(columns=["key", "value"])
    else:
        meta_df = pd.DataFrame(columns=["key", "value"])

    for key, value in updates.items():
        mask = meta_df["key"].astype(str) == key
        if mask.any():
            meta_df.loc[mask, "value"] = value
        else:
            meta_df = pd.concat(
                [meta_df, pd.DataFrame([{"key": key, "value": value}])],
                ignore_index=True,
            )

    _atomic_write_csv(meta_df, metadata_path)


def _aligned_regimen_end(
    regimen_end_input: Optional[float],
    metadata_map: dict[str, str],
    time_min_series: pd.Series,
) -> tuple[float, Optional[float]]:
    """
    Returns (aligned_regimen_end, regimen_end_input_or_none).
    """
    skip_initial = _to_float_or_none(metadata_map.get("skip_initial_minutes"))
    if skip_initial is None:
        skip_initial = 0.0

    applied_shift = _to_float_or_none(metadata_map.get("time_axis_applied_shift_min"))
    if applied_shift is None:
        # Fallback for older exports: infer from shifted timeline itself.
        applied_shift = float(pd.to_numeric(time_min_series, errors="coerce").min())

    if regimen_end_input is not None:
        if regimen_end_input < 0:
            raise ValueError("--regimen-end-minutes must be >= 0")
        aligned = float(regimen_end_input) - float(skip_initial) + float(applied_shift)
        return aligned, float(regimen_end_input)

    existing_aligned = _to_float_or_none(metadata_map.get("regimen_end_minutes_aligned"))
    if existing_aligned is not None:
        return float(existing_aligned), _to_float_or_none(metadata_map.get("regimen_end_minutes_input"))

    existing_input = _to_float_or_none(metadata_map.get("regimen_end_minutes_input"))
    if existing_input is not None:
        aligned = float(existing_input) - float(skip_initial) + float(applied_shift)
        return aligned, float(existing_input)

    raise ValueError(
        "No regimen end time available. Provide --regimen-end-minutes, "
        "or ensure metadata.csv has regimen_end_minutes_input/aligned."
    )


def _apply_columns(df: pd.DataFrame, aligned_regimen_end: float) -> pd.DataFrame:
    if "time_min" not in df.columns:
        raise ValueError("Required column 'time_min' not found.")

    out = df.copy()
    time_vals = pd.to_numeric(out["time_min"], errors="coerce")
    out["regimen_end_time_min"] = float(aligned_regimen_end)
    out["time_from_regimen_end_min"] = time_vals - float(aligned_regimen_end)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add regimen-end-relative time columns to curve_fitter CSV files in an existing export folder."
    )
    parser.add_argument("-help", action="help", help="Show this help message and exit.")
    parser.add_argument("--folder", required=True, help="Folder containing curve_fitter_table.csv and curve_fitter_highlights.csv")
    parser.add_argument(
        "--regimen-end-minutes",
        type=float,
        default=None,
        help="Regimen end time in minutes on the original simulation timeline.",
    )
    args = parser.parse_args()

    folder = Path(args.folder)
    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder}")

    table_path = folder / "curve_fitter_table.csv"
    highlights_path = folder / "curve_fitter_highlights.csv"
    metadata_path = folder / "metadata.csv"

    if not table_path.exists():
        raise FileNotFoundError(f"Missing required file: {table_path}")
    if not highlights_path.exists():
        raise FileNotFoundError(f"Missing required file: {highlights_path}")

    table_df = pd.read_csv(table_path)
    highlights_df = pd.read_csv(highlights_path)
    metadata_map = _read_metadata_map(metadata_path)

    if "time_min" not in table_df.columns:
        raise ValueError(f"'time_min' missing in {table_path.name}")

    aligned_regimen_end, regimen_end_input = _aligned_regimen_end(
        args.regimen_end_minutes,
        metadata_map,
        table_df["time_min"],
    )

    table_df_out = _apply_columns(table_df, aligned_regimen_end)
    highlights_df_out = _apply_columns(highlights_df, aligned_regimen_end)

    _atomic_write_csv(table_df_out, table_path)
    _atomic_write_csv(highlights_df_out, highlights_path)

    updates: dict[str, float] = {
        "regimen_end_minutes_aligned": float(aligned_regimen_end),
    }
    if regimen_end_input is not None:
        updates["regimen_end_minutes_input"] = float(regimen_end_input)

    _update_metadata(metadata_path, updates)

    print(f"Updated: {table_path}")
    print(f"Updated: {highlights_path}")
    print(f"Updated: {metadata_path}")
    print(f"regimen_end_minutes_aligned = {aligned_regimen_end:.10g}")


if __name__ == "__main__":
    main()

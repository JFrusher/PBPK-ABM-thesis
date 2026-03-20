"""Merge chunk-mode Monte Carlo CSV outputs into one combined dataset.

Expected input files (recursive under --root):
- MC_chunk_outputs_*.csv

Usage:
    python merge_MC_chunk_results.py --root MC_results/MC10k_chunks
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge MC chunk output CSV files.")
    parser.add_argument("--root", required=True, help="Root folder containing per-task chunk outputs")
    parser.add_argument(
        "--out",
        default=None,
        help="Optional output CSV path (default: <root>/MC_10k_merged_outputs.csv)",
    )
    args = parser.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        raise FileNotFoundError(f"Root folder not found: {root}")

    csv_files = sorted(root.rglob("MC_chunk_outputs_*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No MC_chunk_outputs_*.csv files found under: {root}")

    frames: list[pd.DataFrame] = []
    for p in csv_files:
        df = pd.read_csv(p)
        df["chunk_file"] = str(p)
        # Try to recover task id from parent folder naming task_00001
        parent = p.parent.name
        if parent.startswith("task_"):
            try:
                df["chunk_task_id"] = int(parent.split("_", 1)[1])
            except ValueError:
                df["chunk_task_id"] = pd.NA
        else:
            df["chunk_task_id"] = pd.NA
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)

    # Stable ordering where task id is available
    if "chunk_task_id" in merged.columns:
        merged = merged.sort_values(["chunk_task_id", "sample_local_index"], na_position="last").reset_index(drop=True)

    out_path = Path(args.out) if args.out else (root / "MC_10k_merged_outputs.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, index=False)

    summary_path = out_path.with_suffix(".summary.txt")
    with summary_path.open("w", encoding="utf-8") as f:
        f.write(f"Merged files: {len(csv_files)}\n")
        f.write(f"Total rows: {len(merged)}\n")
        if "AUC_mg_h_L" in merged.columns:
            f.write(f"Mean AUC_mg_h_L: {pd.to_numeric(merged['AUC_mg_h_L'], errors='coerce').mean()}\n")
        if "Cmax_uM" in merged.columns:
            f.write(f"Mean Cmax_uM: {pd.to_numeric(merged['Cmax_uM'], errors='coerce').mean()}\n")

    print(f"Merged {len(csv_files)} chunk files")
    print(f"Wrote: {out_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()

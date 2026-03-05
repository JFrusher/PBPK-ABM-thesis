import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def load_cells(csv_path: Path):
    by_type = defaultdict(lambda: {"x": [], "y": []})

    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        required_cols = {"x", "y", "type"}
        if reader.fieldnames is None or not required_cols.issubset(set(reader.fieldnames)):
            raise ValueError("CSV must contain headers: x,y,type (z is optional)")

        for row in reader:
            cell_type = row["type"]
            by_type[cell_type]["x"].append(float(row["x"]))
            by_type[cell_type]["y"].append(float(row["y"]))

    return dict(by_type)


def compute_limits(by_type):
    all_x = []
    all_y = []
    for coords in by_type.values():
        all_x.extend(coords["x"])
        all_y.extend(coords["y"])

    if not all_x or not all_y:
        return (-1, 1), (-1, 1)

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)
    span = max(max_x - min_x, max_y - min_y)
    margin = 0.05 * span if span > 0 else 1.0

    center_x = 0.5 * (min_x + max_x)
    center_y = 0.5 * (min_y + max_y)
    half = 0.5 * span + margin

    return (center_x - half, center_x + half), (center_y - half, center_y + half)


def plot_overall(by_type, out_path: Path, point_size: float):
    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)

    for cell_type in sorted(by_type):
        coords = by_type[cell_type]
        ax.scatter(coords["x"], coords["y"], s=point_size, alpha=0.65, label=f"{cell_type} (n={len(coords['x'])})")

    (xmin, xmax), (ymin, ymax) = compute_limits(by_type)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Cells: overall distribution")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right", fontsize=8, frameon=True)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_per_type(by_type, out_path: Path, point_size: float):
    types = sorted(by_type)
    n = len(types)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 4 * rows), dpi=150)
    axes = [axes] if not hasattr(axes, "flatten") else axes.flatten()

    (xmin, xmax), (ymin, ymax) = compute_limits(by_type)

    for i, cell_type in enumerate(types):
        ax = axes[i]
        coords = by_type[cell_type]
        ax.scatter(coords["x"], coords["y"], s=point_size, alpha=0.75)
        ax.set_title(f"{cell_type} (n={len(coords['x'])})", fontsize=10)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.2)

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Cells: per-type distribution", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Quick visualization of cell CSV coordinates (overall + per type).")
    parser.add_argument("--csv", default="cells.csv", help="Path to cell CSV file (default: cells.csv)")
    parser.add_argument("--out-prefix", default="cells_plot", help="Output prefix for PNGs")
    parser.add_argument("--point-size", type=float, default=4.0, help="Scatter point size")
    parser.add_argument("--show", action="store_true", help="Also show figures interactively")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    by_type = load_cells(csv_path)
    total_cells = sum(len(v["x"]) for v in by_type.values())
    if total_cells == 0:
        raise ValueError("No cells found in CSV")

    overall_path = Path(f"{args.out_prefix}_overall.png")
    per_type_path = Path(f"{args.out_prefix}_per_type.png")

    plot_overall(by_type, overall_path, args.point_size)
    plot_per_type(by_type, per_type_path, args.point_size)

    print(f"Loaded {total_cells} cells across {len(by_type)} types from {csv_path}")
    print(f"Saved overall plot: {overall_path}")
    print(f"Saved per-type plot: {per_type_path}")

    if args.show:
        plt.figure(figsize=(6, 1))
        plt.axis("off")
        plt.title("Plots saved. Re-run with your preferred viewer if needed.")
        plt.show()


if __name__ == "__main__":
    main()

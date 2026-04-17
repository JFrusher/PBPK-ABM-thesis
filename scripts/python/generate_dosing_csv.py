import argparse
import csv
import math
from pathlib import Path


HEADER = [
    "start_time_min",
    "end_time_min",
    "dosing_type",
    "dose_amount",
    "infusion_rate",
    "amplitude",
    "frequency_per_min",
    "mean_rate",
    "custom_function",
    "custom_params",
]


def build_eval_context(t: float) -> dict:
    context = {name: getattr(math, name) for name in dir(math) if not name.startswith("_")}
    context["t"] = t
    return context


def evaluate_rate(function_expr: str, t: float) -> float:
    context = build_eval_context(t)
    value = eval(function_expr, {"__builtins__": {}}, context)
    value = float(value)
    return max(value, 0.0)


def generate_rows(function_expr: str, total_time_min: int, chunk_min: int):
    if chunk_min <= 0:
        raise ValueError("chunk_min must be > 0")
    if total_time_min < 0:
        raise ValueError("total_time_min must be >= 0")

    rows = []
    for start_t in range(0, total_time_min, chunk_min):
        end_t = min(start_t + chunk_min, total_time_min)
        infusion_rate = evaluate_rate(function_expr, float(start_t))
        rows.append([
            start_t,
            end_t,
            "continuous",
            0,
            infusion_rate,
            0,
            0,
            0,
            "",
            "",
        ])
    return rows


def write_csv(rows, output_path: Path):
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(HEADER)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate dosing CSV from a function of time t (minutes). "
            "Example function: '0.02*t + 1.5*sin(t/120) + 5'"
        )
    )
    parser.add_argument("--function", required=True, help="Function expression for infusion_rate X using variable t (minutes)")
    parser.add_argument("--total-time-min", type=int, required=True, help="Total simulation time in minutes")
    parser.add_argument("--chunk-min", type=int, required=True, help="Chunk size in minutes (e.g., 30)")
    parser.add_argument("--out", default="generated_dosing.csv", help="Output CSV path")

    args = parser.parse_args()

    rows = generate_rows(args.function, args.total_time_min, args.chunk_min)
    out_path = Path(args.out)
    write_csv(rows, out_path)

    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()

import csv
import math
import random


# (min_r, max_r, CSC, CP, CD, CD8, M2, CAF)
ring_data = [
    (0, 50, 0, 0, 2, 0, 0, 0),
    (50, 100, 4, 0, 17, 0, 1, 0),
    (100, 150, 21, 0, 56, 0, 4, 0),
    (150, 200, 53, 0, 125, 0, 9, 0),
    (200, 250, 88, 0, 231, 0, 16, 0),
    (250, 300, 103, 0, 374, 0, 28, 0),
    (300, 350, 87, 2, 552, 0, 43, 0),
    (350, 400, 54, 6, 756, 0, 61, 0),
    (400, 450, 25, 15, 971, 0, 84, 0),
    (450, 500, 9, 34, 1179, 0, 110, 0),
    (500, 550, 2, 70, 1362, 0, 139, 0),
    (550, 600, 0, 130, 1503, 0, 170, 0),
    (600, 650, 0, 223, 1589, 1, 201, 0),
    (650, 700, 0, 356, 1613, 4, 230, 1),
    (700, 750, 0, 527, 1575, 12, 258, 3),
    (750, 800, 1, 727, 1482, 31, 281, 7),
    (800, 850, 26, 934, 1345, 68, 299, 15),
    (850, 900, 224, 1118, 1178, 135, 311, 30),
    (900, 950, 681, 1249, 997, 235, 317, 55),
    (950, 1000, 756, 1304, 816, 365, 315, 95),
    (1000, 1050, 307, 1272, 646, 503, 308, 154),
    (1050, 1100, 45, 1160, 495, 619, 294, 231),
    (1100, 1150, 2, 989, 367, 678, 275, 325),
    (1150, 1200, 0, 789, 264, 662, 253, 428),
    (1200, 1250, 0, 589, 184, 576, 228, 527),
    (1250, 1300, 0, 412, 124, 447, 201, 608),
    (1300, 1350, 0, 270, 81, 309, 175, 657),
    (1350, 1400, 0, 165, 51, 191, 149, 665),
    (1400, 1450, 0, 95, 31, 105, 124, 630),
    (1450, 1500, 0, 51, 19, 51, 102, 560),
]


cell_type_map = {
    "CSC": "Cancer_stem",
    "CP": "cancer_proliferating",
    "CD": "cancer_differentiated",
    "CD8": "CD8_T_cell",
    "M2": "M2_macrophage",
    "CAF": "CAF",
}


def sample_point_in_annulus(min_radius, max_radius):
    theta = random.uniform(0.0, 2.0 * math.pi)
    radius = math.sqrt(random.uniform(min_radius * min_radius, max_radius * max_radius))
    x = radius * math.cos(theta)
    y = radius * math.sin(theta)
    return x, y


def generate_cells(data):
    cells = []

    for min_r, max_r, csc, cp, cd, cd8, m2, caf in data:
        counts = {
            "CSC": csc,
            "CP": cp,
            "CD": cd,
            "CD8": cd8,
            "M2": m2,
            "CAF": caf,
        }

        for shorthand, count in counts.items():
            full_type = cell_type_map[shorthand]
            for _ in range(count):
                x, y = sample_point_in_annulus(min_r, max_r)
                cells.append((x, y, 0.0, full_type))

    return cells


def write_cells_csv(cells, output_file="cells.csv"):
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "z", "type"])
        writer.writerows(cells)


if __name__ == "__main__":
    generated_cells = generate_cells(ring_data)
    write_cells_csv(generated_cells, "cellsX.csv")
    print(f"Wrote {len(generated_cells)} cells to cellsX.csv")

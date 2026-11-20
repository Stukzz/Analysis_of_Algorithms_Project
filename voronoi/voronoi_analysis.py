import argparse
import csv
import time
from collections import deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon
import numpy as np
from scipy.spatial import QhullError, Voronoi
from shapely.geometry import MultiPoint, Polygon


def sky_to_plane(ra_deg, dec_deg):

    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    ra0 = np.arctan2(np.mean(np.sin(ra)), np.mean(np.cos(ra)))
    dec0 = float(np.mean(dec))

    delta_ra = (ra - ra0 + np.pi) % (2 * np.pi) - np.pi
    cos_c = (
        np.sin(dec0) * np.sin(dec)
        + np.cos(dec0) * np.cos(dec) * np.cos(delta_ra)
    )
    if np.any(cos_c <= 0):
        raise ValueError("Field extends beyond the tangent plane's visible hemisphere")

    x = np.cos(dec) * np.sin(delta_ra) / cos_c
    y = (
        np.cos(dec0) * np.sin(dec)
        - np.sin(dec0) * np.cos(dec) * np.cos(delta_ra)
    ) / cos_c
    return np.column_stack((np.rad2deg(x), np.rad2deg(y)))


def load_coordinates(filename):

    rows = []
    with open(filename, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"alpha_j2000", "delta_j2000"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                "CSV must contain columns named alpha_j2000 and delta_j2000. "
                f"Found: {reader.fieldnames}"
            )
        for row in reader:
            try:
                ra = float(row["alpha_j2000"])
                dec = float(row["delta_j2000"])
            except (ValueError, TypeError):
                continue
            if np.isfinite(ra) and np.isfinite(dec) and 0 <= ra < 360 and -90 <= dec <= 90:
                rows.append((ra, dec))

    if len(rows) < 4:
        raise ValueError("Need at least four valid coordinate rows")
    unique_rows = np.unique(np.asarray(rows, dtype=float), axis=0)
    if len(unique_rows) < 4:
        raise ValueError("Need at least four distinct coordinate rows")
    return unique_rows


def demo_coordinates(n=500, seed=42):
    
    rng = np.random.default_rng(seed)
    n1, n2 = n // 4, n // 4
    background = rng.uniform([-0.85, -0.65], [0.85, 0.65], size=(n - n1 - n2, 2))
    first = rng.normal([-0.35, 0.14], [0.10, 0.09], size=(n1, 2))
    second = rng.normal([0.4, -0.18], [0.12, 0.10], size=(n2, 2))
    plane = np.vstack([background, first, second])
    return np.column_stack([
        49.9 + plane[:, 0] / np.cos(np.deg2rad(41.5)),
        41.5 + plane[:, 1],
    ])


def finite_polygons(vor, radius):

    center = vor.points.mean(axis=0)
    ridge_map = [[] for _ in vor.points]
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        ridge_map[p1].append((p2, v1, v2))
        ridge_map[p2].append((p1, v1, v2))

    vertices = vor.vertices.tolist()
    regions = []

    for point_index, region_index in enumerate(vor.point_region):
        region = vor.regions[region_index]
        if region and all(vertex >= 0 for vertex in region):
            regions.append(region)
            continue

        new_region = [vertex for vertex in region if vertex >= 0]
        for other_index, v1, v2 in ridge_map[point_index]:
            if v1 >= 0 and v2 >= 0:
                continue

            finite_vertex = v1 if v1 >= 0 else v2
            tangent = vor.points[other_index] - vor.points[point_index]
            tangent_length = np.linalg.norm(tangent)
            if tangent_length == 0:
                continue
            tangent = tangent / tangent_length
            normal = np.array([-tangent[1], tangent[0]])
            midpoint = (vor.points[point_index] + vor.points[other_index]) / 2
            direction = np.sign(np.dot(midpoint - center, normal)) * normal
            if np.allclose(direction, 0):
                direction = normal

            vertices.append((vor.vertices[finite_vertex] + direction * radius).tolist())
            new_region.append(len(vertices) - 1)

        region_coords = np.asarray([vertices[index] for index in new_region])
        if len(region_coords) < 3:
            raise ValueError("Could not reconstruct a Voronoi region")
        angles = np.arctan2(
            region_coords[:, 1] - region_coords[:, 1].mean(),
            region_coords[:, 0] - region_coords[:, 0].mean(),
        )
        regions.append([new_region[index] for index in np.argsort(angles)])

    return regions, np.asarray(vertices, dtype=float)


def connected_components(dense_mask, neighbors):
    labels = np.full(len(dense_mask), -1, dtype=int)
    component_sizes = []
    component_id = 0

    for start_index in np.flatnonzero(dense_mask):
        if labels[start_index] != -1:
            continue

        queue = deque([int(start_index)])
        labels[start_index] = component_id
        size = 0

        while queue:
            current = queue.popleft()
            size += 1
            for neighbor in neighbors[current]:
                if dense_mask[neighbor] and labels[neighbor] == -1:
                    labels[neighbor] = component_id
                    queue.append(neighbor)

        component_sizes.append(size)
        component_id += 1

    return labels, component_sizes


def analyze(points, density_quantile=0.85, min_cluster_size=5):

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("Points must be an N x 2 array")
    if len(points) < 4 or len(np.unique(points, axis=0)) != len(points):
        raise ValueError("Need at least four distinct points")
    if not np.all(np.isfinite(points)):
        raise ValueError("Coordinates must all be finite")
    if np.linalg.matrix_rank(points - points.mean(axis=0)) < 2:
        raise ValueError("Points must not be collinear")
    if not 0.5 < density_quantile < 1.0:
        raise ValueError("density_quantile must be between 0.5 and 1.0")
    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be at least 1")

    start_time = time.perf_counter()
    vor = Voronoi(points)

    # Use the point-cloud convex hull as a first-pass study boundary.
    # This is better than a padded rectangle, but is not a substitute for
    # the official Euclid survey mask/footprint.
    study_boundary = MultiPoint(points).convex_hull
    if study_boundary.geom_type != "Polygon" or study_boundary.area <= 0:
        raise ValueError("Could not construct a two-dimensional convex-hull boundary")

    span = np.ptp(points, axis=0)
    extension_radius = 10 * max(float(np.max(span)), 1e-6)
    regions, vertices = finite_polygons(vor, radius=extension_radius)

    polygons = []
    for region in regions:
        polygon = Polygon(vertices[region])
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        clipped = polygon.intersection(study_boundary)
        polygons.append(clipped)

    areas = np.asarray([polygon.area for polygon in polygons], dtype=float)
    if np.any(~np.isfinite(areas)) or np.any(areas <= 0):
        raise ValueError("A Voronoi cell has zero/invalid area; inspect input coordinates")

    density = 1.0 / areas


    edge_tolerance = max(study_boundary.length, 1.0) * 1e-10
    boundary_cells = np.asarray([
        polygon.boundary.intersection(study_boundary.boundary).length > edge_tolerance
        for polygon in polygons
    ], dtype=bool)
    interior_cells = ~boundary_cells

    if np.count_nonzero(interior_cells) < max(10, min_cluster_size):
        raise ValueError(
            "Too few non-boundary cells for reliable density analysis. "
            "Use more source rows or inspect the survey footprint."
        )

    threshold = float(np.quantile(density[interior_cells], density_quantile))
    dense = interior_cells & (density >= threshold)

    neighbors = [set() for _ in points]
    for first, second in vor.ridge_points:
        neighbors[first].add(int(second))
        neighbors[second].add(int(first))

    raw_labels, component_sizes = connected_components(dense, neighbors)
    candidate_labels = np.full(len(points), -1, dtype=int)
    candidate_group_count = 0
    cluster_summaries = []

    for raw_id, size in enumerate(component_sizes):
        if size < min_cluster_size:
            continue
        member_indices = np.flatnonzero(raw_labels == raw_id)
        candidate_labels[member_indices] = candidate_group_count
        member_density = density[member_indices]
        center_x = float(np.mean(points[member_indices, 0]))
        center_y = float(np.mean(points[member_indices, 1]))
        cluster_summaries.append({
            "candidate_group_id": candidate_group_count,
            "object_count": int(size),
            "mean_density": float(np.mean(member_density)),
            "max_density": float(np.max(member_density)),
            "center_x": center_x,
            "center_y": center_y,
        })
        candidate_group_count += 1

    return {
        "vor": vor,
        "polygons": polygons,
        "areas": areas,
        "density": density,
        "boundary_cells": boundary_cells,
        "interior_cells": interior_cells,
        "dense": dense,
        "raw_labels": raw_labels,
        "labels": candidate_labels,
        "component_sizes": component_sizes,
        "cluster_summaries": cluster_summaries,
        "candidate_group_count": candidate_group_count,
        "threshold": threshold,
        "seconds": time.perf_counter() - start_time,
        "study_boundary": study_boundary,
        "density_quantile": density_quantile,
        "min_cluster_size": min_cluster_size,
    }


def save_outputs(points, result, output):
    """Save the density map, candidate-group map, and CSV summaries."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    density_csv = output / "density_results.csv"
    with density_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "x_tangent_deg", "y_tangent_deg", "cell_area_sq_deg",
            "local_density_per_sq_deg", "is_boundary_cell", "is_high_density",
            "connected_component_id", "candidate_group_id",
        ])
        for index, (point, area, density) in enumerate(
            zip(points, result["areas"], result["density"])
        ):
            writer.writerow([
                point[0], point[1], area, density,
                int(result["boundary_cells"][index]),
                int(result["dense"][index]),
                int(result["raw_labels"][index]),
                int(result["labels"][index]),
            ])

    summary_csv = output / "candidate_clusters.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "candidate_group_id", "object_count", "mean_density_per_sq_deg",
            "max_density_per_sq_deg", "center_x_tangent_deg", "center_y_tangent_deg",
        ])
        for summary in result["cluster_summaries"]:
            writer.writerow([
                summary["candidate_group_id"], summary["object_count"],
                summary["mean_density"], summary["max_density"],
                summary["center_x"], summary["center_y"],
            ])

    fig, ax = plt.subplots(figsize=(10, 8))
    patches = [
        MplPolygon(np.asarray(poly.exterior.coords), closed=True)
        for poly in result["polygons"]
    ]
    collection = PatchCollection(
        patches, cmap="viridis", edgecolor="#555555", linewidth=0.12
    )
    collection.set_array(np.log10(result["density"]))
    ax.add_collection(collection)
    fig.colorbar(collection, ax=ax, label="log10(local density / square degree)")

    edge_points = result["boundary_cells"]
    ax.scatter(points[edge_points, 0], points[edge_points, 1],
               s=5, c="lightgray", alpha=0.75, label="Boundary cells (excluded)")
    ax.scatter(points[~edge_points, 0], points[~edge_points, 1],
               s=2, c="black", alpha=0.45)
    ax.scatter(points[result["dense"], 0], points[result["dense"], 1],
               s=13, facecolors="none", edgecolors="red", linewidths=0.7,
               label=f"Top {(1 - result['density_quantile']) * 100:.0f}% interior density")

    hull_x, hull_y = result["study_boundary"].exterior.xy
    ax.plot(hull_x, hull_y, color="white", linewidth=1.1, linestyle="--",
            label="Convex-hull sample boundary")
    ax.set(
        xlabel="Tangent-plane X (degrees)",
        ylabel="Tangent-plane Y (degrees)",
        title="Voronoi local surface-density estimate",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "voronoi_density.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(points[:, 0], points[:, 1], s=4, c="lightgray",
               alpha=0.65, label="Other sources")

    summaries = result["cluster_summaries"]
    color_map = plt.get_cmap("tab20")
    for summary in summaries:
        group_id = summary["candidate_group_id"]
        members = result["labels"] == group_id
        color = color_map(group_id % 20)
        ax.scatter(points[members, 0], points[members, 1], s=18,
                   color=[color], label=f"Group {group_id} (n={summary['object_count']})")
        ax.scatter(summary["center_x"], summary["center_y"], marker="+",
                   s=90, linewidths=1.5, color="black")
        if len(summaries) <= 15:
            ax.annotate(str(group_id), (summary["center_x"], summary["center_y"]),
                        xytext=(4, 4), textcoords="offset points", fontsize=8)

    hull_x, hull_y = result["study_boundary"].exterior.xy
    ax.plot(hull_x, hull_y, color="black", linewidth=1, linestyle="--")
    ax.set(
        xlabel="Tangent-plane X (degrees)",
        ylabel="Tangent-plane Y (degrees)",
        title=f"Candidate connected high-density groups (minimum {result['min_cluster_size']} cells)",
    )
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    if 0 < len(summaries) <= 15:
        ax.legend(loc="best", fontsize=7)
    elif len(summaries) > 15:
        ax.text(0.02, 0.98,
                f"{len(summaries)} candidate groups; colors repeat after 20 groups\nBlack + = group center",
                transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
    else:
        ax.text(0.02, 0.98, "No groups met the minimum-size criterion",
                transform=ax.transAxes, va="top", fontsize=9)
    fig.tight_layout()
    fig.savefig(output / "candidate_clusters.png", dpi=180)
    plt.close(fig)


def benchmark(
    points,
    output,
    sizes=(100, 250, 500, 1000, 2000, 5000, 10000, 25000, 50000),
    repeats=3,
    seed=42,
):

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    sizes = sorted({int(size) for size in sizes if 4 <= int(size) <= len(points)})
    if repeats < 1:
        raise ValueError("Benchmark repeats must be at least 1")

    rng = np.random.default_rng(seed)
    shuffled_points = np.asarray(points)[rng.permutation(len(points))]
    rows = []
    for size in sizes:
        elapsed_times = []
        sample = shuffled_points[:size]
        for _ in range(repeats):
            try:
                elapsed_times.append(analyze(sample)["seconds"])
            except (ValueError, QhullError):
                elapsed_times = []
                break
        if elapsed_times:
            rows.append((
                size,
                float(np.median(elapsed_times)),
                float(np.min(elapsed_times)),
                float(np.max(elapsed_times)),
                repeats,
            ))

    with (output / "benchmark.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "number_of_points", "median_seconds", "min_seconds",
            "max_seconds", "repeats",
        ])
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true",
                        help="Use synthetic data, not actual Euclid observations")
    source.add_argument("--csv", help="CSV with alpha_j2000 and delta_j2000 headers")
    parser.add_argument("--output", default="results", help="Output directory")
    parser.add_argument("--limit", type=int, default=2000,
                        help="Maximum number of coordinate rows to analyze")
    parser.add_argument("--density-quantile", type=float, default=0.85,
                        help="Density quantile for candidate cells (default: 0.85)")
    parser.add_argument("--min-cluster-size", type=int, default=5,
                        help="Minimum connected high-density cells in a candidate group")
    parser.add_argument("--benchmark-repeats", type=int, default=3,
                        help="Number of timing repeats per benchmark size (default: 3)")
    args = parser.parse_args()

    if args.limit < 4:
        parser.error("--limit must be at least 4")

    try:
        coords = demo_coordinates() if args.demo else load_coordinates(args.csv)
        coords = coords[:args.limit]
        if len(coords) < 4:
            raise ValueError("At least four valid coordinates are needed")

        points = sky_to_plane(coords[:, 0], coords[:, 1])
        result = analyze(
            points,
            density_quantile=args.density_quantile,
            min_cluster_size=args.min_cluster_size,
        )
        save_outputs(points, result, args.output)
        measured = benchmark(points, args.output, repeats=args.benchmark_repeats)

    except (ValueError, QhullError, OSError) as error:
        parser.error(str(error))

    print(f"Processed {len(points)} unique celestial positions in {result['seconds']:.3f}s")
    print(f"Boundary cells excluded from candidate selection: {int(result['boundary_cells'].sum())}")
    print(f"Interior cells used for density threshold: {int(result['interior_cells'].sum())}")
    print(f"High-density candidate cells: {int(result['dense'].sum())}")
    print(f"Candidate connected groups (minimum size {args.min_cluster_size}): {result['candidate_group_count']}")
    print(f"Saved outputs to {Path(args.output).resolve()}")
    print("Benchmarks (points, median seconds, min seconds, max seconds, repeats):")
    for row in measured:
        print("  ", row)
    if args.demo:
        print("NOTE: DEMO USES SYNTHETIC DATA, NOT OBSERVED GALAXIES")
    print("NOTE: Candidate groups are projected density concentrations, not confirmed physical galaxy clusters.")
    print("NOTE: The convex hull approximates the sample boundary; use the official survey footprint for scientific analysis.")


if __name__ == "__main__":
    main()

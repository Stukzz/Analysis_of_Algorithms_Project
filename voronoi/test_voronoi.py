import tempfile
import unittest
from pathlib import Path

import numpy as np

from voronoi_analysis import (
    analyze,
    connected_components,
    demo_coordinates,
    save_outputs,
    sky_to_plane,
    benchmark,
)


class VoronoiTests(unittest.TestCase):
    def test_projection_center(self):
        points = sky_to_plane([50, 50, 50], [41, 41, 41])
        self.assertTrue(np.allclose(points, 0, atol=1e-10))

    def test_area_and_density(self):
        points = sky_to_plane(*demo_coordinates(150).T)
        result = analyze(points)
        self.assertTrue(np.all(result["areas"] > 0))
        self.assertTrue(np.allclose(result["density"] * result["areas"], 1))
        self.assertAlmostEqual(
            sum(result["areas"]), result["study_boundary"].area, places=7
        )

    def test_boundary_cells_not_selected_as_dense(self):
        points = sky_to_plane(*demo_coordinates(150).T)
        result = analyze(points)
        self.assertTrue(np.all(~result["dense"][result["boundary_cells"]]))

    def test_candidate_groups_meet_minimum_size(self):
        points = sky_to_plane(*demo_coordinates(180).T)
        result = analyze(points, min_cluster_size=5)
        self.assertTrue(np.all(result["labels"][~result["dense"]] == -1))
        for summary in result["cluster_summaries"]:
            self.assertGreaterEqual(summary["object_count"], 5)
        self.assertEqual(result["candidate_group_count"], len(result["cluster_summaries"]))

    def test_connected_components(self):
        dense = np.array([True, True, False, True, True, True])
        neighbors = [{1}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4}]
        labels, sizes = connected_components(dense, neighbors)
        self.assertEqual(sizes, [2, 3])
        self.assertEqual(labels[2], -1)

    def test_duplicate_points_rejected(self):
        with self.assertRaises(ValueError):
            analyze(np.array([[0, 0], [1, 0], [0, 1], [0, 0]]))


    def test_benchmark_writes_repeated_timing_summary(self):
        points = sky_to_plane(*demo_coordinates(120).T)
        with tempfile.TemporaryDirectory() as folder:
            rows = benchmark(points, folder, sizes=(50, 100), repeats=2)
            self.assertEqual([row[0] for row in rows], [50, 100])
            self.assertTrue((Path(folder) / "benchmark.csv").exists())
            with (Path(folder) / "benchmark.csv").open(encoding="utf-8") as handle:
                header = handle.readline().strip()
            self.assertEqual(
                header,
                "number_of_points,median_seconds,min_seconds,max_seconds,repeats",
            )
            self.assertTrue(all(row[1] >= 0 for row in rows))

    def test_output_files_created(self):
        points = sky_to_plane(*demo_coordinates(150).T)
        result = analyze(points)
        with tempfile.TemporaryDirectory() as folder:
            save_outputs(points, result, folder)
            expected = [
                "density_results.csv",
                "candidate_clusters.csv",
                "voronoi_density.png",
                "candidate_clusters.png",
            ]
            for name in expected:
                self.assertTrue((Path(folder) / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
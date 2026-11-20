cat > README.md <<'EOF'
# Voronoi Tessellation

This folder contains the Python implementation of Voronoi Tessellation for analyzing galaxy distributions in the Euclid Perseus cluster dataset.

## Requirements

- Python 3
- NumPy
- SciPy
- Matplotlib
- Pandas

## 1. Setup

Create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the required packages:

```bash
pip install -r requirements.txt
```

## 2. Run the Program

Run the analysis using the included sample dataset:

```bash
python3 voronoi_analysis.py --csv demo.csv --limit 2000 --output results_2000
```

The program analyzes galaxy coordinates and generates Voronoi tessellation and density-based cluster analysis results.

## 3. Run Tests

```bash
python3 -m unittest -v
```

## 4. Run Performance Benchmark

```bash
python3 voronoi_analysis.py --csv demo.csv --limit 2000 --benchmark-repeats 3 --output results_2000
```

## Notes

- `demo.csv` contains sample Euclid galaxy coordinate data.
- The CSV must contain `alpha_j2000` and `delta_j2000` columns.
- Results are saved in the folder specified by `--output`.
EOF
# 🌌 Spatial Graph Analysis of the Perseus Galaxy Cluster
**Analysis of Algorithms Project**

## 📌 Motivation
Modern astronomical surveys generate massive, high-dimensional datasets that require scalable computational solutions. This project explores the algorithmic challenges of processing large-scale astrophysical data—specifically focusing on spatial partitioning and graph construction across more than 1.5 million celestial records to uncover the underlying macro-structures of the universe.

---

## 🧮 Algorithmic Approach
To analyze the spatial distribution of celestial objects within the **Perseus Galaxy Cluster (Abell 426)**, we construct Euclidean geometric graphs mapped over J2000 celestial coordinates—Right Ascension (`alpha_j2000`) and Declination (`delta_j2000`):

* **Minimum Spanning Tree (MST):** Connects celestial sources into a single, cycle-free network that minimizes total edge distance. In computational astrophysics, MSTs are crucial for identifying **cosmic filamentary networks** and isolating dense galaxy groupings from background noise.
* **Voronoi Tessellation:** Partitions the celestial plane into convex polygons based on source proximity, where each polygon contains a single astronomical object. Calculating the area of these Voronoi cells provides a direct, non-parametric measure of **local environmental density** across the cluster.

---

## 🔭 Dataset Overview
The data is sourced from the **Euclid Space Telescope’s Early Release Observations (ERO-10)** targeting the Perseus Galaxy Cluster. We utilize two primary source extraction catalogs suited for spatial and graph algorithms:

| Catalog Name | Instrument & Spectral Bands | Volume (Rows) | Primary Characteristics |
| :--- | :--- | :---: | :--- |
| `euclid_ero.perseus_nisp` | Near-Infrared Spectrometer and Photometer ($Y$, $J$, $H$ bands) | **1,006,020** | High completeness across multi-band near-infrared wavelengths |
| `euclid_ero.perseus_vis` | Visible Instrument ($I_E$ band) | **546,563** | Superior spatial resolution for precise morphological and positional mapping |

### Key Extracted Features
* **Coordinates:** `alpha_j2000` (Right Ascension), `delta_j2000` (Declination)
* **Spatial Indices:** `nest4096` (HEALPix index), `htm9` (Hierarchical Triangular Mesh)
* **Photometry & Morphology:** `mag_auto`, `magerr_auto`, `flux_auto`, `class_star`, `spread_model`, `fwhm_image`

---

## 🗄️ Data Acquisition (NOIRLab Astro Data Lab)
The raw catalog files can be retrieved directly from the [NOIRLab Astro Data Lab Data Explorer](https://datalab.noirlab.edu/data-explorer?showTable=euclid_ero.perseus_vis) using the following SQL queries. Both queries filter out sentinel error values by constraining `mag_auto < 90`.

### 1. NISP Catalog Query
*(Can be executed as-is to pull all three infrared bands, or filtered to a single band such as `WHERE band = 'J' AND mag_auto < 90`)*

```sql
SELECT 
    number, band, alpha_j2000, delta_j2000, 
    nest4096, htm9, mag_auto, magerr_auto, 
    flux_auto, class_star, spread_model, fwhm_image 
FROM euclid_ero.perseus_nisp 
WHERE mag_auto < 90;



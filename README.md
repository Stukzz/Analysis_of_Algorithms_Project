# Analysis_of_Algorithms_Project

Motivation: To familiarize ourselves with big data problems one such case is processing and making sense of astronomical data

Description of the algorithm: Minimum Spanning Tree (MST) & Voronoi Tessellation: Construct Euclidean geometric graphs on celestial coordinates (ra, dec) to identify cosmic filamentary networks, galaxy groupings, and local environment densities within the Perseus cluster.

description of data: The euclid_ero.perseus_nisp.csv file is a sample source extraction catalog from the Euclid Space Telescope’s Early Release Observations (ERO-10) targeting the Perseus Galaxy Cluster (Abell 426) with the Near-Infrared Spectrometer and Photometer (NISP)
references

Euclid_ero.perseus_nisp + euclid_ero.perseus_vis (Best for Spatial &amp; Graph Algorithms)

1. From Astro Data Lab (https://datalab.noirlab.edu/data-explorer?showTable=euclid_ero.perseus_vis)
   query using the following.
   
Query 1 — NISP Catalog (Filter to J band or pull all 3 bands):

`SELECT number, band, alpha_j2000, delta_j2000, nest4096, htm9,
       mag_auto, magerr_auto, flux_auto, class_star, spread_model, fwhm_image
FROM euclid_ero.perseus_nisp
WHERE mag_auto < 90`

Query 2 — VIS Catalog:

`SELECT number, alpha_j2000, delta_j2000, nest4096, htm9,
       mag_auto, magerr_auto, flux_auto, class_star, spread_model, fwhm_image
FROM euclid_ero.perseus_vis
WHERE mag_auto < 90`

     1.euclid_ero.perseus_nisp: 1,006,020 rows (Near-Infrared Y, J, H bands)
     2.euclid_ero.perseus_vis: 546,563 rows (Visible Ie band, higher spatial resolution)




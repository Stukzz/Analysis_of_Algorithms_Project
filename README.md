# Analysis_of_Algorithms_Project
Euclid_ero.perseus_nisp + euclid_ero.perseus_vis (Best for Spatial &amp; Graph Algorithms)

1. From Astro Data Lab (https://datalab.noirlab.edu/data-explorer?showTable=euclid_ero.perseus_vis)
   query using the following.
__________________________________________________________________________________
Query 1 — NISP Catalog (Filter to J band or pull all 3 bands):

SELECT number, band, alpha_j2000, delta_j2000, nest4096, htm9,
       mag_auto, magerr_auto, flux_auto, class_star, spread_model, fwhm_image
FROM euclid_ero.perseus_nisp
WHERE mag_auto < 90

Query 2 — VIS Catalog:
SELECT number, alpha_j2000, delta_j2000, nest4096, htm9,
       mag_auto, magerr_auto, flux_auto, class_star, spread_model, fwhm_image
FROM euclid_ero.perseus_vis
WHERE mag_auto < 90
__________________________________________________________________________________
     -euclid_ero.perseus_nisp: 1,006,020 rows (Near-Infrared Y, J, H bands)
   
     -euclid_ero.perseus_vis: 546,563 rows (Visible Ie band, higher spatial resolution)


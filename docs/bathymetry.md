# Bathymetry notes

Research notes on real-world bed data. Nothing here is wired into the solver; the solver currently runs on the synthetic archetypes in `crates/wavecore/src/bathymetry.rs`.

## Candidate source: NOAA CUDEM, Hawaii

[NOAA NCEI Continuously Updated DEM (CUDEM)](https://chs.coast.noaa.gov/htdata/raster2/elevation/NCEI_ninth_Topobathy_Hawaii_9428/), 1/9 arc-second (about 3 m) topo-bathymetric tiles. Public, no authentication.

- **Size.** About 8.6 GB for all of Hawaii, in 0.25° tiles named by upper-left corner. The tile covering Pipeline (21.665°N, 158.051°W) is `ncei19_n21x75_w158x25_2021v1.tif`.
- **Partial reads work.** The GeoTIFFs are internally tiled (256 × 256 blocks, deflate), so GDAL's `/vsicurl/` reads only the blocks a window needs. A 2.4 km window around Pipeline came to about 2 MB (713 × 713 cells, roughly 3.2 × 3.4 m each).
- **Datum.** Local mean sea level, float32, nodata `-9999`.

## Provenance

The CUDEM spatial-metadata layer for that tile (`ninth_spatial_meta_hi.zip`, `ncei19_n21x75_w158x25_2021v1_sm.shp`) lists nine source polygons. Overlap with the Pipeline window:

| Source | Overlap | Type |
|---|---|---|
| 2013 USACE NCMP topobathy lidar (Oahu) | 70% | measured |
| SHOALS topobathy lidar (2015, state of Hawaii) | 49% | measured |
| 2013 NOAA lidar (Oahu) | 52% | topographic |
| 2003 NOAA lidar (Oahu coastline) | 36% | topographic |
| 2007 USACE NCMP topobathy lidar | 35% | measured |
| NOAA NOS bathymetric soundings | 3% | point soundings |
| Multibeam, digitised charts, USGS NED | 0% | not used here |

Underwater, the surf zone is covered by three topo-bathymetric lidar campaigns. Multibeam contributes nothing in this window.

## Caveats

- The polygons are coverage extents and overlap, so they do not say which survey decided any single cell.
- Airborne lidar loses returns in breaking whitewater and in murky or deep water. NOAA's lineage text states that unconstrained cells are filled by interpolation.
- 3 m cells cannot resolve features smaller than a few metres.
- Vertical accuracy of the lidar was not checked. A metre of error matters on a shelf a few metres deep.
- Three surveys from different years are blended. A lava reef is stable, so this matters less than it would over sand.

## A sharper source: the 1 m lidar DEM

The 2013 USACE NCMP topobathy lidar of Oahu (CZMIL, collected September to November 2013) is published as a 1 m DEM in the same bucket: `https://noaa-nos-coastal-lidar-pds.s3.amazonaws.com/dem/USACE_Oahu_HI_LMSL_DEM_2013_9365/`, as GeoTIFF tiles with a VRT mosaic in lon/lat (`USACE_Oahu_HI_LMSL_DEM_2013_m9365_EPSG-6322.vrt`) and a tile index. Its metadata says areas without data are masked, not interpolated, so a gap means nothing was measured there.

Along three cross-shore lines of the Pipeline bed (y = 361.5, 511.5 and 601.5 m), sampled every metre from 600 m offshore to the shore:

- it has data on every point in 1 to 10 m of water: the survey was flown in calm conditions, and the surf zone is measured, not filled in;
- it agrees with the 3 m CUDEM bed within 0.22 to 0.25 m on average there;
- its steepest rise over 10 m is 1:3.5, 1:5.2 and 1:5.7, where the 3 m bed reads 1:8.7, 1:6.2 and 1:6.0: short steep steps, a few metres long, that 3 m cells soften.

So at the scale of a wavelength the reef really is as gentle as the 3 m bed says (1:20 to 1:40 where the waves break), and what the 1 m data adds are those steps. They matter for the shape of a break more than for where it is, and they are what a cross-section on 1 m cells should use.

## Other regions

For Portugal, the merged [EMODnet DTM 2024](https://emodnet.ec.europa.eu/en/bathymetry) is about 115 m and too coarse for reef or sandbar geometry. Sentinel-2 satellite-derived bathymetry has errors of roughly 1–5 m to about 10 m depth, which is as large as the signal in a surf zone. Instituto Hidrográfico is mapping Portuguese waters at high resolution under SEAMAP 2030, but no open data for the Lisbon–Peniche coast was found. This search was not exhaustive.

## Attribution

Please credit NOAA National Centers for Environmental Information when publishing results derived from CUDEM data.

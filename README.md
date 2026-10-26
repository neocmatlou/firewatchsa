# FIREWATCH SA — Where South Africa Burns

**Group-selected dashboard: Bright theme with interactive KDE basemap | v5.0 | 23 October 2026**

**Live version previously shared:** https://brudance.shinyapps.io/firewatch-sa/  
**Status:** This ZIP contains the updated KDE-basemap source code. It has **not** been deployed or run against the team's RDS datasets in this environment. The live URL may still show the previously deployed version until you deploy this `app.R`.

## What the project does

Interactive R Shiny spatial statistics dashboard for STAT7006A using the team's **NASA FIRMS VIIRS NOAA-21 2024–2025** South Africa active-fire detections. The group selected the user's dashboard design with a **basemap under the KDE surface**. It retains seven tabs: Start Here; Overview; Seasonal Fire Intelligence; Kernel Density (KDE); Nearest Neighbour; Ripley's K; Methods & Limitations.

### New in this package

- KDE now renders on an **interactive OpenStreetMap basemap**, with a georeferenced raster overlay, legend, pan/zoom and provincial outlines.
- Existing static KDE PNG download is preserved, but that PNG does **not** include the online basemap.
- Existing bright theme, automated interpretations, national/provincial selection and original NOAA-21 datasets are preserved.
- National/provincial analysis samples are limited for responsiveness: KDE 6,000; NN 2,000; Ripley's K 800 detections.

## Folder layout

```
FIREWATCH_SA_FINAL_BASEMAP_SUBMISSION/
├── app.R
├── README.md
├── MANIFEST.txt
├── data/
│   └── ADD_REQUIRED_RDS_FILES.txt
├── docs/
│   ├── 01_Project_and_Script_Guide.docx / .pdf
│   ├── 02_Datasets_and_Sources.docx / .pdf
│   ├── 03_Meeting_and_Presentation_Plan.docx / .pdf
│   └── 04_Final_KDE_Basemap_Handover.docx / .pdf
└── source_archives/
    ├── DL_FIRE_J2V-C2_828479.zip
    └── FIREWATCH_SA_NOAA21_2024_2025_vegetation_detections.zip
```

**Important:** `data/fires_final.rds` and `data/provinces.rds` are NOT in this archive because they were not provided here. Copy the two files from your working `FIREWATCH_SA/data` folder. The source archives are not substitutes for the prepared RDS files.

## Requirements and startup

From the RStudio **Console**, install once:

```r
install.packages(c(
  "shiny", "sf", "dplyr", "ggplot2", "leaflet", "scales",
  "spatstat.geom", "spatstat.explore", "spatstat.random", "terra"
))
```

Ensure the prepared RDS files are in `data/`. Open `app.R` in RStudio and click **Run App**, or run `shiny::runApp(".")` from the project folder.

Expected inputs: `fires_final.rds` is an `sf` point layer with `province`, `year`, `season` and `month` or `acq_date`; `provinces.rds` is an `sf` polygon layer with `province`. The basemap tiles need an internet connection.

## Test the KDE basemap before deployment

1. Run the app locally; verify all seven tabs.
2. Choose a province, year and season with detections; go to **Kernel Density (KDE)** and click **Run KDE Analysis**.
3. Confirm the map shows **OpenStreetMap underneath a colored density surface**, with a legend, provincial outlines, and working pan/zoom.
4. Try **All Provinces** for national KDE. Run NN and Ripley's K on a province, and confirm the interpretation panels update.
5. Confirm the **Download static KDE image** button works. This PNG does not include the basemap.
6. If the KDE overlay errors, send the exact RStudio console error before deploying; raster rendering has not been tested against your local RDS inputs in this environment.

## Publish — Console only

```r
library(rsconnect)
rsconnect::deployApp(
  appDir = "C:/Users/bruda/Documents/FIREWATCH_SA",
  appName = "firewatch-sa",
  appFiles = c("app.R", "data/fires_final.rds", "data/provinces.rds"),
  forceUpdate = TRUE
)
```

The Windows path above is the user's existing project folder. Change it if another group member deploys. Never paste deployment code into `app.R`. Never share or commit deployment tokens/secrets.

## Source links

- NASA FIRMS: https://firms.modaps.eosdis.nasa.gov/
- Download/archive: https://firms.modaps.eosdis.nasa.gov/download/
- FIRMS data guide: https://firms.modaps.eosdis.nasa.gov/download/Readme.txt
- OpenStreetMap copyright and attribution: https://www.openstreetmap.org/copyright
- **Provincial boundaries:** exact original provider/version not confirmed. Confirm before citation. Do not assume GADM or geoBoundaries was used.

## Methodological limitations

FIRMS detections are not distinct fire events or burned area. KDE is descriptive, not a significance test; sampled densities are scaled estimates. NN uses a uniform CSR benchmark and a one-sided Monte Carlo clustering p-value. Ripley's K uses pointwise CSR envelopes, not a global test. Spatial clustering may reflect satellite coverage, land use, repeated observations and the uniform null model. Calendar-year summer groups December with January–February in the same year.

## Presentation note

Start the app a few minutes early and run KDE once. Hosted sessions can disconnect after idle time; refresh if needed. Have the static KDE PNG and screenshots as a backup. The group should validate the basemap locally and online before presenting.

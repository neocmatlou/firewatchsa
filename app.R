# ================================================================
# FIREWATCH SA | WHERE SOUTH AFRICA BURNS
# Interactive spatial statistics + seasonal intelligence dashboard
# NASA FIRMS NOAA-21
# Data: 2024-2025, South Africa | v5.0: bright theme + interactive KDE basemap
#
# Required project files:
#   data/fires_final.rds
#   data/provinces.rds
#
# Install packages ONCE in the RStudio Console (not in this file):
# install.packages(c("shiny","sf","dplyr","ggplot2","leaflet",
#   "scales","spatstat.geom","spatstat.explore","spatstat.random",
#   "terra"))
# ================================================================

library(shiny)
library(sf)
library(dplyr)
library(ggplot2)
library(leaflet)
library(scales)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(terra)

# ------------------------- DATA ---------------------------------
stopifnot(file.exists("data/fires_final.rds"),
          file.exists("data/provinces.rds"))

fires <- st_transform(readRDS("data/fires_final.rds"), 4326)
provinces <- st_transform(readRDS("data/provinces.rds"), 4326)

required_fire_cols <- c("province", "year", "season")
stopifnot(all(required_fire_cols %in% names(fires)))
stopifnot("province" %in% names(provinces))

fires$province <- as.character(fires$province)
fires$season <- as.character(fires$season)
fires$year <- as.integer(as.character(fires$year))
provinces$province <- as.character(provinces$province)

# Coordinate extraction is done only once.
coords <- st_coordinates(fires)
fires$map_lng <- coords[, 1]
fires$map_lat <- coords[, 2]
rm(coords)

province_names <- sort(unique(provinces$province))

# ---------------- SEASONAL INTELLIGENCE DATA ----------------------
# Aggregate ONCE at startup. This avoids repeatedly processing sf geometry
# when the user explores monthly and seasonal comparisons.
# Prefer the archive's numeric month column; derive it from acquisition date
# if the prepared RDS does not contain a month column.
fire_attributes <- sf::st_drop_geometry(fires)
if ("month" %in% names(fire_attributes)) {
  fire_attributes$month_num <- suppressWarnings(
    as.integer(as.character(fire_attributes$month))
  )
} else if ("acq_date" %in% names(fire_attributes)) {
  fire_attributes$month_num <- as.integer(format(
    as.Date(as.character(fire_attributes$acq_date)), "%m"
  ))
} else {
  stop("Seasonal Intelligence needs 'month' or 'acq_date' in fires_final.rds.")
}
stopifnot(all(!is.na(fire_attributes$month_num)),
          all(fire_attributes$month_num %in% 1:12))

season_labels <- c("Summer", "Autumn", "Winter", "Spring")
month_season <- function(month_num) {
  season_labels[c(1L, 1L, 2L, 2L, 2L, 3L,
                  3L, 3L, 4L, 4L, 4L, 1L)[month_num]]
}

monthly_observed <- fire_attributes |>
  dplyr::filter(year %in% c(2024L, 2025L)) |>
  dplyr::count(province, year, month_num, name = "detections")
monthly_grid <- expand.grid(
  province = province_names,
  year = c(2024L, 2025L),
  month_num = 1:12,
  stringsAsFactors = FALSE
)
monthly_counts <- monthly_grid |>
  dplyr::left_join(monthly_observed,
                   by = c("province", "year", "month_num")) |>
  dplyr::mutate(
    detections = dplyr::coalesce(detections, 0L),
    season = month_season(month_num)
  )
rm(fire_attributes, monthly_observed, monthly_grid)

month_abbr <- month.abb
format_change <- function(old, new) {
  if (old == 0) return("N/A (2024 baseline is zero)")
  sprintf("%+.1f%%", 100 * (new - old) / old)
}

# ---------------------- HELPERS ---------------------------------
fmt <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", trim = TRUE,
         scientific = FALSE, nsmall = digits)
}

dark_plot <- function() {
  # Bright, presentation-friendly ggplot styling: warm white panels,
  # deep charcoal text, warm gridlines and vivid fire-colour accents.
  theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "#FFFFFF", colour = NA),
      panel.background = element_rect(fill = "#FFFFFF", colour = NA),
      panel.grid.major = element_line(colour = "#E7EAF0", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      text = element_text(colour = "#1F2937"),
      axis.text = element_text(colour = "#334155", size = 12),
      axis.title = element_text(colour = "#1F2937", size = 13),
      axis.ticks = element_line(colour = "#64748B"),
      plot.title = element_text(face = "bold", size = 16, colour = "#182235"),
      plot.subtitle = element_text(colour = "#475569", size = 11),
      legend.background = element_rect(fill = "#FFFFFF", colour = NA),
      legend.key = element_rect(fill = "#FFFFFF", colour = NA),
      legend.text = element_text(colour = "#334155", size = 11),
      legend.title = element_text(colour = "#1F2937", size = 11),
      strip.text = element_text(colour = "#1F2937"),
      plot.margin = margin(12, 18, 14, 14)
    )
}

# Distances are calculated in projected metres, never longitude/latitude.
# Provincial selections use a local southern-hemisphere UTM zone.
# National selections use a South Africa-centred Lambert azimuthal equal-area
# projection (LAEA). Its distance distortion is modest but non-zero.
# The same source datasets are used for all scopes.
make_point_pattern <- function(fire_sf, province_sf, max_points, seed = 123) {
  if (nrow(fire_sf) < 10) {
    stop("At least 10 fire detections are needed for this analysis.")
  }

  total <- nrow(fire_sf)
  set.seed(seed)
  if (total > max_points) {
    fire_sf <- fire_sf[sample.int(total, max_points), ]
  }

  national <- nrow(province_sf) > 1L
  if (national) {
    # An equal-area national projection is preferable to using one UTM
    # zone for a country spanning several zones. Distances are approximate.
    analysis_crs <- paste(
      "+proj=laea +lat_0=-29 +lon_0=25 +datum=WGS84",
      "+units=m +no_defs"
    )
    projection_label <- "South Africa-centred LAEA (metres)"
  } else {
    bb <- st_bbox(province_sf)
    lon_mid <- (unname(bb["xmin"]) + unname(bb["xmax"])) / 2
    zone <- max(1, min(60, floor((lon_mid + 180) / 6) + 1))
    epsg <- 32700 + zone
    analysis_crs <- epsg
    projection_label <- paste0("EPSG:", epsg, " (UTM, metres)")
  }

  province_proj <- st_transform(province_sf, analysis_crs)
  points_proj <- st_transform(fire_sf, analysis_crs)
  # Dissolve provincial borders so they do not become artificial national
  # edges. st_union preserves coastline and the national exterior.
  boundary_union <- sf::st_union(sf::st_geometry(province_proj))
  win <- spatstat.geom::as.owin(boundary_union)
  xy <- st_coordinates(points_proj)

  # Protect against tiny projection/boundary precision differences.
  keep <- spatstat.geom::inside.owin(xy[, 1], xy[, 2], win)
  xy <- xy[keep, , drop = FALSE]
  if (nrow(xy) < 10) {
    stop("Too few points fall inside the projected study area.")
  }

  duplicate_locations <- sum(duplicated(
    data.frame(easting = xy[, 1], northing = xy[, 2])
  ))

  pp <- spatstat.geom::ppp(
    x = xy[, 1], y = xy[, 2], window = win, checkdup = FALSE
  )

  list(
    pp = pp,
    duplicate_locations = duplicate_locations,
    total = total,
    used = spatstat.geom::npoints(pp),
    projection_label = projection_label,
    analysis_crs = analysis_crs,
    national = national,
    sampled = total > max_points
  )
}

# Convert spatstat's pixel KDE (bottom-to-top rows) to a georeferenced
# terra raster (top-to-bottom rows), then transform to WGS84 for Leaflet.
# Keep the raster resolution modest: only 180 x 180 pixels are calculated.
kde_to_wgs84_raster <- function(result) {
  im <- result$surface
  stopifnot(length(im$xcol) >= 2, length(im$yrow) >= 2)
  dx <- abs(im$xcol[2] - im$xcol[1])
  dy <- abs(im$yrow[2] - im$yrow[1])
  rr <- terra::rast(
    nrows = length(im$yrow), ncols = length(im$xcol),
    xmin = min(im$xcol) - dx / 2,
    xmax = max(im$xcol) + dx / 2,
    ymin = min(im$yrow) - dy / 2,
    ymax = max(im$yrow) + dy / 2,
    crs = sf::st_crs(result$analysis_crs)$wkt
  )
  # spatstat im rows increase in y; terra cells start at the north edge.
  north_first <- im$v[nrow(im$v):1, , drop = FALSE]
  terra::values(rr) <- as.vector(t(north_first))
  terra::project(rr, 'EPSG:4326', method = 'bilinear')
}

draw_kde <- function(x) {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(3, 3, 4, 5), bg = "#FFFFFF",
      fg = "#182235", col.main = "#182235",
      col.axis = "#334155", col.lab = "#334155")
  plot(x$surface, main = paste0(x$province, " | KDE | ",
                               x$bandwidth, " km bandwidth"),
       col = hcl.colors(90, "Inferno"))
  plot(Window(x$pp), add = TRUE, border = "#182235", lwd = 1.2)
}

stat_note <- function(result) {
  paste0(
    "Based on ", fmt(result$used), " ",
    if (result$sampled) "randomly sampled" else "selected",
    " satellite detections from ", fmt(result$total),
    " filtered records. Projection: ", result$projection_label,
    ". Exact-coordinate duplicates in analytical sample: ",
    fmt(result$duplicate_locations),
    ". Repeated satellite observations are not necessarily separate fires."
  )
}

# ------------------------ INTERFACE ------------------------------
ui <- fluidPage(
  tags$head(
    tags$meta(name = "viewport",
              content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      /* FIREWATCH SA: brighter, warm, fire-inspired visual system.
         CSS-only changes; data, filters, and statistical calculations stay intact. */
      :root {
        --fw-ink:#182235; --fw-text:#334155; --fw-muted:#526174;
        --fw-orange:#E85D04; --fw-coral:#E63946; --fw-gold:#FFB703;
        --fw-teal:#087E8B; --fw-canvas:#FFF8F0; --fw-white:#FFFFFF;
        --fw-border:#E9E2DC;
      }
      body {background:var(--fw-canvas);color:var(--fw-ink);
            font-family:'Segoe UI',Arial,sans-serif;}
      .container-fluid {padding:20px;max-width:1800px;}
      .fw-header {
        background:linear-gradient(112deg,#FFF1DA 0%,#FFE3CB 53%,#FFD3BF 100%);
        border-left:7px solid var(--fw-orange);
        border-top:1px solid #F8C99D;
        padding:25px 27px;border-radius:15px;margin-bottom:20px;
        box-shadow:0 7px 24px rgba(116,54,12,.08);
      }
      .fw-header h1 {color:#B93609;font-size:35px;font-weight:850;margin:0;}
      .fw-header p {color:#5F3B2E;margin:8px 0 0;font-weight:500;}
      .fw-sidebar,.fw-panel,.fw-card {
        background:var(--fw-white);border:1px solid var(--fw-border);
        border-radius:14px;padding:19px;margin-bottom:18px;
        box-shadow:0 5px 19px rgba(37,45,62,.055);
      }
      .fw-sidebar h3 {margin-top:3px;font-weight:750;color:var(--fw-ink);}
      .fw-sidebar label {color:var(--fw-ink);font-weight:650;}
      .fw-sidebar p,.fw-muted {color:var(--fw-muted);font-size:13px;}
      .fw-card {border-top:4px solid var(--fw-orange);min-height:106px;}
      .fw-card:nth-child(2n) {border-top-color:var(--fw-coral);}
      .fw-card .fw-label {font-size:12px;color:#526174;
                          text-transform:uppercase;font-weight:700;letter-spacing:.4px;}
      .fw-card .fw-value {color:#C4470B;font-size:29px;
                          font-weight:800;margin-top:8px;}
      .fw-panel h3 {font-size:21px;font-weight:750;margin:3px 0 16px;color:var(--fw-ink);}
      .fw-panel h4 {color:var(--fw-ink);}
      .fw-panel .shiny-plot-output {background:#FFFFFF;
          border:1px solid #F1E8E0;border-radius:9px;overflow:hidden;}
      .fw-panel .shiny-plot-output img {max-width:100%;}
      .form-control,.selectize-input {background:#FFFFFF;color:#182235;
                                     border:1px solid #D4DCE4;border-radius:7px;}
      .selectize-dropdown {color:#182235;background:#FFFFFF;}
      .form-control:focus,.selectize-input.focus {
        border-color:#E85D04;box-shadow:0 0 0 2px rgba(232,93,4,.12);}
      .btn-primary,.btn-default.fw-run {
        background:var(--fw-orange);border-color:var(--fw-orange);
        color:#FFFFFF;font-weight:750;border-radius:8px;
      }
      .btn-primary:hover,.btn-default.fw-run:hover,
      .btn-primary:focus,.btn-default.fw-run:focus {
        background:#C4470B;border-color:#C4470B;color:#FFFFFF;
      }
      .nav-tabs {border-bottom:2px solid #F4D4BB;margin-bottom:18px;}
      .nav-tabs>li>a {color:#3E4B5C;background:#FFF0E2;
                       border:1px solid #F2D9C6;border-radius:9px 9px 0 0;
                       margin-right:5px;font-weight:650;}
      .nav-tabs>li>a:hover {background:#FFD8B7;color:#8F330B;}
      .nav-tabs>li.active>a,
      .nav-tabs>li.active>a:focus,
      .nav-tabs>li.active>a:hover {
        background:var(--fw-orange);color:#FFFFFF;
        border-color:var(--fw-orange);font-weight:750;
      }
      .well {background:#FFF6EC;border:1px solid #F0DFD0;color:#334155;}
      .leaflet-container {border-radius:9px;}
      hr {border-color:#E9E2DC;}
      .fw-note {background:#FFF3E8;border-left:4px solid var(--fw-orange);
                border-radius:7px;padding:13px;color:#563C30;
                margin-top:12px;line-height:1.55;}
      .shiny-output-error-validation {color:#A33019;}
      .fw-actions {display:flex;flex-wrap:wrap;gap:10px;margin:12px 0;}
      .fw-actions .btn-default {background:var(--fw-orange);
                                border-color:var(--fw-orange);
                                color:#FFFFFF;font-weight:750;border-radius:8px;}
      .fw-actions .btn-default:hover {background:#C4470B;
                                     border-color:#C4470B;color:#FFFFFF;}
      .fw-method {background:#FFF5E9;border-radius:9px;padding:16px;
                  margin:12px 0;border-left:4px solid var(--fw-coral);}
      .fw-method h4 {color:#B93609;margin:0 0 9px;font-weight:750;}
      .fw-method p {line-height:1.6;color:#334155;}
      .fw-method a {color:#B93609;text-decoration:underline;}
      .fw-intro {background:#FFF5EB;border:1px solid #F2DCC9;
                 border-radius:10px;padding:16px 19px;margin-bottom:18px;
                 color:#334155;line-height:1.65;}
      .fw-intro strong {color:#B93609;}
      .fw-legend {font-size:13px;color:#526174;margin-top:8px;}
      /* Bright editorial welcome, with strong contrast and warm fire accents. */
      .fw-welcome {
        background:linear-gradient(115deg,#FFF0D9 0%,#FFD8BB 63%,#FFBCA5 100%);
        border:1px solid #F4C49F;border-radius:16px;
        padding:30px 32px;margin-bottom:18px;
        box-shadow:0 7px 24px rgba(116,54,12,.08);
      }
      .fw-eyebrow {color:#B93609;font-weight:800;letter-spacing:1.5px;
                   font-size:12px;text-transform:uppercase;}
      .fw-welcome h2 {font-size:30px;font-weight:850;color:#842F12;
                      margin:12px 0 14px;line-height:1.25;}
      .fw-welcome p {font-size:16px;color:#503B35;line-height:1.75;
                     max-width:900px;}
      .fw-welcome .btn {margin-top:12px;padding:10px 18px;}
      .fw-context-card {background:#FFFFFF;border:1px solid #F1DDCB;
                        border-top:4px solid var(--fw-orange);border-radius:11px;
                        padding:20px;min-height:205px;margin-bottom:15px;
                        box-shadow:0 5px 16px rgba(37,45,62,.045);}
      .fw-context-card h4 {color:#B93609;font-weight:750;font-size:17px;
                           margin:0 0 11px;}
      .fw-context-card p {color:#334155;font-size:14px;line-height:1.65;}
      .fw-step {background:#EAF7F7;border-left:4px solid var(--fw-teal);
                border-radius:8px;padding:15px 18px;margin-bottom:12px;}
      .fw-step strong {color:#14535B;font-size:15px;}
      .fw-step p {color:#334155;line-height:1.55;margin:5px 0 0;}
      .fw-panel a,.fw-welcome a {color:#B93609;text-decoration:underline;}
      .fw-panel a:hover,.fw-welcome a:hover {color:#842F12;}
      @media(max-width:768px) {
        .fw-header h1 {font-size:27px;}
        .fw-card .fw-value {font-size:23px;}
        .fw-welcome {padding:22px;}
        .fw-welcome h2 {font-size:25px;}
      }
    "))
  ),

  div(class = "fw-header",
      h1("FIREWATCH SA"),
      p("WHERE SOUTH AFRICA BURNS"),
      p("Interactive spatial statistics of satellite-detected active fires",
        "| VIIRS NOAA-21 | 2024–2025")),

  fluidRow(
    column(
      width = 3,
      div(class = "fw-sidebar",
          h3("Explore Fire Patterns"),
          selectInput("province", "Province",
                      choices = c("All Provinces", province_names),
                      selected = "All Provinces"),
          selectInput("year", "Year",
                      choices = c("All Years", "2024", "2025"),
                      selected = "All Years"),
          selectInput("season", "Season",
                      choices = c("All Seasons", "Summer", "Autumn",
                                  "Winter", "Spring"),
                      selected = "All Seasons"),
          hr(),
          h4("Dataset"),
          p("Source: NASA FIRMS"),
          p("Satellite: VIIRS NOAA-21"),
          p("Study area: South Africa"),
          p("Period: 2024–2025"),
          p("These are satellite active-fire detections, not unique",
            "wildfire incidents or measured burned area."),
          hr(),
          p("For KDE, Nearest Neighbour and Ripley's K, select",
            "a province or All Provinces (national), then click Run Analysis."),
          p("National analyses use reproducible point samples for speed.",
            "Results are exploratory, not estimates of fire causes.")
      )
    ),

    column(
      width = 9,
      fluidRow(
        column(4, div(class = "fw-card",
                      div(class = "fw-label", "Fire detections"),
                      div(class = "fw-value",
                          textOutput("total_fires", inline = TRUE)))),
        column(4, div(class = "fw-card",
                      div(class = "fw-label", "Provinces with detections"),
                      div(class = "fw-value",
                          textOutput("total_provinces", inline = TRUE)))),
        column(4, div(class = "fw-card",
                      div(class = "fw-label", "Selected period"),
                      div(class = "fw-value",
                          textOutput("selected_period", inline = TRUE))))
      ),

      tabsetPanel(
        id = "main_tab",

        tabPanel(
          title = "Start Here",
          div(class = "fw-welcome",
              div(class = "fw-eyebrow", "South Africa | Fire patterns | 2024–2025"),
              h2("Why study fire patterns in South Africa?"),
              p("Fire is a natural part of many South African grassland and",
                "savanna landscapes, but uncontrolled veld fires can threaten",
                "people, farms, infrastructure and ecosystems. Their timing",
                "and location vary across the country. Understanding where",
                "satellites detect fire activity, and how those detections",
                "are distributed, provides a useful starting point for",
                "spatial investigation."),
              actionButton("start_exploring", "Explore the dashboard →",
                           class = "btn-primary")),
          fluidRow(
            column(4, div(class = "fw-context-card",
                          h4("The South African context"),
                          p("Different vegetation types, seasonal conditions",
                            "and land-use practices shape fire activity.",
                            "Fires can be ecologically important, yet",
                            "uncontrolled fires can have serious social,",
                            "economic and environmental consequences."))),
            column(4, div(class = "fw-context-card",
                          h4("The challenge"),
                          p("Thousands of satellite detections are difficult",
                            "to interpret as a table of coordinates. Raw",
                            "counts do not show whether observations are",
                            "spatially concentrated, how patterns vary",
                            "between provinces, or when activity peaks."))),
            column(4, div(class = "fw-context-card",
                          h4("Our response"),
                          p("FIREWATCH SA turns satellite detections into",
                            "interactive maps, seasonal comparisons and",
                            "three spatial statistical analyses, helping",
                            "users move from seeing points to asking",
                            "testable questions about their distribution.")))
          ),
          div(class = "fw-panel",
              h3("What is FIREWATCH SA trying to answer?"),
              tags$p("Our aim is to explore the spatial and seasonal distribution",
                     "of satellite-detected active-fire observations across",
                     "South Africa during 2024 and 2025."),
              tags$ul(
                tags$li("Where are fire detections concentrated?"),
                tags$li("How do detection counts change across months, seasons and years?"),
                tags$li("Are observed point patterns more clustered than a",
                        "uniform random spatial pattern?"),
                tags$li("At what distances do departures from spatial randomness appear?")
              )),
          div(class = "fw-panel",
              h3("From satellite data to spatial insight"),
              fluidRow(
                column(6, div(class = "fw-step",
                              strong("01 | Collect and prepare"),
                              p("Use NASA FIRMS VIIRS NOAA-21 archive detections",
                                "for 2024–2025; clean the observations and",
                                "assign detections to South African provinces."))),
                column(6, div(class = "fw-step",
                              strong("02 | Explore"),
                              p("Filter the interactive map by province,",
                                "year and season, and inspect detection counts."))),
                column(6, div(class = "fw-step",
                              strong("03 | Compare"),
                              p("Use Seasonal Fire Intelligence to compare",
                                "monthly and seasonal patterns between",
                                "2024 and 2025."))),
                column(6, div(class = "fw-step",
                              strong("04 | Analyze"),
                              p("Apply KDE to visualize intensity, Nearest",
                                "Neighbour to assess local spacing, and",
                                "Ripley's K to examine multiple distances.")))
              )),
          div(class = "fw-panel",
              h3("What exactly are we mapping?"),
              p("The source is NASA's Fire Information for Resource",
                "Management System (FIRMS). The VIIRS instrument on",
                "NOAA-21 detects thermal anomalies at approximately",
                "375 m nominal resolution. The dashboard uses a",
                "historical 2024–2025 dataset, not a live fire feed."),
              div(class = "fw-note",
                  strong("Important distinction: "),
                  "A satellite detection is not necessarily one unique",
                  "wildfire. A fire can be detected repeatedly, some fires",
                  "may be missed, and detections do not measure burned area.",
                  "The dashboard explores detected activity, not confirmed",
                  "fire incidents, causes or forecasts."),
              tags$p(style = "margin-top:14px;",
                     strong("Sources and further reading: "),
                     tags$a(href = "https://firms.modaps.eosdis.nasa.gov/",
                            target = "_blank", rel = "noopener noreferrer",
                            "NASA FIRMS"),
                     " | ",
                     tags$a(href = "https://www.gov.za/news/media-statements/prevent-veld-fires-during-winter-02-jul-2014",
                            target = "_blank", rel = "noopener noreferrer",
                            "South African Government: Veld fire prevention"),
                     " | ",
                     tags$a(href = "https://www.dffe.gov.za/working-fire-0",
                            target = "_blank", rel = "noopener noreferrer",
                            "DFFE: Working on Fire"))),
          div(class = "fw-note",
              "Tip: The sidebar filters and the three summary cards at",
              "the top of the dashboard apply to the selected subset.",
              "This introduction provides national context and does not",
              "change when filters are selected.")
        ),

        tabPanel(
          title = "Overview",
          div(class = "fw-panel",
              h3("Satellite Fire Detection Map"),
              leafletOutput("fire_map", height = 490),
              div(class = "fw-muted",
                  "For map performance, up to 3,000 detections are displayed.",
                  "All filtered detections are used for summary counts."),
              div(class = "fw-actions",
                  downloadButton("download_fires", "Download Filtered Data (CSV)"))),
          div(class = "fw-panel",
              h3("Fire Detections by Province"),
              plotOutput("province_chart", height = 370)),
          div(class = "fw-panel",
              h3("Monthly Fire Activity"),
              p("Compare the monthly distribution of detections across",
                "2024 and 2025. This chart follows the Overview filters."),
              plotOutput("overview_monthly_chart", height = 340)),
          div(class = "fw-panel",
              h3("Understanding the Fire Patterns"),
              textOutput("pattern_summary"),
              div(class = "fw-note",
                  "Counts are not adjusted for provincial area,",
                  "vegetation availability or satellite observation coverage."))
        ),

        tabPanel(
          title = "Seasonal Fire Intelligence",
          div(class = "fw-intro",
              strong("Compare, don't just map."),
              " This page uses the full cleaned dataset to compare",
              " monthly and seasonal detection patterns in 2024 and 2025.",
              " It follows the Province selection in the sidebar,",
              " but intentionally ignores the sidebar's Year and Season",
              " filters so that both years can always be compared."),
          fluidRow(
            column(3, div(class = "fw-card",
                          div(class = "fw-label", "2024 detections"),
                          div(class = "fw-value",
                              textOutput("intel_2024", inline = TRUE)))),
            column(3, div(class = "fw-card",
                          div(class = "fw-label", "2025 detections"),
                          div(class = "fw-value",
                              textOutput("intel_2025", inline = TRUE)))),
            column(3, div(class = "fw-card",
                          div(class = "fw-label", "2025 vs 2024"),
                          div(class = "fw-value",
                              textOutput("intel_change", inline = TRUE)))),
            column(3, div(class = "fw-card",
                          div(class = "fw-label", "Peak month (both years)"),
                          div(class = "fw-value",
                              textOutput("intel_peak", inline = TRUE))))
          ),
          div(class = "fw-panel",
              h3("Monthly Fire Detection Timeline"),
              p("Compare the same calendar months between the two years.",
                " Each line is a different year; zero-detection months",
                " are included."),
              plotOutput("intel_monthly", height = 390),
              div(class = "fw-actions",
                  downloadButton("download_intel_monthly",
                                 "Download Monthly Comparison (CSV)"))),
          div(class = "fw-panel",
              h3("Which Seasons Changed?"),
              p("Compare detections for summer (Dec–Feb), autumn (Mar–May),",
                " winter (Jun–Aug), and spring (Sep–Nov).",
                " Each season is grouped within a calendar year."),
              plotOutput("intel_seasonal", height = 380),
              div(class = "fw-actions",
                  downloadButton("download_intel_seasonal",
                                 "Download Seasonal Comparison (CSV)")),
              uiOutput("intel_story")),
          div(class = "fw-panel",
              h3("South Africa's Seasonal Signatures"),
              p("For each province, this heatmap shows the percentage of",
                " its 2024–2025 detections occurring in each season.",
                " This reveals differences in seasonal composition",
                " rather than simply rewarding provinces with more points."),
              plotOutput("intel_heatmap", height = 410),
              div(class = "fw-legend",
                  "National view: the heatmap includes all nine provinces,",
                  " regardless of the Province filter.")),
          div(class = "fw-note",
              strong("Interpretation limits: "),
              "These are counts of satellite active-fire detections,",
              " not numbers of unique fires or burned area.",
              " Differences may reflect observation coverage, repeated",
              " detections, vegetation and other factors.",
              " Do not interpret a year-to-year change as a measured",
              " change in wildfire incidence or risk.")
        ),

        tabPanel(
          title = "Kernel Density (KDE)",
          div(class = "fw-panel",
              h3("Kernel Density Estimation"),
              p("Explore where detections are most concentrated.",
                "Adjust the smoothing bandwidth and calculate a",
                "density surface for a province or all of South Africa."),
              sliderInput("kde_bandwidth", "Bandwidth (km)",
                          min = 2, max = 30, value = 10, step = 1),
              actionButton("run_kde", "Run KDE Analysis",
                           class = "btn-primary"),
              div(class = "fw-note",
                  "KDE uses up to 6,000 randomly sampled points.",
                  "The density surface is scaled to estimate intensity",
                  "for the full filtered dataset. It is an estimate,",
                  "not a count of independent wildfire events.",
                  "National patterns are smoothed at the selected bandwidth,",
                  "so try a larger bandwidth to see broader patterns."),
              br(),
              p(strong("Interactive KDE basemap: "),
                "Pan and zoom to inspect hotspots over OpenStreetMap. ",
                "The colors show estimated detection density, not burned area."),
              leafletOutput("kde_map", height = 540),
              div(class = "fw-actions",
                  downloadButton("download_kde_png", "Download static KDE image (PNG)")),
              uiOutput("kde_explanation"),
              div(class = "fw-muted",
                  "The density scale reflects detections per square kilometre,",
                  "not burned area, fire intensity or fire risk."))
        ),

        tabPanel(
          title = "Nearest Neighbour",
          div(class = "fw-panel",
              h3("Nearest Neighbour Analysis"),
              p("Compare the observed mean nearest-neighbour distance",
                "with distances expected under Complete Spatial",
                "Randomness (CSR)."),
              actionButton("run_nn", "Run Nearest Neighbour Analysis",
                           class = "btn-primary"),
              div(class = "fw-note",
                  "Uses up to 2,000 randomly sampled detections and",
                  "99 Monte Carlo simulations. The reported p-value",
                  "tests clustering against uniform CSR."),
              br(),
              uiOutput("nn_summary"),
              plotOutput("nn_plot", height = 390),
              div(class = "fw-actions",
                  downloadButton("download_nn", "Download NN Results (CSV)")),
              uiOutput("nn_explanation"),
              div(class = "fw-muted",
                  "The p-value compares the sample with uniform random points",
                  "inside the selected boundary. It does not identify causes."))
        ),

        tabPanel(
          title = "Ripley's K",
          div(class = "fw-panel",
              h3("Ripley's K and L(r) - r"),
              p("Investigate departures from CSR across distances",
                "from 0 to 10 km, using pointwise simulation envelopes."),
              actionButton("run_k", "Run Ripley's K Analysis",
                           class = "btn-primary"),
              div(class = "fw-note",
                  "Uses up to 800 randomly sampled detections and",
                  "19 CSR simulations. Envelopes are pointwise and",
                  "exploratory, not a global significance test."),
              br(),
              plotOutput("k_plot", height = 460),
              div(class = "fw-actions",
                  downloadButton("download_k", "Download K/L Curves (CSV)")),
              uiOutput("k_explanation"),
              div(class = "fw-muted",
                  "Envelopes are pointwise; values outside them are not",
                  "a formal global test of spatial clustering."))
        ),

        tabPanel(
          title = "Methods & Limitations",
          div(class = "fw-panel",
              h3("How FIREWATCH SA Works"),
              p("This dashboard explores the spatial pattern of VIIRS NOAA-21",
                "satellite active-fire detections in South Africa (2024–2025).",
                "All analyses are descriptive or exploratory comparisons against",
                "a defined statistical reference, not proof of fire causes."),
              div(class = "fw-method",
                  h4("Seasonal Fire Intelligence"),
                  p("Monthly and seasonal summaries are computed from all",
                    " cleaned detections rather than the point samples",
                    " used in computationally intensive statistics.",
                    " Summer is December–February, autumn March–May,",
                    " winter June–August, and spring September–November.",
                    " December is grouped with the summer months of its",
                    " calendar year; these are not continuous cross-year",
                    " summer seasons. Percent change uses 2024 as its",
                    " denominator and is undefined if the baseline is zero.",
                    " The seasonal heatmap displays the percentage share",
                    " of each province's total detections across 2024–2025.")),
              div(class = "fw-method",
                  h4("1. Kernel Density Estimation (KDE)"),
                  p("KDE smooths point detections into a continuous estimated",
                    "intensity surface. Brighter colors indicate higher",
                    "estimated detections per square kilometer.",
                    "Bandwidth controls smoothing, so peaks change with",
                    "the selected bandwidth. A random sample of up to 6,000",
                    "detections is used for speed; intensity is scaled to",
                    "the selected full dataset.")),
              div(class = "fw-method",
                  h4("2. Nearest Neighbour Analysis"),
                  p("For each sampled detection, the distance to its closest",
                    "neighbor is calculated. A ratio below 1 means the",
                    "observed mean is smaller than the infinite-plane CSR",
                    "approximation. A one-sided Monte Carlo test compares",
                    "the observed mean with 99 uniformly random patterns",
                    "having the same point count and provincial boundary.",
                    "A small one-sided p-value supports clustering relative",
                    "to this specific uniform-CSR null model.")),
              div(class = "fw-method",
                  h4("3. Ripley's K and L(r) - r"),
                  p("Ripley's K compares point-pair proximity over multiple",
                    "distances. The L(r) - r transformation has a CSR",
                    "reference of zero. The app uses a border correction",
                    "and 19 uniform-CSR simulations for exploratory",
                    "pointwise envelopes. Crossing the envelope at many",
                    "distances is not a formal global significance test.")),
              div(class = "fw-method",
                  h4("Projection and sampling"),
                  p("Distances are calculated in meters after transforming",
                    "each province into a locally chosen southern-hemisphere",
                    "UTM projection. The map uses WGS84. Sampling is",
                    "reproducible (random seed 123) and the app displays",
                    "the number of observations actually analyzed.")),
              div(class = "fw-method",
                  h4("Scientific limitations"),
                  p("A satellite detection is not a unique wildfire event.",
                    "Repeated detections and exact-coordinate duplicates",
                    "may strengthen apparent clustering. The uniform-CSR",
                    "reference does not account for nonuniform vegetation,",
                    "topography, human activity, satellite observation",
                    "coverage or other drivers. Counts are not area-normalized,",
                    "and different years should not be treated as direct",
                    "fire-frequency trends without further checks.",
                    "National analyses combine different environments",
                    "and use a fixed-size random sample; national and",
                    "provincial summary statistics are not directly",
                    "comparable without accounting for these differences.")),
              div(class = "fw-method",
                  h4("Data source"),
                  p("NASA Fire Information for Resource Management System",
                    "(FIRMS), VIIRS NOAA-21 archive. Provincial boundaries",
                    "come from the prepared spatial boundary dataset."),
                  tags$a(href = "https://firms.modaps.eosdis.nasa.gov/",
                         target = "_blank", rel = "noopener noreferrer",
                         "NASA FIRMS website"))
          )
        )
      )
    )
  )
)

# ------------------------- SERVER --------------------------------
server <- function(input, output, session) {

  # Landing-page navigation; no statistical computation on the intro tab.
  observeEvent(input$start_exploring, {
    updateTabsetPanel(session, "main_tab", selected = "Overview")
  })

  filtered_fires <- reactive({
    x <- fires
    if (input$province != "All Provinces") {
      x <- x[x$province == input$province, ]
    }
    if (input$year != "All Years") {
      x <- x[x$year == as.integer(input$year), ]
    }
    if (input$season != "All Seasons") {
      x <- x[x$season == input$season, ]
    }
    x
  })

  selected_boundary <- reactive({
    if (input$province == "All Provinces") {
      provinces
    } else {
      provinces[provinces$province == input$province, ]
    }
  })

  # An existing result must never silently appear as if it reflects
  # newly changed filters. The user must rerun the method.
  filters_match <- function(result, check_bandwidth = FALSE) {
    identical(result$province, input$province) &&
      identical(result$year, input$year) &&
      identical(result$season, input$season) &&
      (!check_bandwidth || identical(result$bandwidth, input$kde_bandwidth))
  }

  # ---- Overview ----
  output$total_fires <- renderText(fmt(nrow(filtered_fires())))
  output$total_provinces <- renderText({
    dplyr::n_distinct(filtered_fires()$province)
  })
  output$selected_period <- renderText({
    if (input$year == "All Years") "2024–2025" else input$year
  })

  output$fire_map <- renderLeaflet({
    x <- filtered_fires()
    boundary <- selected_boundary()
    validate(need(nrow(x) > 0, "No detections for these filters."))

    set.seed(123)
    if (nrow(x) > 3000) {
      x <- x[sample.int(nrow(x), 3000), ]
    }
    bb <- st_bbox(boundary)

    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addTiles() |>
      addPolygons(
        data = boundary, color = "#475569", weight = 2,
        fillColor = "#E85D04", fillOpacity = 0.05,
        label = ~province, group = "Provincial boundaries"
      ) |>
      addCircleMarkers(
        lng = x$map_lng, lat = x$map_lat,
        radius = 2.5, color = "#E63946", fillColor = "#E63946",
        stroke = FALSE, fillOpacity = 0.7, group = "Fire detections"
      ) |>
      addLayersControl(
        overlayGroups = c("Provincial boundaries", "Fire detections"),
        options = layersControlOptions(collapsed = TRUE)
      ) |>
      fitBounds(
        unname(bb["xmin"]), unname(bb["ymin"]),
        unname(bb["xmax"]), unname(bb["ymax"])
      )
  })

  output$province_chart <- renderPlot({
    counts <- filtered_fires() |>
      st_drop_geometry() |>
      count(province, name = "detections")
    validate(need(nrow(counts) > 0, "No data available."))
    ggplot(counts, aes(x = reorder(province, detections),
                       y = detections)) +
      geom_col(fill = "#E85D04", width = 0.72) +
      coord_flip() +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.05))) +
      labs(x = NULL, y = "Satellite fire detections") +
      dark_plot() +
      theme(panel.grid.major.y = element_blank())
  })

  output$pattern_summary <- renderText({
    x <- filtered_fires()
    if (nrow(x) == 0) return("No detections for this selection.")
    if (input$province != "All Provinces") {
      return(paste0(fmt(nrow(x)), " detections in ",
                    input$province, " for the selected period."))
    }
    counts <- x |> st_drop_geometry() |>
      count(province, sort = TRUE)
    paste0(fmt(nrow(x)), " detections. ",
           counts$province[1], " has the highest raw count (",
           fmt(counts$n[1]), ").")
  })

  output$download_fires <- downloadHandler(
    filename = function() "firewatch_filtered_detections.csv",
    content = function(file) {
      out <- sf::st_drop_geometry(filtered_fires())
      utils::write.csv(out, file, row.names = FALSE)
    }
  )

  # ---- Monthly comparison (Overview: follows all sidebar filters) ----
  output$overview_monthly_chart <- renderPlot({
    d <- monthly_counts
    if (input$province != "All Provinces") {
      d <- d[d$province == input$province, , drop = FALSE]
    }
    if (input$year != "All Years") {
      d <- d[d$year == as.integer(input$year), , drop = FALSE]
    }
    if (input$season != "All Seasons") {
      d <- d[d$season == input$season, , drop = FALSE]
    }
    d <- d |>
      dplyr::group_by(year, month_num) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop")
    validate(need(nrow(d) > 0, "No monthly data for this selection."))
    d$year_label <- factor(as.character(d$year), levels = c("2024", "2025"))
    ggplot(d, aes(x = month_num, y = detections,
                  colour = year_label, group = year_label)) +
      geom_line(linewidth = 1.15) +
      geom_point(size = 2.6) +
      scale_colour_manual(values = c("2024" = "#E85D04",
                                     "2025" = "#087E8B"),
                          drop = FALSE) +
      scale_x_continuous(breaks = 1:12, labels = month_abbr,
                         limits = c(1, 12)) +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.08))) +
      labs(x = "Calendar month", y = "Satellite fire detections",
           colour = "Year") +
      dark_plot() +
      theme(legend.position = "top")
  })

  # ---- Seasonal Fire Intelligence: intentionally compares both years ----
  comparison_monthly <- reactive({
    d <- monthly_counts
    if (input$province != "All Provinces") {
      d <- d[d$province == input$province, , drop = FALSE]
    }
    d |>
      dplyr::group_by(year, month_num, season) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop") |>
      dplyr::arrange(year, month_num)
  })

  comparison_seasonal <- reactive({
    comparison_monthly() |>
      dplyr::group_by(year, season) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop") |>
      dplyr::mutate(season = factor(season, levels = season_labels)) |>
      dplyr::arrange(season, year)
  })

  comparison_totals <- reactive({
    d <- comparison_monthly()
    c(`2024` = sum(d$detections[d$year == 2024L]),
      `2025` = sum(d$detections[d$year == 2025L]))
  })

  output$intel_2024 <- renderText(fmt(comparison_totals()[["2024"]]))
  output$intel_2025 <- renderText(fmt(comparison_totals()[["2025"]]))
  output$intel_change <- renderText({
    t <- comparison_totals()
    format_change(t[["2024"]], t[["2025"]])
  })
  output$intel_peak <- renderText({
    d <- comparison_monthly() |>
      dplyr::group_by(month_num) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop")
    month_abbr[d$month_num[which.max(d$detections)]]
  })

  output$intel_monthly <- renderPlot({
    d <- comparison_monthly()
    d$year_label <- factor(as.character(d$year), levels = c("2024", "2025"))
    ggplot(d, aes(x = month_num, y = detections,
                  colour = year_label, group = year_label)) +
      geom_line(linewidth = 1.25) +
      geom_point(size = 2.8) +
      scale_colour_manual(values = c("2024" = "#E85D04",
                                     "2025" = "#087E8B")) +
      scale_x_continuous(breaks = 1:12, labels = month_abbr) +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.08))) +
      labs(x = "Calendar month", y = "Satellite fire detections",
           colour = "Year") +
      dark_plot() +
      theme(legend.position = "top")
  })

  output$intel_seasonal <- renderPlot({
    d <- comparison_seasonal()
    d$year_label <- factor(as.character(d$year), levels = c("2024", "2025"))
    ggplot(d, aes(x = season, y = detections, fill = year_label)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.68) +
      scale_fill_manual(values = c("2024" = "#E85D04",
                                   "2025" = "#087E8B")) +
      scale_y_continuous(labels = scales::comma,
                         expand = expansion(mult = c(0, 0.08))) +
      labs(x = NULL, y = "Satellite fire detections", fill = "Year") +
      dark_plot() +
      theme(legend.position = "top", panel.grid.major.x = element_blank())
  })

  output$intel_story <- renderUI({
    d <- comparison_seasonal()
    t <- comparison_totals()
    peak <- d |>
      dplyr::group_by(season) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(detections))
    peak_season <- as.character(peak$season[1])
    delta <- t[["2025"]] - t[["2024"]]
    direction <- if (delta > 0) "more" else if (delta < 0) "fewer" else "the same number of"
    scope <- if (input$province == "All Provinces") {
      "South Africa"
    } else {
      input$province
    }
    div(class = "fw-note",
        strong("What the data show: "),
        paste0(scope, " recorded ", fmt(abs(delta)), " ", direction,
               " satellite detections in 2025 than in 2024 (",
               format_change(t[["2024"]], t[["2025"]]), "). ",
               peak_season, " accounts for the highest combined",
               " seasonal count across both years. "),
        "These descriptive differences do not establish a change",
        " in the number of unique wildfires.")
  })

  # Heatmap uses ALL provinces so the national seasonal signatures
  # remain comparable even when the sidebar selects a single province.
  output$intel_heatmap <- renderPlot({
    d <- monthly_counts |>
      dplyr::group_by(province, season) |>
      dplyr::summarise(detections = sum(detections), .groups = "drop") |>
      dplyr::group_by(province) |>
      dplyr::mutate(share = 100 * detections / sum(detections)) |>
      dplyr::ungroup()
    d$season <- factor(d$season, levels = season_labels)
    province_order <- d |>
      dplyr::group_by(province) |>
      dplyr::summarise(total = sum(detections), .groups = "drop") |>
      dplyr::arrange(total)
    d$province <- factor(d$province, levels = province_order$province)
    ggplot(d, aes(x = season, y = province, fill = share)) +
      geom_tile(colour = "#FFFFFF", linewidth = 1.1) +
      geom_text(aes(label = sprintf("%.1f%%", share),
                    colour = share >= 35),
                size = 3.7, fontface = "bold") +
      scale_colour_manual(values = c("TRUE" = "#182235",
                                     "FALSE" = "#FFFFFF"),
                          guide = "none") +
      scale_fill_gradient(low = "#FFF1DB", high = "#D94312",
                          limits = c(0, 100),
                          labels = function(z) paste0(z, "%")) +
      labs(x = NULL, y = NULL, fill = "Share of province's detections") +
      dark_plot() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom")
  })

  output$download_intel_monthly <- downloadHandler(
    filename = function() "firewatch_monthly_comparison_2024_2025.csv",
    content = function(file) {
      d <- comparison_monthly()
      d$scope <- if (input$province == "All Provinces") {
        "South Africa"
      } else {
        input$province
      }
      utils::write.csv(d, file, row.names = FALSE)
    }
  )

  output$download_intel_seasonal <- downloadHandler(
    filename = function() "firewatch_seasonal_comparison_2024_2025.csv",
    content = function(file) {
      d <- comparison_seasonal()
      d$season <- as.character(d$season)
      d$scope <- if (input$province == "All Provinces") {
        "South Africa"
      } else {
        input$province
      }
      utils::write.csv(d, file, row.names = FALSE)
    }
  )

  # Both national and provincial analyses use the current selection.
  # Each analysis stores its filter values at the time it is run.
  analysis_ready <- function() {
    if (nrow(filtered_fires()) < 10) {
      showNotification(
        "This selection needs at least 10 detections.",
        type = "warning", duration = 6
      )
      return(FALSE)
    }
    TRUE
  }

  # ---- KDE ----
  kde_result <- eventReactive(input$run_kde, {
    req(analysis_ready())
    withProgress(message = "Calculating KDE...", value = 0, {
      x <- make_point_pattern(
        filtered_fires(), selected_boundary(), max_points = 6000
      )
      incProgress(0.3)
      surface <- spatstat.explore::density.ppp(
        x$pp,
        sigma = input$kde_bandwidth * 1000,
        edge = TRUE,
        at = "pixels",
        dimyx = 180
      )
      # Correct for uniform sampling fraction to estimate full
      # selected-data detection intensity (detections per km^2).
      surface <- surface * (1e6 * x$total / x$used)
      incProgress(0.7)
      x$surface <- surface
      x$bandwidth <- input$kde_bandwidth
      valid_intensity <- as.numeric(surface$v)
      valid_intensity <- valid_intensity[is.finite(valid_intensity)]
      x$peak_density <- if (length(valid_intensity)) max(valid_intensity) else NA_real_
      x$median_density <- if (length(valid_intensity)) median(valid_intensity) else NA_real_
      x$province <- input$province
      x$year <- input$year
      x$season <- input$season
      x
    })
  }, ignoreInit = TRUE)

  output$kde_map <- renderLeaflet({
    x <- kde_result()
    req(x)
    validate(need(filters_match(x, check_bandwidth = TRUE),
                  "Filters or bandwidth changed. Click Run KDE Analysis again."))
    boundary <- selected_boundary()
    bb <- sf::st_bbox(boundary)
    density_raster <- kde_to_wgs84_raster(x)
    vals <- terra::values(density_raster, mat = FALSE)
    vals <- vals[is.finite(vals)]
    validate(need(length(vals) > 0,
                  "No valid KDE cells were available for the map."))
    pal <- leaflet::colorNumeric(
      palette = grDevices::hcl.colors(90, "Inferno"),
      domain = range(vals), na.color = "transparent"
    )
    leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) |>
      leaflet::addTiles(group = "OpenStreetMap") |>
      leaflet::addRasterImage(
        density_raster, colors = pal, opacity = 0.72,
        project = FALSE, group = "KDE density"
      ) |>
      leaflet::addPolygons(
        data = boundary, color = "#182235", weight = 1.5,
        fill = FALSE, group = "Province boundaries"
      ) |>
      leaflet::addLegend(
        position = "bottomright", pal = pal, values = vals,
        title = "Detections / km²", opacity = 0.85
      ) |>
      leaflet::addLayersControl(
        overlayGroups = c("KDE density", "Province boundaries"),
        options = leaflet::layersControlOptions(collapsed = TRUE)
      ) |>
      leaflet::fitBounds(
        unname(bb["xmin"]), unname(bb["ymin"]),
        unname(bb["xmax"]), unname(bb["ymax"])
      )
  })

  output$download_kde_png <- downloadHandler(
    filename = function() "firewatch_kde.png",
    content = function(file) {
      x <- kde_result()
      req(x)
      validate(need(filters_match(x, check_bandwidth = TRUE),
                    "Rerun KDE before downloading."))
      grDevices::png(file, width = 1600, height = 1200, res = 160)
      on.exit(grDevices::dev.off(), add = TRUE)
      draw_kde(x)
    }
  )

  output$kde_explanation <- renderUI({
    x <- kde_result()
    req(x)
    validate(need(filters_match(x, check_bandwidth = TRUE),
                  "Filters changed. Rerun KDE to update this result."))
    scope <- if (x$province == "All Provinces") "South Africa" else x$province
    peak_text <- if (is.finite(x$peak_density)) {
      paste0("The highest modelled grid-cell intensity is approximately ",
             fmt(x$peak_density, 2), " detections per km². ")
    } else ""
    div(class = "fw-note",
        h4("What do these KDE results mean?"),
        p(strong("Finding: "),
          "The brighter parts of the surface are the locations where",
          "satellite detections are most concentrated in ", scope,
          " for the selected year and season. ", peak_text),
        p(strong("How to read the map: "),
          "KDE spreads the contribution of each detection over nearby",
          "space. A bandwidth of ", x$bandwidth,
          " km controls smoothing: smaller values show local peaks,",
          "while larger values show broader concentrations."),
        p(strong("What it does NOT show: "),
          "This is a smoothed estimate of detection intensity, not",
          "a statistical significance test, unique-fire count,",
          "burned-area map or measure of fire severity.",
          "The surface is scaled from a reproducible sample."),
        p(strong("Analysis details: "), stat_note(x)),
        p("Filters at run time: ", x$province, " | ",
          x$year, " | ", x$season, "."))
  })

  # ---- Nearest neighbour and Monte Carlo ----
  nn_result <- eventReactive(input$run_nn, {
    req(analysis_ready())
    withProgress(message = "Running nearest-neighbour simulations...",
                 value = 0, {
      x <- make_point_pattern(
        filtered_fires(), selected_boundary(), max_points = 2000
      )
      pp <- x$pp
      n <- npoints(pp)
      observed <- mean(spatstat.geom::nndist(pp))
      expected <- 1 / (2 * sqrt(n / area.owin(Window(pp))))

      set.seed(123)
      nsim <- 99
      simulated <- numeric(nsim)
      for (i in seq_len(nsim)) {
        random_pp <- spatstat.random::runifpoint(
          n = n, win = Window(pp)
        )
        simulated[i] <- mean(spatstat.geom::nndist(random_pp))
        if (i %% 10 == 0) incProgress(10 / nsim)
      }
      x$observed <- observed
      x$expected <- expected
      x$ratio <- observed / expected
      x$simulated <- simulated
      x$p_value <- (1 + sum(simulated <= observed)) / (nsim + 1)
      x$province <- input$province
      x$year <- input$year
      x$season <- input$season
      x
    })
  }, ignoreInit = TRUE)

  output$nn_summary <- renderUI({
    x <- nn_result()
    req(x)
    validate(need(filters_match(x),
                  "Filters changed. Rerun Nearest Neighbour Analysis."))
    div(class = "fw-note",
        tags$b("Observed mean NN distance: "),
        paste0(fmt(x$observed / 1000, 3), " km"),
        tags$br(),
        tags$b("Expected mean under planar CSR approximation: "),
        paste0(fmt(x$expected / 1000, 3), " km"),
        tags$br(),
        tags$b("Nearest-neighbour ratio: "),
        fmt(x$ratio, 3),
        tags$br(),
        tags$b("Monte Carlo clustering p-value: "),
        fmt(x$p_value, 3))
  })

  output$nn_plot <- renderPlot({
    x <- nn_result()
    req(x)
    validate(need(filters_match(x),
                  "Filters changed. Rerun Nearest Neighbour Analysis."))
    sim <- data.frame(distance_km = x$simulated / 1000)
    limits <- range(c(x$simulated, x$observed)) / 1000
    pad <- max(0.03, diff(limits) * 0.07)
    ggplot(sim, aes(x = distance_km)) +
      geom_histogram(bins = 24, fill = "#FFB703",
                     colour = "#FFFFFF") +
      geom_vline(xintercept = x$observed / 1000,
                 colour = "#E63946", linewidth = 1.2) +
      coord_cartesian(xlim = limits + c(-pad, pad)) +
      labs(title = "Observed vs simulated random patterns",
           subtitle = "Red line: observed mean | Orange: CSR simulations",
           x = "Mean nearest-neighbour distance (km)",
           y = "Number of simulations") +
      dark_plot()
  })

  output$nn_explanation <- renderUI({
    x <- nn_result()
    req(x)
    validate(need(filters_match(x),
                  "Filters changed. Rerun Nearest Neighbour Analysis."))
    ratio_class <- if (x$ratio < 0.95) "clustered" else
      if (x$ratio > 1.05) "more dispersed" else "close to random"
    percent <- 100 * abs(1 - x$ratio)
    ratio_meaning <- if (x$ratio < 1) {
      paste0("The observed mean distance is about ", fmt(percent, 1),
             "% shorter than the simple CSR expectation.")
    } else if (x$ratio > 1) {
      paste0("The observed mean distance is about ", fmt(percent, 1),
             "% longer than the simple CSR expectation.")
    } else "The observed mean distance matches the planar CSR approximation."
    inference <- if (x$p_value < 0.05) {
      paste0("The one-sided Monte Carlo p-value is ", fmt(x$p_value, 3),
             ", below 0.05. This is evidence of clustering relative",
             " to uniform complete spatial randomness (CSR).")
    } else {
      paste0("The one-sided Monte Carlo p-value is ", fmt(x$p_value, 3),
             ". At the 5% level, this analysis does not provide",
             " sufficient evidence of clustering relative to uniform CSR.")
    }
    div(class = "fw-note",
        h4("What do these Nearest Neighbour results mean?"),
        p(strong("Observed pattern: "),
          "The sample is ", ratio_class,
          " relative to the simple nearest-neighbour benchmark. ",
          ratio_meaning),
        p(strong("Actual distances: "),
          "The average observed nearest-neighbour distance is ",
          fmt(x$observed / 1000, 3), " km, compared with ",
          fmt(x$expected / 1000, 3), " km expected under the",
          " infinite-plane CSR approximation."),
        p(strong("Statistical evidence: "), inference),
        p(strong("Important limitation: "),
          "The p-value is a one-sided test for clustering, not a",
          "test for dispersion. The random model assumes points",
          "could occur uniformly within the chosen boundary.",
          "Repeated satellite detections, vegetation distribution,",
          "and observation coverage can also affect clustering."),
        p(strong("Analysis details: "), stat_note(x)),
        p("Filters at run time: ", x$province, " | ",
          x$year, " | ", x$season, "."))
  })

  output$download_nn <- downloadHandler(
    filename = function() "firewatch_nearest_neighbour.csv",
    content = function(file) {
      x <- nn_result()
      req(x)
      validate(need(filters_match(x), "Rerun analysis before downloading."))
      results <- data.frame(
        province = x$province,
        year = x$year,
        season = x$season,
        detections_filtered = x$total,
        detections_analyzed = x$used,
        duplicate_coordinates_in_sample = x$duplicate_locations,
        projected_crs = x$projection_label,
        observed_mean_km = x$observed / 1000,
        expected_mean_km = x$expected / 1000,
        nearest_neighbour_ratio = x$ratio,
        one_sided_csr_monte_carlo_p = x$p_value,
        simulations = length(x$simulated)
      )
      utils::write.csv(results, file, row.names = FALSE)
    }
  )

  # ---- Ripley's K and transformed L with CSR envelopes ----
  k_result <- eventReactive(input$run_k, {
    req(analysis_ready())
    withProgress(message = "Calculating Ripley's K envelopes...",
                 value = 0, {
      x <- make_point_pattern(
        filtered_fires(), selected_boundary(), max_points = 800
      )
      incProgress(0.2)
      set.seed(123)
      distances <- seq(0, 10000, by = 500)
      env <- spatstat.explore::envelope(
        x$pp, fun = spatstat.explore::Kest,
        nsim = 19, r = distances,
        correction = "border", verbose = FALSE
      )
      incProgress(0.8)
      r <- env$r
      transform_l <- function(k) (sqrt(pmax(k, 0) / pi) - r) / 1000
      x$curves <- data.frame(
        distance_km = r / 1000,
        observed = transform_l(env$obs),
        lower = transform_l(env$lo),
        upper = transform_l(env$hi)
      )
      x$province <- input$province
      x$year <- input$year
      x$season <- input$season
      x
    })
  }, ignoreInit = TRUE)

  output$k_plot <- renderPlot({
    x <- k_result()
    req(x)
    validate(need(filters_match(x),
                  "Filters changed. Rerun Ripley's K Analysis."))
    d <- x$curves
    ggplot(d, aes(x = distance_km)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  fill = "#FFE29A", alpha = 0.65) +
      geom_hline(yintercept = 0, colour = "#64748B",
                 linetype = "dashed", linewidth = 0.8) +
      geom_line(aes(y = observed), colour = "#E85D04",
                linewidth = 1.3) +
      labs(title = "Clustering across spatial distances",
           subtitle = "Orange: observed L(r)-r | Yellow: pointwise CSR envelope",
           x = "Distance (km)", y = "L(r) - r (km)") +
      dark_plot()
  })

  output$download_k <- downloadHandler(
    filename = function() "firewatch_ripleys_l_curves.csv",
    content = function(file) {
      x <- k_result()
      req(x)
      validate(need(filters_match(x), "Rerun analysis before downloading."))
      out <- x$curves
      out$province <- x$province
      out$year <- x$year
      out$season <- x$season
      out$detections_analyzed <- x$used
      out$simulations <- 19L
      utils::write.csv(out, file, row.names = FALSE)
    }
  )

  output$k_explanation <- renderUI({
    x <- k_result()
    req(x)
    validate(need(filters_match(x),
                  "Filters changed. Rerun Ripley's K Analysis."))
    d <- x$curves
    positive <- d$distance_km > 0
    above <- sum(d$observed[positive] > d$upper[positive],
                 na.rm = TRUE)
    tested <- sum(positive & is.finite(d$observed) &
                    is.finite(d$upper))
    below <- sum(d$observed[positive] < d$lower[positive], na.rm = TRUE)
    above_zero <- sum(d$observed[positive] > 0, na.rm = TRUE)
    peak_idx <- if (any(positive & is.finite(d$observed))) {
      which.max(replace(d$observed, !positive | !is.finite(d$observed), -Inf))
    } else NA_integer_
    peak_sentence <- if (is.finite(peak_idx) && d$observed[peak_idx] > 0) {
      paste0("The largest positive L(r)-r deviation occurs around ",
             fmt(d$distance_km[peak_idx], 1), " km in the plotted range.")
    } else if (is.finite(peak_idx)) {
      "No positive L(r)-r deviation appears at the evaluated distances."
    } else "The curve could not be summarised at the plotted distances."
    div(class = "fw-note",
        h4("What do these Ripley's K results mean?"),
        p(strong("Overall pattern: "),
          "The observed L(r)-r curve is positive at ", above_zero,
          " of ", tested, " evaluated positive distances.",
          "Positive values mean more nearby pairs than expected",
          "under uniform CSR at that distance; negative values",
          "mean fewer nearby pairs."),
        p(strong("Comparison with random simulations: "),
          "The curve lies above the pointwise CSR envelope at ",
          above, " of ", tested, " evaluated distances and below",
          " it at ", below, " distances. ", peak_sentence),
        p(strong("How to interpret this: "),
          "If the orange curve rises above the yellow band, the",
          "observed clustering at that distance is stronger than",
          "the simulated pointwise CSR range. Values within the",
          "band are compatible with that reference range."),
        p(strong("Important limitation: "),
          "The envelope uses only 19 random simulations and is",
          "pointwise. Counting distances outside the envelope",
          "is descriptive and is NOT a global significance test",
          "or a global p-value. Multiple distances are related,",
          "and uniform CSR ignores environmental heterogeneity."),
        p(strong("Analysis details: "), stat_note(x)),
        p("Filters at run time: ", x$province, " | ",
          x$year, " | ", x$season, "."))
  })
}

shinyApp(ui = ui, server = server)

#!/usr/bin/env Rscript

# ============================================================
# updated_ok_induced_seismicity_regional20150318.R
#
# Data prep and visualization script for Oklahoma induced-seismicity
# using the 2015 regional directive AOI_20150318 (layer 2).  This
# large directive covers a broad area of Oklahoma, so using it
# yields hundreds of events both inside and outside the treated
# region.  Treatment is applied on 18 March 2015 and the post
# window ends on 24 June 2015, just before the next regional
# directive (AOI_20150624) becomes active.  A 300 km buffer around
# the AOI ensures plenty of control observations.  Cells are
# designated as treated if their centroid lies within the AOI.
# A 20 km tessellation grid strikes a balance between spatial
# resolution and sufficient counts per cell.
#
# This script requires the same packages as the other data‑prep
# scripts (httr, jsonlite, sf, data.table, dplyr, lubridate, plotly).
# It writes outputs into the directory
# "oklahoma_induced_seismicity_data_regional20150318".
# ============================================================

req_pkgs <- c("httr", "jsonlite", "sf", "dplyr", "lubridate", "data.table", "plotly")
missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(sf)
  library(dplyr)
  library(lubridate)
  library(data.table)
  library(plotly)
})

sf::sf_use_s2(FALSE)

# Parameters for the 2015 regional directive
params <- list(
  out_dir = "oklahoma_induced_seismicity_data_regional20150318",
  occ_aoi_layer_id = 2,             # AOI_20150318
  pre_start = "2014-01-01T00:00:00Z",
  t_star    = "2015-03-18T00:00:00Z",
  post_end  = "2015-06-24T00:00:00Z",
  min_mag = 2.5,
  keep_reviewed_only = FALSE,
  buffer_km = 300,                  # larger buffer for more controls
  crs_proj = 5070,
  grid_shape = "square",
  grid_cellsize_km = 20,
  treat_rule = "centroid",         # treat cell if centroid inside AOI
  treat_intersect_threshold = 0.0,   # unused for centroid rule
  force_download = FALSE
)

# Helper functions
safe_mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
iso_no_z <- function(x_posix) {
  format(x_posix, "%Y-%m-%dT%H:%M:%S")
}
occ_base_service <- function() {
  "https://gis.occ.ok.gov/server/rest/services/PUBLIC/DIRECTIVE_AOIs/MapServer"
}
occ_layer_query_url <- function(layer_id) {
  paste0(occ_base_service(), "/", layer_id, "/query")
}
download_occ_aoi <- function(layer_id, out_geojson, force = FALSE) {
  if (file.exists(out_geojson) && !force) {
    message("AOI geojson exists, reading: ", out_geojson)
    return(sf::st_read(out_geojson, quiet = TRUE))
  }
  message("Downloading OCC AOI layer ", layer_id, " -> ", out_geojson)
  q <- list(
    where = "1=1",
    outFields = "*",
    returnGeometry = "true",
    outSR = "4326",
    f = "geojson"
  )
  resp <- httr::RETRY(
    "GET",
    occ_layer_query_url(layer_id),
    query = q,
    times = 5,
    pause_base = 2,
    httr::config(connecttimeout = 60, ipresolve = 1),
    httr::timeout(180)
  )
  httr::stop_for_status(resp)
  writeBin(httr::content(resp, "raw"), out_geojson)
  sf::st_read(out_geojson, quiet = TRUE)
}
usgs_base_url <- function() {
  "https://earthquake.usgs.gov/fdsnws/event/1/query"
}
fetch_usgs_page <- function(starttime, endtime, minmag, bbox, limit = 20000, offset = 1, eventtype = "earthquake") {
  q <- list(
    format = "geojson",
    eventtype = eventtype,
    starttime = starttime,
    endtime = endtime,
    minmagnitude = minmag,
    minlatitude  = bbox$minlat,
    maxlatitude  = bbox$maxlat,
    minlongitude = bbox$minlon,
    maxlongitude = bbox$maxlon,
    orderby = "time-asc",
    limit = limit,
    offset = offset
  )
  resp <- httr::GET(usgs_base_url(), query = q, httr::timeout(120))
  httr::stop_for_status(resp)
  jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
}
parse_usgs_geojson <- function(gj) {
  feats <- gj$features
  if (is.null(feats) || length(feats) == 0) {
    return(data.frame())
  }
  get_prop <- function(f, nm) {
    v <- f$properties[[nm]]
    if (is.null(v)) NA else v
  }
  get_coord <- function(f, idx) {
    coords <- f$geometry$coordinates
    if (is.null(coords) || length(coords) < idx) NA else coords[[idx]]
  }
  data.frame(
    event_id = vapply(feats, function(f) f$id, FUN.VALUE = character(1)),
    time_utc = as.POSIXct(vapply(feats, function(f) get_prop(f, "time") / 1000, FUN.VALUE = numeric(1)),
                          origin = "1970-01-01", tz = "UTC"),
    mag      = vapply(feats, function(f) as.numeric(get_prop(f, "mag")), FUN.VALUE = numeric(1)),
    status   = vapply(feats, function(f) as.character(get_prop(f, "status")), FUN.VALUE = character(1)),
    place    = vapply(feats, function(f) as.character(get_prop(f, "place")), FUN.VALUE = character(1)),
    lon      = vapply(feats, function(f) as.numeric(get_coord(f, 1)), FUN.VALUE = numeric(1)),
    lat      = vapply(feats, function(f) as.numeric(get_coord(f, 2)), FUN.VALUE = numeric(1)),
    depth_km = vapply(feats, function(f) as.numeric(get_coord(f, 3)), FUN.VALUE = numeric(1)),
    stringsAsFactors = FALSE
  )
}
fetch_usgs_all <- function(starttime, endtime, minmag, bbox, limit = 20000) {
  offset <- 1
  out <- list()
  i <- 1
  repeat {
    message(sprintf("USGS query page %d (offset=%d, limit=%d)...", i, offset, limit))
    gj <- fetch_usgs_page(starttime = starttime, endtime = endtime, minmag = minmag,
                          bbox = bbox, limit = limit, offset = offset)
    df <- parse_usgs_geojson(gj)
    if (nrow(df) == 0) break
    out[[i]] <- df
    if (nrow(df) < limit) break
    offset <- offset + limit
    i <- i + 1
  }
  dplyr::bind_rows(out)
}

# Parse times
pre_start <- ymd_hms(params$pre_start, tz = "UTC")
t_star    <- ymd_hms(params$t_star, tz = "UTC")
post_end  <- ymd_hms(params$post_end, tz = "UTC")
stopifnot(pre_start < t_star, t_star < post_end)

# Prepare output directory
safe_mkdir(params$out_dir)

aoi_geojson_path <- file.path(params$out_dir, sprintf("occ_aoi_layer_%s.geojson", params$occ_aoi_layer_id))
aoi_gpkg_path    <- file.path(params$out_dir, "aoi.gpkg")
win_gpkg_path    <- file.path(params$out_dir, "analysis_window.gpkg")
grid_gpkg_path   <- file.path(params$out_dir, "grid_cells.gpkg")

events_raw_csv   <- file.path(params$out_dir, "events_raw_usgs.csv")
events_all_csv   <- file.path(params$out_dir, "events_all.csv")
events_pre_csv   <- file.path(params$out_dir, "events_pre.csv")
events_post_csv  <- file.path(params$out_dir, "events_post.csv")
cells_csv        <- file.path(params$out_dir, "cells.csv")
meta_json        <- file.path(params$out_dir, "metadata.json")

# Download AOI
aoi_wgs <- download_occ_aoi(params$occ_aoi_layer_id, aoi_geojson_path, force = params$force_download)
aoi_wgs <- sf::st_make_valid(aoi_wgs)
aoi_union_wgs <- sf::st_union(sf::st_geometry(aoi_wgs))
aoi_union_wgs <- sf::st_sf(aoi_id = params$occ_aoi_layer_id, geometry = aoi_union_wgs, crs = 4326)
aoi_union_proj <- sf::st_transform(aoi_union_wgs, params$crs_proj)

# Spatial window
window_proj <- sf::st_buffer(aoi_union_proj, dist = params$buffer_km * 1000)
window_proj <- sf::st_make_valid(window_proj)
window_wgs <- sf::st_transform(window_proj, 4326)

bb <- sf::st_bbox(window_wgs)
pad_deg <- 0.05
bbox <- list(
  minlon = as.numeric(bb["xmin"]) - pad_deg,
  maxlon = as.numeric(bb["xmax"]) + pad_deg,
  minlat = as.numeric(bb["ymin"]) - pad_deg,
  maxlat = as.numeric(bb["ymax"]) + pad_deg
)
message(sprintf("USGS bbox: lon[%.3f, %.3f], lat[%.3f, %.3f]",
                bbox$minlon, bbox$maxlon, bbox$minlat, bbox$maxlat))

# Save AOI and window
suppressWarnings({
  sf::st_write(sf::st_transform(aoi_union_proj, 4326), aoi_gpkg_path,
               layer = "aoi", delete_layer = TRUE, quiet = TRUE)
  sf::st_write(sf::st_transform(window_proj, 4326), win_gpkg_path,
               layer = "window", delete_layer = TRUE, quiet = TRUE)
})

# Download earthquakes
if (file.exists(events_raw_csv) && !params$force_download) {
  message("Raw events CSV exists, reading: ", events_raw_csv)
  ev_raw <- fread(events_raw_csv)
  ev_raw[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
} else {
  ev_raw <- fetch_usgs_all(
    starttime = iso_no_z(pre_start),
    endtime   = iso_no_z(post_end),
    minmag    = params$min_mag,
    bbox      = bbox,
    limit     = 20000
  )
  if (nrow(ev_raw) == 0) stop("USGS query returned 0 events.")
  fwrite(ev_raw, events_raw_csv)
}

# Clean and filter events

ev_raw <- as.data.table(ev_raw)
ev_raw <- unique(ev_raw, by = "event_id")
ev_raw[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
if (params$keep_reviewed_only) {
  ev_raw <- ev_raw[tolower(status) == "reviewed"]
}

# Filter to buffered window

ev_sf <- sf::st_as_sf(ev_raw, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
ev_proj <- sf::st_transform(ev_sf, params$crs_proj)
inside_idx <- lengths(sf::st_within(ev_proj, window_proj)) > 0
ev_proj <- ev_proj[inside_idx, ]
if (nrow(ev_proj) == 0) stop("No events within buffered window.")

# Build grid
cellsize_m <- params$grid_cellsize_km * 1000
square_flag <- TRUE
if (tolower(params$grid_shape) == "hex") square_flag <- FALSE
grid_sfc <- sf::st_make_grid(
  window_proj,
  cellsize = cellsize_m,
  what = "polygons",
  square = square_flag
)
grid_sf <- sf::st_sf(cell_id = seq_along(grid_sfc), geometry = grid_sfc, crs = params$crs_proj)

hits2 <- lengths(sf::st_intersects(grid_sf, window_proj)) > 0
grid_sf <- grid_sf[hits2, ]

grid_clip <- suppressWarnings(sf::st_intersection(grid_sf, window_proj))

grid_clip <- sf::st_make_valid(grid_clip)
areas2 <- as.numeric(sf::st_area(grid_clip))

# drop tiny geometries

grid_clip <- grid_clip[areas2 > 1, ]

# Assign treatment using centroid rule
centroids <- sf::st_centroid(grid_clip)
inside_aoi <- lengths(sf::st_within(centroids, aoi_union_proj)) > 0
Z <- as.integer(inside_aoi)

grid_clip$Z <- Z

# Join events to grid cells

ev_join <- sf::st_join(ev_proj, grid_clip[, c("cell_id", "Z")], join = sf::st_within, left = FALSE)

if (nrow(ev_join) == 0) stop("No events joined to grid cells.")

coords_join <- sf::st_coordinates(ev_join)
ev_tbl <- as.data.table(sf::st_drop_geometry(ev_join))
ev_tbl[, x_m := coords_join[, 1]]
ev_tbl[, y_m := coords_join[, 2]]
ev_tbl[, is_post := as.integer(time_utc >= t_star & time_utc < post_end)]
ev_tbl[, t_sec_pre0  := as.numeric(difftime(time_utc, pre_start, units = "secs"))]
ev_tbl[, t_sec_post0 := as.numeric(difftime(time_utc, t_star,   units = "secs"))]
ev_tbl <- ev_tbl[time_utc >= pre_start & time_utc < post_end]

# Split pre and post

ev_pre  <- ev_tbl[time_utc <  t_star]
ev_post <- ev_tbl[time_utc >= t_star & time_utc < post_end]

message(sprintf("Events: total=%d, pre=%d, post=%d", nrow(ev_tbl), nrow(ev_pre), nrow(ev_post)))
message(sprintf("Cells: %d (treated: %d, control: %d)", nrow(grid_clip), sum(grid_clip$Z == 1L), sum(grid_clip$Z == 0L)))
if (nrow(ev_post) > 0) {
  post_counts <- ev_post[, .N, by = Z]
  setnames(post_counts, "Z", "CellTreatment")
  message("Post‑treatment event counts by cell treatment:\n", paste0(capture.output(print(post_counts)), collapse = "\n"))
}

# Save outputs
fwrite(ev_tbl,  events_all_csv)
fwrite(ev_pre,  events_pre_csv)
fwrite(ev_post, events_post_csv)

cent_proj <- sf::st_centroid(grid_clip)
cent_xy <- sf::st_coordinates(cent_proj)
cent_ll <- sf::st_transform(cent_proj, 4326)
cent_ll_xy <- sf::st_coordinates(cent_ll)

cells_tbl <- data.table(
  cell_id = grid_clip$cell_id,
  Z = grid_clip$Z,
  area_km2 = as.numeric(sf::st_area(grid_clip)) / 1e6,
  centroid_x_m = cent_xy[, 1],
  centroid_y_m = cent_xy[, 2],
  centroid_lon = cent_ll_xy[, 1],
  centroid_lat = cent_ll_xy[, 2]
)
fwrite(cells_tbl, cells_csv)

suppressWarnings({
  sf::st_write(sf::st_transform(grid_clip, 4326), grid_gpkg_path,
               layer = "cells", delete_layer = TRUE, quiet = TRUE)
})

# Metadata
meta <- list(
  created_utc = format(with_tz(Sys.time(), "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
  sources = list(
    usgs_event_api = usgs_base_url(),
    occ_aoi_service = occ_base_service()
  ),
  design = list(
    occ_aoi_layer_id = params$occ_aoi_layer_id,
    pre_start_utc = params$pre_start,
    t_star_utc = params$t_star,
    post_end_utc = params$post_end,
    min_mag = params$min_mag,
    keep_reviewed_only = params$keep_reviewed_only,
    buffer_km = params$buffer_km,
    crs_proj_epsg = params$crs_proj,
    grid_shape = params$grid_shape,
    grid_cellsize_km = params$grid_cellsize_km,
    treat_rule = params$treat_rule,
    treat_intersect_threshold = params$treat_intersect_threshold
  ),
  bbox_wgs84 = bbox,
  counts = list(
    events_all = nrow(ev_tbl),
    events_pre = nrow(ev_pre),
    events_post = nrow(ev_post),
    cells = nrow(cells_tbl),
    treated_cells = sum(cells_tbl$Z == 1L)
  )
)
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), meta_json)

# Visualization
# Build interactive 3D plot colored by treatment status

# Prepare for plotting

ev_vis <- copy(ev_tbl)
ev_vis[, t_days := as.numeric(difftime(time_utc, t_star, units = "days"))]
ev_vis[, x_km := x_m / 1000]
ev_vis[, y_km := y_m / 1000]
ev_vis[, trt_status := ifelse(Z == 1L, "treated", "control")]
##
## Define a color grouping for plotting:
##  - All events occurring before t_star are considered control, regardless of cell.
##  - Post-treatment events in treated cells (Z == 1) are labelled as treated.
##  - Post-treatment events in control cells (Z == 0) remain labelled as control.
ev_vis[, color_group := ifelse(time_utc < t_star, "control",
                               ifelse(Z == 1L, "treated", "control"))]

# Convert grid and AOI boundaries to plotly-ready coordinates
sf_lines_to_xy <- function(lines_sfc) {
  co <- as.data.frame(st_coordinates(lines_sfc))
  if (!("L2" %in% names(co))) co$L2 <- 1L
  co$grp <- interaction(co$L1, co$L2, drop = TRUE)
  co$ord <- ave(co$X, co$grp, FUN = seq_along)
  co
}
pack_lines <- function(co) {
  co <- co[order(co$grp, co$ord), ]
  spl <- split(co, co$grp)
  x <- unlist(lapply(spl, function(s) c(s$X / 1000, NA)))
  y <- unlist(lapply(spl, function(s) c(s$Y / 1000, NA)))
  list(x = x, y = y)
}

bnd_cells <- st_boundary(grid_clip) |> st_cast("MULTILINESTRING")
co_cells <- sf_lines_to_xy(bnd_cells)
co_cells$Z <- grid_clip$Z[co_cells$L1]
lines_Z0 <- pack_lines(subset(co_cells, Z == 0))
lines_Z1 <- pack_lines(subset(co_cells, Z == 1))
bnd_aoi <- st_boundary(aoi_union_proj) |> st_cast("MULTILINESTRING")
co_aoi <- sf_lines_to_xy(bnd_aoi)
lines_aoi <- pack_lines(co_aoi)

p <- plot_ly() %>%
  # Draw grid boundaries for control and treated cells at z=0
  add_trace(type = "scatter3d", mode = "lines",
            x = lines_Z0$x, y = lines_Z0$y, z = rep(0, length(lines_Z0$x)),
            name = "Tessellation (control cells)",
            line = list(width = 1), opacity = 0.15, showlegend = TRUE) %>%
  add_trace(type = "scatter3d", mode = "lines",
            x = lines_Z1$x, y = lines_Z1$y, z = rep(0, length(lines_Z1$x)),
            name = "Tessellation (treated cells)",
            line = list(width = 3), opacity = 0.5, showlegend = TRUE) %>%
  # AOI boundary at treatment time (z=0)
  add_trace(type = "scatter3d", mode = "lines",
            x = lines_aoi$x, y = lines_aoi$y, z = rep(0, length(lines_aoi$x)),
            name = "AOI boundary at t*", line = list(width = 6), opacity = 0.8, showlegend = TRUE) %>%
  # Control points: either pre-treatment or post-treatment in control cells
  add_markers(data = ev_vis[color_group == "control"], x = ~x_km, y = ~y_km, z = ~t_days,
               name = "Control events", marker = list(size = 2, color = "#1f77b4"), opacity = 0.6) %>%
  # Treated points: post-treatment in treated cells only
  add_markers(data = ev_vis[color_group == "treated"], x = ~x_km, y = ~y_km, z = ~t_days,
               name = "Treated events", marker = list(size = 2, color = "#d62728"), opacity = 0.85) %>%
  layout(scene = list(
             xaxis = list(title = "X (km, projected)"),
             yaxis = list(title = "Y (km, projected)"),
             zaxis = list(title = sprintf("Days since t* (%s)", format(t_star, "%Y-%m-%d")))
         ),
         legend = list(orientation = "h"))

if (interactive()) {
  print(p)
}

message("Processing complete.  Outputs written to: ", normalizePath(params$out_dir))

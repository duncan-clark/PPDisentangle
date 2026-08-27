# Support routines for oklahoma_report.qmd when results RDS lacks pp_data / counties
# (e.g. older pre-bootstrap checkpoints). Mirrors oklahoma_analysis.R Steps 1–2 + pre split.

local({
  # Source shared geometry helpers from the same directory when available.
  this_file <- tryCatch(
    normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
    error = function(e) NA_character_
  )
  candidates <- c(
    if (is.character(this_file) && nzchar(this_file)) {
      file.path(dirname(this_file), "oklahoma_geometry.R")
    } else {
      character(0)
    },
    file.path(getwd(), "inst/oklahoma/oklahoma_geometry.R"),
    file.path(getwd(), "oklahoma_geometry.R")
  )
  if (!exists("oklahoma_sf_to_owin_km", mode = "function", inherits = TRUE)) {
    found <- candidates[file.exists(candidates)]
    if (!length(found)) {
      stop("Could not locate oklahoma_geometry.R; source it before oklahoma_report_support.R.")
    }
    source(found[[1]], local = FALSE)
  }
})

oklahoma_report_rebuild_pp_and_counties <- function(data_dir, etas_m0, crs_proj = 5070L) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf required")
  if (!requireNamespace("tigris", quietly = TRUE)) stop("tigris required")
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) stop("spatstat required")
  if (!exists("oklahoma_sf_to_owin_km", mode = "function", inherits = TRUE)) {
    stop("oklahoma_geometry.R helpers are not loaded.")
  }

  meta <- jsonlite::fromJSON(readLines(file.path(data_dir, "metadata.json")))
  ev_all <- data.table::fread(file.path(data_dir, "events_all.csv"))
  t_star_utc <- as.POSIXct(meta$design$t_star_utc, tz = "UTC")
  ev_all[, time_utc := as.POSIXct(time_utc, tz = "UTC")]
  ev_all[, t_days := as.numeric(difftime(time_utc, t_star_utc, units = "days"))]
  ev_all[, x_km := x_m / 1000]
  ev_all[, y_km := y_m / 1000]
  ev_all <- ev_all[mag >= etas_m0]

  post_end_days <- as.numeric(difftime(
    as.POSIXct(meta$design$post_end_utc, tz = "UTC"), t_star_utc, units = "days"))

  options(tigris_use_cache = TRUE)
  counties_sf <- tigris::counties(state = "OK", cb = TRUE, year = 2022)
  counties_sf <- sf::st_transform(counties_sf, crs_proj)
  counties_sf <- sf::st_make_valid(counties_sf)

  ok_boundary <- sf::st_union(counties_sf)
  ok_boundary <- sf::st_make_valid(ok_boundary)
  win_km <- oklahoma_sf_to_owin_km(ok_boundary)

  county_owins <- oklahoma_sf_features_to_owins_km(counties_sf, name_col = "NAME")
  valid_idx <- !vapply(county_owins, is.null, logical(1))
  county_owins_valid <- county_owins[valid_idx]
  counties_sf_valid <- counties_sf[valid_idx, ]
  partition <- spatstat.geom::tess(tiles = county_owins_valid, window = win_km)

  aoi_sf <- sf::st_read(file.path(data_dir, "occ_aoi_layer_2.geojson"), quiet = TRUE)
  aoi_sf <- sf::st_transform(aoi_sf, crs_proj)
  aoi_sf <- sf::st_make_valid(aoi_sf)
  aoi_union <- sf::st_union(aoi_sf)
  county_centroids <- sf::st_centroid(counties_sf_valid)
  inside_aoi <- lengths(sf::st_within(county_centroids, aoi_union)) > 0
  partition_processes <- ifelse(inside_aoi, "treated", "control")
  names(partition_processes) <- counties_sf_valid$NAME
  treated_idx <- partition_processes == "treated"
  treated_names <- names(partition_processes)[treated_idx]
  control_ss <- spatstat.geom::as.owin(partition[!treated_idx])
  treated_ss <- spatstat.geom::as.owin(partition[treated_idx])

  assign_county <- function(df) {
    df$location_process <- oklahoma_assign_partition_process(
      df$x, df$y, partition, partition_processes
    )
    df$W <- 1.0
    df$n <- nrow(df)
    df$background <- TRUE
    df
  }

  pp_pre <- assign_county(as.data.frame(ev_all[t_days < 0,
    list(x = x_km, y = y_km, t = t_days, mag = mag)]))
  pp_post <- assign_county(as.data.frame(ev_all[t_days >= 0 & t_days <= post_end_days,
    list(x = x_km, y = y_km, t = t_days, mag = mag)]))
  pp_pre <- pp_pre[!is.na(pp_pre$location_process), , drop = FALSE]
  pp_post <- pp_post[!is.na(pp_post$location_process), , drop = FALSE]
  oklahoma_assert_label_support(pp_pre, control_ss, treated_ss, context = "pp_pre")
  oklahoma_assert_label_support(pp_post, control_ss, treated_ss, context = "pp_post")

  pp_pre_all <- pp_pre[order(pp_pre$t), , drop = FALSE]
  n_pre_total <- nrow(pp_pre_all)
  n_pre_holdout <- floor(n_pre_total * 0.5)
  if (n_pre_total > 0 && n_pre_holdout < 1) n_pre_holdout <- 1
  holdout_idx <- if (n_pre_holdout > 0) seq_len(n_pre_holdout) else integer(0)
  keep_idx <- if (n_pre_total > n_pre_holdout) (n_pre_holdout + 1):n_pre_total else integer(0)
  if (length(keep_idx) < 1) {
    stop("Pre-treatment split left zero estimation events; cannot rebuild pp_data for report.")
  }
  pp_pre_holdout <- pp_pre_all[holdout_idx, , drop = FALSE]
  pp_pre <- pp_pre_all[keep_idx, , drop = FALSE]

  pp_pre$process <- "control"
  pp_post$process <- pp_post$location_process
  pp_pre$inferred_process <- "control"
  pp_post$inferred_process <- pp_post$location_process

  counties <- list(
    names = counties_sf_valid$NAME,
    treated_names = treated_names,
    n_counties = partition$n,
    n_treated = sum(treated_idx)
  )
  pp_data <- list(
    pp_pre = pp_pre,
    pp_pre_holdout = pp_pre_holdout,
    pp_post = pp_post
  )
  list(counties = counties, pp_data = pp_data, meta = meta, post_end_days = post_end_days)
}

oklahoma_report_grid_partition <- function(diameter, win, aoi_owin, label, max_tiles = 5000L) {
  if (is.null(diameter) || length(diameter) < 1L || !is.finite(diameter[[1]])) {
    target_tiles <- 100L
    aspect <- diff(win$xrange) / diff(win$yrange)
    nx <- max(2L, round(sqrt(target_tiles * aspect)))
    ny <- max(2L, round(target_tiles / nx))
  } else {
    nx <- max(2L, ceiling(diff(win$xrange) / diameter))
    ny <- max(2L, ceiling(diff(win$yrange) / diameter))
    if ((nx * ny) > max_tiles) {
      shrink <- sqrt((nx * ny) / max_tiles)
      nx <- max(2L, ceiling(nx / shrink))
      ny <- max(2L, ceiling(ny / shrink))
    }
  }
  grid_part <- spatstat.geom::quadrats(win, nx = nx, ny = ny)
  n_tiles <- grid_part$n
  grid_procs <- character(n_tiles)
  tile_list <- spatstat.geom::tiles(grid_part)
  for (i in seq_len(n_tiles)) {
    tile_cent <- c(mean(tile_list[[i]]$xrange), mean(tile_list[[i]]$yrange))
    grid_procs[i] <- if (spatstat.geom::inside.owin(tile_cent[1], tile_cent[2], aoi_owin)) {
      "treated"
    } else {
      "control"
    }
  }
  names(grid_procs) <- spatstat.geom::tilenames(grid_part)
  list(partition = grid_part, processes = grid_procs, label = label)
}

oklahoma_report_tess_rects <- function(tess_obj, processes) {
  tile_list <- spatstat.geom::tiles(tess_obj)
  nms <- names(tile_list)
  do.call(rbind, lapply(seq_along(tile_list), function(i) {
    ow <- tile_list[[i]]
    nm <- if (!is.null(nms) && length(nms) >= i) nms[[i]] else as.character(i)
    data.frame(
      xmin = ow$xrange[1], xmax = ow$xrange[2],
      ymin = ow$yrange[1], ymax = ow$yrange[2],
      treatment = unname(processes[[nm]]),
      stringsAsFactors = FALSE
    )
  }))
}

oklahoma_report_plot_partition_pp <- function(base_map, pp_pre, pp_post) {
  p1 <- base_map +
    ggplot2::geom_point(
      data = pp_pre, ggplot2::aes(x = x, y = y, alpha = t),
      colour = "#2166ac", size = 0.8
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.15, 0.8), guide = "none") +
    ggplot2::labs(title = "Pre-treatment", x = "X (km)", y = "Y (km)")
  p2 <- base_map +
    ggplot2::geom_point(
      data = pp_post,
      ggplot2::aes(x = x, y = y, colour = Process, alpha = t),
      size = 0.8
    ) +
    ggplot2::scale_colour_manual(
      values = c(control = "#2166ac", treated = "#b2182b"),
      guide = "none"
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    ggplot2::labs(title = "Post-treatment", x = "X (km)", y = "Y (km)")
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

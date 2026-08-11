# Shared Oklahoma county / AOI geometry helpers for spatstat.
#
# Two pitfalls caused silent label corruption in earlier runs:
# 1) Flattening sf coordinates via st_coordinates() into a single ring drops
#    holes / multipart structure. Prefer spatstat's as.owin(<sf>) then scale
#    metres -> km.
# 2) spatstat.geom::tileindex(..., close.gaps = TRUE) (the default) can corrupt
#    already-assigned tile labels when any query point falls outside every tile
#    but inside the tessellation window. Always use close.gaps = FALSE here and
#    look up tile names (not positional integers).

oklahoma_sf_to_owin_km <- function(geom) {
  g <- sf::st_geometry(geom)
  if (length(g) != 1L) {
    stop("oklahoma_sf_to_owin_km() expects a single-feature geometry.")
  }
  ow_m <- spatstat.geom::as.owin(g)
  spatstat.geom::affine(ow_m, mat = diag(c(1 / 1000, 1 / 1000)))
}

oklahoma_sf_features_to_owins_km <- function(sf_df, name_col = "NAME") {
  if (!inherits(sf_df, "sf")) stop("sf_df must be an sf object.")
  n <- nrow(sf_df)
  if (n < 1L) stop("sf_df has zero rows.")
  out <- vector("list", n)
  for (i in seq_len(n)) {
    out[[i]] <- tryCatch(
      oklahoma_sf_to_owin_km(sf_df[i, ]),
      error = function(e) NULL
    )
  }
  names(out) <- as.character(sf_df[[name_col]])
  out
}

oklahoma_tile_names <- function(x, y, partition) {
  # close.gaps=FALSE is required; see file header.
  ti <- spatstat.geom::tileindex(x, y, partition, close.gaps = FALSE)
  as.character(ti)
}

oklahoma_assign_partition_process <- function(x, y, partition, partition_processes) {
  tile_nm <- oklahoma_tile_names(x, y, partition)
  if (is.null(names(partition_processes))) {
    stop("partition_processes must be a named character vector (tile name -> process).")
  }
  out <- rep(NA_character_, length(tile_nm))
  ok <- !is.na(tile_nm) & nzchar(tile_nm)
  if (any(ok)) {
    mapped <- unname(partition_processes[tile_nm[ok]])
    out[ok] <- mapped
  }
  out
}

oklahoma_assert_label_support <- function(df,
                                          control_ss,
                                          treated_ss,
                                          label_col = "location_process",
                                          context = "events") {
  if (is.null(df) || nrow(df) < 1L) return(invisible(0L))
  if (!label_col %in% names(df)) {
    stop("Missing label column '", label_col, "' while checking ", context, ".")
  }
  lab <- as.character(df[[label_col]])
  ic <- spatstat.geom::inside.owin(df$x, df$y, control_ss)
  it <- spatstat.geom::inside.owin(df$x, df$y, treated_ss)
  uniq <- xor(ic, it)
  support <- ifelse(ic & !it, "control", ifelse(it & !ic, "treated", NA_character_))
  n_unique <- sum(uniq)
  n_mismatch <- if (n_unique > 0L) {
    sum(uniq & !is.na(lab) & lab != support, na.rm = TRUE)
  } else {
    0L
  }
  if (n_mismatch > 0L) {
    stop(
      sprintf(
        "Geometry/label inconsistency for %s: %d / %d uniquely-supported points have labels opposite their control/treated background mask.",
        context, n_mismatch, n_unique
      )
    )
  }
  invisible(as.integer(n_mismatch))
}

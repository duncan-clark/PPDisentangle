#' @importFrom ggplot2 ggplot aes geom_point geom_sf scale_fill_manual
#'   scale_alpha_continuous scale_shape_manual labs theme_minimal theme
#'   element_text guide_legend
#' @importFrom sf st_sf st_sfc st_polygon st_union st_make_valid
NULL

#' Convert a spatstat tessellation to an sf POLYGON layer
#'
#' @param x A spatstat tess object
#' @return An sf data frame with tile_id and geometry columns
#' @export
tess_to_sf <- function(x) {
  tile_list <- tiles(x)

  owin_to_sfc <- function(w) {
    w_poly <- as.polygonal(w)
    rings <- w_poly$bdry
    closed_rings <- lapply(rings, function(r) {
      mat <- cbind(r$x, r$y)
      if (!isTRUE(all.equal(mat[1, ], mat[nrow(mat), ], check.attributes = FALSE))) {
        mat <- rbind(mat, mat[1, ])
      }
      return(list(mat))
    })
    polys <- lapply(closed_rings, sf::st_polygon)
    geom_col <- sf::st_sfc(polys)
    return(sf::st_union(geom_col))
  }

  sfc_list <- lapply(tile_list, owin_to_sfc)
  combined_sfc <- do.call(c, sfc_list)

  sf_df <- sf::st_sf(
    tile_id = names(tile_list),
    geometry = combined_sfc,
    crs = 2263
  )
  sf_df <- sf::st_make_valid(sf_df)
  return(sf_df)
}

#' Plot a superposed point process with true labels
#'
#' @param obs_data Data frame with columns x, y, t, process, background
#' @param partition A spatstat tess object
#' @param title Plot title
#' @return Prints a ggplot and returns invisibly
#' @export
plot_pp <- function(obs_data, partition, title = "Point Process Plot") {
  df <- obs_data
  df$Event_Type <- factor(df$background * 1, levels = c(0, 1), labels = c("Triggered", "Background"))
  df$Process <- factor(df$process, levels = c("control", "treated"))
  df$Time <- df$t

  full_pal <- c(control = "#377eb8", treated = "#e41a1c")
  levs <- levels(df$Process)[levels(df$Process) %in% df$Process]
  vals <- full_pal[levs]

  part_sf <- try(tess_to_sf(partition), silent = TRUE)
  if (inherits(part_sf, "try-error")) part_sf <- NA

  plot <- ggplot2::ggplot()
  if (!is.null(dim(part_sf))) {
    plot <- plot + ggplot2::geom_sf(data = part_sf, fill = NA, colour = "grey70")
  }
  plot <- plot +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Process, shape = .data$Event_Type, alpha = .data$Time),
      size = 2.5, stroke = 0.6, na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(name = "Process", values = vals,
      guide = ggplot2::guide_legend(order = 1, override.aes = list(fill = vals, colour = vals, alpha = 1))) +
    ggplot2::scale_alpha_continuous(name = "Time", range = c(0.3, 1)) +
    ggplot2::scale_shape_manual(name = "Event Type", values = c(Background = 21, Triggered = 23),
      guide = ggplot2::guide_legend(order = 3)) +
    ggplot2::labs(title = title, x = "x", y = "y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"))
  print(plot)
}

#' Plot a point process with inferred labels
#'
#' @param obs_data Data frame with columns x, y, t, process, inferred_process
#' @param partition A spatstat tess object
#' @param title Plot title
#' @return Prints a ggplot and returns invisibly
#' @export
plot_pp_inferred <- function(obs_data, partition, title = "Point Process Plot") {
  df <- obs_data
  df$Time <- df$t

  full_pal <- c(control = "#377eb8", treated = "#e41a1c")
  if (is.factor(df$process)) df$process <- as.character(df$process)
  levs <- unique(df$process)
  vals <- full_pal[levs]
  vals <- vals[which(names(vals) %in% df$inferred_process)]

  part_sf <- try(tess_to_sf(partition), silent = TRUE)
  if (inherits(part_sf, "try-error")) part_sf <- NA

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$inferred_process, shape = .data$process, alpha = .data$Time),
      size = 2.5, stroke = 0.6, na.rm = TRUE
    )
  if (!is.null(dim(part_sf))) {
    plot <- plot + ggplot2::geom_sf(data = part_sf, fill = NA, colour = "black", lwd = 0.8)
  }
  plot <- plot +
    ggplot2::scale_fill_manual(name = "inferred_process", values = vals,
      guide = ggplot2::guide_legend(order = 1, override.aes = list(fill = vals, colour = vals, alpha = 1))) +
    ggplot2::scale_alpha_continuous(name = "Time", range = c(0.3, 1)) +
    ggplot2::scale_shape_manual(name = "process", values = c(control = 21, treated = 23),
      guide = ggplot2::guide_legend(order = 3)) +
    ggplot2::labs(title = title, x = "x", y = "y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"))
  print(plot)
}

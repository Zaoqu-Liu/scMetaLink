#' @title scMetaLink Spatial Visualization Functions
#' @description Visualization functions for spatial metabolite communication analysis
#' @name spatial_visualization
#' @keywords internal
NULL

# =============================================================================
# Spatial Communication Plots
# =============================================================================

#' @title Plot Spatial Communication Heatmap
#' @description Visualize metabolite communication on spatial coordinates
#'
#' @param object A scMetaLink object with spatial information
#' @param metabolite Character. Metabolite ID or name to visualize
#' @param type Character. What to plot: "production", "sensing", or "communication"
#' @param cell_type Character. For communication, specify sender or receiver cell type
#' @param point_size Numeric. Size of spot points
#' @param alpha Numeric. Transparency (0-1)
#' @param low_color Character. Color for low values
#' @param high_color Character. Color for high values
#' @param title Character. Plot title (NULL for auto)
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#' \donttest{
#' data(st_expr)
#' data(st_meta)
#' data(st_scalefactors)
#' obj <- createScMetaLinkFromSpatial(st_expr, st_meta[,c("x","y")],
#'                                    st_meta, "cell_type", st_scalefactors)
#' obj <- inferProduction(obj)
#'
#' # Plot lactate production
#' plotSpatialFeature(obj, metabolite = "HMDB0000190", type = "production")
#' }
plotSpatialFeature <- function(object,
                               metabolite,
                               type = "production",
                               cell_type = NULL,
                               point_size = 1.5,
                               alpha = 0.8,
                               low_color = "#FFFFCC",
                               high_color = "#E31A1C",
                               title = NULL) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  coords <- object@parameters$spatial_coords
  cell_meta <- object@cell_meta

  # Get the appropriate scores
  if (type == "production") {
    if (is.null(object@production_scores)) {
      stop("Run inferProduction() first")
    }
    scores <- object@production_scores
    score_name <- "Production"
  } else if (type == "sensing") {
    if (is.null(object@sensing_scores)) {
      stop("Run inferSensing() first")
    }
    scores <- object@sensing_scores
    score_name <- "Sensing"
  } else {
    stop("type must be 'production' or 'sensing'")
  }

  # Find metabolite
  if (metabolite %in% rownames(scores)) {
    met_id <- metabolite
  } else {
    db <- object@database
    met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite, ignore.case = TRUE)]
    if (length(met_match) == 0 || !met_match[1] %in% rownames(scores)) {
      stop("Metabolite not found in scores")
    }
    met_id <- met_match[1]
  }

  # Get metabolite name
  met_name <- object@database$metabolites$metabolite[object@database$metabolites$hmdb == met_id][1]

  # Get cell type scores
  cell_type_col <- object@cell_type_column
  cell_types <- cell_meta[[cell_type_col]]
  met_scores <- scores[met_id, ]

  # Map cell type scores to spots
  spot_scores <- met_scores[cell_types]
  names(spot_scores) <- rownames(cell_meta)

  # Create plot data
  plot_data <- data.frame(
    x = coords$x,
    y = coords$y,
    score = spot_scores,
    cell_type = cell_types,
    stringsAsFactors = FALSE
  )

  # Auto title
  if (is.null(title)) {
    title <- sprintf("%s: %s", met_name, score_name)
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = score)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_gradient(low = low_color, high = high_color, name = score_name) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = title, x = "X", y = "Y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )

  p
}

#' @title Plot Spatial Cell Types
#' @description Visualize cell type distribution on spatial coordinates
#'
#' @param object A scMetaLink object with spatial information
#' @param point_size Numeric. Size of spot points
#' @param alpha Numeric. Transparency (0-1)
#' @param colors Character vector. Colors for cell types (NULL for default)
#' @param show_legend Logical. Show legend
#'
#' @return A ggplot2 object
#' @export
plotSpatialCellTypes <- function(object,
                                 point_size = 1.5,
                                 alpha = 0.8,
                                 colors = NULL,
                                 show_legend = TRUE) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  coords <- object@parameters$spatial_coords
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column

  # Create plot data
  plot_data <- data.frame(
    x = coords$x,
    y = coords$y,
    cell_type = cell_meta[[cell_type_col]],
    stringsAsFactors = FALSE
  )

  # Get colors
  n_types <- length(unique(plot_data$cell_type))
  if (is.null(colors)) {
    colors <- .get_colors(n_types)
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = cell_type)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::scale_color_manual(values = colors, name = "Cell Type") +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = "Spatial Cell Type Distribution", x = "X", y = "Y") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )

  p
}

#' @title Plot Spatial Communication Network
#' @description Visualize cell-cell communication with spatial context
#'
#' @param object A scMetaLink object with spatial communication results
#' @param metabolite Character. Metabolite to visualize (NULL for aggregate)
#' @param sender_type Character. Sender cell type (NULL for all)
#' @param receiver_type Character. Receiver cell type (NULL for all)
#' @param top_n Integer. Number of top interactions to show arrows for
#' @param arrow_scale Numeric. Scale factor for arrow thickness
#' @param point_size Numeric. Size of cell type center points
#' @param show_labels Logical. Show cell type labels
#'
#' @return A ggplot2 object
#' @export
plotSpatialCommunicationNetwork <- function(object,
                                            metabolite = NULL,
                                            sender_type = NULL,
                                            receiver_type = NULL,
                                            top_n = 20,
                                            arrow_scale = 1,
                                            point_size = 5,
                                            show_labels = TRUE) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  if (is.null(object@communication_scores)) {
    stop("Run computeSpatialCommunication() first")
  }

  coords <- object@parameters$spatial_coords
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  comm_scores <- object@communication_scores

  # Calculate cell type centers
  cell_types <- unique(cell_meta[[cell_type_col]])
  type_centers <- data.frame(
    cell_type = cell_types,
    x = sapply(cell_types, function(ct) {
      mean(coords$x[cell_meta[[cell_type_col]] == ct])
    }),
    y = sapply(cell_types, function(ct) {
      mean(coords$y[cell_meta[[cell_type_col]] == ct])
    }),
    stringsAsFactors = FALSE
  )

  # Get communication data
  if (is.null(metabolite)) {
    # Sum across metabolites
    comm_mat <- apply(comm_scores, c(1, 2), sum)
  } else {
    # Find metabolite
    if (metabolite %in% dimnames(comm_scores)[[3]]) {
      met_idx <- metabolite
    } else {
      db <- object@database
      met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite, ignore.case = TRUE)]
      if (length(met_match) == 0 || !met_match[1] %in% dimnames(comm_scores)[[3]]) {
        stop("Metabolite not found")
      }
      met_idx <- met_match[1]
    }
    comm_mat <- comm_scores[, , met_idx]
  }

  # Filter by cell types if specified
  if (!is.null(sender_type)) {
    if (!sender_type %in% rownames(comm_mat)) stop("Sender type not found")
    comm_mat <- comm_mat[sender_type, , drop = FALSE]
  }
  if (!is.null(receiver_type)) {
    if (!receiver_type %in% colnames(comm_mat)) stop("Receiver type not found")
    comm_mat <- comm_mat[, receiver_type, drop = FALSE]
  }

  # Create edge data
  edges <- data.frame(
    sender = rep(rownames(comm_mat), ncol(comm_mat)),
    receiver = rep(colnames(comm_mat), each = nrow(comm_mat)),
    strength = as.vector(comm_mat),
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$strength > 0, ]
  edges <- edges[order(-edges$strength), ]

  if (nrow(edges) > top_n) {
    edges <- edges[1:top_n, ]
  }

  # Add coordinates
  edges$x_start <- type_centers$x[match(edges$sender, type_centers$cell_type)]
  edges$y_start <- type_centers$y[match(edges$sender, type_centers$cell_type)]
  edges$x_end <- type_centers$x[match(edges$receiver, type_centers$cell_type)]
  edges$y_end <- type_centers$y[match(edges$receiver, type_centers$cell_type)]

  # Normalize strength for visualization
  edges$strength_norm <- edges$strength / max(edges$strength)

  # Get colors
  n_types <- nrow(type_centers)
  type_colors <- .get_colors(n_types)
  names(type_colors) <- type_centers$cell_type

  # Create plot
  p <- ggplot2::ggplot() +
    # Background spots (optional)
    ggplot2::geom_point(
      data = data.frame(x = coords$x, y = coords$y, ct = cell_meta[[cell_type_col]]),
      ggplot2::aes(x = x, y = y, color = ct),
      size = 0.5, alpha = 0.2
    ) +
    # Communication arrows
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
                   linewidth = strength_norm),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"), type = "closed"),
      alpha = 0.6, color = "gray30"
    ) +
    # Cell type centers
    ggplot2::geom_point(
      data = type_centers,
      ggplot2::aes(x = x, y = y, fill = cell_type),
      shape = 21, size = point_size, color = "black", stroke = 1
    ) +
    ggplot2::scale_color_manual(values = type_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = type_colors, name = "Cell Type") +
    ggplot2::scale_linewidth_continuous(range = c(0.5, 3) * arrow_scale, guide = "none") +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = "Spatial Communication Network",
      subtitle = if (!is.null(metabolite)) metabolite else "All metabolites (aggregated)",
      x = "X", y = "Y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray50"),
      panel.grid = ggplot2::element_blank()
    )

  # Add labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = type_centers,
      ggplot2::aes(x = x, y = y, label = cell_type),
      vjust = -1.5, size = 3, fontface = "bold"
    )
  }

  p
}

#' @title Plot Spatial Distance Distribution
#' @description Visualize the distribution of spatial distances between spots
#'
#' @param object A scMetaLink object with spatial information
#' @param by_celltype Logical. Show distances stratified by cell type pairs
#' @param max_distance Numeric. Maximum distance to show (NULL for all)
#' @param bins Integer. Number of histogram bins
#'
#' @return A ggplot2 object
#' @export
plotSpatialDistanceDistribution <- function(object,
                                            by_celltype = FALSE,
                                            max_distance = NULL,
                                            bins = 50) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  coords <- object@parameters$spatial_coords
  scale_factors <- object@parameters$scale_factors
  pixels_per_um <- scale_factors$pixels_per_um %||% 1

  # Convert to micrometers
  coords_um <- coords / pixels_per_um

  # Compute distances
  dist_mat <- as.matrix(dist(coords_um))

  if (!by_celltype) {
    # Simple histogram of all distances
    dist_vec <- dist_mat[lower.tri(dist_mat)]

    if (!is.null(max_distance)) {
      dist_vec <- dist_vec[dist_vec <= max_distance]
    }

    plot_data <- data.frame(distance = dist_vec)

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = distance)) +
      ggplot2::geom_histogram(bins = bins, fill = "#4DBBD5", color = "white", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = c(100, 200, 500), linetype = "dashed", color = "red", alpha = 0.5) +
      ggplot2::labs(
        title = "Spatial Distance Distribution",
        x = "Distance (um)",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

  } else {
    # Stratified by cell type pairs
    cell_meta <- object@cell_meta
    cell_type_col <- object@cell_type_column
    cell_types <- cell_meta[[cell_type_col]]
    unique_types <- unique(cell_types)

    # Sample for efficiency
    dist_samples <- list()

    for (ct1 in unique_types) {
      for (ct2 in unique_types) {
        idx1 <- which(cell_types == ct1)
        idx2 <- which(cell_types == ct2)

        if (length(idx1) > 0 && length(idx2) > 0) {
          dists <- as.vector(dist_mat[idx1, idx2])
          dists <- dists[dists > 0]  # Remove self-distances

          if (length(dists) > 1000) {
            dists <- sample(dists, 1000)
          }

          if (length(dists) > 0) {
            dist_samples[[paste(ct1, ct2, sep = " -> ")]] <- data.frame(
              pair = paste(ct1, ct2, sep = " -> "),
              distance = dists,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    plot_data <- do.call(rbind, dist_samples)

    if (!is.null(max_distance)) {
      plot_data <- plot_data[plot_data$distance <= max_distance, ]
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = distance, fill = pair)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::labs(
        title = "Spatial Distance by Cell Type Pairs",
        x = "Distance (um)",
        y = "Density",
        fill = "Cell Type Pair"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  }

  p
}

#' @title Plot Spatial Hotspots
#' @description Visualize communication hotspots on spatial coordinates
#'
#' @param object A scMetaLink object with spatial communication results
#' @param metabolite Character. Metabolite to analyze (NULL for aggregate)
#' @param type Character. "sender" or "receiver"
#' @param point_size Numeric. Size of spot points
#' @param hotspot_size Numeric. Size of hotspot markers
#' @param n_hotspots Integer. Number of hotspots to highlight
#'
#' @return A ggplot2 object
#' @export
plotSpatialHotspots <- function(object,
                                metabolite = NULL,
                                type = "sender",
                                point_size = 1,
                                hotspot_size = 5,
                                n_hotspots = 5) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  coords <- object@parameters$spatial_coords

  # Try to get hotspots
  tryCatch({
    hotspots <- identifyCommunicationHotspots(
      object, metabolite = metabolite, type = type, n_hotspots = n_hotspots
    )
  }, error = function(e) {
    stop("Could not identify hotspots. Make sure you've run computeSpatialCommunication() with analysis_level='spot'")
  })

  # Create plot
  p <- ggplot2::ggplot() +
    # Background spots
    ggplot2::geom_point(
      data = data.frame(x = coords$x, y = coords$y),
      ggplot2::aes(x = x, y = y),
      size = point_size, alpha = 0.3, color = "gray70"
    ) +
    # Hotspot centers
    ggplot2::geom_point(
      data = hotspots,
      ggplot2::aes(x = center_x, y = center_y, size = mean_score),
      shape = 21, fill = "#E64B35", color = "black", stroke = 1, alpha = 0.8
    ) +
    # Hotspot labels
    ggplot2::geom_text(
      data = hotspots,
      ggplot2::aes(x = center_x, y = center_y, label = hotspot_id),
      vjust = -1.5, size = 4, fontface = "bold"
    ) +
    ggplot2::scale_size_continuous(range = c(3, 10), name = "Communication\nScore") +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = sprintf("Communication Hotspots (%s)", type),
      subtitle = if (!is.null(metabolite)) metabolite else "Aggregated",
      x = "X", y = "Y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray50"),
      panel.grid = ggplot2::element_blank()
    )

  p
}

#' @title Plot Spatial Communication Comparison
#' @description Compare communication patterns for multiple metabolites
#'
#' @param object A scMetaLink object with spatial communication results
#' @param metabolites Character vector. Metabolites to compare
#' @param type Character. "production" or "sensing"
#' @param ncol Integer. Number of columns in facet grid
#' @param point_size Numeric. Size of spot points
#'
#' @return A ggplot2 object
#' @export
plotSpatialComparison <- function(object,
                                  metabolites,
                                  type = "production",
                                  ncol = 2,
                                  point_size = 1) {
  .check_dependencies("ggplot2")

  if (!.is_spatial(object)) {
    stop("Object does not contain spatial information")
  }

  if (type == "production") {
    scores <- object@production_scores
  } else {
    scores <- object@sensing_scores
  }

  if (is.null(scores)) {
    stop(sprintf("Run infer%s() first", tools::toTitleCase(type)))
  }

  coords <- object@parameters$spatial_coords
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  db <- object@database

  # Process metabolites
  plot_list <- list()

  for (met in metabolites) {
    # Find metabolite ID
    if (met %in% rownames(scores)) {
      met_id <- met
    } else {
      met_match <- db$metabolites$hmdb[grepl(met, db$metabolites$metabolite, ignore.case = TRUE)]
      if (length(met_match) == 0 || !met_match[1] %in% rownames(scores)) {
        warning(sprintf("Metabolite '%s' not found, skipping", met))
        next
      }
      met_id <- met_match[1]
    }

    met_name <- db$metabolites$metabolite[db$metabolites$hmdb == met_id][1]

    # Get scores per spot
    cell_types <- cell_meta[[cell_type_col]]
    met_scores <- scores[met_id, ]
    spot_scores <- met_scores[cell_types]

    plot_list[[met_name]] <- data.frame(
      x = coords$x,
      y = coords$y,
      score = spot_scores,
      metabolite = met_name,
      stringsAsFactors = FALSE
    )
  }

  if (length(plot_list) == 0) {
    stop("No valid metabolites found")
  }

  plot_data <- do.call(rbind, plot_list)

  # Create faceted plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = score)) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_gradient(low = "#FFFFCC", high = "#E31A1C", name = tools::toTitleCase(type)) +
    ggplot2::facet_wrap(~metabolite, ncol = ncol) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(
      title = sprintf("Spatial %s Comparison", tools::toTitleCase(type)),
      x = "X", y = "Y"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  p
}

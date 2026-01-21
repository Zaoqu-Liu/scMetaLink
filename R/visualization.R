#' @title Plot Communication Heatmap
#' @description Create a heatmap of cell-cell communication
#'
#' @param object scMetaLink object
#' @param metabolite Character. Specific metabolite to show (NULL for aggregated)
#' @param cluster_rows Logical. Cluster rows
#' @param cluster_cols Logical. Cluster columns
#' @param show_values Logical. Show values in cells
#' @param colors Character vector. Color palette
#' @param title Character. Plot title
#'
#' @return A ComplexHeatmap object
#' @export
#'
#' @examples
#' \donttest{
#' # Plot aggregated communication heatmap
#' plotCommunicationHeatmap(obj)
#'
#' # Plot specific metabolite
#' plotCommunicationHeatmap(obj, metabolite = "HMDB0000148")
#' }
plotCommunicationHeatmap <- function(object,
                                     metabolite = NULL,
                                     cluster_rows = TRUE,
                                     cluster_cols = TRUE,
                                     show_values = FALSE,
                                     colors = NULL,
                                     title = NULL) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required")
  }
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(metabolite)) {
    # Aggregated heatmap
    mat <- getCommunicationMatrix(object)
    if (is.null(title)) title <- "Aggregated Cell Communication"
  } else {
    # Specific metabolite
    if (is.null(object@communication_scores)) {
      stop("Communication scores not calculated")
    }

    # Find metabolite in communication array
    met_idx <- which(dimnames(object@communication_scores)[[3]] == metabolite)
    if (length(met_idx) == 0) {
      stop("Metabolite not found in communication scores")
    }

    mat <- object@communication_scores[, , met_idx]

    # Get metabolite name
    db <- object@database
    met_name <- db$metabolites$metabolite[db$metabolites$hmdb == metabolite]
    if (length(met_name) == 0) met_name <- metabolite
    if (is.null(title)) title <- paste("Communication via", met_name)
  }

  # Set colors
  if (is.null(colors)) {
    colors <- viridis::viridis(100)
  }

  col_fun <- circlize::colorRamp2(
    seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(colors)),
    colors
  )

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Score",
    col = col_fun,
    column_title = title,
    row_title = "Sender",
    column_title_side = "top",
    row_names_side = "left",
    column_names_side = "bottom",
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_dend = cluster_rows,
    show_column_dend = cluster_cols,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    column_names_rot = 45,
    cell_fun = if (show_values) {
      function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", mat[i, j]), x, y,
          gp = grid::gpar(fontsize = 8)
        )
      }
    } else {
      NULL
    },
    heatmap_legend_param = list(
      title = "Communication\nScore",
      title_position = "topcenter"
    )
  )

  ht
}

#' @title Plot Communication Circle
#' @description Create a chord diagram of cell-cell communication
#'
#' @param object scMetaLink object
#' @param top_n Integer. Number of top interactions to show
#' @param metabolite Character. Specific metabolite (NULL for aggregated)
#' @param colors Named vector. Colors for cell types
#' @param transparency Numeric. Link transparency (0-1)
#' @param title Character. Plot title
#'
#' @return Invisibly returns NULL, plots chord diagram
#' @export
plotCommunicationCircle <- function(object,
                                    top_n = 50,
                                    metabolite = NULL,
                                    colors = NULL,
                                    transparency = 0.5,
                                    title = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required")
  }
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (top_n < 1) {
    stop("top_n must be at least 1")
  }
  if (transparency < 0 || transparency > 1) {
    stop("transparency must be between 0 and 1")
  }

  if (is.null(metabolite)) {
    summary_df <- summarizeCommunicationPairs(object, "sum")
  } else {
    sig <- object@significant_interactions
    sig <- sig[sig$metabolite_id == metabolite, ]
    if (nrow(sig) == 0) {
      stop("No interactions found for this metabolite")
    }
    summary_df <- aggregate(communication_score ~ sender + receiver, data = sig, FUN = sum)
  }

  # Select top interactions
  summary_df <- summary_df[order(-summary_df$communication_score), ]
  summary_df <- head(summary_df, top_n)

  # Prepare for chord diagram
  cell_types <- unique(c(summary_df$sender, summary_df$receiver))

  # Set colors
  if (is.null(colors)) {
    n_types <- length(cell_types)
    if (n_types <= 12) {
      colors <- setNames(
        RColorBrewer::brewer.pal(max(3, n_types), "Set3")[1:n_types],
        cell_types
      )
    } else {
      colors <- setNames(
        c(
          RColorBrewer::brewer.pal(12, "Set3"),
          viridis::viridis(n_types - 12)
        ),
        cell_types
      )
    }
  }

  # Create chord diagram
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)

  circlize::chordDiagram(
    summary_df[, c("sender", "receiver", "communication_score")],
    grid.col = colors,
    transparency = transparency,
    annotationTrack = c("name", "grid"),
    annotationTrackHeight = c(0.03, 0.01),
    preAllocateTracks = list(
      track.height = 0.1,
      track.margin = c(0.01, 0)
    )
  )

  # Add title
  if (!is.null(title)) {
    graphics::title(main = title)
  }

  circlize::circos.clear()
  invisible(NULL)
}

#' @title Plot Communication Network
#' @description Create a network visualization of cell-cell communication
#'
#' @param object scMetaLink object
#' @param metabolite Character. Specific metabolite (NULL for aggregated)
#' @param min_score Numeric. Minimum score threshold
#' @param layout Character. Network layout algorithm
#' @param node_size_by Character. Size nodes by "degree" or "centrality"
#' @param edge_width_scale Numeric. Scale factor for edge widths
#' @param colors Named vector. Colors for cell types
#'
#' @return A ggplot object
#' @export
plotCommunicationNetwork <- function(object,
                                     metabolite = NULL,
                                     min_score = 0,
                                     layout = "fr",
                                     node_size_by = "degree",
                                     edge_width_scale = 2,
                                     colors = NULL) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Package 'ggraph' is required")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required")
  }
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (!node_size_by %in% c("degree", "centrality")) {
    stop("node_size_by must be 'degree' or 'centrality'")
  }

  # Get communication data
  if (is.null(metabolite)) {
    summary_df <- summarizeCommunicationPairs(object, "sum")
  } else {
    sig <- object@significant_interactions
    sig <- sig[sig$metabolite_id == metabolite, ]
    if (nrow(sig) == 0) {
      stop("No interactions found for this metabolite")
    }
    summary_df <- aggregate(communication_score ~ sender + receiver, data = sig, FUN = sum)
  }

  # Filter by minimum score
  summary_df <- summary_df[summary_df$communication_score >= min_score, ]

  if (nrow(summary_df) == 0) {
    stop("No interactions above minimum score threshold")
  }

  # Create igraph object
  g <- igraph::graph_from_data_frame(
    summary_df[, c("sender", "receiver", "communication_score")],
    directed = TRUE
  )

  # Calculate node metrics
  if (node_size_by == "degree") {
    node_size <- igraph::degree(g, mode = "all")
  } else {
    node_size <- igraph::betweenness(g)
  }

  # Set colors
  cell_types <- igraph::V(g)$name
  if (is.null(colors)) {
    n_types <- length(cell_types)
    colors <- setNames(
      scales::hue_pal()(n_types),
      cell_types
    )
  }

  # Create plot
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_arc(
      ggplot2::aes(
        width = communication_score,
        alpha = communication_score
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed"),
      end_cap = ggraph::circle(3, "mm"),
      strength = 0.2
    ) +
    ggraph::scale_edge_width_continuous(
      range = c(0.5, 3) * edge_width_scale,
      guide = "none"
    ) +
    ggraph::scale_edge_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
    ggraph::geom_node_point(
      ggplot2::aes(size = node_size, fill = name),
      shape = 21,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_size_continuous(range = c(5, 15)) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      size = 4
    ) +
    ggraph::theme_graph() +
    ggplot2::labs(
      fill = "Cell Type",
      size = if (node_size_by == "degree") "Degree" else "Centrality"
    ) +
    ggplot2::theme(legend.position = "right")

  p
}

#' @title Plot Metabolite Profile
#' @description Visualize production and sensing profiles for a metabolite
#'
#' @param object scMetaLink object
#' @param metabolite Character. Metabolite HMDB ID or name
#' @param show_genes Logical. Show contributing genes
#'
#' @return A ggplot object
#' @export
plotMetaboliteProfile <- function(object, metabolite, show_genes = FALSE) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  # Find metabolite
  db <- object@database
  if (!metabolite %in% db$metabolites$hmdb) {
    met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite,
      ignore.case = TRUE
    )]
    if (length(met_match) == 0) stop("Metabolite not found")
    metabolite <- met_match[1]
  }

  met_name <- db$metabolites$metabolite[db$metabolites$hmdb == metabolite]

  # Get scores
  prod_scores <- NULL
  sens_scores <- NULL

  if (!is.null(object@production_scores) && metabolite %in% rownames(object@production_scores)) {
    prod_scores <- object@production_scores[metabolite, ]
  }
  if (!is.null(object@sensing_scores) && metabolite %in% rownames(object@sensing_scores)) {
    sens_scores <- object@sensing_scores[metabolite, ]
  }

  if (is.null(prod_scores) && is.null(sens_scores)) {
    stop("No scores available for this metabolite")
  }

  # Prepare data
  plot_data <- data.frame()

  if (!is.null(prod_scores)) {
    plot_data <- rbind(plot_data, data.frame(
      cell_type = names(prod_scores),
      score = as.numeric(prod_scores),
      type = "Production",
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(sens_scores)) {
    plot_data <- rbind(plot_data, data.frame(
      cell_type = names(sens_scores),
      score = as.numeric(sens_scores),
      type = "Sensing",
      stringsAsFactors = FALSE
    ))
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cell_type, y = score, fill = type)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::scale_fill_manual(values = c("Production" = "#E64B35", "Sensing" = "#4DBBD5")) +
    ggplot2::labs(
      title = met_name,
      x = "Cell Type",
      y = "Score",
      fill = "Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  p
}

#' @title Plot Top Interactions
#' @description Dot plot of top metabolite-mediated interactions
#'
#' @param object scMetaLink object
#' @param top_n Integer. Number of top interactions
#' @param group_by Character. Group by "sender", "receiver", or "metabolite"
#'
#' @return A ggplot object
#' @export
plotTopInteractions <- function(object, top_n = 30, group_by = "metabolite") {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (!group_by %in% c("sender", "receiver", "metabolite")) {
    stop("group_by must be 'sender', 'receiver', or 'metabolite'")
  }

  sig <- object@significant_interactions

  if (nrow(sig) == 0) {
    stop("No significant interactions")
  }

  # Get top interactions
  sig <- sig[order(-sig$communication_score), ]
  sig <- head(sig, top_n)

  # Create interaction label
  sig$interaction <- paste(sig$sender, "->", sig$receiver)

  # Set group and color mapping
  if (group_by == "metabolite") {
    sig$group <- sig$metabolite_name
    x_var <- "interaction"
    y_var <- "group"
  } else if (group_by == "sender") {
    sig$group <- sig$sender
    x_var <- "metabolite_name"
    y_var <- "interaction"
  } else {
    sig$group <- sig$receiver
    x_var <- "metabolite_name"
    y_var <- "interaction"
  }

  p <- ggplot2::ggplot(sig, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]])) +
    ggplot2::geom_point(ggplot2::aes(
      size = communication_score,
      color = -log10(pvalue_adjusted + 1e-10)
    )) +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(
      x = if (group_by == "metabolite") "Cell Type Pair" else "Metabolite",
      y = if (group_by == "metabolite") "Metabolite" else "Interaction",
      size = "Communication\nScore",
      color = "-log10(p-value)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}

#' @title Plot Pathway Communication
#' @description Visualize pathway-level communication
#'
#' @param object scMetaLink object
#' @param top_pathways Integer. Number of top pathways to show
#' @param type Character. "heatmap" or "bar"
#'
#' @return A ggplot or ComplexHeatmap object
#' @export
plotPathwayCommunication <- function(object, top_pathways = 20, type = "bar") {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (!type %in% c("bar", "heatmap")) {
    stop("type must be 'bar' or 'heatmap'")
  }
  if (nrow(object@pathway_aggregated) == 0) {
    stop("Pathway aggregation not done. Run aggregateByPathway() first.")
  }

  pw_summary <- summarizePathwayActivity(object)
  pw_summary <- head(pw_summary, top_pathways)

  if (type == "bar") {
    p <- ggplot2::ggplot(
      pw_summary,
      ggplot2::aes(
        x = stats::reorder(pathway, communication_score),
        y = communication_score,
        fill = n_metabolites
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Total Communication Score",
        fill = "# Metabolites"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 9)
      )

    return(p)
  }

  if (type == "heatmap") {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("Package 'ComplexHeatmap' is required for heatmap")
    }

    # Create pathway x cell type matrix
    pw_agg <- object@pathway_aggregated
    pw_agg <- pw_agg[pw_agg$pathway %in% pw_summary$pathway, ]

    # Sender activity
    sender_mat <- tapply(
      pw_agg$communication_score,
      list(pw_agg$pathway, pw_agg$sender),
      sum
    )
    sender_mat[is.na(sender_mat)] <- 0

    return(ComplexHeatmap::Heatmap(
      sender_mat,
      name = "Score",
      column_title = "Sender Cell Type",
      row_title = "Pathway",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      col = viridis::viridis(100)
    ))
  }
}

#' @title Plot Differential Communication
#' @description Visualize differential communication between conditions
#'
#' @param diff_results data.frame from compareCommunication()
#' @param top_n Integer. Number of top changes to show
#' @param type Character. "bar" or "volcano"
#'
#' @return A ggplot object
#' @export
plotDifferentialCommunication <- function(diff_results, top_n = 30, type = "bar") {
  if (!is.data.frame(diff_results)) {
    stop("diff_results must be a data.frame from compareCommunication()")
  }
  if (!"change" %in% names(diff_results)) {
    stop("diff_results must contain a 'change' column")
  }
  if (!type %in% c("bar", "volcano")) {
    stop("type must be 'bar' or 'volcano'")
  }

  if (type == "bar") {
    # Get top positive and negative changes
    diff_sorted <- diff_results[order(diff_results$change), ]
    top_neg <- head(diff_sorted, ceiling(top_n / 2))
    top_pos <- tail(diff_sorted, floor(top_n / 2))
    plot_data <- rbind(top_neg, top_pos)

    plot_data$label <- paste(
      plot_data$sender, "->", plot_data$receiver,
      "\n(", plot_data$metabolite, ")"
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = stats::reorder(label, change),
        y = change,
        fill = change > 0
      )
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
        guide = "none"
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Log2 Fold Change"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

    return(p)
  }

  if (type == "volcano") {
    # Need score columns for volcano plot
    score_cols <- grep("^score_", names(diff_results), value = TRUE)
    if (length(score_cols) < 2) {
      stop("Volcano plot requires score columns from compareCommunication()")
    }

    # Calculate mean score for point size
    diff_results$mean_score <- (diff_results[[score_cols[1]]] + diff_results[[score_cols[2]]]) / 2

    # Significance threshold
    diff_results$significant <- abs(diff_results$change) > 1

    p <- ggplot2::ggplot(
      diff_results,
      ggplot2::aes(x = change, y = mean_score, color = significant)
    ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#E64B35", "FALSE" = "grey60"),
        guide = "none"
      ) +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
      ggplot2::labs(
        x = "Log2 Fold Change",
        y = "Mean Communication Score"
      ) +
      ggplot2::theme_minimal()

    return(p)
  }
}

#' @title Plot Enriched Pathways
#' @description Visualize enriched pathways from enrichPathways()
#'
#' @param enrichment_results data.frame from enrichPathways()
#' @param top_n Integer. Number of top pathways to show
#' @param show_overlap Logical. Show overlap count
#'
#' @return A ggplot object
#' @export
plotEnrichedPathways <- function(enrichment_results, top_n = 20, show_overlap = TRUE) {
  if (!is.data.frame(enrichment_results)) {
    stop("enrichment_results must be a data.frame from enrichPathways()")
  }
  if (nrow(enrichment_results) == 0) {
    stop("No enriched pathways to plot")
  }

  plot_data <- head(enrichment_results, top_n)

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = stats::reorder(pathway, fold_enrichment),
      y = fold_enrichment
    )
  ) +
    ggplot2::geom_segment(ggplot2::aes(xend = pathway, y = 0, yend = fold_enrichment),
      color = "grey60"
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n_overlap, color = -log10(pvalue_adjusted)),
      alpha = 0.8
    ) +
    ggplot2::scale_color_viridis_c(option = "plasma") +
    ggplot2::scale_size_continuous(range = c(3, 10)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Fold Enrichment",
      color = "-log10(p-adj)",
      size = "Overlap"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9)
    )

  p
}

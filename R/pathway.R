#' @title Aggregate Communication by Pathway
#' @description Aggregate metabolite-mediated communication at the pathway level.
#'   Due to the large number of pathway associations, this function uses a 
#'   simplified approach focusing on top pathways.
#'
#' @param object A scMetaLink object with significant interactions
#' @param top_pathways Integer. Number of top pathways to analyze (default: 50)
#' @param min_metabolites Integer. Minimum metabolites per pathway (default: 3)
#'
#' @return Updated scMetaLink object with pathway_aggregated slot filled
#' @export
aggregateByPathway <- function(object,
                               top_pathways = 50,
                               min_metabolites = 3) {

  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@significant_interactions) == 0) {
    stop("No significant interactions. Run filterSignificantInteractions() first.")
  }
  if (top_pathways < 1) {
    stop("top_pathways must be at least 1")
  }
  if (min_metabolites < 1) {
    stop("min_metabolites must be at least 1")
  }

  db <- object@database
  sig <- object@significant_interactions

  # Get metabolites from significant interactions
  met_ids <- unique(sig$metabolite_id)
  
  # Get pathway info for these metabolites
  pathway_info <- db$pathway[db$pathway$hmdb %in% met_ids, ]
  
  if (nrow(pathway_info) == 0) {
    message("No pathway information found")
    object@pathway_aggregated <- data.frame()
    return(object)
  }

  # Count metabolites per pathway
  pw_met_counts <- table(pathway_info$pathway)
  
  # Filter pathways with minimum metabolites
  valid_pathways <- names(pw_met_counts)[pw_met_counts >= min_metabolites]
  
  if (length(valid_pathways) == 0) {
    message("No pathways with enough metabolites")
    object@pathway_aggregated <- data.frame()
    return(object)
  }

  # Calculate pathway activity score
  # For each pathway: sum of communication scores of metabolites in that pathway
  pathway_scores <- sapply(valid_pathways, function(pw) {
    pw_mets <- pathway_info$hmdb[pathway_info$pathway == pw]
    pw_sig <- sig[sig$metabolite_id %in% pw_mets, ]
    sum(pw_sig$communication_score)
  })
  
  # Sort and get top pathways
  pathway_scores <- sort(pathway_scores, decreasing = TRUE)
  top_pw_names <- names(pathway_scores)[1:min(top_pathways, length(pathway_scores))]
  
  # Build summary for top pathways
  pathway_agg <- lapply(top_pw_names, function(pw) {
    pw_mets <- pathway_info$hmdb[pathway_info$pathway == pw]
    pw_sig <- sig[sig$metabolite_id %in% pw_mets, ]
    
    # Summarize by cell type pairs
    if (nrow(pw_sig) > 0) {
      pair_summary <- aggregate(
        communication_score ~ sender + receiver,
        data = pw_sig,
        FUN = sum
      )
      pair_summary$pathway <- pw
      pair_summary$n_metabolites <- length(unique(pw_sig$metabolite_id))
      pair_summary
    }
  })
  
  pathway_agg <- do.call(rbind, pathway_agg)
  
  if (is.null(pathway_agg) || nrow(pathway_agg) == 0) {
    message("No pathway associations found")
    object@pathway_aggregated <- data.frame()
    return(object)
  }

  pathway_agg <- pathway_agg[order(-pathway_agg$communication_score), ]
  rownames(pathway_agg) <- NULL
  
  # Reorder columns
  pathway_agg <- pathway_agg[, c("pathway", "sender", "receiver", "communication_score", "n_metabolites")]

  object@pathway_aggregated <- pathway_agg
  object@parameters$pathway <- list(
    top_pathways = top_pathways,
    min_metabolites = min_metabolites
  )

  message(sprintf("Aggregated %d pathway-cell type pair interactions", nrow(pathway_agg)))
  object
}

#' Get Pathway Communication Matrix
#' @param object scMetaLink object
#' @param pathway Character. Pathway name
#' @return Matrix of pathway-specific communication
#' @export
getPathwayCommunicationMatrix <- function(object, pathway) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@pathway_aggregated) == 0) {
    stop("Pathway aggregation not done. Run aggregateByPathway() first.")
  }

  pw_data <- object@pathway_aggregated[object@pathway_aggregated$pathway == pathway, ]

  if (nrow(pw_data) == 0) {
    stop("Pathway not found in aggregated results")
  }

  cell_types <- unique(c(pw_data$sender, pw_data$receiver))
  mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
                dimnames = list(cell_types, cell_types))

  for (i in seq_len(nrow(pw_data))) {
    mat[pw_data$sender[i], pw_data$receiver[i]] <- pw_data$communication_score[i]
  }

  mat
}

#' Summarize Pathway Activity
#' @param object scMetaLink object
#' @return data.frame with pathway activity summary
#' @export
summarizePathwayActivity <- function(object) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@pathway_aggregated) == 0) {
    stop("Pathway aggregation not done. Run aggregateByPathway() first.")
  }

  pw_agg <- object@pathway_aggregated

  # Aggregate by pathway
  result <- aggregate(
    communication_score ~ pathway,
    data = pw_agg, FUN = sum
  )
  
  # Add n_metabolites (use max since it's the same for each pathway)
  n_mets <- aggregate(n_metabolites ~ pathway, data = pw_agg, FUN = max)
  result <- merge(result, n_mets, by = "pathway")
  
  # Add n_cell_pairs
  result$n_cell_pairs <- as.numeric(table(pw_agg$pathway)[result$pathway])
  result <- result[order(-result$communication_score), ]
  rownames(result) <- NULL
  
  result
}

#' Get Metabolites in Pathway
#' @param pathway Character. Pathway name (can be partial match)
#' @param only_signaling Logical. Only return metabolites with receptors
#' @return data.frame with metabolites in the pathway
#' @export
getPathwayMetabolites <- function(pathway, only_signaling = FALSE) {
  db <- .load_metalinksdb()

  matching_pathways <- unique(db$pathway$pathway[grepl(pathway, db$pathway$pathway, ignore.case = TRUE)])

  if (length(matching_pathways) == 0) {
    stop("No pathways found matching the query")
  }

  if (length(matching_pathways) > 1) {
    message(sprintf("Found %d matching pathways, using first match: %s", 
                    length(matching_pathways), matching_pathways[1]))
  }

  pw <- matching_pathways[1]
  pw_mets <- db$pathway$hmdb[db$pathway$pathway == pw]
  mets <- db$metabolites[db$metabolites$hmdb %in% pw_mets, ]

  if (only_signaling) {
    lr_mets <- unique(db$edges$hmdb[db$edges$type == "lr"])
    mets <- mets[mets$hmdb %in% lr_mets, ]
  }

  mets
}

#' List Top Pathways
#' @param object scMetaLink object  
#' @param n Integer. Number of top pathways
#' @return data.frame with top pathways
#' @export
listTopPathways <- function(object, n = 20) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@pathway_aggregated) == 0) {
    stop("Pathway aggregation not done. Run aggregateByPathway() first.")
  }
  if (n < 1) {
    stop("n must be at least 1")
  }
  
  summary <- summarizePathwayActivity(object)
  head(summary, n)
}

#' @title Pathway Enrichment Analysis
#' @description Perform hypergeometric test to identify enriched pathways
#'   among significant metabolite-mediated interactions
#'
#' @param object scMetaLink object with significant interactions
#' @param pvalue_threshold Numeric. P-value cutoff for enrichment (default: 0.05)
#' @param min_overlap Integer. Minimum overlap between significant metabolites and pathway (default: 2)
#' @param adjust_method Character. P-value adjustment method (default: "BH")
#'
#' @return data.frame with enriched pathways and statistics
#' @export
#'
#' @examples
#' \donttest{
#' # Perform pathway enrichment
#' enriched <- enrichPathways(result)
#' head(enriched)
#' }
enrichPathways <- function(object,
                           pvalue_threshold = 0.05,
                           min_overlap = 2,
                           adjust_method = "BH") {
  
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@significant_interactions) == 0) {
    stop("No significant interactions. Run filterSignificantInteractions() first.")
  }
  if (pvalue_threshold <= 0 || pvalue_threshold > 1) {
    stop("pvalue_threshold must be between 0 and 1")
  }
  if (min_overlap < 1) {
    stop("min_overlap must be at least 1")
  }
  
  db <- object@database
  sig <- object@significant_interactions
  
  # Get significant metabolites
  sig_mets <- unique(sig$metabolite_id)
  n_sig <- length(sig_mets)
  
  # Get all metabolites in the database (background)
  all_mets <- unique(db$metabolites$hmdb)
  n_total <- length(all_mets)
  
  # Get pathway-metabolite associations
  pathway_info <- db$pathway
  
  # Get all unique pathways
  all_pathways <- unique(pathway_info$pathway)
  
  # Calculate enrichment for each pathway using hypergeometric test
  results <- lapply(all_pathways, function(pw) {
    # Metabolites in this pathway
    pw_mets <- pathway_info$hmdb[pathway_info$pathway == pw]
    pw_mets <- intersect(pw_mets, all_mets)  # Ensure they're in background
    n_pathway <- length(pw_mets)
    
    # Overlap with significant metabolites
    overlap <- intersect(sig_mets, pw_mets)
    n_overlap <- length(overlap)
    
    if (n_overlap < min_overlap) {
      return(NULL)
    }
    
    # Hypergeometric test (one-tailed, enrichment)
    # phyper(q, m, n, k, lower.tail = FALSE)
    # q = overlap - 1 (for strictly greater)
    # m = pathway size
    # n = non-pathway size
    # k = number of significant
    pvalue <- phyper(
      n_overlap - 1,
      n_pathway,
      n_total - n_pathway,
      n_sig,
      lower.tail = FALSE
    )
    
    # Calculate fold enrichment
    expected <- n_sig * n_pathway / n_total
    fold_enrichment <- n_overlap / expected
    
    data.frame(
      pathway = pw,
      n_pathway_mets = n_pathway,
      n_sig_mets = n_sig,
      n_overlap = n_overlap,
      overlap_mets = paste(overlap, collapse = ";"),
      expected = round(expected, 2),
      fold_enrichment = round(fold_enrichment, 2),
      pvalue = pvalue,
      stringsAsFactors = FALSE
    )
  })
  
  results <- do.call(rbind, results)
  
  if (is.null(results) || nrow(results) == 0) {
    message("No pathways with sufficient overlap found")
    return(data.frame())
  }
  
  # Adjust p-values
  results$pvalue_adjusted <- p.adjust(results$pvalue, method = adjust_method)
  
  # Filter by adjusted p-value
  results <- results[results$pvalue_adjusted < pvalue_threshold, ]
  
  if (nrow(results) == 0) {
    message("No significantly enriched pathways found")
    return(data.frame())
  }
  
  # Sort by fold enrichment
  results <- results[order(-results$fold_enrichment), ]
  rownames(results) <- NULL
  
  message(sprintf("Found %d enriched pathways (p < %g)", nrow(results), pvalue_threshold))
  
  results
}

#' @title Get Pathway-Metabolite Network
#' @description Extract pathway-metabolite relationships for network visualization
#'
#' @param object scMetaLink object
#' @param pathways Character vector. Specific pathways to include (NULL for top pathways)
#' @param top_n Integer. Number of top pathways if pathways is NULL
#'
#' @return data.frame with pathway-metabolite edges
#' @export
getPathwayMetaboliteNetwork <- function(object, pathways = NULL, top_n = 10) {
  
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (nrow(object@significant_interactions) == 0) {
    stop("No significant interactions")
  }
  
  db <- object@database
  sig <- object@significant_interactions
  sig_mets <- unique(sig$metabolite_id)
  
  # Determine which pathways to include
  if (is.null(pathways)) {
    if (nrow(object@pathway_aggregated) > 0) {
      pw_summary <- summarizePathwayActivity(object)
      pathways <- head(pw_summary$pathway, top_n)
    } else {
      # Use pathways of significant metabolites
      pathway_info <- db$pathway[db$pathway$hmdb %in% sig_mets, ]
      pw_counts <- table(pathway_info$pathway)
      pathways <- names(sort(pw_counts, decreasing = TRUE))[1:min(top_n, length(pw_counts))]
    }
  }
  
  # Get pathway-metabolite relationships
  pathway_info <- db$pathway[db$pathway$pathway %in% pathways & 
                             db$pathway$hmdb %in% sig_mets, ]
  
  if (nrow(pathway_info) == 0) {
    message("No pathway-metabolite relationships found")
    return(data.frame())
  }
  
  # Add metabolite names
  pathway_info <- merge(pathway_info, db$metabolites[, c("hmdb", "metabolite")],
                        by = "hmdb", all.x = TRUE)
  
  # Add communication score for the metabolite
  met_scores <- aggregate(communication_score ~ metabolite_id, data = sig, FUN = sum)
  pathway_info <- merge(pathway_info, met_scores, 
                        by.x = "hmdb", by.y = "metabolite_id", all.x = TRUE)
  pathway_info$communication_score[is.na(pathway_info$communication_score)] <- 0
  
  pathway_info[, c("pathway", "hmdb", "metabolite", "communication_score")]
}

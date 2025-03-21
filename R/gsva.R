#' @title Custom GSVA Function with AUCell Integration
#' @description Extends GSVA with AUCell for single-sample gene set enrichment analysis.
#' @param expr Gene expression matrix (genes × samples).
#' @param gset.idx.list List of gene sets (e.g., pathways).
#' @param method Enrichment method: "gsva", "ssgsea", "zscore", "plage", or "aucell".
#' @param kcdf Kernel function for CDF estimation (default: "Gaussian").
#' @param aucMaxRank Top-ranked genes to consider for AUCell (default: top 5%).
#' @param ... Additional arguments passed to GSVA or AUCell.
#' @return A matrix of enrichment scores (gene sets × samples).
#' @importFrom GSVA gsvaParam ssgseaParam zscoreParam plageParam gsva
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC
#' @importFrom SummarizedExperiment assay
#' @export
gsva <- function(expr, gset.idx.list, method = c("gsva", "ssgsea", "zscore", "plage", "aucell"),
                 kcdf = c("Gaussian", "Poisson", "none"), aucMaxRank = NULL, ...) {
  
  # Match method and kcdf arguments
  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  
  # Handle AUCell separately
  if (method == "aucell") {
    if (!requireNamespace("AUCell", quietly = TRUE))
      stop("Install the 'AUCell' package to use method='aucell'.")
    
    # Build gene rankings
    rankings <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
    
    # Calculate AUC scores
    auc_scores <- AUCell::AUCell_calcAUC(
      geneSets = gset.idx.list,
      rankings = rankings,
      aucMaxRank = if (is.null(aucMaxRank)) ceiling(0.05 * nrow(expr)) else aucMaxRank,
      ...
    )
    
    # Extract AUC scores
    es <- SummarizedExperiment::assay(auc_scores)
  } else {
    # GSVA 2.0.6+ requires method-specific parameter objects
    param <- switch(
      method,
      "gsva"    = GSVA::gsvaParam(expr, gset.idx.list, kcdf = kcdf, ...),
      "ssgsea"  = GSVA::ssgseaParam(expr, gset.idx.list, ...),
      "zscore"  = GSVA::zscoreParam(expr, gset.idx.list, ...),  # Remove kcdf for zscore
      "plage"   = GSVA::plageParam(expr, gset.idx.list, ...)    # Remove kcdf for plage
    )
    
    # Run GSVA
    es <- GSVA::gsva(param)
  }
  
  return(es)
}

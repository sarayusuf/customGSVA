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
#' @examples
#' # Example expression matrix
#' expr <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' rownames(expr) <- paste0("Gene", 1:100)
#' 
#' # Example gene set list
#' gset.idx.list <- list(
#'   Pathway1 = c("Gene1", "Gene2", "Gene3"),
#'   Pathway2 = c("Gene4", "Gene5", "Gene6"),
#'   Pathway3 = c("Gene7", "Gene8", "Gene9")
#' )
#' 
#' # Run GSVA
#' es_gsva <- gsva(expr, gset.idx.list, method = "gsva")
#' 
#' # Run AUCell
#' es_aucell <- gsva(expr, gset.idx.list, method = "aucell")
gsva <- function(expr, gset.idx.list, method = c("gsva", "ssgsea", "zscore", "plage", "aucell"),
                 kcdf = c("Gaussian", "Poisson", "none"), aucMaxRank = NULL, ...) {
  
  # Validate inputs
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("expr must be a numeric matrix.")
  }
  
  if (!is.list(gset.idx.list) || is.null(names(gset.idx.list))) {
    stop("gset.idx.list must be a named list.")
  }
  
  # Match method and kcdf arguments
  method <- match.arg(method)
  kcdf <- match.arg(kcdf)
  
  # Handle AUCell separately
  if (method == "aucell") {
    if (!requireNamespace("AUCell", quietly = TRUE)) {
      stop(
        "The 'AUCell' package is required for method='aucell'. ",
        "Install it using: BiocManager::install('AUCell')"
      )
    }
    
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
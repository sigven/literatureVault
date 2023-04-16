#' Vault of literature references for precision oncology
#'
#' A small dataset with identifiers and links to literature references
#' of relevance for precision cancer medicine (biomarkers etc.)
#'
#' @format \bold{litvault} - A data frame with 8,371 rows and 3 columns:
#' \itemize{
#'   \item{source_id} - {source identifier (e.g. PMID)}
#'   \item{source} - {type of source (PubMed, clinicaltrials.gov, FDA etc)}
#'   \item{source_entity} - {gwas or biomarker}
#'   \item{name} - {name of reference}
#'   \item{link} - {HTML link of reference}
#' }
#'
"vault"

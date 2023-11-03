chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

#' A function that returns a citation with first author, journal and year for a PubMed ID
#'
#' @param pmid An array of Pubmed IDs
#' @param chunk_size Size of PMID chunks
#'
#' @export
#'
#' @return citation PubMed citation, with first author, journal and year
#'
get_citations_pubmed <- function(
    pmid,
    chunk_size = 100){


  ## set logging layout
  lgr::lgr$appenders$console$set_layout(
    lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T"))

  ## make chunk of maximal 400 PMIDs from input array (limit by EUtils)
  pmid_chunks <- chunk(
    pmid, ceiling(length(pmid)/chunk_size))
  j <- 0
  all_citations <- data.frame()
  lgr::lgr$info(
    paste0('Retrieving PubMed citations for PMID list, total length: ',
           length(pmid)))
  while (j < length(pmid_chunks)) {
    pmid_chunk <- pmid_chunks[[as.character(j)]]
    if(!is.null(pmid_chunk)){
      lgr::lgr$info(
        paste0('Processing chunk ',j,' with ',length(pmid_chunk),' PMIDS'))
      pmid_string <- paste(pmid_chunk,collapse = " ")
      res <- RISmed::EUtilsGet(
        RISmed::EUtilsSummary(
          pmid_string, type = "esearch", db = "pubmed", retmax = 5000)
      )


      year <- RISmed::YearPubmed(res)
      authorlist <- RISmed::Author(res)
      pmid_list <- RISmed::PMID(res)
      keywords <- RISmed::Keywords(res)
      mesh <- RISmed::Mesh(res)
      mesh_descriptors <- c()
      keyword_elements <- c()


      k <- 1
      while(k <= length(mesh)){
        if(identical(typeof(mesh[[k]]), "logical")){
          mesh_descriptors <- c(mesh_descriptors,NA)
        }else{
          descriptors <- mesh[[k]] |> dplyr::filter(Type == "Descriptor")
          descriptors_all <- paste(descriptors$Heading, collapse=";")
          mesh_descriptors <- c(mesh_descriptors,descriptors_all)
        }
        k <- k + 1
      }

      m <- 1
      while(m <= length(keywords)){
        if(identical(typeof(keywords[[m]]), "logical")){
          keyword_elements <- c(keyword_elements, NA)
        }else{
          keyword_elements <- c(keyword_elements,
                                paste(unlist(keywords[[m]]), collapse=";"))
        }
        m <- m + 1
      }

      i <- 1
      first_author <- c()
      while (i <= length(authorlist)) {

        if (length(authorlist[[i]]) == 5) {
          first_author <- c(
            first_author,
            paste(authorlist[[i]][1,]$LastName," et al.",sep = ""))
        } else{
          first_author <- c(
            first_author, as.character("Unknown et al.")
          )
        }
        i <- i + 1
      }
      journal <- RISmed::ISOAbbreviation(res)
      citations <- data.frame(
        'pmid' = as.integer(pmid_list),
        'keywords' = keyword_elements,
        'mesh_descriptors' = mesh_descriptors,
        'citation' = paste(
          first_author, year, journal, sep = ", "),
        stringsAsFactors = F)
      citations$link <- paste0(
        '<a href=\'https://www.ncbi.nlm.nih.gov/pubmed/',
        citations$pmid,'\' target=\'_blank\'>',
        citations$citation,'</a>')
      all_citations <- dplyr::bind_rows(
        all_citations, citations)
    }
    j <- j + 1
  }

  return(all_citations)

}

get_literature <- function(
    literature_df, chunk_size = 100) {

  stopifnot(!is.na(literature_df))
  stopifnot(is.data.frame(literature_df))
  assertable::assert_colnames(
    literature_df, c("source","source_id", "source_entity"),
    quiet = T)

  pubmed_source_ids <- as.data.frame(literature_df |>
    dplyr::filter(
      source == "PubMed") |>
    dplyr::mutate(
      source_id = as.character(source_id)) |>
    dplyr::group_by(
      source_id, source
    ) |>
    dplyr::summarise(
      source_entity = paste(source_entity, collapse=";"),
      .groups = "drop"
    ) |>
    dplyr::arrange(source_id) |>
    dplyr::distinct())

  other_source_ids <- literature_df |>
    dplyr::filter(source != "PubMed")


  if (NROW(pubmed_source_ids) > 0) {
    pubmed_source_ids <- get_citations_pubmed(
      pmid = unique(pubmed_source_ids$source_id), chunk_size = chunk_size) |>
      dplyr::rename(
        name = citation,
      ) |>
      dplyr::mutate(source_id = as.character(pmid),
                    source_keywords = keywords,
                    source_mesh_descriptors = mesh_descriptors) |>
      dplyr::mutate(
        source_keywords = stringr::str_replace_all(
          source_keywords, "\\t", " "
        )
      ) |>
      dplyr::mutate(
        source_keywords = stringr::str_replace_all(
          source_keywords, "\\n", ""
        )
      ) |>
      dplyr::mutate(
        source_mesh_descriptors = stringr::str_replace_all(
          source_mesh_descriptors, "\\t", " "
        )
      ) |>
      dplyr::mutate(
        source_mesh_descriptors = stringr::str_replace_all(
          source_mesh_descriptors, "\\n", ""
        )
      ) |>
      dplyr::select(-pmid) |>
      dplyr::left_join(
        pubmed_source_ids, by = "source_id")
  }

  if (NROW(other_source_ids) > 0) {
    other_source_ids <- other_source_ids |>
      dplyr::filter(source != "PubMed") |>
      dplyr::mutate(name = source_id,
                    source_keywords = NA,
                    source_mesh_descriptors = NA) |>
      dplyr::mutate(link = dplyr::case_when(
        source == "FDA" ~
          "<a href='https://www.fda.gov/drugs/resources-information-approved-drugs/oncology-cancer-hematologic-malignancies-approval-notifications' target='_blank'>FDA approvals</a>",
        source == "NCCN" ~
          "<a href='https://www.nccn.org/guidelines/category_1' target='_blank'>NCCN guidelines</a>",
        source == "EMA" ~ "<a href='https://www.ema.europa.eu/en/medicines/field_ema_web_categories%253Aname_field/Human' target='_blank'>European Medicines Agency (EMA)</a>",
        source == "clinicaltrials.gov" ~ paste0(
          "<a href='https://clinicaltrials.gov/ct2/show/",
          source_id, "' target='_blank'>",
          source_id, "</a>")
      ))
  }

  all_literature <- pubmed_source_ids |>
    dplyr::select(source,
                  source_id,
                  source_entity,
                  source_keywords,
                  source_mesh_descriptors,
                  name, link)

  if(NROW(other_source_ids) > 0){
    all_literature <-
      other_source_ids |>
      dplyr::bind_rows(pubmed_source_ids)
  }

  return(all_literature)


}



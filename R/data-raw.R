source_ids <- as.data.frame(
 readr::read_csv("data-raw/source_data.csv",
                 show_col_types = F)
)

missing <- source_ids |>
  dplyr::anti_join(literatureVault::vault, by = "source_id")

add_to_vault <- literatureVault:::get_literature(
 literature_df = missing, chunk_size = 240
)
#
vault <- literatureVault::vault |>
  dplyr::bind_rows(add_to_vault)

usethis::use_data(vault, overwrite = T)

library(tidyverse)
library(magrittr)

data_path <- "./data/parsed"
files <- dir(data_path)

load_cdr <- function(file_names) {
  df <- NULL
  for (i in file_names) {
    print(i)
    if (!is_tibble(df)) {
      df <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
      tmp <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
    } else {
      df %<>% full_join(tmp)
    }
  }
  df
}

files <- files[str_detect(files, "(?i)isaura")]
# match only lines that does not contains "falho" or "teste"
files <- files[str_detect(files, "^(?:(?!falho|teste).)*$")]
# merge file name and file path
files <- sapply(files, function(x) {
  file.path(data_path, x)
})

cdr_df  <- read_rds("./data/binary/isaura_compressed.rds")

cdr_df %<>% mutate(type = case_when(
                             str_detect(file, "Final") ~ "final",
                             str_detect(file, "Initial") ~ "initial"
)) %>%
            select(cdr3, type, everything())

cdr_df %<>% group_by(type, file)

cdr_df %>% summarise(quantity = sum(quantity)) %>% arrange(-quantity)

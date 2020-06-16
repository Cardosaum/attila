library(tidyverse)
library(magrittr)

data_path <- "./data/parsed"
files <- dir(data_path)

load_cdr <- function(file_names){
    df <- NULL
    for (i in file_names){
        print(i)
        if (!is_tibble(df)){
            df <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
        } else {
            tmp <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
            df %<>% full_join(tmp)
        }
    }
    df
}

files <- files[str_detect(files, "(?i)isaura")]
# match only lines that does not contains "falho" or "teste"
files <- files[str_detect(files, "^(?:(?!falho|teste).)*$")]
# merge file name and file path
files <- sapply(files, function(x){file.path(data_path, x)})

cdr_df <- load_cdr(files)

cdr_df %<>% group_by(file)

cdr_df %>% mutate(cdrp = quantity/sum(quantity)) %>%
            select(cdrp, cdr3, quantity, file) %>% 
            filter(cdrp == max(cdrp)) -> ham

# write_rds(cdr_df, file.path(data_path, "isaura_compressed.rds"))

# cdr_df <- read_rds(file.path(data_path, "isaura_compressed.rds"))

files_duplicated <- ham$file[duplicated(ham$file)] %>% unique()

ham %>% select(cdrp, cdr3, quantity, file) %>% 
        filter(file %in% files_duplicated) %>% 
        filter(cdrp == max(cdrp)) -> bagar

wh <- read_csv(file.path(data_path, "mariajac_isaura_2_2_isaura_H_0eX4g_L_0bX4c_VH_FinalRound_VLR42c_S10_L001_R1_001aafreq.csv"), col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")

cdr_df %>% mutate(cdrp = quantity/sum(quantity)) %>% 
            slice_head(n = 5) -> cdr

cdr %<>% select(cdr3, cdrp, everything()) %>% 
            arrange(-cdrp, .by_group = TRUE)

cdr %>% slice_head(n = 1) %>% 
ggplot() +
    geom_point(aes(cdrp, log(quantity), color = file)) +
    scale_color_discrete(labels = NULL)

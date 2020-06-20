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
            df %<>% mutate(cdrp = quantity/sum(quantity))
        } else {
            tmp <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
            tmp %<>% mutate(cdrp = quantity/sum(quantity))
            df %<>% full_join(tmp)
        }
    }
    df %<>%
            select(cdr3, cdrp, quantity, everything()) %>% 
            arrange(-cdrp, -quantity)
}

######################################################################
#
# load_files_cdr.R
# 
# This file is a helper script containing code to handle
# loading of cdr3 raw files into tidy dataframes
#
# Copyright: (c) 2019-2020 Matheus Cardoso <github.com/cardosaum>
# License: Apache 2.0 <apache.org/licenses/LICENSE-2.0>
#
######################################################################
    

##############################
# User configuration
# 
# default file paths
# note that this paths are also loaded into the session of all other
# scripts that `source` this file. Please, consider with caution if
# you should to change this paths.

data_path        <- "./data"
data_path_raw    <- file.path(data_path, "raw")
data_path_parsed <- file.path(data_path, "parsed")
data_path_binary <- file.path(data_path, "binary")
files_parsed     <- dir(data_path_parsed)
files_binary     <- dir(data_path_binary)
##############################

library(tidyverse, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(knitr, quietly = TRUE)
set.seed(42)

load_and_merge_cdr_data <- function(file_names){
    # load all passed files and merge them in a tidy
    # dataframe format
    
    # we need to ensure if we are in the first or 
    # subsequent loops; use df to verify.
    # Reason: we need to have a tibble datafreme 
    # in order to merge the subsequent ones.
    df <- NULL
    
    for (i in file_names){
        
        # if this is our first loop, initialize the
        # dataframe
        if (!is_tibble(df)){
            df <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
            df %<>% mutate(cdrp = quantity/sum(quantity))
        } else {
            tmp <- read_csv(i, col_types = "ciiddddddddiiiiiiiiiiiiiiiiiiiiiiiiiic")
            tmp %<>% mutate(cdrp = quantity/sum(quantity))
            df %<>% full_join(tmp)
        }
    }
    
    # we reorder the collumns. 
    # cdr3, cdrp and quantity should appear firts
    # (easier to read in the data analysis part)
    df %<>%
            select(cdr3, cdrp, quantity, everything()) %>% 
            arrange(-cdrp, -quantity)
}

check_sha1_cdr <- function(file_names){
    # check if the same list of files already exists as a binary.
    # doing this, we reduce the need to always load the csv file
    # and merge them all over again. 
    #
    # (depending on the list of passed files, this process of
    # calling the `load_and_merge_cdr_data` function can be
    # really time consuming)
    
    files_sha1 <- get_sha1_cdr(file_names)
    if (files_sha1 %in% files_binary){
        return(TRUE)
    } else{
        return(FALSE)
    }
}

# in current implementation we rely sollely in the file names to asses if they 
# have already been parsed previously.
# TODO: use some more enduring way. we need to assess the binary content of the 
# file, not their names.
get_sha1_cdr <- function(file_names){digest::sha1(paste(file_names, collapse = ""))}

load_cdr <- function(file_names){
    # with all set, now we actualy read the data
    
    files_sha1 <- get_sha1_cdr(file_names)
    files_sha1_path <- file.path(data_path_binary, files_sha1)
    
    if (check_sha1_cdr(file_names)){
        df <- read_rds(files_sha1_path)
    } else{
        df <- load_and_merge_cdr_data(file_names)
        write_rds(df, file.path(data_path_binary, files_sha1))
    }
    df
}

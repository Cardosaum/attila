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
library(tidymodels, quietly = TRUE)
library(broom.mixed, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(plotly, quietly = TRUE)
library(ggforce, quietly = TRUE)
library(knitr, quietly = TRUE)
library(ggsci, quietly = TRUE)
library(AMR, quietly = TRUE)

# set custom theme for all plots
# credit: https://stackoverflow.com/a/16437625/10719703
ggplot <- function(...) ggplot2::ggplot(...) +
  scale_color_d3() +
  scale_fill_d3()

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

get_files_experiment <-
  function(regexPattern){
    # This function search for the expecified regexPattern and return the file names

    files_experiment <-
      list.files(
        pattern = regexPattern,
        recursive = TRUE,
        ignore.case = TRUE,
      )

    files_experiment
  }

get_files_experiment_rafael <-
  function(){
    # set default search for rafael's experiments
    files_experiment <-
      get_files_experiment(".*rafael.*resultado.*csv$")

    files_experiment
  }

get_files_experiment_carol <-
  function(){
    # set default search for carol's experiments
    files_experiment <-
      get_files_experiment(".*carolnovoattila.*csv$")

    files_experiment
  }

get_files_experiment_thais <-
  function(){
    # set default search for thais's experiments
    files_experiment <-
      get_files_experiment(".*thais_.*csv$")

    files_experiment
  }


load_cdr_carol <-
  function(){
    # set default load for carol's experiment

    cdr <-
      load_cdr(get_files_experiment_carol())
    if (cdr_exist_store(cdr)) {
      cdr <- cdr_retrieve(cdr)
    } else {
      cdr %>%
        cdr_preprocess(experiment = "carol")
    }

    cdr
  }

load_cdr_rafael <-
  function(){
    # set default load for rafael's experiment

    cdr <-
      load_cdr(get_files_experiment_rafael())

    if (cdr_exist_store(cdr)) {
      cdr <- cdr_retrieve(cdr)
    } else {
      cdr %<>%
        cdr_preprocess(experiment = "rafael")
    }

    cdr
  }

load_cdr_thais <-
  function(){
    # set default load for thais's experiment

    cdr <-
      load_cdr(get_files_experiment_thais())

    if (cdr_exist_store(cdr)) {
      cdr <- cdr_retrieve(cdr)
    } else {
      cdr %>%
        cdr_preprocess(experiment = "thais")
    }

    cdr
  }

cdr_preprocess <-
  function(cdr, experiment, store_cdr = TRUE, ...){
    # do basic preprocessing in data passed
    # this step is common to all experiments;
    # 1) calculate fold-changes
    # 2) calculate rich cdrs
    #
    # cdr: dataframe to preprocess.
    # experiment: character, from which experiment is this dataset?
    #             (current available experiments:
    #             Thais;
    #             Carol;
    #             Rafael.)
    # store_cdr: logical, store preprocessed dataframe in "binary" folder?
    # ...: argumets passed to other functions.

    # check if already exists a preproced dataframe in binary format

    if (cdr_exist_store(cdr)) {
      cdr <- cdr_retrieve(cdr)
    } else{
      if (store_cdr) {
        cdr_sha1 <- digest::sha1(cdr)
      }
      if (experiment == "rafael") {
        cdr %<>%
          mutate(
            expgroup = case_when(
              str_detect(file, "rafael") ~ "rafael",
              TRUE ~ "unknown"),
            cycle = case_when(
              str_detect(file, "R0") ~ "R0",
              str_detect(file, "R4") ~ "R4",
              TRUE ~ "unknown"),
            time = case_when(
                str_detect(file, "Initial") ~ "initial",
                str_detect(file, "Final") ~ "final",
                TRUE ~ "unknown")) %>%
          select(cdr3, cycle, time, expgroup, everything())
      } else if (experiment == "thais") {
        cdr %<>%
          mutate(
            expgroup = case_when(
              str_detect(file, "thaisnovoheader") ~ "nh",
              str_detect(file, "thais_29") ~ "thais_29",
              str_detect(file, "thais_66") ~ "thais_66",
              TRUE ~ "unknown"),
            cycle = case_when(
              str_detect(file, "R0_R2") ~ "R0_R2",
              str_detect(file, "R0_R3") ~ "R0_R3",
              str_detect(file, "R0_R4") ~ "R0_R4",
              str_detect(file, "R2_R3") ~ "R2_R3",
              str_detect(file, "R2_R4") ~ "R2_R4",
              str_detect(file, "R3_R4") ~ "R3_R4",
              TRUE ~ "unknown"),
            time = case_when(
                str_detect(file, "Initial") ~ "initial",
                str_detect(file, "Final") ~ "final",
                TRUE ~ "unknown")) %>%
          select(cdr3, cycle, time, expgroup, everything())
      } else if (experiment == "carol") {
        cdr %<>%
          mutate(
            expgroup = case_when(
              str_detect(file, "carol") ~ "carol",
              TRUE ~ "unknown"),
            cycle = case_when(
              str_detect(file, "R3") ~ "R3",
              str_detect(file, "R4") ~ "R4",
              TRUE ~ "unknown"),
            time = case_when(
                str_detect(file, "Initial") ~ "initial",
                str_detect(file, "Final") ~ "final",
                TRUE ~ "unknown")) %>%
          select(cdr3, cycle, time, expgroup, everything())
      } else {
        stop("This experiment is not yet implemented.")
      }
      cdr %<>%
        group_by(cdr3, expgroup) %>%
        arrange(cycle, desc(time), .by_group = TRUE) %>%
        mutate(
          fcp = cdrp / lag(cdrp, default = first(cdrp)),
          fcq = quantity / lag(quantity, default = first(quantity))
        ) %>%
        select(cdr3:quantity, fcp, fcq, everything())

      cdr %<>%
        group_by(expgroup, cycle, time) %>%
        arrange(desc(fcp)) %>%
        slice_head(prop = .1) %>%
        mutate(
          threshold = mean(log10(fcp))) %>%
        mutate(
          rich = if_else(
              (log10(fcp) >= threshold) &
                # (time == "final") &
                # (str_detect(cycle, "R0")),
                (time == "final"),
            "rich",
            "medium")) %>%
        full_join(cdr) %>%
        mutate(
          rich = if_else(
            is.na(rich),
            "poor",
            rich)) %>%
        mutate(
          rich = factor(rich,
                        levels = c("rich", "medium", "poor"))) %>%
        mutate(
          threshold = if_else(
            is.na(threshold),
            0,
            threshold))

      if (store_cdr) {
        write_rds(cdr, file.path(data_path_binary, cdr_sha1))
      }
    }

    cdr
  }

cdr_exist_store <- function(cdr){
  # check if there already exists this same dataframe as binary format
  # this avoid desnecessary preprocess steps

  cdr_sha1 <- digest::sha1(cdr)
  cdr_path <- file.path(data_path_binary, cdr_sha1)

  if (file.exists(cdr_path)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

cdr_store <- function(cdr, overwrite = TRUE){
  # this function save the preproced dataframe as a RDS binary.
  # preprocess steps can take a considerable amount of time,
  # and, since they do not change (or, at least, it's not common
  # to change), makes sense to store de dataframe as a binary.
  #
  # this funcition is by default active in the cdr_preprocess
  # function.

  cdr_sha1 <-
    digest::sha1(cdr)

  # store the preproced dataframe in "binary" folder
  if (overwrite) {
    write_rds(cdr, file.path(data_path_binary, cdr_sha1))
  }

}

cdr_retrieve <- function(cdr) {
  # this function loads an already existing cdr dataframe

  cdr_path <-
    file.path(
      data_path_binary,
      digest::sha1(cdr))

  cdr <-
    read_rds(cdr_path)

  cdr
}

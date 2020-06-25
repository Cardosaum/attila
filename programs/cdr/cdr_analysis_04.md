---
title: "Análise exploratória de sequências CDR3"
subtitle: "Quarta Iteração"
author: "Matheus Cardoso"
date: "Jun 20, 2020"
output:
  html_document:
    fig_caption: yes
    highlight: zenburn
    keep_md: yes
    theme: simplex
    toc: yes
    toc_float: yes
    number_sections: yes
    code_folding: hide
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: console
---


```r
knitr::opts_chunk$set(echo = TRUE, include = TRUE, cache = TRUE, warning = FALSE, message = FALSE)
```

# Introdução

Essa análise será baseada na segunda iteração, porém mudando os dados analisados.
Dessa vez, realizarei as análises nos dados de sequenciamento da Thais.


```r
# this R script does all the setup for us
source("./load_files.R")


data_files_selected  <-
    list.files(pattern = ".*carolnovoattila.*csv$", recursive = TRUE, ignore.case = TRUE)

data_files_and_size <- sapply(data_files_selected, file.size)

files_to_include_in_dataframe <- tibble(
                                    "Files" = names(data_files_and_size),
                                    "Size (in MB)" = (data_files_and_size/1E6))

kable(
     files_to_include_in_dataframe,
     caption = "List of files to be analysed")
```



Table: List of files to be analysed

Files                                                                                                      Size (in MB)
--------------------------------------------------------------------------------------------------------  -------------
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_carolnovoattila_VH_FinalRound_R4_S4_L001_R1_001aafreq.csv          14.51756
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_carolnovoattila_VH_InitialRound_R3_S3_L001_R1_001aafreq.csv        65.19080

# métodos

## Processamento dos dados


```r
cdr_df <- load_cdr(names(data_files_and_size))
```


```r
cdr_df %<>% 
          mutate(expgroup = case_when(
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

cdr_df %<>% 
          group_by(cdr3, expgroup) %>% 
          arrange(cycle, desc(time), .by_group = TRUE) %>% 
          mutate(
            fcp = cdrp / lag(cdrp, default = first(cdrp)),
            fcq = quantity / lag(quantity, default = first(quantity))) %>% 
          select(cdr3:quantity, fcp, fcq, everything())

glimpse(cdr_df)
```

```
## Rows: 358,265
## Columns: 44
## Groups: cdr3, expgroup [337,156]
## $ cdr3      <chr> "A", "A", "AA", "AAAAADTRPPLHLHVDY", "AAAAAGGGNWFDP", "AAAA…
## $ cycle     <chr> "R3", "R4", "R3", "R3", "R3", "R3", "R3", "R3", "R3", "R3",…
## $ time      <chr> "initial", "final", "initial", "initial", "initial", "initi…
## $ expgroup  <chr> "carol", "carol", "carol", "carol", "carol", "carol", "caro…
## $ cdrp      <dbl> 4.651968e-05, 3.725103e-05, 8.893468e-06, 6.841129e-07, 1.3…
## $ quantity  <int> 68, 56, 13, 1, 2, 1, 2, 3, 2, 1, 9, 1, 1, 1, 1, 10, 14, 1, …
## $ fcp       <dbl> 1.0000000, 0.8007586, 1.0000000, 1.0000000, 1.0000000, 1.00…
## $ fcq       <dbl> 1.0000000, 0.8235294, 1.0000000, 1.0000000, 1.0000000, 1.00…
## $ length    <int> 1, 1, 2, 17, 13, 13, 13, 11, 9, 12, 12, 13, 9, 15, 15, 7, 8…
## $ MW        <dbl> 89.0932, 89.0932, 160.1711, 1817.9976, 1204.2476, 1156.2048…
## $ AV        <dbl> 0.0000, 0.0000, 0.0000, 0.0588, 0.1538, 0.0769, 0.1538, 0.1…
## $ IP        <dbl> 5.5700, 5.5700, 5.5700, 5.9823, 4.0500, 4.0500, 5.0781, 8.6…
## $ flex      <dbl> 0.7040, 0.7040, 0.7040, 0.7363, 0.7498, 0.7514, 0.7037, 0.7…
## $ gravy     <dbl> 1.8000, 1.8000, 1.8000, -0.1353, 0.0846, 0.1923, 0.9231, -0…
## $ SSF_Helix <dbl> 0.0000, 0.0000, 0.0000, 0.2353, 0.1538, 0.1538, 0.3077, 0.1…
## $ SSF_Turn  <dbl> 0.0000, 0.0000, 0.0000, 0.1176, 0.3846, 0.3846, 0.1538, 0.0…
## $ SSF_Sheet <dbl> 1.0000, 1.0000, 1.0000, 0.4118, 0.3846, 0.3846, 0.3846, 0.4…
## $ n_A       <int> 1, 1, 2, 5, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4,…
## $ n_C       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_D       <int> 0, 0, 0, 2, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 1, 1, 0, 2, 1,…
## $ n_E       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,…
## $ n_F       <int> 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1,…
## $ n_G       <int> 0, 0, 0, 0, 3, 3, 2, 1, 1, 1, 1, 1, 0, 3, 3, 0, 0, 1, 1, 1,…
## $ n_H       <int> 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,…
## $ n_I       <int> 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ n_K       <int> 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_L       <int> 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0,…
## $ n_M       <int> 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_N       <int> 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_P       <int> 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_Q       <int> 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_R       <int> 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ n_S       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ n_T       <int> 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_V       <int> 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 2, 0, 0, 0, 0, 0,…
## $ n_W       <int> 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_Y       <int> 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 2, 2, 0, 0, 1, 1, 1, 1, 1,…
## $ aliphatic <int> 1, 1, 2, 10, 9, 10, 9, 6, 8, 8, 8, 9, 5, 10, 10, 4, 5, 6, 5…
## $ aromatic  <int> 0, 0, 0, 1, 2, 1, 2, 2, 0, 1, 2, 3, 3, 1, 1, 1, 2, 2, 2, 2,…
## $ neutral   <int> 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,…
## $ positive  <int> 0, 0, 0, 3, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1,…
## $ negative  <int> 0, 0, 0, 2, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 2, 1, 0, 2, 1,…
## $ invalid   <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ file      <chr> "heidi_ATTILA_ANALISYS_TO_BACKUP_carolnovoattila_VH_Initial…
```

## Análise exploratória


```r
cdr_df %>% 
          group_by(expgroup, cycle) %>% 
          summarise(
            "No. distinct CDR3 sequences" = n(),
            "Sum of all CDR3 sequences" = sum(quantity)) -> cdr_df_summary_per_experiment

ggplot(cdr_df_summary_per_experiment) +
  geom_bar(aes(
            cycle,
            `No. distinct CDR3 sequences`,
            fill = cycle),
          stat = "identity") +
  facet_grid(. ~ expgroup) +
  labs(
    title = "Difference in number of distinct CDR3 sequences"
  ) +
  theme(axis.text.x = element_text(angle = 45))
```

![](cdr_analysis_04_files/figure-html/data_exploration_1-1.png)<!-- -->

```r
ggplot(cdr_df_summary_per_experiment) +
  geom_bar(aes(
            cycle,
            `Sum of all CDR3 sequences`,
            fill = cycle),
          stat = "identity") +
  facet_grid(. ~ expgroup) +
  labs(
    title = "Difference in number of total CDR3 sequences"
  ) +
  theme(axis.text.x = element_text(angle = 45))
```

![](cdr_analysis_04_files/figure-html/data_exploration_1-2.png)<!-- -->


```r
cdr_df %>% 
            ggplot() +
              geom_density(aes(cdrp, color = cycle, fill = cycle), alpha = 0.2) +
              labs(
                title = "Density plot of CDR3 Percentage of Prevalence",
                x = "CDR3 Percentage of Prevalence (CDRP)"
              ) +
              xlim(0.01, 0.33) +
              facet_grid(cycle ~ .)
```

![](cdr_analysis_04_files/figure-html/data_exploration_2-1.png)<!-- -->

```r
cdr_df %>% 
            ggplot() +
              geom_density(
                aes(
                  x = quantity,
                  color = cycle,
                  fill = cycle),
                alpha = 0.2) +
              labs(
                title = "Density plot of Quantity Distribuition",
                x = "Quantity"
              ) + 
            facet_zoom(xlim = c(0, 500))
```

![](cdr_analysis_04_files/figure-html/data_exploration_2-2.png)<!-- -->


```r
ggplot(cdr_df) +
  geom_point(aes(cdrp, quantity, color = cycle)) +
  facet_grid(expgroup ~ .) +
  labs(
    title = "Correlation between CDR3 percentage of prevalence and Quantity",
    x = "CDR3 Percentage of Prevalence (CDRP)",
    y = "Quantity"
  )
```

![](cdr_analysis_04_files/figure-html/data_exploration_3-1.png)<!-- -->

```r
ggplot(cdr_df) +
  geom_point(aes(fcp, fcq, color = cycle), alpha = 0.3) +
  facet_grid(expgroup ~ .) + 
  labs(
    title = "Correlation between Fold-Changes",
    caption = "Note how different the correlation is between experiment groups.",
    x = "Fold-Change by Percentage",
    y = "Fold-Change by Quantity"
  )
```

![](cdr_analysis_04_files/figure-html/data_exploration_3-2.png)<!-- -->


```r
# search for enriched CDR3s
for (i in c("fcp", "fcq")){
cdr_df %>% 
          filter(cycle  == "R4") %>% 
          group_by(cycle) %>% 
          arrange(desc(.[[i]]), .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          select(cdr3:fcq) -> tmp_plot_1
  
print(
  ggplot(tmp_plot_1) +
    geom_density(aes(.data[[i]], color = cycle, fill = cycle), alpha = 0.3)
)

print(
  ggplot(tmp_plot_1, aes(cycle, .data[[i]])) +
    geom_violin(aes(fill = cycle)) +
    geom_jitter(aes(shape = cycle), alpha = 0.7)
)

}
```

![](cdr_analysis_04_files/figure-html/data_exploration_4-1.png)<!-- -->![](cdr_analysis_04_files/figure-html/data_exploration_4-2.png)<!-- -->![](cdr_analysis_04_files/figure-html/data_exploration_4-3.png)<!-- -->![](cdr_analysis_04_files/figure-html/data_exploration_4-4.png)<!-- -->


### PCA for the two experiments


```r
# PCA for all groups
cdr_df %>%
          group_by(cycle) %>% 
          arrange(-quantity, .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          ungroup() %>% 
          select(where(is.numeric), cycle) -> cdr_rich_df_all

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_all %<>% select(!c(which(apply(cdr_rich_df_all, 2, var)==0)))

cdr_rich_pca_all_groups <- cdr_rich_df_all %>% pca()

cdr_rich_df_all %>%
          group_by(cycle) %>% 
          arrange(-quantity, .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          .$cycle -> cdr_rich_pca_color_all_groups

ggplot_pca(
  cdr_rich_pca_all_groups,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_all_groups
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences"
  )
```

![](cdr_analysis_04_files/figure-html/pca_1-1.png)<!-- -->

```r
summary(cdr_rich_pca_all_groups)
```

```
## Importance of components:
##                           PC1    PC2     PC3     PC4     PC5    PC6     PC7
## Standard deviation     2.3540 1.9951 1.82334 1.66805 1.57822 1.4718 1.39500
## Proportion of Variance 0.1458 0.1047 0.08749 0.07322 0.06555 0.0570 0.05121
## Cumulative Proportion  0.1458 0.2506 0.33807 0.41129 0.47684 0.5338 0.58505
##                            PC8     PC9    PC10    PC11    PC12   PC13    PC14
## Standard deviation     1.35521 1.30097 1.27632 1.21780 1.12437 1.0479 1.01313
## Proportion of Variance 0.04833 0.04454 0.04287 0.03903 0.03327 0.0289 0.02701
## Cumulative Proportion  0.63338 0.67792 0.72079 0.75982 0.79308 0.8220 0.84899
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     0.98092 0.93018 0.84571 0.81564 0.77335 0.73108 0.67746
## Proportion of Variance 0.02532 0.02277 0.01882 0.01751 0.01574 0.01407 0.01208
## Cumulative Proportion  0.87432 0.89709 0.91591 0.93341 0.94915 0.96322 0.97530
##                           PC22    PC23    PC24    PC25    PC26    PC27    PC28
## Standard deviation     0.58942 0.56206 0.37125 0.27487 0.17582 0.13508 0.09635
## Proportion of Variance 0.00914 0.00831 0.00363 0.00199 0.00081 0.00048 0.00024
## Cumulative Proportion  0.98444 0.99275 0.99638 0.99837 0.99918 0.99966 0.99991
##                           PC29    PC30      PC31      PC32      PC33      PC34
## Standard deviation     0.05915 0.01011 0.0002584 1.049e-15 5.183e-16 4.472e-16
## Proportion of Variance 0.00009 0.00000 0.0000000 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.00000 1.00000 1.0000000 1.000e+00 1.000e+00 1.000e+00
##                             PC35      PC36      PC37      PC38
## Standard deviation     3.423e-16 3.237e-16 2.883e-16 1.152e-16
## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00
```

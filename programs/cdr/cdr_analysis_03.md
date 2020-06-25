---
title: "Análise exploratória de sequências CDR3"
subtitle: "Terceira iteração - Dados Thais"
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
    list.files(pattern = ".*thais_.*csv$", recursive = TRUE, ignore.case = TRUE)

data_files_and_size <- sapply(data_files_selected, file.size)

files_to_include_in_dataframe <- tibble(
                                    "Files" = names(data_files_and_size),
                                    "Size (in MB)" = (data_files_and_size/1E6))

kable(
     files_to_include_in_dataframe,
     caption = "List of files to be analysed")
```



Table: List of files to be analysed

Files                                                                                                                                 Size (in MB)
-----------------------------------------------------------------------------------------------------------------------------------  -------------
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R2_thais_29_R0_R2_VH_FinalRound_VCL29VHR2_S2_L001_R1_001aafreq.csv          2.077444
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R2_thais_29_R0_R2_VH_InitialRound_VCL29VHR0_S1_L001_R1_001aafreq.csv        0.316471
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R3_thais_29_R0_R3_VH_FinalRound_VCL29VHR3_S3_L001_R1_001aafreq.csv          0.343624
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R3_thais_29_R0_R3_VH_InitialRound_VCL29VHR0_S1_L001_R1_001aafreq.csv        0.316471
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R4_thais_29_R0_R4_VH_FinalRound_VCL29VHR4_S4_L001_R1_001aafreq.csv          0.395904
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R4_thais_29_R0_R4_VH_InitialRound_VCL29VHR0_S1_L001_R1_001aafreq.csv        0.316471
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R2_R3_thais_29_R2_R3_VH_FinalRound_VCL29VHR3_S3_L001_R1_001aafreq.csv          0.343624
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R2_R3_thais_29_R2_R3_VH_InitialRound_VCL29VHR2_S2_L001_R1_001aafreq.csv        2.094224
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R2_R4_thais_29_R2_R4_VH_FinalRound_VCL29VHR4_S4_L001_R1_001aafreq.csv          0.395904
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R2_R4_thais_29_R2_R4_VH_InitialRound_VCL29VHR2_S2_L001_R1_001aafreq.csv        2.094224
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R3_R4_thais_29_R3_R4_VH_FinalRound_VCL29VHR4_S4_L001_R1_001aafreq.csv          0.395904
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R3_R4_thais_29_R3_R4_VH_InitialRound_VCL29VHR3_S3_L001_R1_001aafreq.csv        0.346386
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R2_thais_66_R0_R2_VH_FinalRound_VCL66VHR2_S6_L001_R1_001aafreq.csv          0.719382
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R2_thais_66_R0_R2_VH_InitialRound_VCL66VHR0_S5_L001_R1_001aafreq.csv       10.994050
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R3_thais_66_R0_R3_VH_FinalRound_VCL66VHR3_S7_L001_R1_001aafreq.csv          0.484975
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R3_thais_66_R0_R3_VH_InitialRound_VCL66VHR0_S5_L001_R1_001aafreq.csv       10.994050
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R4_thais_66_R0_R4_VH_FinalRound_VCL66VHR4_S8_L001_R1_001aafreq.csv          2.115096
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R0_R4_thais_66_R0_R4_VH_InitialRound_VCL66VHR0_S5_L001_R1_001aafreq.csv       10.994050
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R2_R3_thais_66_R2_R3_VH_FinalRound_VCL66VHR3_S7_L001_R1_001aafreq.csv          0.484975
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R2_R3_thais_66_R2_R3_VH_InitialRound_VCL66VHR2_S6_L001_R1_001aafreq.csv        0.725188
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R2_R4_thais_66_R2_R4_VH_FinalRound_VCL66VHR4_S8_L001_R1_001aafreq.csv          2.115096
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R2_R4_thais_66_R2_R4_VH_InitialRound_VCL66VHR2_S6_L001_R1_001aafreq.csv        0.725188
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R3_R4_thais_66_R3_R4_VH_FinalRound_VCL66VHR4_S8_L001_R1_001aafreq.csv          2.115096
data/parsed/heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_66_R3_R4_thais_66_R3_R4_VH_InitialRound_VCL66VHR3_S7_L001_R1_001aafreq.csv        0.488909

# métodos

## Processamento dos dados


```r
cdr_df <- load_cdr(names(data_files_and_size))
```


```r
cdr_df %<>% 
          mutate(expgroup = case_when(
                                  str_detect(file, "thaisnovoheader") ~ "nh",
                                  str_detect(file, "thais_29") ~ "29",
                                  str_detect(file, "thais_66") ~ "66",
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

cdr_df %<>% 
          group_by(cdr3, expgroup, cycle) %>% 
          arrange(desc(time), .by_group = TRUE) %>% 
          mutate(
            fcp = cdrp / lag(cdrp, default = first(cdrp)),
            fcq = quantity / lag(quantity, default = first(quantity))) %>% 
          select(cdr3:quantity, fcp, fcq, everything())

glimpse(cdr_df)
```

```
## Rows: 210,432
## Columns: 44
## Groups: cdr3, expgroup, cycle [205,669]
## $ cdr3      <chr> "A", "A", "A", "A", "A", "A", "AA", "AA", "AA", "AAAAGHFDY"…
## $ cycle     <chr> "R0_R2", "R2_R3", "R2_R4", "R0_R2", "R0_R3", "R0_R4", "R0_R…
## $ time      <chr> "final", "initial", "initial", "initial", "initial", "initi…
## $ expgroup  <chr> "29", "29", "29", "66", "66", "66", "66", "66", "66", "66",…
## $ cdrp      <dbl> 3.974199e-06, 3.974199e-06, 3.974199e-06, 3.210531e-05, 3.2…
## $ quantity  <int> 1, 1, 1, 4, 4, 4, 2, 2, 2, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ fcp       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ fcq       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ length    <int> 1, 1, 1, 1, 1, 1, 2, 2, 2, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8,…
## $ MW        <dbl> 89.0932, 89.0932, 89.0932, 89.0932, 89.0932, 89.0932, 160.1…
## $ AV        <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0…
## $ IP        <dbl> 5.5700, 5.5700, 5.5700, 5.5700, 5.5700, 5.5700, 5.5700, 5.5…
## $ flex      <dbl> 0.7040, 0.7040, 0.7040, 0.7040, 0.7040, 0.7040, 0.7040, 0.7…
## $ gravy     <dbl> 1.8000, 1.8000, 1.8000, 1.8000, 1.8000, 1.8000, 1.8000, 1.8…
## $ SSF_Helix <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0…
## $ SSF_Turn  <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0…
## $ SSF_Sheet <dbl> 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0…
## $ n_A       <int> 1, 1, 1, 1, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,…
## $ n_C       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,…
## $ n_D       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ n_E       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_F       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_G       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ n_H       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_I       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_K       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_L       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_M       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_N       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_P       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1,…
## $ n_Q       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_R       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_S       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_T       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ n_V       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_W       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_Y       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ aliphatic <int> 1, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5, 5, 6, 6, 6, 5, 5, 5, 6, 6,…
## $ aromatic  <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ neutral   <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 1, 1,…
## $ positive  <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ negative  <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ invalid   <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ file      <chr> "heidi_ATTILA_ANALISYS_TO_BACKUP_thais_thais_29_R0_R2_thais…
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

![](cdr_analysis_03_files/figure-html/data_exploration_1-1.png)<!-- -->

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

![](cdr_analysis_03_files/figure-html/data_exploration_1-2.png)<!-- -->


```r
cdr_df %>% 
            ggplot() +
              geom_density(aes(cdrp, color = expgroup, fill = expgroup), alpha = 0.2) +
              labs(
                title = "Density plot of CDR3 Percentage of Prevalence",
                x = "CDR3 Percentage of Prevalence (CDRP)"
              ) +
              xlim(0.01, 0.33) +
              facet_grid(expgroup ~ .)
```

![](cdr_analysis_03_files/figure-html/data_exploration_2-1.png)<!-- -->

```r
cdr_df %>% 
            ggplot() +
              geom_density(
                aes(
                  x = quantity,
                  color = expgroup,
                  fill = expgroup),
                alpha = 0.2) +
              labs(
                title = "Density plot of Quantity Distribuition",
                x = "Quantity"
              ) + 
            facet_zoom(xlim = c(0, 100))
```

![](cdr_analysis_03_files/figure-html/data_exploration_2-2.png)<!-- -->


```r
ggplot(cdr_df) +
  geom_point(aes(cdrp, quantity, color = cycle)) +
  facet_grid(expgroup ~ .) +
  labs(
    title = "Correlation between CDR3 percentage of prevalence and Quantity",
    subtitle = "With outlier",
    x = "CDR3 Percentage of Prevalence (CDRP)",
    y = "Quantity"
  )
```

![](cdr_analysis_03_files/figure-html/data_exploration_3-1.png)<!-- -->

```r
ggplot(cdr_df) +
  geom_point(aes(fcp, fcq, color = cycle), alpha = 0.3) +
  facet_grid(expgroup ~ .) + 
  labs(
    title = "Correlation between Fold-Changes",
    subtitle = "With outlier",
    caption = "Note how different the correlation is between experiment groups.",
    x = "Fold-Change by Percentage",
    y = "Fold-Change by Quantity"
  )
```

![](cdr_analysis_03_files/figure-html/data_exploration_3-2.png)<!-- -->


```r
# search for enriched CDR3s
for (i in c("fcp", "fcq")){
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(.[[i]]), .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          select(cdr3:fcq) -> tmp_plot_1
  
print(
  ggplot(tmp_plot_1) +
    geom_density(aes(.data[[i]], color = expgroup, fill = expgroup), alpha = 0.3)
)

print(
  ggplot(tmp_plot_1, aes(expgroup, log10(.data[[i]]))) +
    geom_violin(aes(fill = expgroup)) +
    geom_jitter(aes(shape = expgroup), alpha = 0.7)
)

}
```

![](cdr_analysis_03_files/figure-html/data_exploration_4-1.png)<!-- -->![](cdr_analysis_03_files/figure-html/data_exploration_4-2.png)<!-- -->![](cdr_analysis_03_files/figure-html/data_exploration_4-3.png)<!-- -->![](cdr_analysis_03_files/figure-html/data_exploration_4-4.png)<!-- -->


### PCA for the two experiments


```r
# PCA for all groups
cdr_df %>%
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_all

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_all %<>% select(!c(which(apply(cdr_rich_df_all, 2, var)==0)))

cdr_rich_pca_all_groups <- cdr_rich_df_all %>% pca()

cdr_rich_df_all %>%
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 100) %>% 
          .$expgroup -> cdr_rich_pca_color_all_groups

ggplot_pca(
  cdr_rich_pca_all_groups,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_all_groups
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences",
    subtitle = "All 3 experiments used to perform the PCA"
  )
```

![](cdr_analysis_03_files/figure-html/pca_1-1.png)<!-- -->

```r
summary(cdr_rich_pca_all_groups)
```

```
## Importance of components:
##                           PC1    PC2    PC3     PC4    PC5    PC6     PC7
## Standard deviation     2.4195 2.1329 1.9712 1.79589 1.6227 1.4522 1.32594
## Proportion of Variance 0.1541 0.1197 0.1022 0.08487 0.0693 0.0555 0.04627
## Cumulative Proportion  0.1541 0.2738 0.3760 0.46090 0.5302 0.5857 0.63196
##                            PC8     PC9    PC10    PC11    PC12    PC13   PC14
## Standard deviation     1.21802 1.16329 1.09564 1.08576 0.98651 0.97260 0.9628
## Proportion of Variance 0.03904 0.03561 0.03159 0.03102 0.02561 0.02489 0.0244
## Cumulative Proportion  0.67100 0.70661 0.73820 0.76922 0.79484 0.81973 0.8441
##                           PC15   PC16    PC17    PC18   PC19    PC20    PC21
## Standard deviation     0.93236 0.8911 0.86395 0.80100 0.7797 0.77105 0.70424
## Proportion of Variance 0.02288 0.0209 0.01964 0.01688 0.0160 0.01565 0.01305
## Cumulative Proportion  0.86700 0.8879 0.90754 0.92442 0.9404 0.95607 0.96912
##                           PC22    PC23    PC24    PC25    PC26    PC27    PC28
## Standard deviation     0.64015 0.44208 0.41301 0.37296 0.30010 0.27833 0.21170
## Proportion of Variance 0.01078 0.00514 0.00449 0.00366 0.00237 0.00204 0.00118
## Cumulative Proportion  0.97990 0.98505 0.98954 0.99320 0.99557 0.99760 0.99878
##                           PC29    PC30    PC31      PC32     PC33      PC34
## Standard deviation     0.17537 0.10125 0.07207 1.074e-15 4.42e-16 4.415e-16
## Proportion of Variance 0.00081 0.00027 0.00014 0.000e+00 0.00e+00 0.000e+00
## Cumulative Proportion  0.99959 0.99986 1.00000 1.000e+00 1.00e+00 1.000e+00
##                             PC35      PC36      PC37      PC38
## Standard deviation     3.127e-16 3.111e-16 2.646e-16 2.218e-16
## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00
```

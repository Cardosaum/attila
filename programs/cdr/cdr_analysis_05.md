---
title: "Análise exploratória de sequências CDR3"
subtitle: "Quinta Iteração - Dados Rafael"
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
Dessa vez, realizarei as análises nos dados de sequenciamento do Rafael.


```r
# this R script does all the setup for us
source("./load_files.R")


data_files_selected  <-
    list.files(pattern = ".*rafael.*resultado.*csv$", recursive = TRUE, ignore.case = TRUE)

data_files_and_size <- sapply(data_files_selected, file.size)

files_to_include_in_dataframe <- tibble(
                                    "Files" = names(data_files_and_size),
                                    "Size (in MB)" = (data_files_and_size/1E6))

kable(
     files_to_include_in_dataframe,
     caption = "List of files to be analysed")
```



Table: List of files to be analysed

Files                                                                                             Size (in MB)
-----------------------------------------------------------------------------------------------  -------------
data/parsed/brigido_rafael_Nestor_resultado_nestor_VH_FinalRound_rafaCD20_Vh_R4_R1aafreq.csv          1.214258
data/parsed/brigido_rafael_Nestor_resultado_nestor_VH_InitialRound_rafaCD20_Vh_R0_R1aafreq.csv        0.518894

# métodos

## Processamento dos dados


```r
cdr_df <- load_cdr(names(data_files_and_size))
```


```r
cdr_df %<>% 
          mutate(expgroup = case_when(
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
## Rows: 8,197
## Columns: 44
## Groups: cdr3, expgroup [7,924]
## $ cdr3      <chr> "A", "AAAGPIFDS", "AAARREFDI", "AADWFGAPFFDY", "AAEYYDSSVNF…
## $ cycle     <chr> "R0", "R0", "R4", "R4", "R4", "R4", "R4", "R4", "R4", "R4",…
## $ time      <chr> "initial", "initial", "final", "final", "final", "final", "…
## $ expgroup  <chr> "rafael", "rafael", "rafael", "rafael", "rafael", "rafael",…
## $ cdrp      <dbl> 4.450477e-05, 2.225239e-05, 1.809300e-05, 1.809300e-05, 1.8…
## $ quantity  <int> 2, 1, 1, 1, 1, 23, 1, 1, 1, 1, 1, 1, 1, 16, 2, 1, 1, 1, 1, …
## $ fcp       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ fcq       <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
## $ length    <int> 1, 9, 9, 12, 13, 13, 13, 15, 16, 16, 10, 17, 15, 15, 13, 13…
## $ MW        <dbl> 89.0932, 847.9116, 1048.1532, 1406.4948, 1467.4456, 1543.54…
## $ AV        <dbl> 0.0000, 0.1111, 0.1111, 0.4167, 0.2308, 0.3077, 0.3077, 0.2…
## $ IP        <dbl> 5.5700, 4.0500, 6.1150, 4.0500, 4.0500, 4.0500, 4.0500, 5.0…
## $ flex      <dbl> 0.7040, 0.7470, 0.7472, 0.7074, 0.7540, 0.7325, 0.7265, 0.7…
## $ gravy     <dbl> 1.8000, 0.7111, -0.3667, 0.2167, -0.6462, -0.6846, -0.4462,…
## $ SSF_Helix <dbl> 0.0000, 0.2222, 0.2222, 0.4167, 0.3077, 0.3846, 0.3846, 0.3…
## $ SSF_Turn  <dbl> 0.0000, 0.3333, 0.0000, 0.1667, 0.3077, 0.2308, 0.3077, 0.2…
## $ SSF_Sheet <dbl> 1.0000, 0.3333, 0.4444, 0.2500, 0.2308, 0.2308, 0.2308, 0.1…
## $ n_A       <int> 1, 3, 3, 3, 2, 2, 2, 2, 3, 2, 4, 4, 2, 2, 3, 3, 3, 2, 2, 2,…
## $ n_C       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_D       <int> 0, 1, 1, 2, 2, 2, 1, 2, 2, 1, 0, 1, 2, 2, 3, 3, 3, 1, 1, 1,…
## $ n_E       <int> 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_F       <int> 0, 1, 1, 3, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0,…
## $ n_G       <int> 0, 1, 0, 1, 0, 0, 1, 1, 3, 3, 2, 1, 1, 1, 0, 0, 0, 1, 1, 1,…
## $ n_H       <int> 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0,…
## $ n_I       <int> 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,…
## $ n_K       <int> 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_L       <int> 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_M       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_N       <int> 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ n_P       <int> 0, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0,…
## $ n_Q       <int> 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0,…
## $ n_R       <int> 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,…
## $ n_S       <int> 0, 1, 0, 0, 3, 2, 2, 0, 1, 2, 0, 2, 0, 0, 1, 1, 1, 1, 1, 1,…
## $ n_T       <int> 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,…
## $ n_V       <int> 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,…
## $ n_W       <int> 0, 0, 0, 1, 0, 0, 0, 2, 1, 1, 0, 1, 2, 2, 1, 1, 1, 0, 1, 1,…
## $ n_Y       <int> 0, 0, 0, 1, 2, 3, 3, 1, 1, 2, 0, 1, 1, 1, 1, 1, 2, 1, 3, 4,…
## $ aliphatic <int> 1, 6, 4, 5, 3, 3, 4, 6, 7, 8, 8, 7, 6, 7, 3, 4, 4, 3, 3, 3,…
## $ aromatic  <int> 0, 1, 1, 5, 3, 4, 4, 4, 2, 3, 0, 3, 4, 4, 4, 3, 3, 1, 5, 5,…
## $ neutral   <int> 0, 1, 0, 0, 4, 3, 3, 1, 2, 3, 1, 2, 1, 1, 3, 3, 3, 2, 1, 1,…
## $ positive  <int> 0, 0, 2, 0, 0, 0, 0, 2, 3, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 1,…
## $ negative  <int> 0, 1, 2, 2, 3, 3, 2, 2, 2, 1, 1, 4, 2, 2, 3, 3, 3, 1, 1, 1,…
## $ invalid   <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
## $ file      <chr> "brigido_rafael_Nestor_resultado_nestor_VH_InitialRound_raf…
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

![](cdr_analysis_05_files/figure-html/data_exploration_1-1.png)<!-- -->

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

![](cdr_analysis_05_files/figure-html/data_exploration_1-2.png)<!-- -->


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

![](cdr_analysis_05_files/figure-html/data_exploration_2-1.png)<!-- -->

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

![](cdr_analysis_05_files/figure-html/data_exploration_2-2.png)<!-- -->


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

![](cdr_analysis_05_files/figure-html/data_exploration_3-1.png)<!-- -->

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

![](cdr_analysis_05_files/figure-html/data_exploration_3-2.png)<!-- -->


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

![](cdr_analysis_05_files/figure-html/data_exploration_4-1.png)<!-- -->![](cdr_analysis_05_files/figure-html/data_exploration_4-2.png)<!-- -->![](cdr_analysis_05_files/figure-html/data_exploration_4-3.png)<!-- -->![](cdr_analysis_05_files/figure-html/data_exploration_4-4.png)<!-- -->


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

![](cdr_analysis_05_files/figure-html/pca_1-1.png)<!-- -->

```r
summary(cdr_rich_pca_all_groups)
```

```
## Importance of components:
##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
## Standard deviation     2.4290 1.9651 1.78283 1.69324 1.55049 1.47371 1.40788
## Proportion of Variance 0.1553 0.1016 0.08364 0.07545 0.06326 0.05715 0.05216
## Cumulative Proportion  0.1553 0.2569 0.34054 0.41599 0.47925 0.53640 0.58856
##                           PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     1.3505 1.25641 1.21050 1.16938 1.06165 1.02210 0.97900
## Proportion of Variance 0.0480 0.04154 0.03856 0.03599 0.02966 0.02749 0.02522
## Cumulative Proportion  0.6366 0.67810 0.71666 0.75265 0.78231 0.80980 0.83502
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     0.96813 0.91945 0.87865 0.85308 0.84700 0.78037 0.70196
## Proportion of Variance 0.02466 0.02225 0.02032 0.01915 0.01888 0.01603 0.01297
## Cumulative Proportion  0.85969 0.88194 0.90225 0.92140 0.94028 0.95631 0.96928
##                           PC22   PC23    PC24    PC25   PC26    PC27    PC28
## Standard deviation     0.69359 0.5375 0.38382 0.32315 0.2823 0.19929 0.11536
## Proportion of Variance 0.01266 0.0076 0.00388 0.00275 0.0021 0.00105 0.00035
## Cumulative Proportion  0.98194 0.9895 0.99342 0.99616 0.9983 0.99931 0.99966
##                           PC29    PC30      PC31      PC32      PC33      PC34
## Standard deviation     0.09628 0.06147 0.0006186 1.079e-15 6.774e-16 4.377e-16
## Proportion of Variance 0.00024 0.00010 0.0000000 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  0.99990 1.00000 1.0000000 1.000e+00 1.000e+00 1.000e+00
##                             PC35      PC36      PC37      PC38
## Standard deviation     4.303e-16 4.247e-16 2.764e-16 2.103e-16
## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00
```

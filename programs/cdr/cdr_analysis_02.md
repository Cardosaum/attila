---
title: "Análise exploratória de sequências CDR3"
subtitle: "Segunda iteração"
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

Essa análise será baseada na primeira iteração, porém mudando os dados analisados.
Dessa vez, usarei os seguintes dados:


```r
# this R script does all the setup for us
source("./load_files.R")

data_files_selected  <- readLines(file("cdr_analysis_02_selected_files.txt"))


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
data/parsed/brigido_renato_zika_ago18_ac_phage_zika_acid_VH_FinalRound_R4ac_VH_R1aafreq.csv           0.365703
data/parsed/brigido_renato_zika_ago18_ac_phage_zika_acid_VH_InitialRound_R0_VH_R1aafreq.csv           6.587344
data/parsed/brigido_renato_zika_ago18_phage_zika_VH_FinalRound_R4pep_VH_R1aafreq.csv                  0.567803
data/parsed/brigido_renato_zika_ago18_phage_zika_VH_InitialRound_R0_VH_R1aafreq.csv                   6.336536

# métodos

## Processamento dos dados


```r
cdr_df <- load_cdr(names(data_files_and_size))
```


```r
cdr_df %<>% 
          mutate(experiment = case_when(
                                  str_detect(file, "rafael.*R4") ~ "rafael_R4",
                                  str_detect(file, "rafael.*R0") ~ "rafael_R0",
                                  str_detect(file, "renato.*R4ac") ~ "renato_ac_R4",
                                  str_detect(file, "renato.*acid.*R0") ~ "renato_ac_R0",
                                  str_detect(file, "renato.*R4pep") ~ "renato_pep_R4",
                                  str_detect(file, "renato.*R0") ~ "renato_pep_R0"),
                cycle = case_when(
                                  str_detect(file, "R4") ~ "R4",
                                  str_detect(file, "R0") ~ "R0"),
                expgroup = case_when(
                                  str_detect(file, "rafael") ~ "rafael",
                                  str_detect(file, "renato.*acid") ~ "renato_acid",
                                  TRUE ~ "renato_peptide")) %>% 
          select(cdr3, cycle, expgroup, experiment, everything())

cdr_df %<>% 
          # select(cdr3:quantity) %>% 
          # filter(quantity > 400) %>% 
          group_by(cdr3, expgroup) %>% 
          # filter(n()>1) %>% 
          arrange(cycle, .by_group = TRUE) %>% 
          mutate(
            fcp = cdrp / lag(cdrp, default = first(cdrp)),
            fcq = quantity / lag(quantity, default = first(quantity))) %>% 
          select(cdr3:quantity, fcp, fcq, everything())

glimpse(cdr_df)
```

```
## Rows: 75,459
## Columns: 44
## Groups: cdr3, expgroup [74,244]
## $ cdr3       <chr> "A", "A", "A", "A", "A", "AAAAAGGGNWFDP", "AAAAAGGGNWFDP",…
## $ cycle      <chr> "R0", "R0", "R4", "R0", "R4", "R0", "R0", "R0", "R0", "R0"…
## $ expgroup   <chr> "rafael", "renato_acid", "renato_acid", "renato_peptide", …
## $ experiment <chr> "rafael_R0", "renato_ac_R0", "renato_ac_R4", "renato_pep_R…
## $ cdrp       <dbl> 4.450477e-05, 9.368705e-05, 1.267748e-04, 9.368705e-05, 4.…
## $ quantity   <int> 2, 6, 2, 6, 3, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 1…
## $ fcp        <dbl> 1.0000000, 1.0000000, 1.3531736, 1.0000000, 0.4409217, 1.0…
## $ fcq        <dbl> 1.0000000, 1.0000000, 0.3333333, 1.0000000, 0.5000000, 1.0…
## $ length     <int> 1, 1, 1, 1, 1, 13, 13, 9, 9, 9, 9, 8, 8, 9, 9, 9, 9, 8, 8,…
## $ MW         <dbl> 89.0932, 89.0932, 89.0932, 89.0932, 89.0932, 1204.2476, 12…
## $ AV         <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1538, 0.1538, 0.…
## $ IP         <dbl> 5.5700, 5.5700, 5.5700, 5.5700, 5.5700, 4.0500, 4.0500, 5.…
## $ flex       <dbl> 0.7040, 0.7040, 0.7040, 0.7040, 0.7040, 0.7498, 0.7498, 0.…
## $ gravy      <dbl> 1.8000, 1.8000, 1.8000, 1.8000, 1.8000, 0.0846, 0.0846, 0.…
## $ SSF_Helix  <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.1538, 0.1538, 0.…
## $ SSF_Turn   <dbl> 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.3846, 0.3846, 0.…
## $ SSF_Sheet  <dbl> 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.3846, 0.3846, 0.…
## $ n_A        <int> 1, 1, 1, 1, 1, 5, 5, 4, 4, 4, 4, 4, 4, 5, 5, 4, 4, 3, 3, 3…
## $ n_C        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_D        <int> 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3…
## $ n_E        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_F        <int> 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0…
## $ n_G        <int> 0, 0, 0, 0, 0, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1…
## $ n_H        <int> 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_I        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_K        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_L        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_M        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_N        <int> 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_P        <int> 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1…
## $ n_Q        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_R        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_S        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 2…
## $ n_T        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ n_V        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0…
## $ n_W        <int> 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0…
## $ n_Y        <int> 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 2…
## $ aliphatic  <int> 1, 1, 1, 1, 1, 9, 9, 5, 5, 5, 5, 4, 4, 6, 6, 4, 4, 3, 3, 5…
## $ aromatic   <int> 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 2, 2, 2…
## $ neutral    <int> 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 2…
## $ positive   <int> 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ negative   <int> 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3…
## $ invalid    <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
## $ file       <chr> "brigido_rafael_Nestor_resultado_nestor_VH_InitialRound_ra…
```

## Análise exploratória


```r
cdr_df %>% 
          group_by(experiment, expgroup, cycle) %>% 
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
  )
```

![](cdr_analysis_02_files/figure-html/data_exploration_1-1.png)<!-- -->

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
  )
```

![](cdr_analysis_02_files/figure-html/data_exploration_1-2.png)<!-- -->


```r
cdr_df %>% 
            group_by(experiment) %>% 
            ggplot() +
              geom_density(aes(cdrp, color = experiment, fill = experiment), alpha = 0.2) +
              labs(
                title = "Density plot of CDR3 Percentage of Prevalence",
                x = "CDR3 Percentage of Prevalence (CDRP)"
              ) +
              xlim(0.01, 0.33) +
              facet_grid(expgroup ~ .)
```

![](cdr_analysis_02_files/figure-html/data_exploration_2-1.png)<!-- -->

```r
cdr_df %>% 
            group_by(experiment) %>% 
            ggplot() +
              geom_density(
                aes(
                  x = quantity,
                  color = experiment,
                  fill = experiment),
                alpha = 0.2) +
              labs(
                title = "Density plot of Quantity Distribuition",
                x = "Quantity"
              ) + 
            facet_zoom(xlim = c(0, 50))
```

![](cdr_analysis_02_files/figure-html/data_exploration_2-2.png)<!-- -->


```r
ggplot(cdr_df) +
  geom_point(aes(cdrp, quantity, color = cycle)) +
  facet_grid(expgroup ~ .) +
  labs(
    title = "Correlation between CDR3 percentage of prevalence and Quantity",
    caption = "Note the interesting different trend between R0 and R4 in \"Renato Acid\" experiment.",
    x = "CDR3 Percentage of Prevalence (CDRP)",
    y = "Quantity"
  )
```

![](cdr_analysis_02_files/figure-html/data_exploration_3-1.png)<!-- -->

```r
ggplot(cdr_df) +
  geom_point(aes(fcp, fcq, color = expgroup)) +
  facet_grid(expgroup ~ .) + 
  labs(
    title = "Correlation between Fold-Changes",
    subtitle = "With outlier",
    caption = "Note how different the correlation is between experiment groups.",
    x = "Fold-Change by Percentage",
    y = "Fold-Change by Quantity"
  )
```

![](cdr_analysis_02_files/figure-html/data_exploration_3-2.png)<!-- -->

```r
ggplot(filter(cdr_df, fcp < 2E3)) +
  geom_point(aes(fcp, fcq, color = expgroup)) +
  facet_grid(expgroup ~ .) + 
  labs(
    title = "Correlation between Fold-Changes",
    subtitle = "Without outlier",
    caption = "Note how different the correlation is between experiment groups.",
    x = "Fold-Change by Percentage",
    y = "Fold-Change by Quantity"
  )
```

![](cdr_analysis_02_files/figure-html/data_exploration_3-3.png)<!-- -->

```r
# filtering the outlier.
cdr_df %>% filter(fcp > 1000) %>% select(cdr3:fcq)
```

```
## # A tibble: 1 x 8
## # Groups:   cdr3, expgroup [1]
##   cdr3         cycle expgroup    experiment     cdrp quantity   fcp   fcq
##   <chr>        <chr> <chr>       <chr>         <dbl>    <int> <dbl> <dbl>
## 1 AHIAAEYNWFDP R4    renato_acid renato_ac_R4 0.0938     1479 3002.  740.
```

```r
cdr_df %>% filter(cdr3 == "AHIAAEYNWFDP") %>% select(cdr3:fcq)
```

```
## # A tibble: 4 x 8
## # Groups:   cdr3, expgroup [2]
##   cdr3         cycle expgroup       experiment         cdrp quantity   fcp   fcq
##   <chr>        <chr> <chr>          <chr>             <dbl>    <int> <dbl> <dbl>
## 1 AHIAAEYNWFDP R0    renato_acid    renato_ac_R0  0.0000312        2    1     1 
## 2 AHIAAEYNWFDP R4    renato_acid    renato_ac_R4  0.0938        1479 3002.  740.
## 3 AHIAAEYNWFDP R0    renato_peptide renato_pep_R0 0.0000312        2    1     1 
## 4 AHIAAEYNWFDP R4    renato_peptide renato_pep_R4 0.00390        283  125.  142.
```


```r
# search for enriched CDR3s
for (i in c("fcp", "fcq")){
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(.[[i]]), .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          select(cdr3:fcq) -> tmp_plot_1
  
print(
  ggplot(tmp_plot_1) +
    geom_density(aes(.data[[i]], color = expgroup, fill = expgroup), alpha = 0.3)
)

print(
  ggplot(tmp_plot_1, aes(expgroup, .data[[i]])) +
    geom_violin(aes(fill = expgroup)) +
    geom_jitter(aes(shape = expgroup), alpha = 0.7)
)

}
```

![](cdr_analysis_02_files/figure-html/data_exploration_4-1.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_exploration_4-2.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_exploration_4-3.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_exploration_4-4.png)<!-- -->


### PCA for all groups


```r
# PCA for all groups
cdr_df %>%
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_all

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_all %<>% select(!c(which(apply(cdr_rich_df_all, 2, var)==0)))

cdr_rich_pca_all_groups <- cdr_rich_df_all %>% pca()

cdr_rich_df_all %>%
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
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

![](cdr_analysis_02_files/figure-html/pca_1-1.png)<!-- -->

```r
summary(cdr_rich_pca_all_groups)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5    PC6     PC7
## Standard deviation     2.7519 2.3224 2.2950 1.9197 1.60109 1.5073 1.32699
## Proportion of Variance 0.2047 0.1458 0.1424 0.0996 0.06928 0.0614 0.04759
## Cumulative Proportion  0.2047 0.3504 0.4928 0.5924 0.66167 0.7231 0.77066
##                           PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     1.2975 1.09811 1.03683 0.98567 0.96154 0.84203 0.74689
## Proportion of Variance 0.0455 0.03259 0.02905 0.02626 0.02499 0.01916 0.01508
## Cumulative Proportion  0.8162 0.84875 0.87781 0.90407 0.92905 0.94822 0.96329
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     0.64239 0.60184 0.48230 0.33095 0.29088 0.25643 0.25130
## Proportion of Variance 0.01115 0.00979 0.00629 0.00296 0.00229 0.00178 0.00171
## Cumulative Proportion  0.97445 0.98424 0.99052 0.99348 0.99577 0.99755 0.99925
##                           PC22    PC23    PC24      PC25      PC26      PC27
## Standard deviation     0.16190 0.03461 0.01458 7.802e-16 2.104e-16 1.998e-16
## Proportion of Variance 0.00071 0.00003 0.00001 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  0.99996 0.99999 1.00000 1.000e+00 1.000e+00 1.000e+00
##                             PC28      PC29      PC30
## Standard deviation     1.998e-16 1.998e-16 1.998e-16
## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00
```

### PCA Renato's groups


```r
# PCA for renato's groups
cdr_df %>%
          filter(str_detect(expgroup, "renato")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_renato_groups

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_renato_groups %<>% select(!c(which(apply(cdr_rich_df_renato_groups, 2, var)==0)))

cdr_rich_pca_renato_groups <- cdr_rich_df_renato_groups %>% 
                              pca()

cdr_rich_df_all %>%
          filter(str_detect(expgroup, "renato")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          .$expgroup -> cdr_rich_pca_color_renato_groups

ggplot_pca(
  cdr_rich_pca_renato_groups,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_renato_groups
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences",
    subtitle = "Only for Renato's experiments"
  )
```

![](cdr_analysis_02_files/figure-html/pca_2-1.png)<!-- -->

```r
summary(cdr_rich_pca_renato_groups)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
## Standard deviation     2.9959 2.5982 2.4532 2.1606 1.74641 1.66164 1.40577
## Proportion of Variance 0.2426 0.1825 0.1627 0.1262 0.08243 0.07462 0.05341
## Cumulative Proportion  0.2426 0.4250 0.5877 0.7139 0.79629 0.87091 0.92433
##                            PC8     PC9    PC10    PC11    PC12    PC13    PC14
## Standard deviation     1.04492 0.92190 0.79227 0.36597 0.29748 0.08359 0.03343
## Proportion of Variance 0.02951 0.02297 0.01696 0.00362 0.00239 0.00019 0.00003
## Cumulative Proportion  0.95383 0.97680 0.99377 0.99739 0.99978 0.99997 1.00000
##                             PC15      PC16      PC17      PC18      PC19
## Standard deviation     1.188e-15 3.584e-16 2.169e-16 1.521e-16 8.159e-17
## Proportion of Variance 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
## Cumulative Proportion  1.000e+00 1.000e+00 1.000e+00 1.000e+00 1.000e+00
##                             PC20
## Standard deviation     6.349e-17
## Proportion of Variance 0.000e+00
## Cumulative Proportion  1.000e+00
```

### PCA Renato Acid


```r
# PCA for renato acid
cdr_df %>%
          filter(str_detect(expgroup, "renato_acid")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_renato_acid

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_renato_acid %<>% select(!c(which(apply(cdr_rich_df_renato_acid, 2, var)==0)))

cdr_rich_pca_renato_acid <- cdr_rich_df_renato_acid %>%
                              pca()

cdr_df %>%
          filter(str_detect(expgroup, "renato_acid")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          .$expgroup -> cdr_rich_pca_color_renato_acid

ggplot_pca(
  cdr_rich_pca_renato_acid,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_renato_acid
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences",
    subtitle = "Only for Renato Acid experiment"
  )
```

![](cdr_analysis_02_files/figure-html/pca_3-1.png)<!-- -->

```r
summary(cdr_rich_pca_renato_acid)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
## Standard deviation     3.0728 2.7559 2.4962 2.2001 1.65639 1.39806 1.08458
## Proportion of Variance 0.2698 0.2170 0.1780 0.1383 0.07839 0.05584 0.03361
## Cumulative Proportion  0.2698 0.4868 0.6648 0.8031 0.88148 0.93732 0.97093
##                            PC8     PC9      PC10
## Standard deviation     0.90124 0.45291 1.149e-15
## Proportion of Variance 0.02321 0.00586 0.000e+00
## Cumulative Proportion  0.99414 1.00000 1.000e+00
```

### PCA Renato Peptide


```r
# PCA for renato peptide
cdr_df %>%
          filter(str_detect(expgroup, "renato_pep")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_renato_pep

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_renato_pep %<>% select(!c(which(apply(cdr_rich_df_renato_pep, 2, var)==0)))

cdr_rich_pca_renato_pep <- cdr_rich_df_renato_pep %>%
                              pca()

cdr_df %>%
          filter(str_detect(expgroup, "renato_pep")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          .$expgroup -> cdr_rich_pca_color_renato_pep

ggplot_pca(
  cdr_rich_pca_renato_pep,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_renato_pep
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences",
    subtitle = "Only for Renato Peptide experiment"
  )
```

![](cdr_analysis_02_files/figure-html/pca_4-1.png)<!-- -->

```r
summary(cdr_rich_pca_renato_pep)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
## Standard deviation     3.0927 2.6323 2.5744 2.1816 1.69511 1.46515 1.22007
## Proportion of Variance 0.2657 0.1925 0.1841 0.1322 0.07982 0.05963 0.04135
## Cumulative Proportion  0.2657 0.4582 0.6423 0.7745 0.85428 0.91391 0.95526
##                            PC8     PC9      PC10
## Standard deviation     0.98524 0.79992 1.071e-15
## Proportion of Variance 0.02696 0.01777 0.000e+00
## Cumulative Proportion  0.98223 1.00000 1.000e+00
```

### PCA Rafael


```r
# PCA for rafael
cdr_df %>%
          filter(str_detect(expgroup, "rafael")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          ungroup() %>% 
          select(where(is.numeric), expgroup) -> cdr_rich_df_rafael

# remove collumns with variance equal to 0
# credit: https://stackoverflow.com/a/40317343
cdr_rich_df_rafael %<>% select(!c(which(apply(cdr_rich_df_rafael, 2, var)==0)))

cdr_rich_pca_rafael <- cdr_rich_df_rafael %>%
                              pca()

cdr_df %>%
          filter(str_detect(expgroup, "rafael")) %>% 
          group_by(expgroup) %>% 
          arrange(-fcp, .by_group = TRUE) %>% 
          slice_head(n = 10) %>% 
          .$expgroup -> cdr_rich_pca_color_rafael

ggplot_pca(
  cdr_rich_pca_rafael,
  arrows_size = 0.2,
  labels = cdr_rich_pca_color_rafael
  ) +
  labs(
    title = "PCA for enriched CDR3 sequences",
    subtitle = "Only for Rafael's experiment"
  )
```

![](cdr_analysis_02_files/figure-html/pca_5-1.png)<!-- -->

```r
summary(cdr_rich_pca_renato_pep)
```

```
## Importance of components:
##                           PC1    PC2    PC3    PC4     PC5     PC6     PC7
## Standard deviation     3.0927 2.6323 2.5744 2.1816 1.69511 1.46515 1.22007
## Proportion of Variance 0.2657 0.1925 0.1841 0.1322 0.07982 0.05963 0.04135
## Cumulative Proportion  0.2657 0.4582 0.6423 0.7745 0.85428 0.91391 0.95526
##                            PC8     PC9      PC10
## Standard deviation     0.98524 0.79992 1.071e-15
## Proportion of Variance 0.02696 0.01777 0.000e+00
## Cumulative Proportion  0.98223 1.00000 1.000e+00
```

##  Draft section


```r
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(fcp)) %>% 
          slice_head(n = 10) %>% 
          select(cdr3:fcq) -> tmp_plot_1
    
    
          print(ggplot(tmp_plot_1) +
            geom_density(aes(quantity))) +
            facet_grid(expgroup ~ .)
```

![](cdr_analysis_02_files/figure-html/drafts1-1.png)<!-- -->![](cdr_analysis_02_files/figure-html/drafts1-2.png)<!-- -->

```r
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(fcp)) %>% 
          slice_head(n = 10) -> tmp_plot_2

          print(ggplot(tmp_plot_2, aes(expgroup, log10(quantity))) +
            geom_violin(aes(fill = expgroup)) +
            geom_jitter(aes(shape = expgroup)))
```

![](cdr_analysis_02_files/figure-html/drafts1-3.png)<!-- -->

```r
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(fcq)) %>% 
          slice_head(n = 10) -> tmp_plot_3

          print(ggplot(tmp_plot_3, aes(expgroup, log10(quantity))) +
            geom_violin(aes(fill = expgroup)) +
            geom_jitter(aes(shape = expgroup)))
```

![](cdr_analysis_02_files/figure-html/drafts1-4.png)<!-- -->

```r
cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(fcq), .by_group=T) %>% 
          slice_head(n = 10) %>% 
          ggplot() +
            geom_point(aes(fcp, fcq)) +
            facet_grid(expgroup ~ .)
```

![](cdr_analysis_02_files/figure-html/drafts1-5.png)<!-- -->


```r
cdr_df %>% 
          filter(cycle == "R4", fcp == 1, fcq == 1) %>%
          select(cdr3, cdrp, fcp, quantity, fcq) %>%
          ggplot() +
            geom_density(aes(quantity)) +
            facet_grid(expgroup ~ .)
```

![](cdr_analysis_02_files/figure-html/drafts2-1.png)<!-- -->

```r
cdr_df %>% group_by(experiment) -> df_cdr_exp

df_cdr_exp %>% summarise(quantiles = quantile(quantity))
```

```
## # A tibble: 30 x 2
## # Groups:   experiment [6]
##    experiment quantiles
##    <chr>          <dbl>
##  1 rafael_R0          1
##  2 rafael_R0          1
##  3 rafael_R0          1
##  4 rafael_R0          2
##  5 rafael_R0       5294
##  6 rafael_R4          1
##  7 rafael_R4          1
##  8 rafael_R4          1
##  9 rafael_R4          2
## 10 rafael_R4       3036
## # … with 20 more rows
```

```r
df_cdr_exp %>% 
              filter(experiment == "renato_pep_R4") %>% 
              mutate(ham = if_else(quantity > 20, 1, 0)) %>% 
              select(cdr3, quantity, ham) %>% 
              mutate(per = (sum(ham)/length(ham)) * 100) %>% 
              group_by(ham) %>% sample_n(3)
```

```
## # A tibble: 6 x 5
## # Groups:   ham [2]
##   experiment    cdr3       quantity   ham   per
##   <chr>         <chr>         <int> <dbl> <dbl>
## 1 renato_pep_R4 DGVAVAGLDN        3     0  2.60
## 2 renato_pep_R4 DRDHRFDS          2     0  2.60
## 3 renato_pep_R4 GLYSSGRIDV        1     0  2.60
## 4 renato_pep_R4 GRWGSF           39     1  2.60
## 5 renato_pep_R4 DLGIPDDY         52     1  2.60
## 6 renato_pep_R4 PLTGLHY        3440     1  2.60
```

```r
cdr_df %>% 
          group_by(expgroup, cdr3, cycle) %>% 
          filter(cdr3 == "FIVESK") %>% 
          select(cdr3:quantity) %>% 
          arrange(expgroup)
```

```
## # A tibble: 6 x 6
## # Groups:   expgroup, cdr3, cycle [6]
##   cdr3   cycle expgroup       experiment       cdrp quantity
##   <chr>  <chr> <chr>          <chr>           <dbl>    <int>
## 1 FIVESK R0    rafael         rafael_R0     0.118       5294
## 2 FIVESK R4    rafael         rafael_R4     0.0549      3036
## 3 FIVESK R0    renato_acid    renato_ac_R0  0.104       6655
## 4 FIVESK R4    renato_acid    renato_ac_R4  0.0126       199
## 5 FIVESK R0    renato_peptide renato_pep_R0 0.104       6655
## 6 FIVESK R4    renato_peptide renato_pep_R4 0.00617      448
```

```r
cdr_df %>% 
          group_by(expgroup, cycle) %>% 
          select(cdr3:quantity) %>% 
          sample_n(1)
```

```
## # A tibble: 6 x 6
## # Groups:   expgroup, cycle [6]
##   cdr3             cycle expgroup       experiment         cdrp quantity
##   <chr>            <chr> <chr>          <chr>             <dbl>    <int>
## 1 FYSNNWNEVFCDY    R0    rafael         rafael_R0     0.0000223        1
## 2 EEELTGTGYYYYGMDV R4    rafael         rafael_R4     0.0000362        2
## 3 APYDILTGYSN      R0    renato_acid    renato_ac_R0  0.0000156        1
## 4 GTMYGLVPSDS      R4    renato_acid    renato_ac_R4  0.0000634        1
## 5 DPYDLWSGNSIVY    R0    renato_peptide renato_pep_R0 0.0000156        1
## 6 VGGRRALDY        R4    renato_peptide renato_pep_R4 0.0000275        2
```

```r
cdr_df %>%
          group_by(experiment) %>% 
          select(cdr3:quantity)
```

```
## # A tibble: 75,459 x 6
## # Groups:   experiment [6]
##    cdr3          cycle expgroup       experiment         cdrp quantity
##    <chr>         <chr> <chr>          <chr>             <dbl>    <int>
##  1 A             R0    rafael         rafael_R0     0.0000445        2
##  2 A             R0    renato_acid    renato_ac_R0  0.0000937        6
##  3 A             R4    renato_acid    renato_ac_R4  0.000127         2
##  4 A             R0    renato_peptide renato_pep_R0 0.0000937        6
##  5 A             R4    renato_peptide renato_pep_R4 0.0000413        3
##  6 AAAAAGGGNWFDP R0    renato_acid    renato_ac_R0  0.0000156        1
##  7 AAAAAGGGNWFDP R0    renato_peptide renato_pep_R0 0.0000156        1
##  8 AAAAGHFDY     R0    renato_acid    renato_ac_R0  0.0000312        2
##  9 AAAAGHFDY     R0    renato_peptide renato_pep_R0 0.0000312        2
## 10 AAAAGQFDY     R0    renato_acid    renato_ac_R0  0.0000156        1
## # … with 75,449 more rows
```

```r
# analysing cdrp

cdr_df %>% 
          select(cdr3:quantity) %>% 
          # filter(quantity > 400) %>% 
          group_by(cdr3, expgroup) %>% 
          # filter(n()>1) %>% 
          arrange(cycle, .by_group = TRUE) %>% 
          mutate(
            fcp = cdrp / lag(cdrp, default = first(cdrp)),
            fcq = quantity / lag(quantity, default = first(quantity))) %>% 
          select(cdr3:quantity, fcp, fcq, everything()) -> b

b %>% 
      group_by(cdr3, expgroup) %>% 
      filter(n()>1) %>% 
      group_by(cycle, .add=T) %>% 
      summarise(n = n(), s = sum(fcp), q = sum(quantity), c = sum(cdrp)) %>% print()
```

```
## # A tibble: 2,430 x 7
## # Groups:   cdr3, expgroup [1,215]
##    cdr3       expgroup       cycle     n     s     q         c
##    <chr>      <chr>          <chr> <int> <dbl> <int>     <dbl>
##  1 A          renato_acid    R0        1 1         6 0.0000937
##  2 A          renato_acid    R4        1 1.35      2 0.000127 
##  3 A          renato_peptide R0        1 1         6 0.0000937
##  4 A          renato_peptide R4        1 0.441     3 0.0000413
##  5 AAGGRPFDY  renato_peptide R0        1 1         1 0.0000156
##  6 AAGGRPFDY  renato_peptide R4        1 0.882     1 0.0000138
##  7 AARQPDY    renato_acid    R0        1 1         5 0.0000781
##  8 AARQPDY    renato_acid    R4        1 0.812     1 0.0000634
##  9 ADGYNNEVDY renato_acid    R0        1 1        13 0.000203 
## 10 ADGYNNEVDY renato_acid    R4        1 0.312     1 0.0000634
## # … with 2,420 more rows
```

```r
b %>% filter(fcp > 200) %>% arrange(-fcp)
```

```
## # A tibble: 18 x 8
## # Groups:   cdr3, expgroup [18]
##    cdr3            cycle expgroup      experiment      cdrp quantity   fcp   fcq
##    <chr>           <chr> <chr>         <chr>          <dbl>    <int> <dbl> <dbl>
##  1 AHIAAEYNWFDP    R4    renato_acid   renato_ac_R4 0.0938      1479 3002. 740. 
##  2 EPSS            R4    renato_acid   renato_ac_R4 0.0138       218  885. 218  
##  3 VGGGRALDY       R4    renato_acid   renato_ac_R4 0.0401       633  642. 158. 
##  4 VGGGRALDY       R4    renato_pepti… renato_pep_… 0.0354      2574  567. 644. 
##  5 GRWGSY          R4    renato_pepti… renato_pep_… 0.292      21192  425. 482. 
##  6 EPSS            R4    renato_pepti… renato_pep_… 0.00642      466  411. 466  
##  7 DGVAVAGLDY      R4    renato_pepti… renato_pep_… 0.252      18304  329. 374. 
##  8 GGIVGAPDY       R4    renato_pepti… renato_pep_… 0.0252      1829  323. 366. 
##  9 GRGYSGYDRPFDY   R4    renato_acid   renato_ac_R4 0.00488       77  313.  77  
## 10 GRGYSGYDRPFDY   R4    renato_pepti… renato_pep_… 0.00450      327  288. 327  
## 11 GRWGSY          R4    renato_acid   renato_ac_R4 0.196       3094  285.  70.3
## 12 DAHRKGYYGMDV    R4    renato_pepti… renato_pep_… 0.00836      607  268. 304. 
## 13 GGIVGAPDY       R4    renato_acid   renato_ac_R4 0.0209       329  267.  65.8
## 14 GGVNWNDQ        R4    rafael        rafael_R4    0.0112       620  252. 310  
## 15 DAHRKGYYGMDV    R4    renato_acid   renato_ac_R4 0.00697      110  223.  55  
## 16 GRWGGY          R4    renato_pepti… renato_pep_… 0.00343      249  220. 249  
## 17 PQQWLAWTGAEGYF… R4    renato_pepti… renato_pep_… 0.0273      1985  219. 248. 
## 18 APAGREFDY       R4    rafael        rafael_R4    0.00456      252  205. 252
```

```r
b %>%
      filter(cycle == "R4", fcp == 1, fcq == 1) %>%
      select(cdr3, cdrp, fcp, quantity, fcq) %>% 
      ggplot() + 
        geom_density(aes(quantity))
```

![](cdr_analysis_02_files/figure-html/drafts2-2.png)<!-- -->

```r
ggplot(filter(b, cycle == "R4")) +
  geom_point(aes(fcp, fcq))
```

![](cdr_analysis_02_files/figure-html/drafts2-3.png)<!-- -->

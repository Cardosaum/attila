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
          arrange(desc(i)) %>% 
          slice_head(n = 10) %>% 
          select(cdr3:fcq) -> tmp_plot_1
    
    
          print(ggplot(tmp_plot_1) +
            geom_density(aes(quantity)))

cdr_df %>% 
          group_by(expgroup) %>% 
          arrange(desc(i)) %>% 
          slice_head(n = 10) -> tmp_plot_2

          print(ggplot(tmp_plot_2, aes(expgroup, log10(quantity))) +
            geom_violin(aes(fill = expgroup)) +
            geom_jitter(aes(shape = expgroup)))
}
```

![](cdr_analysis_02_files/figure-html/data_processing_2-1.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_processing_2-2.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_processing_2-3.png)<!-- -->![](cdr_analysis_02_files/figure-html/data_processing_2-4.png)<!-- -->

##  Draft section


```r
cdr_df %>% 
          filter(cycle == "R4", fcp == 1, fcq == 1) %>%
          select(cdr3, cdrp, fcp, quantity, fcq) %>%
          ggplot() +
            geom_density(aes(quantity)) +
            facet_grid(expgroup ~ .)
```

![](cdr_analysis_02_files/figure-html/drafts-1.png)<!-- -->

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
##   experiment    cdr3             quantity   ham   per
##   <chr>         <chr>               <int> <dbl> <dbl>
## 1 renato_pep_R4 RRDNSGNTPFDD            1     0  2.60
## 2 renato_pep_R4 DYGGPRGARYYYGMDV        5     0  2.60
## 3 renato_pep_R4 GRWRSF                  2     0  2.60
## 4 renato_pep_R4 PLAGLHY                22     1  2.60
## 5 renato_pep_R4 EMWGPEY                65     1  2.60
## 6 renato_pep_R4 GRGYSGYDRPFDY         327     1  2.60
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
## 1 EYSGDPRRIDY      R0    rafael         rafael_R0     0.0000223        1
## 2 RISMMGSQH        R4    rafael         rafael_R4     0.0000181        1
## 3 DQGYYYDSSDY      R0    renato_acid    renato_ac_R0  0.0000156        1
## 4 GGSSSPGLVGFHSMDV R4    renato_acid    renato_ac_R4  0.000254         4
## 5 DRGMTTVTTVDY     R0    renato_peptide renato_pep_R0 0.0000468        3
## 6 GGWGSS           R4    renato_peptide renato_pep_R4 0.0000275        2
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

![](cdr_analysis_02_files/figure-html/drafts-2.png)<!-- -->

```r
ggplot(filter(b, cycle == "R4")) +
  geom_point(aes(fcp, fcq))
```

![](cdr_analysis_02_files/figure-html/drafts-3.png)<!-- -->
<!-- ### Isolando apenas as sequências CDR3 enriquecidas -->

<!-- Como é possível perceber pelos dados acima mostrados, temos muitas reads no dataframe. -->
<!-- Entretanto, nosso interesse por agora é nas sequências que foram enriquecidas após várias etapas de seleção. -->
<!-- Para isso, nós precisaremos criar um subset do dataframe inicial, contendo apenas CDR3s que apresentam alto percentual de predominância em seu respectivo arquivo de leitura. -->

<!-- Vou mostrar um exemplo do que quero dizer: -->

<!-- ```{r data_processing_1} -->
<!-- cdr %>% -->
<!--       select(cdr3, type, cdrp, quantity, file) %>%  -->
<!--       head() -> exemplo_unico_cdr -->

<!-- knitr::kable(exemplo_unico_cdr) -->
<!-- ``` -->

<!-- Como é possível observar, nas duas primeras linhas temos uma mesma sequência, que apresenta um percentual de 100% predôminancia em seu respectivo arquivo de leiura. -->
<!-- (coluna `cdrp` - cdr percentage, variando de 0 a 1). -->
<!-- Porém, observamos também que a mesma sequência aparece nesse arquivo somente uma vez. -->
<!-- Ou seja, esses dois primeiros arquivos contém só uma leitura, e, portanto, seu percentual de predominância será de 100%. -->
<!-- Isso, por outro lado, não reflete enriquecimento de CDR3, e, portanto, nós precisamos remover esses casos. -->


<!-- Pensando em como fazer a seleção dessas sequências enriquecidas, fiz algumas análises: -->

<!-- ```{r data_exploration_1} -->
<!-- ggplot(filter(cdr, type == "final")) + -->
<!--   geom_histogram(aes(quantity)) -->

<!-- ggplot(filter(cdr, type == "final")) + -->
<!--   geom_histogram(aes(quantity)) + -->
<!--   xlim(0, 300) -->

<!-- cdr %>% -->
<!--         filter(type == "final") %>%  -->
<!--         mutate(level = case_when( -->
<!--                   quantity <= 300 ~ "quantity <= 300", -->
<!--                   TRUE ~ "quantity > 300")) %>% -->
<!--         group_by(level) %>%  -->
<!--         summarise("Number of CDR3 sequences" = n()) -> cdr_quantity_comparison_1 -->

<!-- knitr::kable(cdr_quantity_comparison_1) -->
<!-- ``` -->


<!-- ```{r data_exploration_2} -->
<!-- ggplot(filter(cdr, type == "final")) + -->
<!--   geom_histogram(aes(cdrp)) -->

<!-- ggplot(filter(cdr, type == "final")) + -->
<!--   geom_histogram(aes(cdrp)) + -->
<!--   xlim(0.5, 1) -->

<!-- cdr %>% -->
<!--         filter(type == "final") %>%  -->
<!--         mutate(level = case_when( -->
<!--                   cdrp <= 0.3 ~ "cdrp <= 0.3", -->
<!--                   TRUE ~ "cdrp > 0.3")) %>% -->
<!--         group_by(level) %>%  -->
<!--         summarise("Percentage" = n()) -> cdr_cdrp_comparison_1 -->

<!-- knitr::kable(cdr_cdrp_comparison_1, caption = "Percentage of prevalence of CDR3 sequence") -->


<!-- cdr %>% -->
<!--         filter(type == "final") %>%  -->
<!--         mutate(level = case_when( -->
<!--                   cdrp < 0.5 ~ "cdrp < 0.5", -->
<!--                   TRUE ~ "cdrp > 0.5")) %>% -->
<!--         group_by(level) %>%  -->
<!--         summarise("Percentage" = n()) -> cdr_cdrp_comparison_2 -->

<!-- knitr::kable(cdr_cdrp_comparison_2, caption = "Percentage of prevalence of CDR3 sequence") -->
<!-- ``` -->

<!-- Como é possível notar, temos 23 sequências de CDR3 que apresentam prevalência maior que 30% em arquivos de leitura individual, e 22 se considerarmos 50% de prevalência. -->

<!-- Para termos noção do que isso significa, vejamos o seguinte: -->

<!-- ```{r data_exploration_3} -->
<!-- cdr$file %>% unique() %>% length() -> total_arquivos_leitura -->

<!-- filter(cdr, type == "final")$file %>% unique() %>% length() -> total_arquivos_leitura_final_read -->

<!-- tibble( -->
<!--   "Arquivo de leitura" = c("Todos (Inicial + Final)", "Apenas Final", "Final com CDR3 prevalência >= 50%"), -->
<!--   "Quantidade de Arquivos" = c(total_arquivos_leitura, total_arquivos_leitura_final_read, cdr_cdrp_comparison_2$Percentage[2]) -->
<!-- ) %>% knitr::kable() -->
<!-- ``` -->

<!-- E, para mostrar todos os arquivos com prevalência maior que 50%: -->

<!-- ```{r data_exploration_4} -->
<!-- cdr %>%  -->
<!--         filter(type == "final" & cdrp >= 0.5) %>%  -->
<!--         select(cdr3, cdrp, quantity, file) %>%  -->
<!--         knitr::kable() -->
<!-- ``` -->


<!-- Portanto, eu resolvi salvar esse dataframe como aquele contendo as sequências enriquecidas. -->

<!-- ```{r data_processing_2} -->
<!-- cdr_rich <- cdr %>% filter(type == "final" & cdrp >= 0.5)  -->
<!-- ``` -->

<!-- **Todo o código feito a partir daqui é um rascunho** -->

<!-- Peço perdão pela bagunça nos próximos blocos. -->
<!-- Eu escrevi isso para me ajudar a entender os dados, sem a intenção de apresentar isso para ninguém. -->

<!-- ## Análise Exploratória -->

<!-- ```{r eda1} -->
<!-- cdr %>% -->
<!--         ungroup() %>%  -->
<!--         arrange(-cdrp, type, file) %>%  -->
<!--         filter(quantity > 1) %>%  -->
<!--         filter(type == "final") -> cdr_final -->

<!-- cdr_final %>%  -->
<!--                 filter(quantity > 1) %>%  -->
<!--                 group_by(file) %>%  -->
<!--                 slice_head(n = 1) -> cdr_enriched -->

<!-- library(GGally) -->
<!-- cdr_enriched %<>% -->
<!--                 select(cdr3:SSF_Sheet, aromatic:file) -->

<!-- cdr_enriched %>% -->
<!--                 ungroup() %>% -->
<!--                 select(-file) %>% -->
<!--                 ggpairs(aes(alpha = 0.4)) -->
<!-- ``` -->



<!-- ```{r eda2} -->
<!-- cdr_final %>%  -->
<!--              ungroup() %>%  -->
<!--              select(!c(cdr3, type, file, invalid)) -> cdr_final_pca -->

<!-- pca_result <- prcomp(cdr_final_pca, center = T, scale. = T) -->
<!-- summary(pca_result) -->

<!-- plot(pca_result$x[,1], pca_result$x[,2]) -->
<!-- cdr_final_pca -->

<!-- cdr_final %>%  -->
<!--               group_by(file) %>%  -->
<!--               arrange(-cdrp) %>%  -->
<!--               slice_head(n = 1) %>%  -->
<!--               ungroup() %>%  -->
<!--               select(!c(cdr3, type, file, invalid)) %>%  -->
<!--               arrange(-cdrp) -> a -->

<!-- # in this line we remove all collumns that have variance equal to 0 -->
<!-- # Doing this, we can apply a pca to the dataframe without erros -->
<!-- # credit goes to: https://stackoverflow.com/a/40317343 -->
<!-- a <- select(a, !c(which(apply(a, 2, var)==0))) -->
<!-- pca_a <- prcomp(a, center = T, scale. = T) -->
<!-- summary(pca_a) -->
<!-- plot(pca_a$x[,1], pca_a$x[,2]) -->

<!-- ggplot(as_tibble(pca_a$x)) + -->
<!--   geom_point(aes(PC1, PC2)) -->

<!-- pca_a$x -->
<!-- str(pca_a) -->

<!-- pca_cdr_result <- cdr %>% -->
<!--                         select(!c(cdr3, type, file, invalid)) %>% -->
<!--                         prcomp(center = T, scale. = T) -->
<!-- summary(pca_cdr_result) -->
<!-- ``` -->

<!-- ```{r eda3} -->
<!-- cdr_final %>%  -->
<!--               group_by(file) %>%  -->
<!--               arrange(-cdrp) %>%  -->
<!--               slice_head(n = 10) %>%  -->
<!--               ungroup() %>%  -->
<!--               select(!c(cdr3, type, file, invalid)) %>%  -->
<!--               arrange(-cdrp) -> b -->

<!-- b <- select(b, !c(which(apply(b, 2, var)==0))) -->
<!-- b -->
<!-- pca_b <- prcomp(b, center = T, scale. = T) -->
<!-- summary(pca_b) -->
<!-- plot(pca_b$x[,1], pca_b$x[,2]) -->

<!-- ggplot(as_tibble(pca_b$x)) + -->
<!--   geom_point(aes(PC1, PC2)) -->

<!-- summary(pca_b) -->
<!-- ``` -->

<!-- ```{r eda4} -->
<!-- summary(cdr$quantity) -->
<!-- cdr %>%  -->
<!--         filter(quantity >= 100) -> a -->
<!-- a -->

<!-- ggplot(a) + -->
<!--   geom_density(aes(quantity)) -->

<!-- ggplot(cdr) + -->
<!--   geom_bar(aes(quantity)) + -->
<!--   xlim(0, 30) -->

<!-- ggplot(cdr) + -->
<!--   geom_histogram(aes(quantity)) + -->
<!--   xlim(0, 300) -->

<!-- ggplot(cdr) + -->
<!--   geom_density(aes(quantity), fill = "lightblue") + -->
<!--   xlim(0, 300) -->

<!-- quantile(cdr$quantity) -->

<!-- dim(cdr) -->
<!-- cdr %>% filter(quantity >= 1E3) %>% dim() -->
<!-- cdr %>% filter(quantity >= 1E4) %>% dim() -->
<!-- cdr %>% filter(quantity >= 1E5) %>% dim() -->

<!-- cdr %>% filter(quantity >= 1E3) -> b -->
<!-- b %>% group_by(type) %>% summarise(total = n()) -->
<!-- b %>% group_by(type) %>% summarise(quantile = quantile(cdrp)) -> b_quantiles -->
<!-- b_quantiles <- add_column(b_quantiles, quantiles = rep(attr(quantile(b$quantity), "names"), 2)) -->
<!-- knitr::kable(b_quantiles) -->

<!-- ggplot(b) + -->
<!--   geom_density(aes(quantity)) -->

<!-- b %>%  -->
<!--       group_by(cdr3, type) %>%  -->
<!--       arrange(-cdrp) -->

<!-- b %>%  -->
<!--       group_by(cdr3, type) %>%  -->
<!--       select(cdr3, type, cdrp, quantity) %>%  -->
<!--       arrange(-cdrp, -quantity) %>%  -->
<!--       slice_head(n = 1) %>%  -->
<!--       arrange(-cdrp, -quantity)  -->

<!-- b %>%  -->
<!--       group_by(type, cdr3) %>%  -->
<!--       summarise(total = n()) %>%  -->
<!--       arrange(-total) -->
<!-- b %>%  -->
<!--       group_by(cdr3, type) %>%  -->
<!--       select(cdr3, type, cdrp, quantity) %>%  -->
<!--       arrange(-cdrp, -quantity) -> c -->

<!-- c %>% filter(type == "initial") %>% slice_head(n = 1) -->

<!-- ggplot(c) + -->
<!--   geom_density(aes(cdrp, color = type), alpha = .4) -->
<!-- ``` -->

<!-- ```{r eda5} -->
<!-- b %>% -->
<!--     group_by(cdr3, type) %>%  -->
<!--     summarise( -->
<!--       quantity = sum(quantity), -->
<!--       reads    = n()) %>%  -->
<!--     arrange(-quantity, -reads) -> d -->

<!-- d -->

<!-- d %>% group_by(type) %>% summarise(n = n()) -->

<!-- ggplot(d) + -->
<!--   geom_density(aes(reads, color = type)) -->

<!-- ggplot(d) + -->
<!--   geom_boxplot(aes(type, log10(quantity), fill = type)) + -->
<!--   geom_jitter(aes(type, log10(quantity), fill = type)) -->

<!-- d -->
<!-- ``` -->


<!-- # Resultados -->

<!-- # Conclusão -->

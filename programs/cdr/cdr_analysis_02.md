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
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: console
---



# Introdução

Essa análise será baseada na primeira iteração, porém mudando os dados analisados.
Dessa vez, usarei os seguintes dados:


```r
library(tidyverse, quietly = TRUE)
library(magrittr, quietly = TRUE)

data_path            <- "./data"
data_path_raw        <- file.path(data_path, "raw")
data_path_csv        <- file.path(data_path, "parsed")
data_files_selected  <- readLines(file("cdr_analysis_02_selected_files.txt"))


data_files_and_size <- sapply(data_files_selected, file.size)
files_to_include_in_dataframe <- tibble("Files" = names(data_files_and_size), "Size (in MB)" = (data_files_and_size/1E6))
knitr::kable(files_to_include_in_dataframe, caption = "List of files to be analysed")
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
source("load_files.R")

cdr_df <- load_cdr(names(data_files_and_size))
```

```
## [1] "data/parsed/brigido_rafael_Nestor_resultado_nestor_VH_FinalRound_rafaCD20_Vh_R4_R1aafreq.csv"
## [1] "data/parsed/brigido_rafael_Nestor_resultado_nestor_VH_InitialRound_rafaCD20_Vh_R0_R1aafreq.csv"
## [1] "data/parsed/brigido_renato_zika_ago18_ac_phage_zika_acid_VH_FinalRound_R4ac_VH_R1aafreq.csv"
## [1] "data/parsed/brigido_renato_zika_ago18_ac_phage_zika_acid_VH_InitialRound_R0_VH_R1aafreq.csv"
## [1] "data/parsed/brigido_renato_zika_ago18_phage_zika_VH_FinalRound_R4pep_VH_R1aafreq.csv"
## [1] "data/parsed/brigido_renato_zika_ago18_phage_zika_VH_InitialRound_R0_VH_R1aafreq.csv"
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
                                  str_detect(file, "rafael") ~ "rafael_CD20",
                                  str_detect(file, "renato.*acid") ~ "renato_acid",
                                  TRUE ~ "renato_peptide")) %>% 
          select(cdr3, cycle, expgroup, experiment, everything())

names(cdr_df)
```

```
##  [1] "cdr3"       "cycle"      "expgroup"   "experiment" "cdrp"      
##  [6] "quantity"   "length"     "MW"         "AV"         "IP"        
## [11] "flex"       "gravy"      "SSF_Helix"  "SSF_Turn"   "SSF_Sheet" 
## [16] "n_A"        "n_C"        "n_D"        "n_E"        "n_F"       
## [21] "n_G"        "n_H"        "n_I"        "n_K"        "n_L"       
## [26] "n_M"        "n_N"        "n_P"        "n_Q"        "n_R"       
## [31] "n_S"        "n_T"        "n_V"        "n_W"        "n_Y"       
## [36] "aliphatic"  "aromatic"   "neutral"    "positive"   "negative"  
## [41] "invalid"    "file"
```

```r
dim(cdr_df)
```

```
## [1] 75459    42
```


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
library(ggsci)

# search for enriched CDR3s

# analysing cdrp
cdr_df %>% 
            group_by(experiment) %>% 
            ggplot() +
              geom_density(aes(cdrp, color = experiment, fill = experiment), alpha = 0.2) +
              xlim(0.01, 0.3) +
              scale_color_d3() + scale_fill_d3()
```

![](cdr_analysis_02_files/figure-html/data_exploration_2-1.png)<!-- -->

```r
# analysing fold change
```

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

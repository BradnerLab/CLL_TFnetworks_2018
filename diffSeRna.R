## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = F,
                      warning = F,
                      message = F,
                      tidy = F)

## ----message=FALSE, warning=FALSE----------------------------------------
library(magrittr)
library(tidyverse)
library(readr)
library(readxl)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(iheatmapr)

set.seed(12345)
genome <- BSgenome.Hsapiens.UCSC.hg19
final_metadata <- read_tsv("results/primary_sample_metadata.txt")

## ----results='hide'------------------------------------------------------
read_tsv("tables/CLL_WGS_RNA-Seq.cufflinks.txt", skip = 2) %>%
  set_tidy_names(syntactic = T) %>%
  select(-Name) %>%
  gather(sample, expression, -Description) %>%
  group_by(Description, sample) %>%
  summarize(expression = sum(expression)) %>% # sum all transcripts for each gene
  ungroup %>%
  mutate(Condition = if_else(str_detect(sample, "CLL"), "pCLL", "CD19")) %>%
  group_by(Description, Condition) %>%
  summarize(expression = list(expression)) %>%
  ungroup() %>%
  spread(Condition, expression) %>%
  rowwise() %>%
  mutate(lfc = log2(mean(unlist(pCLL))+1) - log2(mean(unlist(CD19))+1)) %>%
  ungroup() -> res

## ----results='hide'------------------------------------------------------
cn <- c("name","chrom","start","stop",
        "n_constit","width","bam","rank1",
        "super1","overlap","prox","closest",
        "rank","isSuper")

gains <- read_tsv("results/gained_ENHANCER_TO_GENE.txt", 
                  col_names = cn, skip = 1) %>% 
  select(name,prox, closest) %>%
  mutate(prox = str_split(prox,",")) %>% 
  .$closest %>% 
  unlist() %>% 
  table() %>% 
  as.tibble() %>% #changed to closest
  mutate(type = "gain") %>% 
  set_colnames(c("Description","ct","type"))

losses <- read_tsv("results/lost_ENHANCER_TO_GENE.txt", 
                  col_names = cn, skip = 1) %>% 
  select(name,prox, closest) %>%
  mutate(prox = str_split(prox,",")) %>% 
  .$closest %>% 
  unlist %>% 
  table %>% 
  as.tibble() %>% 
  mutate(type = "loss") %>% 
  set_colnames(c("Description","ct","type"))

non_diff_ses <- read_tsv("results/searchspace_ENHANCER_TO_GENE.txt", 
                  col_names = cn, skip = 1) %>% 
  select(name,prox, closest) %>%
  mutate(prox = str_split(prox,",")) %>% 
  .$closest %>% 
  unlist %>% 
  table %>% 
  as.tibble() %>% #
  mutate(type = "se.non.diff") %>% 
  set_colnames(c("Description","ct","type")) %>%
  filter(!(Description %in% c(gains$Description, losses$Description)))

gene_to_diff_se <- bind_rows(gains,losses, non_diff_ses)

## ------------------------------------------------------------------------
expr_to_se_gains <- left_join(res, gene_to_diff_se, 
                              by = "Description") %>%
  mutate(type = if_else(is.na(type),"non.se",type))

expr_to_se_gains %>% 
  mutate(CD19 = map_chr(CD19, paste, collapse = ","), 
         pCLL = map_chr(pCLL, paste, collapse = ",")) %>% 
  write_tsv("results/gene_expr_lfc_to_supers.txt")

ggplot(filter(expr_to_se_gains, type != "non.se"), aes(x = type, y = lfc)) + 
  geom_jitter(width = 0.2, aes(color = type)) +
  geom_boxplot(outlier.shape = NA, fill = NA) +
  ylim(-2,2) +
  theme_classic() +
  scale_fill_manual(values = c(gain = "red", loss = "blue", se.non.diff = "gray")) +
  scale_color_manual(values = c(gain = "red", loss = "blue", se.non.diff = "gray"))
  

## ------------------------------------------------------------------------
wilcox.test(subset(expr_to_se_gains, type =="gain")[["lfc"]],
            subset(expr_to_se_gains, type =="se.non.diff")[["lfc"]])

wilcox.test(subset(expr_to_se_gains, type =="loss")[["lfc"]],
            subset(expr_to_se_gains, type =="se.non.diff")[["lfc"]])

## ----results='hide'------------------------------------------------------
read_tsv("tables/CLL_WGS_RNA-Seq.cufflinks.txt", skip = 2) %>% # inefficient
  set_tidy_names(syntactic = T) %>%
  group_by(Description) %>%
  mutate_at(vars(CLL.CW114.Tumor.SM.4CQKE:No9_Normal), sum) %>%
  ungroup() %>%
  select(-Name) %>% 
  left_join(., gene_to_diff_se) %>% 
  distinct() -> heat1

unambiguous_genes <- heat1 %>%
  group_by(Description, type) %>% 
  summarise(ct = n()) %>% ungroup() %>%
  group_by(Description) %>%
  filter(n() == 1 & !is.na(type)) %>% 
  .$Description # exclude non-diff-supers + +/- SE genes
  
heat2 <- filter(heat1, Description %in% unambiguous_genes)
  
se_gene_grps <- heat2$type
sample_grps <-  ifelse(str_detect(colnames(heat2)[2:39],
                                  "Normal"),"CD19","pCLL")

row_z_mat <- heat2 %>% select(-type, -ct) %>%
  as.data.frame() %>% 
  column_to_rownames("Description") %>% 
  as.matrix() %>% 
  t() %>% 
  scale %>% 
  t()

## ------------------------------------------------------------------------
se_expr_heat <- main_heatmap(row_z_mat, zmin = -3, zmax = 3, 
             name = "Row Z-score Expression") %>% 
  add_row_clusters(se_gene_grps ) %>% 
  add_col_clusters(sample_grps) # %>%
  # add_row_labels(tickvals = 1:length(is_hilite), 
  #               ticktext = is_hilite, side = "right", 
  #               font = list(size = 6))

se_expr_heat

## ----save-se-expr-heat, echo=F, include=F--------------------------------
save_iheatmap(se_expr_heat, "results/diff-supers-expression.heat.pdf")


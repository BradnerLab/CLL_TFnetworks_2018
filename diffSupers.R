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
library(Rsubread)
library(EDASeq)
library(DESeq2)
library(apeglm)
library(BSgenome.Hsapiens.UCSC.hg19)
library(iheatmapr)

set.seed(12345)
genome <- BSgenome.Hsapiens.UCSC.hg19
final_metadata <- read_tsv("results/primary_sample_metadata.txt")

## ----results='hide'------------------------------------------------------
rose_se_table_colnames <- c("REGION_ID","CHROM","START", "STOP", "NUM_LOCI", "CONSTITUENT_SIZE", 
                            "chip.signal", "input.signal", "enhancerRank", "isSuper")

metadata_for_diff_se <- final_metadata %>% 
  filter(INCL_K27)

supers <- metadata_for_diff_se$K27.SUPERS %>%
  lapply(read_tsv, skip = 6, col_names = rose_se_table_colnames) %>%
  bind_rows() %>%
  select(CHROM, START, STOP) %>%
  mutate(START = START + 1) %>% # GRanges are 1-indexed
  GRanges()

## ----quantify-all-primary-se, cache=T, message=F, warning=F, results='hide'----
quantifyReads <- function(gr, bamlist, nthreads = 8, paired_end = T) {
    GenomicRanges::strand(gr) <- "*"
    saf <- data.frame(GeneID = as.character(gr), 
                      Chr = GenomicRanges::seqnames(gr),
                      Start = GenomicRanges::start(gr),
                      End = GenomicRanges::end(gr), 
                      Strand = GenomicRanges::strand(gr))
    cts <- Rsubread::featureCounts(bamlist, annot.ext = saf, nthreads = nthreads, 
                                   isPairedEnd = paired_end,
                                   allowMultiOverlap = F, 
                                   largestOverlap = T, 
                                   requireBothEndsMapped = F)
    cts$counts
}

cts <- quantifyReads(supers, metadata_for_diff_se$K27.BAM, 
                     paired_end = T, # These are paired k27 chips
                     nthreads = 22)
counts <- cts %>% set_colnames(metadata_for_diff_se$ID)

## ----se-deseq------------------------------------------------------------
saveRDS(counts, file = "results/k27_counts.rds")

k27_rse <- SummarizedExperiment(assays = list(counts = counts),
                            rowRanges = GRanges(rownames(counts)))

k27_rse %<>% sort %>% chromVAR::filterPeaks(non_overlapping =T) # remove overlaps

gcview <- Biostrings::Views(genome, rowRanges(k27_rse))
gcFrequency <- Biostrings::letterFrequency(gcview, 
                                           letters = "GC", 
                                           as.prob = TRUE) %>% 
  set_colnames("GC")

mcols(k27_rse) <- cbind(mcols(k27_rse), gcFrequency)

coldata <- metadata_for_diff_se[c("ID","K27.bench.batch","K27.seq.batch","CONDITION")] %>%
  as.data.frame() %>%
  tibble::column_to_rownames("ID")

coldata[c("K27.bench.batch","K27.seq.batch","CONDITION")] %<>% 
  lapply(as.factor)

colData(k27_rse) <- DataFrame(coldata)

eda_data <- newSeqExpressionSet(counts = as.matrix(counts(k27_rse)), 
                                featureData = as.data.frame(mcols(k27_rse, 
                                                                  use.names = T)), 
                                phenoData = coldata["K27.bench.batch"]) 
# for color coding corrected signal plots

dataOffset <- EDASeq::withinLaneNormalization(eda_data, "GC", 
        which = "full", offset = T)
dataOffset <- EDASeq::betweenLaneNormalization(dataOffset, 
        which = "full", offset = T)

EDASeq::biasPlot(eda_data, "GC", log = TRUE, ylim =c(0,10))

EDASeq::biasPlot(dataOffset, "GC", log = TRUE, ylim = c(0, 10))


EDASeqNormFactors <- exp(-1 * EDASeq::offst(dataOffset))
EDASeqNormFactors <- EDASeqNormFactors/exp(rowMeans(log(EDASeqNormFactors)))

counts(k27_rse) <- as.matrix(counts(k27_rse)) # deseq2 wants this to a vanilla matrix

dds <- DESeqDataSet(k27_rse, design = ~ CONDITION)

dds$CONDITION <- relevel(dds$CONDITION, ref = "CD19")

normalizationFactors(dds) <- EDASeqNormFactors

dds <- DESeq(dds,quiet = T)

saveRDS(dds, file = "results/pCLL_se_dds.rds")

## ----se-deseq-res--------------------------------------------------------
res <- lfcShrink(dds, coef = "CONDITION_pCLL_vs_CD19", type = "apeglm") %>% 
  as.data.frame %>% 
  rownames_to_column("Locus") %>%
  as_tibble()

sig <- res %>% filter(padj < 0.1) %>% arrange(desc(log2FoldChange))

## ----export-diff-se-table------------------------------------------------
export_enhancer_table <- function(gr, file) {
  region_id <- as.character(gr)
  chrom <- seqnames(gr) %>% as.vector
  start <- start(gr) %>% as.vector
  stop <- end(gr) %>% as.vector
  num_loci <- 1
  constituent_size <- width(gr)
  bam <- "meta"
  enhancer_rank <- 1
  is_super <- 1
  
  rangedata <- tibble(region_id,chrom,start,stop,num_loci,
                      constituent_size,bam,enhancer_rank,is_super) %>%
    set_names(paste0("V",1:9))
  
  exp <- matrix(nrow = 6, ncol = 9)
  exp[1,] <- c("# Differential Results",rep("",8))
  exp[2,] <- c("# DESeq2",rep("",8))
  exp[3,] <- c("# Multiple bams",rep("",8))
  exp[4,] <- c(paste0("# ", Sys.Date()),rep("",8))
  exp[5,] <- c("# ",rep("",8))
  
  cn <- c("REGION_ID", "CHROM","START", "STOP", 
          "NUM_LOCI", "CONSTITUENT_SIZE", "[bam]",
          "enhancerRank","isSuper")
  
  exp[6,] <- cn
  exp %<>% as.tibble()
  rbind(exp, rangedata) %>% write_tsv(file, col_names = F)
}

sig %>% subset(log2FoldChange > 0) %>% 
  .$Locus %>% GRanges %>% 
  export_enhancer_table("results/gained.SuperEnhancers.txt")

sig %>% subset(log2FoldChange < 0) %>% 
  .$Locus %>% GRanges %>% 
  export_enhancer_table("results/lost.SuperEnhancers.txt")

rownames(counts) %>%
  GRanges() %>%
  export_enhancer_table("results/searchspace.SuperEnhancers.txt")

## python2.7 ~/pipeline/ROSE2_geneMapper.py -i ./results/gained.SuperEnhancers.txt -g HG19 -o ./results/

## ----se-norm-ct-z, results='hide'----------------------------------------
gained_gene_calls <- read_tsv("results/gained_ENHANCER_TO_GENE.txt", skip = 1, col_names = F) %>%
  select(X1, X12) %>% 
  set_colnames(c("Locus", "Closest.Gene"))

lost_gene_calls <- read_tsv("results/lost_ENHANCER_TO_GENE.txt", skip = 1, col_names = F) %>%
  select(X1, X12) %>% 
  set_colnames(c("Locus", "Closest.Gene"))

to_hilite <- read_excel("tables/DifSEs_to_highlight.xlsx", 
                        col_names = c("Locus", "Gene"))

sig_cts <- counts(dds, normalized = T) %>%
  subset(rownames(.) %in% sig$Locus)

row_z <- t(scale(t(sig_cts)))

is_hilite <- rownames(row_z) %>% 
  lapply(FUN = function(x ) { 
    ifelse(x %in% to_hilite$Locus, subset(to_hilite, Locus == x)[["Gene"]], "")
  } ) %>% unlist

## ----se-diff-heat--------------------------------------------------------
se_heat <- main_heatmap(row_z, name = "Row Z-score Norm Cts") %>% 
  add_col_annotation(colData(dds)["CONDITION"]) %>%
  add_row_clustering() %>% 
  add_col_clustering() %>%
  add_row_labels(tickvals = 1:length(is_hilite), 
                 ticktext = is_hilite, side = "right", 
                 font = list(size = 6))
se_heat

## ----save-se-heat, include=F, echo=F-------------------------------------
  save_iheatmap(se_heat, "results/diff-supers.heat.pdf")


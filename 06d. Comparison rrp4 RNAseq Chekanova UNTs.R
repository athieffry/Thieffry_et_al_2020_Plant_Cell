#### Chekanova UNTs: look at rrp4 RNAseq
#### Axel Thieffry
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(rtracklayer)
library(GenomicRanges)
library(CAGEfightR)
library(TeMPO)
library(ggplot2)
library(RColorBrewer)
library(BiocParallel)
register(MulticoreParam(workers=3))
library(TxDb.Athaliana.BioMart.plantsmart28)
options(scipen=10) # disable scientific notation
'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))
setwd('~/masked_path/6d. Chekanova UNTs')


# 1. GET INPUT DATA ####
# ----------------------
# seqinfo
myseqinfo <- read_rds('~/masked_path/myseqinfo.rds')
# Chekanova UNTs AGIs
chekanova_genes <- readxl::read_xlsx('Chekanova 2007 UNTs in rrp4est.xlsx')$UNTs_AGI
# 73 CAGE TCs only up-regulated in hen2-4 (not in rrp4-2)
TCs <- readRDS('~/masked_path/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')

seventythree <- rowRanges(TCs) %>%
  subset(genotypehen2==1 & genotyperrp4 != 1)

# Chekanova RRP4i induced exosome transcripts
chekanova_rrp4i_up <- readxl::read_xls('~/masked_path/mmc2.xls', col_names=c('geneid', 'type'))
chekanova_rrp4i_up$geneid %<>% toupper()
n_distinct(chekanova_rrp4i_up$geneid)

table(seventythree$geneID_anti %in% chekanova_rrp4i_up$geneid)


# 2. export UNTs as bed for metaplot with deeptools (computeMatrix) ####
# ----------------------------------------------------------------------
genes <- genes(TxDb.Athaliana.BioMart.plantsmart28)
  seqlevelsStyle(genes) <- seqlevelsStyle(myseqinfo)
  seqlevels(genes) <- seqlevels(myseqinfo)
  seqinfo(genes) <- myseqinfo
  
  export.bed(subset(genes, gene_id %in% chekanova_genes), '~/masked_path/chekanova_UNTs_genes.bed')

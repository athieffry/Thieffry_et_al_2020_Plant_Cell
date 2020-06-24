#### Arabidopsis TSS story: comparison with Liu 2016 & gene loops
#### Axel Thieffry - October 2019
set.seed(42)
library(patchwork)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(CAGEfightR)
library(readxl)
library(tidyverse)
library(tidylog)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
options(scipen=999)
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { GR[-idx]}
                                     else {GR}}

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6f. Li 2016 - HiC')



# 1. READ DATA ####
# -----------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
# Liu gene loops
geneloops <- read.table('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6f. Li 2016 - HiC/Supplemental_Table_S5.txt', h=T, sep='\t') %>% as_tibble()
# TCs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')
# 3'UTRs
threeUTRs_with_aTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_with_aTSS_gr.rds')
threeUTRs_wout_aTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_without_aTSS_gr.rds')


# 2. OVERLAPS ####
# ----------------
# 2a. get unique gene IDs which loops
geneloops_ids <- geneloops$Gene_ID %>% unique() %>% as.character()

# 2b. get genes with 3'UTR antisense TC
genes_with_3UTR_DE_aTSSs <- rowRanges(TCs) %>%
  subset(txType_TAIR10extended=='antisense_threeUTR') %>%
  subset(genotyperrp4 ==1 | genotypehen2 == 1) %$%
  .$geneID_anti %>%
  unique()

mean(genes_with_3UTR_aTSSs %in% geneloops_ids) * 100





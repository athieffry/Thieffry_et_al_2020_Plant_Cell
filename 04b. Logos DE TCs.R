#### Arabidopsis flg22 : Logos of DE TCs
#### Axel Thieffry - February 2019
set.seed(42)
library(tidyverse)
library(magrittr)
library(reshape2)
library(CAGEfightR)
library(ggplot2)
library(fancycut)
library(RColorBrewer)
library(pheatmap)
library(BiocParallel)
register(MulticoreParam(workers=4))
library(TxDb.Athaliana.BioMart.plantsmart28)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(ggseqlogo)
library(motifStack)
library(seqPattern)
library(TFBSTools)
library(JASPAR2016)
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}
setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE/')



# 1. GET ALL INPUT DATA ####
# --------------------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
      # remove non-canonical chromosomes
      seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
# TCs from DE
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_for_seqlogos.rds')
      # remove non-canonical chromosomes
      seqlevels(TCs, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
# read TAIR10 genome
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack  
# make promoter sequence regions for (shift everything by -1 to account for the safe-trimming)
upstream=40
downstream=15



# 2. SEQLOGO FOR SENSE DE TCs ####
# --------------------------------
# 2.1 get DE TCs and split by sense annotation
TCs_DE_sense_byTxType <- rowRanges(TCs) %>%
  subset(genotypehen2 == -1 | genotyperrp4 == -1) %>%
  split(.$txType_TAIR10)

# 2.2 get sequences after extending range
sense_DE_TCs_seq <- lapply(TCs_DE_sense_byTxType, function(x)
                           swapRanges(x) %>%
                           flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
                           promoters(upstream=upstream, downstream=downstream) %>%
                           remove_out_of_bound() %>%
                           getSeq(genome, .))

(names_sense_DE_TCs <- paste0(names(sense_DE_TCs_seq), ' (N=', lapply(sense_DE_TCs_seq, length), ')'))

sense_DE_TCs_pfm <- mapply(function(x, y) consensusMatrix(x, as.prob=T) %>%
                                          subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
                                          as.matrix() %>%
                                          new('pfm', mat=., name=y),
                           sense_DE_TCs_seq, names_sense_DE_TCs)

sense_DE_TCs_pfm <- mapply(function(x, y) { x@background <- colSums(letterFrequency(y, DNA_BASES)) / sum(colSums(letterFrequency(y, DNA_BASES))) ; x },
                           sense_DE_TCs_pfm, sense_DE_TCs_seq)

# 2.3 plot logos
plusone <- seq(-upstream, downstream-1, 10)
plusone[plusone==0] <- '+1'

mylist_sense <- lapply(sense_DE_TCs_pfm, function(x) x$mat)
names(mylist_sense) <- names_sense_DE_TCs

ggseqlogo(data=mylist_sense, ncol=2, facet='wrap') +
  scale_x_continuous(breaks=seq(1, upstream + downstream, 10), labels=plusone) +
  theme(aspect.ratio=.5) +
  labs(title='Logos at DOWNREGULATED TCs',
       subtitle='sense to gene, DE (down) in any mutant',
       x='Position relative to TSSs (bp)')




# 3. SEQLOGO FOR ANTISENSE DE TCs ####
# --------------------------------
# 2.1 get DE TCs and split by sense annotation
TCs_DE_anti_byTxTypextended <- rowRanges(TCs) %>%
  subset(txType_TAIR10=='antisense') %>%
  subset(genotypehen2 != 0 | genotyperrp4 != 0) %>%
  split(.$txType_TAIR10extended)
      # remove the empty GRs
      TCs_DE_anti_byTxTypextended <- TCs_DE_anti_byTxTypextended[sapply(TCs_DE_anti_byTxTypextended, function(x) length(x)!=0)]
      # sanity check
      sapply(TCs_DE_anti_byTxTypextended, length)

# 2.2 get sequences after extending range
anti_DE_TCs_seq <- lapply(TCs_DE_anti_byTxTypextended, function(x)
                          swapRanges(x) %>%
                          flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
                          promoters(upstream=upstream, downstream=downstream) %>%
                          remove_out_of_bound() %>%
                          getSeq(genome, .))

(names_anti_DE_TCs <- paste0(names(anti_DE_TCs_seq), ' (N=', lapply(anti_DE_TCs_seq, length), ')')) %>% str_replace('_', ' ')

anti_DE_TCs_pfm <- mapply(function(x, y) consensusMatrix(x, as.prob=T) %>%
                          subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
                          as.matrix() %>%
                          new('pfm', mat=., name=y),
                          anti_DE_TCs_seq, names_anti_DE_TCs)

anti_DE_TCs_pfm <- mapply(function(x, y) { x@background <- colSums(letterFrequency(y, DNA_BASES)) / sum(colSums(letterFrequency(y, DNA_BASES))) ; x },
                          anti_DE_TCs_pfm, anti_DE_TCs_seq)

# 2.3 plot logos
plusone <- seq(-upstream, downstream-1, 10)
plusone[plusone==0] <- '+1'

mylist_anti <- lapply(anti_DE_TCs_pfm, function(x) x$mat)
names(mylist_anti) <- names_anti_DE_TCs

ggseqlogo(data=mylist_anti, ncol=2, facet='wrap') +
  scale_x_continuous(breaks=seq(1, upstream + downstream, 10), labels=plusone) +
  theme(aspect.ratio=.5) +
  labs(title='Logos at REGULATED antisense TCs',
       subtitle='antisense to gene, DE (up or down) in any mutant',
       x='Position relative to TSSs (bp)')

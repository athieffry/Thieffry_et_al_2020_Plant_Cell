#### 03d_bis. Seqlogos at TSSs: sanity check with strand splitted
#### Axel Thieffry - July 2018
set.seed(42)
library(tidyverse)
library(tidylog)
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
library(patchwork)
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

setwd('~/masked_path/03 - TSS analysis')


# 1. GET ALL INPUT DATA ####
# --------------------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
# remove non-canonical chromosomes
seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))

# TSS dataset (comparable ones after +/-100bp extension, quantification with CAGE wt T=0 and >= 1TPM in >= 2 libs)
cage_tss <- readRDS('~/masked_path/SE_CAGE_wt_comparable_TSSs_redone_17June2019.rds')
tair_tss <- readRDS('~/masked_path/SE_TAIR10_comparable_TSSs_redone_17June2019.rds')
araport_tss <- readRDS('~/masked_path/SE_ARAPORT11_comparable_TSSs_redone_17June2019.rds')
peat_tss <- readRDS('~/masked_path/SE_PEAT_comparable_TSSs_redone_17June2019.rds')
nanopare_tss <- readRDS('~/masked_path/SE_nanoPARE_comparable_TSSs_redone_17June2019.rds')
rowRanges(nanopare_tss)


# 2. SEQLOGO AT COMPARABLE TSSs ####
# ----------------------------------
#  Read TAIR10 genome
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack  

# make promoter sequence regions for (shift everything by -1 to account for the safe-trimming)
upstream=40
downstream=40

# CAGE WT ctrl (splitted by strand)
cage_seq_p <- cage_tss %>%
  rowRanges() %>%
  subset(strand=='+') %>%
  flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  getSeq(genome, .) # fetch DNAStringSet

cage_seq_m <- cage_tss %>%
  rowRanges() %>%
  subset(strand=='-') %>%
  flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  getSeq(genome, .) # fetch DNAStringSet

cage_pfm_p <- cage_seq_p %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='CAGE wt')

cage_pfm_m <- cage_seq_m %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='CAGE wt')

cage_pfm_p@background <- colSums(letterFrequency(cage_seq_p, DNA_BASES)) / sum(colSums(letterFrequency(cage_seq_p, DNA_BASES)))
cage_pfm_m@background <- colSums(letterFrequency(cage_seq_m, DNA_BASES)) / sum(colSums(letterFrequency(cage_seq_m, DNA_BASES)))


# PEAT
peat_seq_p <- peat_tss %>%
  rowRanges() %>%
  subset(strand=='+') %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

peat_seq_m <- peat_tss %>%
  rowRanges() %>%
  subset(strand=='-') %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

peat_pfm_p <- peat_seq_p %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='PEAT')

peat_pfm_m <- peat_seq_m %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='PEAT')

peat_pfm_p@background <- colSums(letterFrequency(peat_seq_p, DNA_BASES)) / sum(colSums(letterFrequency(peat_seq_p, DNA_BASES)))
peat_pfm_m@background <- colSums(letterFrequency(peat_seq_m, DNA_BASES)) / sum(colSums(letterFrequency(peat_seq_m, DNA_BASES)))

# nanoPARE
nanopare_seq_p <- nanopare_tss %>%
  rowRanges() %>%
  subset(strand=='+') %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

nanopare_seq_m <- nanopare_tss %>%
  rowRanges() %>%
  subset(strand=='-') %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

nanopare_pfm_p <- nanopare_seq_p %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='nanoPARE')

nanopare_pfm_m <- nanopare_seq_m %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='nanoPARE')

nanopare_pfm_p@background <- colSums(letterFrequency(nanopare_seq_p, DNA_BASES)) / sum(colSums(letterFrequency(nanopare_seq_p, DNA_BASES)))
nanopare_pfm_m@background <- colSums(letterFrequency(nanopare_seq_m, DNA_BASES)) / sum(colSums(letterFrequency(nanopare_seq_m, DNA_BASES)))

# TAIR10
tair_seq_p <- tair_tss %>%
  rowRanges() %>%
  subset(strand=='+') %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

tair_seq_m <- tair_tss %>%
  rowRanges() %>%
  subset(strand=='-') %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

tair_pfm_p <- tair_seq_p %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='TAIR10')

tair_pfm_m <- tair_seq_m %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='TAIR10')

tair_pfm_p@background <- colSums(letterFrequency(tair_seq_p, DNA_BASES)) / sum(colSums(letterFrequency(tair_seq_p, DNA_BASES)))
tair_pfm_m@background <- colSums(letterFrequency(tair_seq_m, DNA_BASES)) / sum(colSums(letterFrequency(tair_seq_m, DNA_BASES)))

# ARAPORT11
araport_seq_p <- araport_tss %>%
  rowRanges() %>%
  subset(strand=='+') %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

araport_seq_m <- araport_tss %>%
  rowRanges() %>%
  subset(strand=='-') %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

araport_pfm_p <- araport_seq_p %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='ARAPORT11')

araport_pfm_m <- araport_seq_m %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='ARAPORT11')

araport_pfm_p@background <- colSums(letterFrequency(araport_seq_p, DNA_BASES)) / sum(colSums(letterFrequency(araport_seq_p, DNA_BASES)))
araport_pfm_m@background <- colSums(letterFrequency(araport_seq_m, DNA_BASES)) / sum(colSums(letterFrequency(araport_seq_m, DNA_BASES)))

# plot logos
mylist <- list('CAGE wt (plus)'=cage_pfm_p$mat,
               'CAGE wt (minus)'=cage_pfm_p$mat,
               'PEAT (plus)'=peat_pfm_p$mat,
               'PEAT (minus)'=peat_pfm_m$mat,
               'nanoPARE (plus)'=nanopare_pfm_p$mat,
               'nanoPARE (minus)'=nanopare_pfm_m$mat,
               'TAIR10 (plus)'=tair_pfm_p$mat,
               'TAIR10 (minus)'=tair_pfm_m$mat,
               'ARAPORT11 (plus)'=araport_pfm_p$mat,
               'ARAPORT11 (minus)'=araport_pfm_m$mat)

plusone <- seq(-upstream, downstream-1, 10)
plusone[plusone==0] <- '+1'

list_plus <- mylist[seq(1, 10, 2)]
list_minus <- mylist[seq(2, 10, 2)]

ggseqlogo(data=mylist, ncol=2) +
  scale_x_continuous(breaks=seq(1, upstream + downstream, 10), labels=plusone) +
  theme(aspect.ratio=.25) +
  labs(x='Distance to TSSs (bp)') +
  geom_vline(xintercept=41, lty=2)

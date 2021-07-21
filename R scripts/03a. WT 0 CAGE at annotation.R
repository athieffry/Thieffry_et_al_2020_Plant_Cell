#### 03a. TCs IQR and CAGE at annotations
#### Axel Thieffry - July 2018
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(reshape2)
library(CAGEfightR)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(TeMPO)
library(ggalluvial)
library(RColorBrewer)
library(BiocParallel)
register(MulticoreParam(workers=4))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
genome <- BSgenome.Athaliana.TAIR.TAIR9
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}

setwd('~/masked_path/03 - TSS analysis/')


# 1. CAGE AROUND ANNOTATION TSSs ####
# -----------------------------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
    # remove non-canonical chromosomes
    seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
    # relevel txdb
    seqlevels(txdb) <- seqlevels(myseqinfo)

# TAIR10
txdb_tair <- TxDb.Athaliana.BioMart.plantsmart28
    # remove non-canonical chromosomes
    seqlevelsStyle(txdb_tair) <- seqlevelsStyle(myseqinfo)[1]
    seqlevels(txdb_tair, pruning.mode='coarse') <- seqlevels(myseqinfo)

# ARAPORT11
araport <- import.gff3('~/masked_path/Araport11_GFF3_genes_transposons.201606.gff3')
    # fix seqinfo
    seqlevels(araport, pruning.mode='coarse') <- seqlevels(myseqinfo)
    seqinfo(araport) <- myseqinfo
    # make TxDb object
    txdb_araport <- suppressWarnings( makeTxDbFromGRanges(araport) )

# make promoter regions
tair_tss <- txdb_tair %>% promoters(upstream=0, downstream=1) %>% remove_out_of_bound()
araport_tss <- txdb_araport %>% promoters(upstream=0, downstream=1) %>% remove_out_of_bound()

# combine promoter regions in one GR object
bla <- GRangesList('TAIR10'=tair_tss, 'ARAPORT11'=araport_tss)

# CAGE wt T=0
wt_ctrl_p <- BigWigFile('~/masked_path/wt_0_R123.plus.tpm.bw')
wt_ctrl_m <- BigWigFile('~/masked_path/wt_0_R123.minus.tpm.bw')

# compute coverage around promoters
footprint <- tidyMetaProfile(sites=bla, forward=wt_ctrl_p, reverse=wt_ctrl_m, upstream=500, downstream=501, trimLower=0.01, trimUpper=0.99)

footprint %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
  ggplot(aes(x=pos0, y=score, col=direction, fill=direction)) +
         geom_line() + geom_area() + facet_grid(~sites) +
         cowplot::theme_cowplot() + theme(legend.position=c(.1, .9), aspect.ratio=.8) +
         scale_color_brewer(palette='Set1', direction=-1, name='') +
         scale_fill_brewer(palette='Set1', direction=-1, name='') +
         labs(x='Position relative to TSS (bp)', y='CAGE WT TPM (0.1-0.99%tile)')

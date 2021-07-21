#### comparison with Tokizawa
#### Axel Thieffry
library(tidyverse)
library(magrittr)
library(reshape2)
library(rtracklayer)
library(CAGEfightR)
library(ggpubr)
library(TeMPO)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(RColorBrewer)
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(workers=3))
library(readxl)
library(WriteXLS)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}

setwd('~/masked_path/6b. comparison Tokizawa/')

# 1. READ DATA ####
# -----------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
# CAGE wildtype TSSs (t=0)
TCs_ctrl <- readRDS('~/masked_path/SE_TCs_ctrl.rds') # 21,221
# TAIR10 TSSs (comparable to CAGE WT ctrl)
tair_tss <- readRDS('~/masked_path/SE_TAIR10_comparable_TSSs_with_WTctrl0.rds')
seqinfo(tair_tss) <- myseqinfo
# ARAPORT11 TSSs (comparable to CAGE WT ctrl)
araport_tss <- readRDS('~/masked_path/SE_ARAPORT11_comparable_TSSs_with_WTctrl0.rds')
seqinfo(araport_tss) <- myseqinfo
# Tokizawa
toki_tair <- read.table('Tokizawa data/Data S1 - Information about maximum 5 utr regions - TAIR10.txt', col.names=c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'parent')) %>% as_tibble()
toki_arap <- read.table('Tokizawa data/Data S2 - Information about maximum 5 utr regions - ARAPORT11.txt', col.names=c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'parent')) %>% as_tibble()
      # fix the start & end because those incredible stupid motherfuckers inverted end/start in case TSS on minus strand
      toki_tair %<>%
        mutate('tmp'=start, score=NULL, phase=NULL) %>%
        mutate('start'=ifelse(strand=='-', end, start)) %>%
        mutate('end'=ifelse(strand=='-', tmp, end)) %>%
        select(-tmp)
      
      with(toki_tair, table(start <= end))
      
      toki_arap %<>%
        mutate('tmp'=start, score=NULL, phase=NULL) %>%
        mutate('start'=ifelse(strand=='-', end, start)) %>%
        mutate('end'=ifelse(strand=='-', tmp, end)) %>%
        select(-tmp)
      
      with(toki_arap, table(start <= end))
      
      # make GenomicRanges objects
      toki_tair_gr <- makeGRangesFromDataFrame(df=toki_tair, keep.extra.columns=T, seqnames.field='chr', start.field='start', end.field='end', strand.field='strand')
      seqinfo(toki_tair_gr) <- myseqinfo
      toki_arap_gr <- makeGRangesFromDataFrame(df=toki_arap, keep.extra.columns=T, seqnames.field='chr', start.field='start', end.field='end', strand.field='strand')
      seqinfo(toki_arap_gr) <- myseqinfo
      # look at widths
      rbind(as.data.frame(toki_tair_gr),
            as.data.frame(toki_arap_gr)) %>%
        ggplot(aes(x=width, col=source)) +
               geom_density() +
               scale_x_log10() +
               theme(aspect.ratio=1)

# remove all non-canonical chromosomes
seqlevels(TCs_ctrl, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
seqlevels(tair_tss, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
seqlevels(araport_tss, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
seqlevels(toki_tair_gr, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
seqlevels(toki_arap_gr, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))



# 1. QUICK FOOTPRINT ####
# -----------------------
# MNaseI PlantDHS.org
mnase <- BigWigFile('~/masked_path/plantDHS_mnase_leaf_TPM99pc.bw')

# GRO-seq Nature Plants 2018 (Jacobsen)
# GRO-seq PNAS 2016 (Hetzel)
# 5'GRO-seq PNAS 2016 (Hetzel)
gro_p <- BigWigFileList(list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*plus', full.names=T))
gro_m <- BigWigFileList(list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*minus', full.names=T))
names(gro_p) <- names(gro_m) <- list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*plus') %>% str_remove('_TPM99pc.plus.bw')

# make 1bp GRanges
cage_1bp_tss <- swapRanges(rowRanges(TCs_ctrl))
toki_tair_tss <- resize(toki_tair_gr, width=1, fix='center')
toki_arap_tss <- resize(toki_arap_gr, width=1, fix='center')

# make TSS GRangesList
tss_list <-list("CAGE"=cage_1bp_tss,
                "TAIR10"=tair_tss,
                "ARAPORT11"=araport_tss,
                "TOKI_TAIR10"=toki_tair_tss,
                "TOKI_ARAPORT11"=toki_arap_tss)

# anticipate TeMPO and remove out-of-bound ranges
tss_list_filtered_OoB_400bp <- lapply(tss_list, function(x) promoters(x, upstream=400, downstream=400))
tss_list_filtered_OoB_400bp %<>% lapply(remove_out_of_bound)
tss_list_filtered_OoB_400bp %<>% lapply(function(x) resize(x, width=1, fix='center', use.names=T))

gro_signal <- lapply(tss_list_filtered_OoB_400bp, function(x) tidyMetaProfile(sites=x, forward=gro_p, reverse=gro_m, upstream=400, downstream=400))
mnase_signal <- lapply(tss_list_filtered_OoB_400bp, function(x) tidyMetaProfile(sites=x, forward=mnase, reverse=NULL, upstream=400, downstream=400))

gro_signal %<>% plyr::ldply() %>% as_tibble()
mnase_signal %<>% plyr::ldply() %>% as_tibble()

gro_signal %<>% gather(key="direction", value="score", sense, anti, factor_key=T)

rbind(gro_signal,
      mnase_signal %>% mutate('signal'='MNase', 'direction'='sense', 'score'=sense, 'sense'=NULL)) %>%
  ggplot(aes(x=pos0, y=score, lty=direction, col=.id)) +
         geom_vline(xintercept=0, lty=2) + geom_line(lwd=.8) +
         facet_wrap(~signal, scales='free_y', nrow=2) +
         cowplot::theme_cowplot() + theme(aspect.ratio=.75, legend.position='bottom', legend.direction=1) +
         scale_color_brewer(palette='Set1', name='TSS sources') +
         labs(title='GRO-Seq at TSSs', x='TSSs (bp)', y='Normalized signal')

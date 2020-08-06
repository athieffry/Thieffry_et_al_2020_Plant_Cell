#### Arabidopsis TSS story: comparison with Chekanova 2019 potential enhancers
#### Axel Thieffry
set.seed(42)
library(patchwork)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(TeMPO)
library(rtracklayer)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(CAGEfightR)
library(readxl)
library(tidyverse)
library(tidyquant)
library(tidylog)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
options(scipen=999)
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { GR[-idx]}
                                     else {GR}}

setwd('~/masked_path/6e. Chekanova 2019')



# 1. READ DATA ####
# -----------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
# Chekanova PEs
pes <- readxl::read_xlsx('~/masked_path/Supplemental_Tables_S1_to_S14.xlsx', sheet='PEs', col_names=T, trim_ws=T)
# Chekanova rrp41 sensitive PEs
tmp <- readxl::read_xlsx('~/masked_path/Supplemental_Tables_S1_to_S14.xlsx', sheet='exosensitive_PEs', col_names=T, trim_ws=T) %>% select(PEid, PE_group)
pes %<>% left_join(tmp, by='PEid')
pes %<>% mutate('PE_group'=ifelse(is.na(PE_group), 'Group 0', PE_group)) %>%
                  mutate('rrp41_sensitive'=ifelse(PE_group=='Group 0', FALSE, TRUE))
      # make as GR
      pes <- makeGRangesFromDataFrame(df=pes, keep.extra.columns=T, seqinfo=myseqinfo, seqnames.field='chr', start.field='start', end.field='end', ignore.strand=T)
      # make PE midpoint as IRange
      pes$PE_midpoint <-IRanges(start=pes$PE_midpoint, width=1)
      # sanity check & cleaning
      table(pes$PE_group)
      table(pes$rrp41_sensitive)
      rm(tmp)
# TCs
TCs <- readRDS('~/masked_path/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')
# 3'UTRs
threeUTRs_with_aTSSs <- readRDS('~/masked_path/threeUTRstory_3utrs_with_aTSS_gr.rds')
threeUTRs_wout_aTSSs <- readRDS('~/masked_path/threeUTRstory_3utrs_without_aTSS_gr.rds')
# CAGE data
    cage_p <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*plus'))
    cage_m <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*minus'))
    cage_names <- list.files('~/masked_path/bw_files_R123', pattern='_0.*plus') %>% str_remove('_0_R123.plus.tpm.bw')
    names(cage_p) <- names(cage_m) <- cage_names
# RNA-Seq BIGWIG FILES
    rnaseq_forward <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Forward', full.names=T))
    rnaseq_reverse <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Reverse', full.names=T))
    rnaseq_names <- list.files(path='~/masked_path/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
    names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names
# ALL GROSEQ
    gro_p <- BigWigFileList(c('~/masked_path/GROseq_col_mergedAveraged_R123456.plus.bw',
                              '~/masked_path/5GROcap.aln.unique.plus.bw',
                              '~/masked_path/GROseq_mergedAveraged_R12.plus.bw'))
    gro_m <- BigWigFileList(c('~/masked_path/GROseq_col_mergedAveraged_R123456.minus.bw',
                              '~/masked_path/5GROcap.aln.unique.minus.bw',
                              '~/masked_path/GROseq_mergedAveraged_R12.minus.bw'))
    names(gro_p) <- names(gro_m) <- c('GROseq_Hetzel', 'GROcap_Duttke', 'GROseq_Duttke')



# 2. FOOTPRINTS AT ALL PEs ####
# -----------------------------
# 2a) swap ranges
pes_swapped <- swapRanges(pes, inputColumn='PE_midpoint')
pes_list <- split(pes_swapped, pes_swapped$rrp41_sensitive)
names(pes_list) <- c('insensitive', 'sensitive')
    
# cage
cage_at_pes <- tidyMetaProfile(sites=pes_list, forward=cage_p, reverse=cage_m, upstream=500, downstream=500, trimLower=0.1, trimUpper=0.9) 

cage_at_pes %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  select(-pos1) %>%
  ggplot(aes(x=pos0, y=score, col=direction, group=interaction(sites, signal, direction))) +
         geom_hline(yintercept=0, col='black') +
         geom_ma(n=20, lwd=.6, lty=1) +
         facet_grid(signal~sites) +
         geom_vline(xintercept=0, lty=2) +
         scale_color_brewer(palette='Set1', direction=-1) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold')) +
         labs(title='CAGE footprints', x="Position relative PE midpoints (bp)",
              y='Average CAGE TPM (trimmed 10%, ma 20 bp)')

# rnaseq
rnaseq_at_pes <- tidyMetaProfile(sites=pes_list, forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=1500, downstream=1500, trimLower=0.1, trimUpper=0.9)

rnaseq_at_pes %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  select(-pos1) %>%
  ggplot(aes(x=pos0, y=score, fill=direction, col=direction, group=interaction(sites, direction, signal))) +
         geom_area(alpha=.5) + geom_line() +
         geom_hline(yintercept=0, lty=1) +
         geom_vline(xintercept=0, lty=2) +
         facet_grid(signal~sites) +
         scale_color_brewer(palette='Set1', direction=-1) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold'), legend.position='none') +
         labs(title='RNA-seq footprints', x="Position relative PE midpoints (bp)", y='Average RPM (trimmed 10%, ma 10 bp)')

# groseq
groseq_at_pes <- tidyMetaProfile(sites=pes_list, forward=gro_p, reverse=gro_m, upstream=500, downstream=500, trimLower=0.1, trimUpper=0.9)

groseq_at_pes %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  select(-pos1) %>%
  ggplot(aes(x=pos0, y=score, fill=direction, col=direction, group=interaction(sites, direction, signal))) +
         geom_area(alpha=.5) + geom_line() +
         geom_hline(yintercept=0, lty=1) +
         geom_vline(xintercept=0, lty=2) +
         facet_grid(signal~sites, scales='free_y') +
         scale_fill_brewer(palette='Set1', direction=-1) +
         scale_color_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(legend.position='none', strip.text=element_text(face='bold')) +
         labs(title='GRO-seq footprints', x='Position relative PE midpoints (bp)', y='Normalized signal (trimmed 10%)')



# 3. STATS COMPARED TO US ####
# ----------------------------
# 3a. PEs are overlapping with a 3'UTR with an antisense CAGE TC (N=323)?
subsetByOverlaps(threeUTRs_with_aTSSs, pes, ignore.strand=T)
    # -> 26/323 of 3'UTR with aTC overlap with a PE (8%)

# 3b. CAGE TCs in Chekanova PE's
TCs_in_pes <- subsetByOverlaps(rowRanges(TCs), pes, ignore.strand=T)

gg_a <- TCs_in_pes %>%
  as.data.frame() %>%
  select(txType_TAIR10extended) %>%
  separate(txType_TAIR10extended, c('direction', 'annotation'), sep='_', fill='left') %>%
  mutate('direction'=ifelse(is.na(direction), 'sense', direction)) %>%
  mutate('direction'=ifelse(annotation=='intergenic', 'strandless', direction)) %>%
  ggplot(aes(x=annotation, fill=direction)) +
         geom_bar(position=position_dodge(), lwd=.3, col='black', alpha=.7) +
         geom_text(stat='count', aes(label=..count.., col=direction), position=position_dodge(width=.9), vjust=-.5) +
         cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
         labs(x='CAGE TC annotation', title='CAGE TCs within a Chekanova PE', subtitle='PE region=3000bp centered on DHS midpoint',
              y='Number of TCs')

# 3c. CAGE TCs in Chekanova DHS that are PE's
pe_dhs <- pes
start(pe_dhs) <- pe_dhs$DHS_start
end(pe_dhs) <- pe_dhs$DHS_end

TCs_in_pe_dhs <- subsetByOverlaps(rowRanges(TCs), pe_dhs, ignore.strand=T)

gg_b <- TCs_in_pe_dhs %>%
  as.data.frame() %>%
  select(txType_TAIR10extended) %>%
  separate(txType_TAIR10extended, c('direction', 'annotation'), sep='_', fill='left') %>%
  mutate('direction'=ifelse(is.na(direction), 'sense', direction)) %>%
  mutate('direction'=ifelse(annotation=='intergenic', 'strandless', direction)) %>%
  ggplot(aes(x=annotation, fill=direction)) +
         geom_bar(position=position_dodge(), lwd=.3, col='black', alpha=.7) +
         geom_text(stat='count', aes(label=..count.., col=direction), position=position_dodge(width=.9), vjust=-.5) +
         cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
         labs(x='CAGE TC annotation', title='CAGE TCs within a Chekanova PE DHS', subtitle='PE region = actual DHS region',
              y='Number of TCs')

gg_a + gg_b + plot_layout(ncol=2)

# 3d. how many of our 113 Bidirectional TCs overlap with Chekanova PE's?
BCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_Enhancers.rds')
subsetByOverlaps(rowRanges(BCs), pes, ignore.strand=T) # 5/113 => 4%





# 4. FOOTPRINT OF HISTONE MODIFICATIONS ####
# ------------------------------------------
# ChIP-Seq
chipseq <- c(list.files('~/masked_path/bigwigs', pattern='H3K4me1', full.names=T),
  list.files('~/masked_path/bigwigs', pattern='H3K4me2', full.names=T),
  list.files('~/masked_path/bigwigs', pattern='H3K4me3', full.names=T),
  list.files('~/masked_path/bigwigs', pattern='H3K27me3', full.names=T),
  list.files('~/masked_path/bigwigs', pattern='H3K27ac', full.names=T)) %>% BigWigFileList()

names(chipseq) <- c(list.files('~/masked_path/bigwigs', pattern='H3K4me1'),
  list.files('~/masked_path/bigwigs', pattern='H3K4me2'),
  list.files('~/masked_path/bigwigs', pattern='H3K4me3'),
  list.files('~/masked_path/bigwigs', pattern='H3K27me3'),
  list.files('~/masked_path/bigwigs', pattern='H3K27ac')) %>% str_remove('\\.bw')


# Footprint and plot
pes_swapped_4chip <- pes_swapped
seqlevelsStyle(pes_swapped_4chip) <- seqlevelsStyle(chipseq$Inagaki2017_H3K4me1)
seqlevels(pes_swapped_4chip) <- seqlevels(chipseq$Inagaki2017_H3K4me1)
seqinfo(pes_swapped_4chip) <- seqinfo(chipseq$Inagaki2017_H3K4me1)

chipseq_at_pes <- tidyMetaProfile(sites=pes_swapped_4chip, forward=chipseq,
                                  upstream=1500, downstream=1500, trimLower=0.01, trimUpper=0.99)

chipseq_at_pes %>%
  separate(signal, c('author', 'histone'), sep='_') %>%
  ggplot(aes(x=pos0, y=sense, col=histone, group=interaction(author, histone))) +
         geom_vline(xintercept=0, lty=2) +
         geom_line(lwd=1) +
         facet_wrap(~histone, scales='free_y') +
         cowplot::theme_cowplot() + theme(legend.position='none') +
         labs(title='Histone marks', x='Position relative PE midpoints (bp)',
              y='Average signal (trimmed 1%)')



# 5. FOOTPRINT OF  MNaseI & DNaseI from PlantDHS.org ####
# -------------------------------------------------------
dmnase <- BigWigFileList(c('~/masked_path/plantDHS_mnase_leaf_TPM99pc.bw',
                           '~/masked_path/plantDHS_dnase_leaf_TPM99pc.bw'))
names(dmnase) <- c('MNase_leaf', 'DNase_leaf')

# Footprint and plot
dmnase_at_pes <- tidyMetaProfile(sites=pes_swapped, forward=dmnase, upstream=1000, downstream=1000, trimLower=0.01, trimUpper=0.99)

ggplot(dmnase_at_pes, aes(x=pos0, y=sense, col=signal)) +
       geom_vline(xintercept=0, lty=2) +
       geom_ma(n=15, lty=1) +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       labs(title='MNase & DNaseI Footprints', x='Position relative PE midpoints (bp)',
            y='Average signal (trimmed 1% ma 15 bp)', caption='leaf samples')

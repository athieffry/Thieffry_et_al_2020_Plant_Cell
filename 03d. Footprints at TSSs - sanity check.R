#### 03c. All sort of footprints around TSSs: sanity check
#### Axel Thieffry - July 2018
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(reshape2)
library(CAGEfightR)
library(patchwork)
library(rtracklayer)
library(GenomicRanges)
library(TeMPO)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(Hmisc)
library(BiocParallel)
register(MulticoreParam(workers=3))
library(TxDb.Athaliana.BioMart.plantsmart28)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(ggseqlogo)
library(motifStack)
library(seqPattern)
options(scipen=10) # disable scientific notation
'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}
setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/03 - TSS analysis')




# 1. GET INPUT DATA AND FIX SEQINFO ####
# --------------------------------------
# myseqinfo (remove non-canonical chromosomes)
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
    seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
# TSS dataset (comparable ones after +/-100bp extension, quantification with CAGE wt T=0 and >= 1TPM in >= 2 libs)
cage_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CAGE_wt_comparable_TSSs_redone_17June2019.rds')
tair_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TAIR10_comparable_TSSs_redone_17June2019.rds')
araport_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_ARAPORT11_comparable_TSSs_redone_17June2019.rds')
peat_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_PEAT_comparable_TSSs_redone_17June2019.rds')
nanopare_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_nanoPARE_comparable_TSSs_redone_17June2019.rds')




# 2. NORMALIZE INTO TPM AND KEEP 99 PERCENTILE ####
# -------------------------------------------------
# ALL Histone marks (as provided by mivanov)
# bigwig files
chipseq <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ChIP-seq/bigwigs', pattern='.bw', full.names=T) %>% BigWigFileList()
names(chipseq) <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ChIP-seq/bigwigs', pattern='.bw') %>% str_remove('.bw')
# bigwig normalization factors (to RPM)
chipseq_rpm_factor <- read.table('~/Dropbox/Axel_Arabidopsis_Flagellin/ChIP-seq/total_bedgraph_reads.txt', h=T, sep='\t') %>%
  as_tibble() %>%
  mutate('Sample'=str_remove(Sample, '.bedgraph'), 
         'Author'=str_split(Sample, '_', simplify=T)[,1],
         'Histone_mark'=str_split(Sample, '_', simplify=T)[,2],
         'RPM_factor'=total_reads_in_bedgraph / 1000000)
# MNaseI PlantDHS.org
mnase <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/03 - TSS analysis/plantDHS_mnase_leaf_TPM99pc.bw')
# DNaseI PlantDHS.org
dnase <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/03 - TSS analysis/plantDHS_dnase_leaf_TPM99pc.bw')
# GRO-seq
gro_p <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/02 - Enhancers analyses/', pattern='GRO.*plus', full.names=T))
gro_m <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/02 - Enhancers analyses/', pattern='GRO.*minus', full.names=T))
names(gro_p) <- names(gro_m) <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/02 - Enhancers analyses/', pattern='GRO.*plus') %>% str_remove('_TPM99pc.plus.bw')
# make TSS GRangesList
tss_list_plus <-list("CAGE"=rowRanges(cage_tss) %>% subset(strand=='+'),
                     "nanoPARE"=rowRanges(nanopare_tss) %>% subset(strand=='+'),
                     "PEAT"=rowRanges(peat_tss) %>% subset(strand=='+'),
                     "TAIR10"=rowRanges(tair_tss) %>% subset(strand=='+'),
                     "ARAPORT11"=rowRanges(araport_tss) %>% subset(strand=='+'))

tss_list_minus <-list("CAGE"=rowRanges(cage_tss) %>% subset(strand=='-'),
                      "nanoPARE"=rowRanges(nanopare_tss) %>% subset(strand=='-'),
                      "PEAT"=rowRanges(peat_tss) %>% subset(strand=='-'),
                      "TAIR10"=rowRanges(tair_tss) %>% subset(strand=='-'),
                      "ARAPORT11"=rowRanges(araport_tss) %>% subset(strand=='-'))



# 3. MNASE & GROSEQ FOOTPRINTS ####
# ---------------------------------
gro_signal_p <- lapply(tss_list_plus, function(x) tidyMetaProfile(sites=x, forward=gro_p, reverse=gro_m, upstream=400, downstream=400))
gro_signal_m <- lapply(tss_list_minus, function(x) tidyMetaProfile(sites=x, forward=gro_p, reverse=gro_m, upstream=400, downstream=400))
mnase_signal_p <- lapply(tss_list_plus, function(x) tidyMetaProfile(sites=x, forward=mnase, reverse=NULL, upstream=400, downstream=400))
mnase_signal_m <- lapply(tss_list_minus, function(x) tidyMetaProfile(sites=x, forward=mnase, reverse=NULL, upstream=400, downstream=400))

gro_signal_p %<>% plyr::ldply() %>% as_tibble()
gro_signal_m %<>% plyr::ldply() %>% as_tibble()
mnase_signal_p %<>% plyr::ldply() %>% as_tibble()
mnase_signal_m %<>% plyr::ldply() %>% as_tibble()

# finding respective maxima
gro_signal_p %>%
  subset(direction=='sense') %>%
  select(-pos1) %>%
  group_by(.id, signal, direction) %>%
  top_n(wt=score, n=1)

gro_signal_m %>%
  select(-pos1) %>%
  group_by(.id, signal, direction) %>%
  top_n(wt=score, n=1) %>%
  pull(pos0)

mnase_signal_p %>%
  select(-pos1) %>%
  group_by(.id) %>%
  top_n(wt=sense, n=1)

mnase_signal_m %>%
  select(-pos1) %>%
  group_by(.id) %>%
  top_n(wt=sense, n=1)

# plot footprint
gro_signal_p %<>% gather(key="direction", value="score", sense, anti, factor_key=T)
gro_signal_m %<>% gather(key="direction", value="score", sense, anti, factor_key=T)

rbind(gro_signal_p %>% mutate('strand'='plus'),
      gro_signal_m %>% mutate('strand'='minus'),
      mnase_signal_p %>% mutate('signal'='MNase', 'direction'='sense', 'score'=sense, 'sense'=NULL, 'strand'='plus'),
      mnase_signal_m %>% mutate('signal'='MNase', 'direction'='sense', 'score'=sense, 'sense'=NULL, 'strand'='minus')) %>%
  ggplot(aes(x=pos0, y=score, lty=direction, col=.id, group=interaction(.id, strand, direction))) +
         geom_vline(xintercept=0, lty=2) +
         geom_line(lwd=.2) +
         facet_wrap(~signal, scales='free_y') +
         cowplot::theme_cowplot() + theme(aspect.ratio=.75, legend.position='right') +
         scale_color_brewer(palette='Set1', name='TSS sources') +
         labs(title='Footprints at TSSs (plus)', x='Distance to TC peaks/TSSs (bp)', y='Normalized signal')

# re-finding respective maxima
gro_signal %>%
  subset(direction=='sense') %>%
  select(-pos1, -pos1, -direction) %>%
  group_by(.id, signal) %>%
  top_n(wt=score, n=1) %>%
  subset(signal=='GROseq_jacobsen')

# distribution of distance between ARAPORT11 and CAGE TSSs
distanceToNearest(tss_list$ARAPORT11, tss_list$CAGE) %>%
  as.data.frame() %$%
  summary(.$distance)



# 4. QUANTIFY WITH CAGE ####
# --------------------------
# Extend TSSs +/- 100bp
tair_extended_tss <- promoters(rowRanges(tair_tss), upstream=100, downstream=101)
araport_extended_tss <- promoters(rowRanges(araport_tss), upstream=100, downstream=101)
peat_extended_tss <- promoters(rowRanges(peat_tss), upstream=100, downstream=101)
nanopare_extended_tss <- promoters(rowRanges(nanopare_tss), upstream=100, downstream=101)

# function for quantifying any features (GR) across all samples
quantifyAcrossSample <- function(features, scoring_SE, inputAssay='TPM', sample){
  hits <- as(findOverlaps(query=features, subject=rowRanges(scoring_SE)), "List")
  newSE <- suppressWarnings( swapScores(object=scoring_SE, outputColumn='score', inputAssay=inputAssay, sample=sample) )
  sum(extractList(rowRanges(newSE)$score, hits))
}
# quantify TSSs across PEAT extended regions
CTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CTSSs_1count_min3lib_TSSstory.rds')
CTSSs_wt <- subset(CTSSs, select=genotype=='wt')
sampleNames <- CTSSs_wt$Name %>% as.character()

tair_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=tair_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
araport_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=araport_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
peat_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=peat_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
nanopare_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=nanopare_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))

# average expression across WT replicates
peat_meanTPM <- peat_extended_tss_TPM %>%
  melt() %>%
  separate(Var2, c('wt', 'zero', 'rep'), sep='_') %>%
  select(Var1, value) %>%
  group_by(Var1) %>%
  summarise('meanTPM_PEAT'=mean(value)) %>%
  ungroup() %>%
  select(-Var1) %>%
  as_tibble()

cage_meanTPM <- assay(cage_tss, 'TPM') %>%
  melt() %>%
  separate(Var2, c('wt', 'zero', 'rep'), sep='_') %>%
  select(-wt:-rep) %>%
  group_by(Var1) %>%
  summarise('meanTPM_CAGE'=mean(value)) %>%
  ungroup() %>%
  select(-Var1) %>%
  as_tibble()

tair_meanTPM <- tair_extended_tss_TPM %>%
  melt() %>%
  separate(Var2, c('wt', 'zero', 'rep'), sep='_') %>%
  select(Var1, value) %>%
  group_by(Var1) %>%
  summarise('meanTPM_TAIR'=mean(value)) %>%
  ungroup() %>%
  select(-Var1) %>%
  as_tibble()

araport_meanTPM <- araport_extended_tss_TPM %>%
  melt() %>%
  separate(Var2, c('wt', 'zero', 'rep'), sep='_') %>%
  select(Var1, value) %>%
  group_by(Var1) %>%
  summarise('meanTPM_ARAPORT'=mean(value)) %>%
  ungroup() %>%
  select(-Var1) %>%
  as_tibble()

nanopare_meanTPM <- nanopare_extended_tss_TPM %>%
  melt() %>%
  separate(Var2, c('wt', 'zero', 'rep'), sep='_') %>%
  select(Var1, value) %>%
  group_by(Var1) %>%
  summarise('meanTPM_NanoPARE'=mean(value)) %>%
  ungroup() %>%
  select(-Var1) %>%
  as_tibble()

ggplot() +
  geom_density(data=peat_meanTPM, aes(x=meanTPM_PEAT, col='PEAT'), size=1) +
  geom_density(data=cage_meanTPM, aes(x=meanTPM_CAGE, col='CAGE wt'), size=1) +
  geom_density(data=tair_meanTPM, aes(x=meanTPM_TAIR, col='TAIR'), size=1) +
  geom_density(data=araport_meanTPM, aes(x=meanTPM_ARAPORT, col='ARAPORT'), size=1) +
  geom_density(data=nanopare_meanTPM, aes(x=meanTPM_NanoPARE, col='NanoPARE'), size=1) +
  scale_x_log10() +
  labs(title='Expression density of comparable TSSs',
       subtitle='PEAT, NanoPARE, TAIR & ARAPORT TSSs were\nextended +/-100bp and quantified with CAGE',
       x='TPM expression (CAGE wt, average rep.)') +
  scale_color_brewer(palette='Set1', name='TSS sources') + cowplot::theme_cowplot() +
  theme(aspect.ratio=1)

# subsample PEAT to match CAGE ctrl expression distribution
# add PEAT mean CAGE ctrl TPM to PEAT GR
peat_gr$meanTPM <- peat_meanTPM$meanTPM_PEAT
peat_gr$logMeanTPM <- log(peat_meanTPM$meanTPM_PEAT+1)
peat_gr
# do the same for CAGE WT
rowRanges(TCs_ctrl)$meanTPM <- CAGE_meanTPM$meanTPM_CAGE
rowRanges(TCs_ctrl)$logMeanTPM <- log(CAGE_meanTPM$meanTPM_CAGE+1)
rowRanges(TCs_ctrl)
# a. break down PEAT distribution into bins
bins_peat <- seq(from=min(peat_gr$logMeanTPM), to=max(peat_gr$logMeanTPM), length.out=30)
# b. associate PEAT & CAGE TSSs to PEAT bins
peat_gr$bin <- cut2(x=peat_gr$logMeanTPM, cuts=bins_peat) %>% as.character()
rowRanges(TCs_ctrl)$bin <- cut2(x=rowRanges(TCs_ctrl)$logMeanTPM, cuts=bins_peat) %>% as.character()
# c. make dataframes
cage_df <- rowRanges(TCs_ctrl) %>% as.data.frame()
peat_df <- as.data.frame(peat_gr)
# sanity check of number of bins
left_join(data.frame(table(peat_df$bin)) %>% set_colnames(c('bins', 'n_PEAT')),
          data.frame(table(cage_df$bin)) %>% set_colnames(c('bins', 'n_CAGE')), by='bins') %>%
  mutate('ok'=ifelse(n_CAGE >= n_PEAT, TRUE, FALSE))
# d. loop
# initialize empty final DFs
final_cage_df <- final_peat_df <- data.frame()

for (i in 1:length(peat_gr)) {
  # 1. sample random PEAT
  index_peat <- sample(1:nrow(peat_df), size=1, replace=F)
  peat_rnd_sample <- peat_df[sample(1:nrow(peat_df), size=1, replace=F), ]
  # 2. remove sampled PEAT from dataset
  peat_df <- peat_df[-index_peat, ]
  # 2. get all CAGE with same bin as the sampled PEAT
  cage_with_same_bin <- subset(cage_df, bin == peat_rnd_sample$bin)
  # check if there are CAGE left with in that bin
  if(nrow(cage_with_same_bin)==0) {message('TSS skipped') ; next }
  # 3. sample randomly
  index_cage <- sample(1:nrow(cage_with_same_bin), size=1, replace=F)
  cage_rnd_sample <- cage_with_same_bin[index_cage, ]
  # 5. remove sampled CAGE from dataset
  cage_df %<>% subset(rownames(cage_df)!=rownames(cage_rnd_sample))
  # 4. add both to their final df
  final_peat_df <- rbind(peat_rnd_sample, final_peat_df)
  final_cage_df <- rbind(cage_rnd_sample, final_cage_df)
}

# sanity check
table(final_cage_df$seqnames, useNA='ifany')
dim(final_peat_df) ; dim(final_cage_df)

# make density plots
gg_sub_density <- rbind(data.frame('source'='PEAT',
                                   'logMeanTPM'=final_peat_df$logMeanTPM,
                                   'type'='subsampled'),
                        data.frame('source'='CAGE',
                                   'logMeanTPM'=final_cage_df$logMeanTPM,
                                   'type'='subsampled'),
                        data.frame('source'='PEAT',
                                   'logMeanTPM'=log(peat_meanTPM$meanTPM_PEAT+1),
                                   'type'='original'),
                        data.frame('source'='CAGE',
                                   'logMeanTPM'=log(CAGE_meanTPM$meanTPM_CAGE+1),
                                   'type'='original')) %>%
  ggplot(aes(x=logMeanTPM, col=source, lty=type)) +
  geom_density(size=1) +
  cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position=c(0.6, 0.7)) +
  labs(title='TPM expression\nCAGE & PEAT TSSs after subsampling') +
  scale_color_brewer(palette='Set1')
gg_sub_density

# re-transform into GR
subsampled_peat_gr <- makeGRangesFromDataFrame(df=final_peat_df, keep.extra.columns=T, seqnames.field='seqnames', start.field='start', end.field='end')
subsampled_cage_gr <- makeGRangesFromDataFrame(df=final_cage_df, keep.extra.columns=T, seqnames.field='seqnames', start.field='start', end.field='end')
# fix internal IRanges for TSS peak (thick.xxx)
subsampled_cage_gr$thick <- IRanges(start=subsampled_cage_gr$thick.start, end=subsampled_cage_gr$thick.end)
subsampled_cage_gr$thick.start <- subsampled_cage_gr$thick.end <- NULL

# make footprint plots
tss_subsampled_list <-list("CAGE"=swapRanges(subsampled_cage_gr),
                           "PEAT"=subsampled_peat_gr)

tss_subsampled_list_filtered_OoB_400bp <- lapply(tss_subsampled_list, function(x) promoters(x, upstream=400, downstream=400))
tss_subsampled_list_filtered_OoB_400bp %<>% lapply(remove_out_of_bound)
tss_subsampled_list_filtered_OoB_400bp %<>% lapply(function(x) resize(x, width=1, fix='center', use.names=T))

gro_signal <- lapply(tss_subsampled_list_filtered_OoB_400bp, function(x) tidyMetaProfile(sites=x, forward=gro_p, reverse=gro_m, upstream=400, downstream=400))
mnase_signal <- lapply(tss_subsampled_list_filtered_OoB_400bp, function(x) tidyMetaProfile(sites=x, forward=mnase, reverse=NULL, upstream=400, downstream=400))

gro_signal %<>% plyr::ldply() %>% as_tibble()
mnase_signal %<>% plyr::ldply() %>% as_tibble()

gro_signal %<>% gather(key="direction", value="score", sense, anti, factor_key=T)

gg_sub_footprints <- rbind(gro_signal,
                           mnase_signal %>% mutate('signal'='MNase', 'direction'='sense', 'score'=sense, 'sense'=NULL)) %>%
  ggplot(aes(x=pos0, y=score, lty=direction, col=.id)) +
  geom_vline(xintercept=0, lty=2) + geom_line(lwd=.8) +
  facet_wrap(~signal, scales='free_y', nrow=1, ncol=4) +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75, legend.position='bottom', legend.direction=1) +
  scale_color_brewer(palette='Set1', name='TSS sources', direction=-1) +
  labs(title='Footprint at subsampled TSSs',
       caption=paste0('N(CAGE TSSs)=', length(subsampled_cage_gr),' - N(PEAT TSSs)=', length(subsampled_peat_gr)),
       x='TSSs (bp)', y='Normalized signal')

# make Logo
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack  

upstream=40
downstream=10

# CAGE WT subsampled
cage_sub_seq <- subsampled_cage_gr %>%
  swapRanges() %>% # take TC summit positions
  flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  remove_out_of_bound() %>% # remove out of Chr
  getSeq(genome, .) # fetch DNAStringSet

cage_sub_pfm <- cage_sub_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='CAGE ctrl subsampled')

cage_sub_pfm@background <- colSums(letterFrequency(cage_sub_seq, DNA_BASES)) / sum(colSums(letterFrequency(cage_sub_seq, DNA_BASES)))

# PEAT
peat_sub_seq <- subsampled_peat_gr %>%
  flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  remove_out_of_bound() %>% # remove out of Chr
  getSeq(genome, .) # fetch DNAStringSet

peat_sub_pfm <- peat_sub_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='PEAT subsampled')

peat_sub_pfm@background <- colSums(letterFrequency(peat_sub_seq, DNA_BASES)) / sum(colSums(letterFrequency(peat_sub_seq, DNA_BASES)))

# plot logos
mylist <- list('CAGE subsampled (9,125)'=cage_sub_pfm$mat,
               'PEAT subsampled (Morton et al.) (9,125)'=peat_sub_pfm$mat)
plusone <- seq(-upstream, downstream-1, 10)
plusone[plusone==0] <- '+1'

gg_logo <- ggseqlogo(data=mylist, nrow=2) +
  scale_x_continuous(breaks=seq(1, upstream + downstream, 10), labels=plusone) +
  cowplot::theme_cowplot() + theme(aspect.ratio=.5) +
  labs(title='Logos', x='Position relative to TSSs (bp)')


Rmisc::multiplot(gg_sub_density, gg_logo, gg_sub_footprints, layout=matrix(c(1,1,2,2,3,3,3,3), byrow=T, nrow=2, ncol=4))

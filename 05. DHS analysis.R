#### DHS analysis
#### Axel Thieffry
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(rtracklayer)
library(CAGEfightR)
library(ggpubr)
library(TeMPO)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(RColorBrewer)
library(VennDiagram)
library(BiocParallel)
library(patchwork)
library(tidyquant)
library(ggExtra)
register(SerialParam())
register(MulticoreParam(workers=3))
library(WriteXLS)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}

setwd('~/masked_path/05 - DHS analysis')



# 1. INPUT DATA ####
# ------------------
# seqinfo
myseqinfo <- read_rds('~/masked_path/myseqinfo.rds')

# RNA-Seq BIGWIG FILES
rnaseq_forward <- list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Forward_R123_fixed', full.names=T)
rnaseq_reverse <- list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Reverse_R123_fixed', full.names=T)
rnaseq_forward %<>% BigWigFileList()
rnaseq_reverse %<>% BigWigFileList()
rnaseq_names <- list.files(path='~/masked_path/bigwigs_RPM_R123_fixed', pattern='Forward_R123_fixed', full.names=F) %>% str_remove('_Forward_R123_fixed.bw')
names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names

# CAGE BIGWIG FILES
cage_p <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*plus'))
cage_m <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*minus'))
cage_names <- list.files('~/masked_path/bw_files_R123', pattern='_0.*plus') %>% str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names

# DHS sites
dhss <- readRDS('~/masked_path/new_dhssumits_annotated_plantDHSonly_Summits_with_flower.rds')

# CAGE WT TSSs
CTSSs <- readRDS('~/masked_path/SE_CTSSs_1count_min3lib_TSSstory.rds')
wt_CTSSs <- subset(CTSSs, select=genotype=='wt')
wt_CTSSs <- subsetBySupport(wt_CTSSs, inputAssay='counts', unexpressed=0, minSamples=2)
    # percentage of CAGE signal captured by DHSs
    total_signal_WT_CTSSs <- sum(Matrix::rowSums(assay(wt_CTSSs, 'counts')))
    wt_signal_in_DHS <- subsetByOverlaps(wt_CTSSs, promoters(dhss, upstream=400, downstream=400), ignore.strand=T)
    total_signal_WT_CTSSs_in_DHS <- sum(Matrix::rowSums(assay(wt_signal_in_DHS, 'counts')))
    total_signal_WT_CTSSs_in_DHS / total_signal_WT_CTSSs * 100
    # percentage of CAGE CTSSs captured by DHSs
    total_WT_CTSSs <- length(wt_CTSSs)
    wt_CTSSs_in_DHS <- subsetByOverlaps(wt_CTSSs, promoters(dhss, upstream=400, downstream=400), ignore.strand=T) %>% length()
    wt_CTSSs_in_DHS / total_WT_CTSSs * 100

# Histone marks (H3K27me3a is from PlantDHS.org - H3K27me3b is from V. Colot 2011 paper)
histones <- list.files('~/masked_path/03 - TSS analysis', pattern='H3', full.names=T) %>% BigWigFileList()
names(histones) <- list.files('~/masked_path/03 - TSS analysis', pattern='H3') %>%
              str_remove('plantDHS_leaf_') %>% str_remove('_TPM99pc.bw') %>% str_remove('GSE50636_At_') %>%
              str_remove('GSE74841_coverageNormSize_') %>% str_remove('_coVcomp.sorted.bw') %>% str_remove('_Col_wt_R1.bw')
# DNaseI
dnase <- list.files('~/masked_path/DNase', pattern='flower', full.names=T) %>% BigWigFileList()
names(dnase) <- 'DNaseI'

# GROcap Hetzel
grocap_p <- list.files('~/masked_path/02 - Enhancers analyses', pattern='GROcap.*plus', full.names=T) %>% BigWigFileList()
grocap_m <- list.files('~/masked_path/02 - Enhancers analyses', pattern='GROcap.*minus', full.names=T) %>% BigWigFileList()
names(grocap_p) <- names(grocap_m) <- "5GRO-cap"

# GROseq Hetzel
grohet_p <- list.files('~/masked_path/02 - Enhancers analyses', pattern='GROseq_hetzel.*plus', full.names=T) %>% BigWigFileList()
grohet_m <- list.files('~/masked_path/02 - Enhancers analyses', pattern='GROseq_hetzel.*minus', full.names=T) %>% BigWigFileList()
names(grohet_p) <- names(grohet_m) <- 'GROseq_Hetzel'




# 2. CAGE profile at DHSs (intergenic and promoters) ####
# -------------------------------------------------------
dhss_prom_inter <- subset(dhss, txType %in% c('promoter', 'intergenic'))
dhss_prom_inter$txType %<>% droplevels()
dhss_prom_inter_byTxType <- split(dhss_prom_inter, dhss_prom_inter$txType)

cage_profile_promoter_dhss <- tidyMetaProfile(sites=dhss_prom_inter_byTxType, forward=cage_p, reverse=cage_m, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)

# average ratio between sense and antisense signals at promoter-DHSs
cage_profile_promoter_dhss %>%
  select(-pos1) %>%
  group_by(signal, sites) %>%
  summarise('mean_sense'=mean(sense), 
            'mean_antisense'=mean(anti)) %>%
  mutate('sense_ratio'=mean_sense/mean_antisense)

# plot
cage_profile_promoter_dhss %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('sites'=ifelse(sites=='promoter', paste0('promoter (N=', length(dhss_prom_inter_byTxType$promoter), ')'), paste0('intergenic (N=', length(dhss_prom_inter_byTxType$intergenic), ')'))) %>%
  mutate('sites'=factor(sites, levels=unique(sites))) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
         geom_line(alpha=.4, lty=1) +
         geom_ma(lwd=1, n=15, lty=1) +
         facet_grid(sites~signal, scales='free_y') +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         scale_color_brewer(palette='Set1', direction=-1) +
         labs(title='CAGE at promoter and intergenic DHSs', y='CAGE TPM (15bp moving average)', x='DHS Summits (bp)')



# 3. FOOTPRINTS AT INTERGENIC DHSs ####
# -------------------------------------
# 3a. get intergenic DHSs only
dhss_intergenic <- subset(dhss, txType=='intergenic') # 9130

# 3b. compute average coverage for all histone marks + DNaseI + CAGE
hist_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=histones, reverse=NULL, upstream=1000, downstream=1000, trimLower=0.01, trimUpper=0.99)
dnase_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=dnase, reverse=NULL, upstream=1000, downstream=1000, trimLower=0.01, trimUpper=0.99)
cage_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=cage_p, reverse=cage_m, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
rnaseq_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
groseqhet_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=grohet_p, reverse=grohet_m, upstream=400, downstream=400, trimLower=0.1, trimUpper=0.9)
grocap_intergenic_DHS <- tidyMetaProfile(sites=dhss_intergenic, forward=grocap_p, reverse=grocap_m, upstream=400, downstream=400, trimLower=0.1, trimUpper=0.9)

# 3c. plot histone marks (selected by Albin) + DNaseI
gg_hist_intergenic_DHS <- rbind(hist_intergenic_DHS, dnase_intergenic_DHS) %>%
  subset(signal %in% c('DNaseI', 'H3K27ac', 'H3K4me1', 'H3K4me3')) %>%
  subset(signal=='DNaseI') %>%
  ggplot(aes(x=pos0, y=sense, col=signal)) +
         geom_line() +
         geom_vline(xintercept=0, lty=2) +
         cowplot::theme_cowplot() + theme(legend.position='none', aspect.ratio=1) +
         facet_wrap(~signal, scales='free_y',ncol=1) +
         labs(title='Histones marks', x='Distance form DHS summits (bp)', y='Average normalized signal') +
         scale_color_brewer(palette='Dark2', name='')
gg_hist_intergenic_DHS + xlim(-800, 800)

# 3d. look if there's a TSS bias in one direction, which would explain the systematic histone mark shift to the right
gg_cage_intergenic_DHS <- cage_intergenic_DHS %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) + 
         geom_line(alpha=.3) + geom_ma(n=15, lty=1, lwd=1) +
         geom_vline(xintercept=0, lty=2) +
         facet_wrap(~signal, nrow=3) +
         cowplot::theme_cowplot() + theme(legend.position='none') +
         labs(title='CAGE', x='Distance from DHS summits (bp)', y='Average CAGE TPM') +
         scale_color_brewer(palette='Set1', direction=-1)

# 3e. investigate RNA-seq profiles as well
gg_rnaseq_intergenic_DHS <- rnaseq_intergenic_DHS %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=ifelse(signal=='DM', 'lsm8-2.rrp4', 
                         ifelse(signal=='WT', 'wt',
                                ifelse(signal=='RRP4', 'rrp4', 'lsm8-2')))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) + 
         geom_line(lwd=1) +
         geom_vline(xintercept=0, lty=2) +
         facet_wrap(~signal, nrow=4) +
         cowplot::theme_cowplot() + theme(legend.position='none') +
         labs(title='RNA-seq', x='Distance from DHS summits (bp)', y='Average RNAseq RPM') +
         scale_color_brewer(palette='Set1', direction=-1)

# 3f. GROseq footprints as well
gg_groseqhet_intergenic_DHS <- groseqhet_intergenic_DHS %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('direction'=ifelse(direction=='sense', 'forward strand', 'reverse strand')) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) + 
         geom_line(lwd=1) +
         geom_vline(xintercept=0, lty=2) +
         cowplot::theme_cowplot() +  theme(legend.position='none') +
         labs(title='GRO-seq (Hetzel)', x='Distance from DHS summits (bp)', y='Average normalized signal') +
         scale_color_brewer(palette='Set1', direction=-1)

gg_grocap_intergenic_DHS <- grocap_intergenic_DHS %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('direction'=ifelse(direction=='sense', 'forward strand', 'reverse strand')) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) + 
         geom_line(lwd=1) +
         geom_vline(xintercept=0, lty=2) +
         cowplot::theme_cowplot() +  theme(legend.position='bottom', legend.direction=2) +
         labs(title="5'GRO-cap (Hetzel)", x='Distance from DHS summits (bp)', y='Average normalized signal') +
         scale_color_brewer(palette='Set1', direction=-1)

gg_hist_intergenic_DHS + gg_cage_intergenic_DHS + gg_rnaseq_intergenic_DHS + (gg_groseqhet_intergenic_DHS + gg_grocap_intergenic_DHS + plot_layout(ncol=1)) + plot_layout(ncol=4)




# 4. CAGE profile at DHSs (all other types) ####
# ----------------------------------------------
dhss_others <- subset(dhss, txType %!in% c('promoter', 'intergenic'))
dhss_others$txType %<>% droplevels()

dhss_others_byTxType <- split(dhss_others, dhss_others$txType)

cage_profile_other_dhss <- tidyMetaProfile(sites=dhss_others_byTxType, forward=cage_p, reverse=cage_m, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)

gg_5UTR <- cage_profile_other_dhss %>%
  subset(sites=='fiveUTR') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='', x='')

gg_3UTR <- cage_profile_other_dhss %>%
  subset(sites=='threeUTR') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='', x='DHS Summits (bp)')

gg_CDS <- cage_profile_other_dhss %>%
  subset(sites=='CDS') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='CAGE TPM (15bp moving average)')

gg_exon <- cage_profile_other_dhss %>%
  subset(sites=='exon') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='', x='')

gg_intron <- cage_profile_other_dhss %>%
  subset(sites=='intron') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='', x='')

gg_proximal <- cage_profile_other_dhss %>%
  subset(sites=='proximal') %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
  geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
  geom_line(alpha=.4, lty=1) +
  geom_ma(lwd=1, n=15, lty=1) +
  facet_grid(sites~signal, scales='free_y') +
  cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
  scale_color_brewer(palette='Set1', direction=-1) +
  labs(y='', x='')

table(dhss_others$txType)

ggarrange(gg_proximal, gg_5UTR, gg_CDS, gg_exon, gg_intron, gg_3UTR, align='hv', common.legend=T, nrow=6, legend='bottom')



# 5. RNA-seq profile at DHSs (intergenic and promoters) ####
# ----------------------------------------------------------
# fix seqinfo and re-output new BigWigs as 'fixed'
if(FALSE){
  # fix the fucking seqinfo
  rnaseq_forward <- list.files('~/masked_path/bigwigs_RPM_R123', pattern='Forward', full.names=T)
  rnaseq_reverse <- list.files('~/masked_path/bigwigs_RPM_R123', pattern='Reverse', full.names=T)
  rnaseq_forward %<>% lapply(import.bw)
  rnaseq_reverse %<>% lapply(import.bw)
  rnaseq_names <- list.files(path='~/masked_path/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
  names(rnaseq_forward) <- rnaseq_names
  names(rnaseq_reverse) <- rnaseq_names
  
  
  rnaseq_forward <- lapply(rnaseq_forward, function(x) {seqlevels(x) <- seqlevels(myseqinfo); x})
  rnaseq_reverse <- lapply(rnaseq_reverse, function(x) {seqlevels(x) <- seqlevels(myseqinfo); x})
  
  rnaseq_forward <- lapply(rnaseq_forward, function(x) {seqinfo(x) <- myseqinfo; x})
  rnaseq_reverse <- lapply(rnaseq_reverse, function(x) {seqinfo(x) <- myseqinfo; x})
  
  rnaseq_forward <- lapply(rnaseq_forward, trim)
  rnaseq_reverse <- lapply(rnaseq_reverse, trim)
  
  setwd('~/masked_path/bigwigs_RPM_R123_fixed')
  mapply(function(x, y) export.bw(x, y), rnaseq_forward, paste0(names(rnaseq_forward), '_Forward_R123_fixed.bw'))
  mapply(function(x, y) export.bw(x, y), rnaseq_reverse, paste0(names(rnaseq_forward), '_Reverse_R123_fixed.bw'))}

# update seqlengths of DHSs to match those of RNA-seq
dhss_prom_inter_forRNAseq <- dhss_prom_inter
seqlevels(dhss_prom_inter_forRNAseq) <- seqlevels(rnaseq_forward$DM)
seqlengths(dhss_prom_inter_forRNAseq) <- seqlengths(rnaseq_forward$DM)
dhss_prom_inter_forRNAseq_byTxType <- split(dhss_prom_inter_forRNAseq, dhss_prom_inter_forRNAseq$txType)

# compute signal
rnaseq_profile_promoter_dhss <- tidyMetaProfile(sites=dhss_prom_inter_forRNAseq_byTxType, forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)

# average ratio between sense and antisense signals at promoter-DHSs
rnaseq_profile_promoter_dhss %>%
  select(-pos1) %>%
  group_by(signal, sites) %>%
  summarise('mean_sense'=mean(sense), 
            'mean_antisense'=mean(anti)) %>%
  mutate('sense_ratio'=mean_sense/mean_antisense)

# plot
rnaseq_profile_promoter_dhss %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  mutate('sites'=ifelse(sites=='promoter', paste0('promoter (N=', length(dhss_prom_inter_forRNAseq_byTxType$promoter), ')'), paste0('intergenic (N=', length(dhss_prom_inter_forRNAseq_byTxType$intergenic), ')'))) %>%
  mutate('sites'=factor(sites, levels=unique(sites))) %>%
  mutate('signal'=ifelse(signal=='DM', 'rrp4 lsm8-2',
                         ifelse(signal=='WT', 'wt',
                                ifelse(signal=='RRP4', 'rrp4', 'lsm8-2')))) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'lsm8-2', 'rrp4', 'rrp4 lsm8-2'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_vline(xintercept=0, lty=2, col='grey50') + geom_hline(yintercept=0, col='grey50') +
         geom_line(lwd=1) +
         facet_grid(sites~signal, scales='free_y') +
         cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
         scale_color_brewer(palette='Set1', direction=-1) +
         labs(title='RNA-seq at promoter and intergenic DHSs', y='Average RPM', x='DHS Summits (bp)') 

# plot signals by sites rather than by sequencing technology, and with two different Y-axes
  # cage DF
  cage_df <- cage_profile_promoter_dhss %>%
    gather(key="direction", value="score", sense, anti, factor_key=T) %>%
    mutate('sites'=ifelse(sites=='promoter', paste0('promoter (N=', length(dhss_prom_inter_byTxType$promoter), ')'), paste0('intergenic (N=', length(dhss_prom_inter_byTxType$intergenic), ')'))) %>%
    mutate('sites'=factor(sites, levels=unique(sites))) %>%
    mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
    mutate('tech'='CAGE')
  
  rnaseq_df <- rnaseq_profile_promoter_dhss %>%
    gather(key="direction", value="score", sense, anti, factor_key=T) %>%
    mutate('sites'=ifelse(sites=='promoter', paste0('promoter (N=', length(dhss_prom_inter_forRNAseq_byTxType$promoter), ')'), paste0('intergenic (N=', length(dhss_prom_inter_forRNAseq_byTxType$intergenic), ')'))) %>%
    mutate('sites'=factor(sites, levels=unique(sites))) %>%
    mutate('signal'=ifelse(signal=='DM', 'rrp4 lsm8-2',
                           ifelse(signal=='WT', 'wt',
                                  ifelse(signal=='RRP4', 'rrp4', 'lsm8-2')))) %>%
    mutate('signal'=factor(signal, levels=c('wt', 'lsm8-2', 'rrp4', 'rrp4 lsm8-2'))) %>%
    mutate('tech'='RNA-seq')
  
# plot both
  all_levels <- factor(c('wt', 'hen2', 'lsm8-2', 'rrp4', 'rrp4 lsm8-2'), levels=c('wt', 'hen2', 'lsm8-2', 'rrp4', 'rrp4 lsm8-2'))
  
  cage_df %>%
  mutate('signal'=factor(signal, levels=all_levels)) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_vline(xintercept=0, lty=2, col='grey50') +
         geom_hline(yintercept=0, col='grey50') +
         geom_line(alpha=.4, lty=1) +
         geom_ma(lwd=1, n=15, lty=1) +
         facet_grid(sites~signal, scales='free_y', drop=F) +
         scale_color_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() +
         theme(aspect.ratio=1, legend.position='bottom', legend.direction=1) +
         labs(title='CAGE at promoter and intergenic DHSs',
              y='CAGE TPM (15bp moving average)', x='DHS Summits (bp)')
  
  rnaseq_df %>%
  mutate('signal'=factor(signal, levels=all_levels)) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
    geom_vline(xintercept=0, lty=2, col='grey50') +
    geom_hline(yintercept=0, col='grey50') +
    geom_line(lwd=1, lty=2) +
    facet_grid(sites~signal, scales='free_y', drop=F) +
    cowplot::theme_cowplot() +
    theme(aspect.ratio=1, legend.position='bottom', legend.direction=1) +
    scale_color_brewer(palette='Set1', direction=-1) +
    labs(title='RNA-seq at promoter and intergenic DHSs',
         y='CAGE TPM (15bp false abcdefg)', x='DHS Summits (bp)') 



# 6. INTERGENIC DHS DIRECTIONALITY (CAGE) ####
# --------------------------------------------
# 6a. export intergenic DHSSs +/- 800bp as BED file
export.bed(dhss_intergenic %>% resize(width=600, fix='center', use.names=T), 'dhss_intergenic_600bp.bed.txt')

# 6b. prepare left and right DHS regions
dhss_intergenic_left_gr  <- dhss_intergenic
dhss_intergenic_right_gr <- dhss_intergenic

strand(dhss_intergenic_left_gr) <- '-'
strand(dhss_intergenic_right_gr) <- '+'

export.bed(promoters(dhss_intergenic_left_gr, upstream=0, downstream=300), 'dhss_intergenic_LEFT_300bp.bed.txt')
export.bed(promoters(dhss_intergenic_right_gr, upstream=0, downstream=300), 'dhss_intergenic_RIGHT_300bp.bed.txt')

# 6c. compute CAGE directionality at intergenic DHSs
cage_left_arm <- lapply(cage_m, function(x) wideMetaProfile(sites=dhss_intergenic_left_gr, forward=x, reverse=NULL, upstream=1, downstream=300))
cage_right_arm <- lapply(cage_p, function(x) wideMetaProfile(sites=dhss_intergenic_right_gr, forward=x, reverse=NULL, upstream=1, downstream=300))
      # sum arm DHS signal at each arm and add name
      cage_left_arm_sum <- cage_left_arm %>% sapply(rowSums) %>% as_tibble() %>% mutate('dhs_id'=dhss_intergenic_left_gr$name)
      cage_right_arm_sum <- cage_right_arm %>% sapply(rowSums) %>% as_tibble() %>% mutate('dhs_id'=dhss_intergenic_right_gr$name)
      # merge
      cage_left_arm_sum  %<>% melt(id.vars='dhs_id', variable.name='genotype', value.name='leftSum') %>% as_tibble()
      cage_right_arm_sum %<>% melt(id.vars='dhs_id', variable.name='genotype', value.name='rightSum') %>% as_tibble()
      
      cage_arm <- left_join(cage_left_arm_sum, cage_right_arm_sum, by=c('dhs_id', 'genotype')) %>% as_tibble()
      # compute directionality
      cage_arm %<>% mutate('total' = leftSum + rightSum,
                           'dir' = rightSum / total)
      # make colors
      display.brewer.all()
      cage_cols <- brewer.pal(name='Set1', n=9)[c(3, 4, 1)]
      rnaseq_cols <- brewer.pal(name='Set1', n=9)[c(3, 5, 1, 7)]
      # plot directionality: barplot
      gg_dirbar_cage <- cage_arm %>%
        mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
        subset(!is.na(dir)) %>%
        ggplot(aes(x=dir, fill=genotype)) +
               geom_histogram(stat='bin', binwidth=.1, position=position_dodge(preserve='total'), col='black', lwd=.2, alpha=.7) +
               cowplot::theme_cowplot() + theme(legend.position=c(0.25, 0.75), axis.text=element_text(size=8), axis.title=element_text(size=9)) +
               scale_fill_manual(values=cage_cols, name='') +
               labs(title='Intergenic DHSSs: directionality',
                    x='CAGE directionality', y='Nb. DHSs / 9130')
      
      # plot directionality and expression
      gg_direxp_cage <- cage_arm %>%
        mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
        subset(!is.na(dir)) %>%
        ggplot(aes(x=dir, y=total, col=genotype)) +
               geom_point(size=0.5) +
               scale_y_log10() +
               cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none', axis.text=element_text(size=8), axis.title=element_text(size=9)) +
               facet_grid(~genotype) +
               labs(title='Intergenic DHSs: expression vs directionality', x='', y='DHS expression (CAGE TPM)') +
               scale_color_manual(values=cage_cols, name='')



# 7. INTERGENIC DHS DIRECTIONALITY (RNA-SEQ) ####
# -----------------------------------------------
rnaseq_left_arm <- lapply(rnaseq_forward, function(x) wideMetaProfile(sites=dhss_intergenic_left_gr, forward=x, reverse=NULL, upstream=1, downstream=300))
rnaseq_right_arm <- lapply(rnaseq_reverse, function(x) wideMetaProfile(sites=dhss_intergenic_right_gr, forward=x, reverse=NULL, upstream=1, downstream=300))
# sum arm DHS signal at each arm and add name
rnaseq_left_arm_sum <- rnaseq_left_arm %>% sapply(rowSums) %>% as_tibble() %>% mutate('dhs_id'=dhss_intergenic_left_gr$name)
rnaseq_right_arm_sum <- rnaseq_right_arm %>% sapply(rowSums) %>% as_tibble() %>% mutate('dhs_id'=dhss_intergenic_right_gr$name)
# merge
rnaseq_left_arm_sum  %<>% melt(id.vars='dhs_id', variable.name='genotype', value.name='leftSum') %>% as_tibble()
rnaseq_right_arm_sum %<>% melt(id.vars='dhs_id', variable.name='genotype', value.name='rightSum') %>% as_tibble()

rnaseq_arm <- left_join(rnaseq_left_arm_sum, rnaseq_right_arm_sum, by=c('dhs_id', 'genotype')) %>% as_tibble()
# compute directionality
rnaseq_arm %<>% mutate('total' = leftSum + rightSum,
                       'dir' = rightSum / total)
# rename genotypes
rnaseq_arm %<>% mutate('genotype'=ifelse(genotype=='WT', 'wt',
                                         ifelse(genotype=='RRP4', 'rrp4',
                                                ifelse(genotype=='LSM8', 'lsm8-2', 'lsm8-2.rrp4'))))
# plot directionality: barplot
gg_dirbar_rnaseq <- rnaseq_arm %>%
  mutate('genotype'=factor(genotype, levels=c('wt', 'lsm8-2', 'rrp4', 'lsm8-2.rrp4'))) %>%
  subset(!is.na(dir)) %>%
  ggplot(aes(x=dir, fill=genotype)) +
         geom_histogram(stat='bin', binwidth=.1, position=position_dodge(preserve='total'), col='black', lwd=.2, alpha=.7) +
         cowplot::theme_cowplot() + theme(legend.position=c(0.25, 0.75), axis.text=element_text(size=8), axis.title=element_text(size=9)) +
         scale_fill_manual(values=rnaseq_cols, name='') +
         labs(x='RNA-seq directionality', y='')

# plot directionality and expression
gg_direxp_rnaseq <- rnaseq_arm %>%
  mutate('genotype'=factor(genotype, levels=c('wt', 'lsm8-2', 'rrp4', 'lsm8-2.rrp4'))) %>%
  subset(!is.na(dir)) %>%
  ggplot(aes(x=dir, y=total, col=genotype)) +
         geom_point(size=0.5) +
         scale_y_log10() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none', axis.text=element_text(size=8), axis.title=element_text(size=9)) +
         facet_grid(~genotype) +
         labs(x='RNA-seq directionality', y='DHS expression (RNAseq RPM)') +
         scale_color_manual(values=rnaseq_cols, name='')

# plot directionality of CAGE -vs- RNA-seq
cage_arm_wt <- subset(cage_arm, genotype=='wt') %>% select(dhs_id, 'cage_wt_dir'=dir)
cage_arm_rrp4 <- subset(cage_arm, genotype=='rrp4') %>% select(dhs_id, 'cage_rrp4_dir'=dir)
rnaseq_arm_wt <- subset(rnaseq_arm, genotype=='wt') %>% select(dhs_id, 'rnaseq_wt_dir'=dir)
rnaseq_arm_rrp4 <- subset(rnaseq_arm, genotype=='rrp4') %>% select(dhs_id, 'rnaseq_rrp4_dir'=dir)

rrp4_dirByTech <-  plyr::join_all(dfs=list(cage_arm_wt, cage_arm_rrp4, rnaseq_arm_wt, rnaseq_arm_rrp4), type='left', by='dhs_id') %>% as_tibble()

gg_dirdir_wt <- ggplot(rrp4_dirByTech, aes(x=cage_wt_dir, y=rnaseq_wt_dir)) +
      geom_point(col=cage_cols[1]) +
      cowplot::theme_cowplot() + theme(aspect.ratio=1) +
      labs(title='Intergenic DHS directionality', x='CAGE directionality', y='RNA-seq directionality')

gg_dirdir_rrp4 <- ggplot(rrp4_dirByTech, aes(x=cage_rrp4_dir, y=rnaseq_rrp4_dir)) +
       geom_point(col=cage_cols[3], size=0.5) +
       geom_rug(data=subset(rrp4_dirByTech, is.na(cage_rrp4_dir)), alpha=.5) + # to plot rnaseq dir in margin when cage dir is NA
       geom_rug(data=subset(rrp4_dirByTech, is.na(rnaseq_rrp4_dir)), alpha=.5) + # to plot cage dir in margin when rnaseq dir is NA
       cowplot::theme_cowplot() + theme(aspect.ratio=1, axis.text=element_text(size=8), axis.title=element_text(size=9)) +
       labs(title='rrp4 comparison', x='CAGE dir.', y='RNA-seq dir.')

((gg_dirbar_cage + gg_dirbar_rnaseq + plot_layout(ncol=2, nrow=1)) / (gg_direxp_cage + gg_dirdir_rrp4 + plot_layout(ncol=2, nrow=1, widths=c(3,1))) + gg_direxp_rnaseq)



# get CAGE DHS-enhancer candidates
# --------------------------------
# make expression DFs
cage_expression_df <- cage_arm %>%
  select(dhs_id, genotype, total) %>%
  spread(genotype, total) %>%
  rename('hen2_tpm'='hen2', 'rrp4_tpm'='rrp4', 'wt_tpm'='wt')

rnaseq_expression_df <- rnaseq_arm %>%
  select(dhs_id, genotype, total) %>%
  spread(genotype, total) %>%
  rename('lsm8_rpm'=`lsm8-2`, 'dm_rpm'=`lsm8-2.rrp4`, 'rrp4_rpm'='rrp4', 'wt_rpm'='wt')

# reshape dataframes
# remove unidirectional WT's
# add total expression
cage_dir_df <- cage_arm %>%
  select(dhs_id, genotype, dir) %>%
  spread(genotype, dir) %>%
  subset(wt >= 0.65 | wt <= 0.35 | is.na(wt)) %>%
  left_join(cage_expression_df, by='dhs_id')

rnaseq_dir_df <- rnaseq_arm %>%
  select(dhs_id, genotype, dir) %>%
  spread(genotype, dir) %>%
  subset(wt >= 0.65 | wt <= 0.35 | is.na(wt)) %>%
  rename('lsm8'=`lsm8-2`, 'dm'=`lsm8-2.rrp4`) %>%
  left_join(rnaseq_expression_df, by='dhs_id')

# CAGE candidates (on directionality & mutant expression higher than in wt:
# hen2 (76)
cage_hen2_candidates <- cage_dir_df %>%
  subset(hen2 <= 0.65 & hen2 >= 0.35) %>%
  subset(hen2_tpm > wt_tpm)
# rrp4 (90)
cage_rrp4_candidates <- cage_dir_df %>%
  subset(rrp4 <= 0.65 & rrp4 >= 0.35) %>%
  subset(rrp4_tpm > wt_tpm)

# RNA-seq candidates (on directionality & mutant expression higher than in wt:
# lsm8 (5)
rnaseq_lsm8_candidates <- rnaseq_dir_df %>%
  subset(lsm8 <= 0.64 & lsm8 >= 0.35) %>%
  subset(lsm8_rpm > wt_rpm)
# rrp4 (72)
rnaseq_rrp4_candidates <- rnaseq_dir_df %>%
  subset(rrp4 <= 0.64 & rrp4 >= 0.35) %>%
  subset(rrp4_rpm > wt_rpm)
# dm (61)
rnaseq_dm_candidates <- rnaseq_dir_df %>%
  subset(dm <= 0.64 & dm >= 0.35) %>%
  subset(dm_rpm > wt_rpm)

# make Venn diagrams
data.frame('all_candidates'=unique(c(cage_hen2_candidates$dhs_id, cage_rrp4_candidates$dhs_id, rnaseq_lsm8_candidates$dhs_id, rnaseq_rrp4_candidates$dhs_id, rnaseq_dm_candidates$dhs_id))) %>%
  mutate(#'CAGE_hen2'=ifelse(all_candidates %in% cage_hen2_candidates$dhs_id, T, F),
         #'CAGE_rrp4'=ifelse(all_candidates %in% cage_rrp4_candidates$dhs_id, T, F))
         'RNAseq_lsm8'=ifelse(all_candidates %in% rnaseq_lsm8_candidates$dhs_id, T, F),
         'RNAseq_rrp4'=ifelse(all_candidates %in% rnaseq_rrp4_candidates$dhs_id, T, F),
         'RNAseq_lsm8.rrp4'=ifelse(all_candidates %in% rnaseq_dm_candidates$dhs_id, T, F)) %>%
  column_to_rownames('all_candidates') %>%
  eulerr::euler() %>% plot(quantities=T)

# graphical table
data.frame('genotype'=c('hen2-4', 'rrp4', 'lsm8-2', 'rrp4', 'lsm8-2 rrp4'),
           'sequencing_platform'=c(rep('CAGE', 2), rep('RNAseq', 3)),
           'bidir_intergenic_DHSs'=c(82, 95, 7, 76, 65),
           '..._with_expMUT_above_expWT'=c(76, 90, 5, 72, 61)) %>%
  gt::gt()

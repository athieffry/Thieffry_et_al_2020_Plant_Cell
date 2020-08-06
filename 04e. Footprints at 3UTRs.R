#### Arabidopsis flg22 : footprints at 3'UTRs with and without aTSSs
#### Axel Thieffry
set.seed(42)
library(WriteXLS)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggalluvial)
library(ggridges)
library(pheatmap)
library(patchwork)
library(venn)
library(GGally)
library(RColorBrewer)
library(edgeR)
library(limma)
library(DESeq2)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(TeMPO)
library(rtracklayer)
library(gridExtra)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(CAGEfightR)
library(tidyverse)
library(gProfileR)
library(tidylog)
library(viridis)
library(GeneOverlap)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                    if(length(idx) != 0) { GR[-idx]}
                                    else {GR}}
                                    get_density <- function(x, y, ...) {dens <- MASS::kde2d(x, y, ...)
                                    ix <- findInterval(x, dens$x)
                                    iy <- findInterval(y, dens$y)
                                    ii <- cbind(ix, iy)
                                    return(dens$z[ii])}
setwd('~/masked_path/04 - TSS_Level DE/')



# 1. LOAD ALL INPUT FILES ####
# ----------------------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
# 3'UTRs
threeUTRs_with_aTSSs <- readRDS('~/masked_path/threeUTRstory_3utrs_with_aTSS_gr.rds')
threeUTRs_wout_aTSSs <- readRDS('~/masked_path/threeUTRstory_3utrs_without_aTSS_gr.rds')
# genes
genes <- genes(TxDb.Athaliana.BioMart.plantsmart28)
      seqlevelsStyle(genes) <- seqlevelsStyle(myseqinfo)
      seqinfo(genes) <- myseqinfo
# TCs
TCs <- readRDS('~/masked_path/SE_TCs_TPM1_min3lib.rds')
# CAGE BIGWIG FILES
cage_p <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*plus'))
cage_m <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*minus'))
cage_names <- list.files('~/masked_path/bw_files_R123', pattern='_0.*plus') %>% str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names
# RNA-Seq BIGWIG FILES
rnaseq_forward <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Forward', full.names=T))
rnaseq_reverse <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Reverse', full.names=T))
rnaseq_names <- list.files(path='~/masked_path/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names
# all GROseq at one
gro_p <- BigWigFileList(c('~/masked_path/GROseq_col_mergedAveraged_R123456.plus.bw',
                          '~/masked_path/5GROcap.aln.unique.plus.bw',
                          '~/masked_path/GROseq_mergedAveraged_R12.plus.bw'))
gro_m <- BigWigFileList(c('~/masked_path/GROseq_col_mergedAveraged_R123456.minus.bw',
                          '~/masked_path/5GROcap.aln.unique.minus.bw',
                          '~/masked_path/GROseq_mergedAveraged_R12.minus.bw'))
names(gro_p) <- names(gro_m) <- c('GROseq_Hetzel', 'GROcap_Duttke', 'GROseq_Duttke')
# ChIP-Seq PlantDHS.org
chipseq <- BigWigFileList(c(list.files('~/masked_path/03 - TSS analysis', pattern='leaf_H', full.names=T),
             list.files('~/masked_path/03 - TSS analysis', pattern='GSE', full.names=T)))

names_chipseq <- c(list.files('~/masked_path/03 - TSS analysis', pattern='leaf_H'),
                   list.files('~/masked_path/03 - TSS analysis', pattern='GSE')) %>%
  str_remove('plantDHS_leaf_') %>%
  str_remove('_TPM99pc.bw') %>%
  str_remove('GSE74841_coverageNormSize_') %>%
  str_remove('_Col_wt_R1.bw') %>%
  str_remove('GSE50636_At_') %>%
  str_remove('_coVcomp.sorted.bw')

names(chipseq) <- names_chipseq
# MNaseI & DNaseI from PlantDHS.org
dmnase <- BigWigFileList(c('~/masked_path/plantDHS_mnase_leaf_TPM99pc.bw',
                           '~/masked_path/plantDHS_dnase_leaf_TPM99pc.bw'))
names(dmnase) <- c('MNase_leaf', 'DNase_leaf')



# 2. DISCARD 3'UTRs OVERLAPPING ANYTHING ELSE ####
# ------------------------------------------------
threeUTRs_with_aTSSs_noOverlap <- threeUTRs_with_aTSSs[countOverlaps(threeUTRs_with_aTSSs + 250, genes, ignore.strand=T) == 1] # 248 out of 323 don't overlap anything else
threeUTRs_wout_aTSSs_noOverlap <- threeUTRs_wout_aTSSs[countOverlaps(threeUTRs_wout_aTSSs + 250, genes, ignore.strand=T) == 1] # 19,060 out of 22,366 don't overlap anyting else


# 3. MAKE 3'UTR REGIONS ####
# --------------------------
# 5'end
threeUTRs_with_aTSSs_5end <- resize(threeUTRs_with_aTSSs_noOverlap, width=1, fix='start')
threeUTRs_wout_aTSSs_5end <- resize(threeUTRs_wout_aTSSs_noOverlap, width=1, fix='start')
threeUTRs_5ends <- GRangesList(threeUTRs_with_aTSSs_5end, threeUTRs_wout_aTSSs_5end)
names(threeUTRs_5ends) <- c('with aTSSs', 'witout aTSSs')
# center
threeUTRs_with_aTSSs_center <- resize(threeUTRs_with_aTSSs_noOverlap, width=1, fix='center')
threeUTRs_wout_aTSSs_center <- resize(threeUTRs_wout_aTSSs_noOverlap, width=1, fix='center')
threeUTRs_centers <- GRangesList(threeUTRs_with_aTSSs_center, threeUTRs_wout_aTSSs_center)
names(threeUTRs_centers) <- c('with aTSSs', 'witout aTSSs')
# 3'end
threeUTRs_with_aTSSs_3end <- resize(threeUTRs_with_aTSSs_noOverlap, width=1, fix='end')
threeUTRs_wout_aTSSs_3end <- resize(threeUTRs_wout_aTSSs_noOverlap, width=1, fix='end')
threeUTRs_3ends <- GRangesList(threeUTRs_with_aTSSs_3end, threeUTRs_wout_aTSSs_3end)
names(threeUTRs_3ends) <- c('with aTSSs', 'witout aTSSs')



# 3. COMPUTE FOOTPRINTS ####
# --------------------------
# 3a. CAGE
cage_at_threeUTRs_5ends   <- tidyMetaProfile(sites=threeUTRs_5ends,   forward=cage_p, reverse=cage_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
cage_at_threeUTRs_centers <- tidyMetaProfile(sites=threeUTRs_centers, forward=cage_p, reverse=cage_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
cage_at_threeUTRs_3ends   <- tidyMetaProfile(sites=threeUTRs_3ends,   forward=cage_p, reverse=cage_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)

gg_cage_a <- cage_at_threeUTRs_5ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction, col=direction)) +
         geom_area() + facet_grid(signal ~ sites) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         scale_color_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold'), legend.position='none') +
         labs(title='', x="3'UTR 5'ends (bp)", y='Average CAGE TPM', subtitle="5'end")

gg_cage_b <- cage_at_threeUTRs_centers %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction, col=direction)) +
         geom_area() + facet_grid(signal ~ sites) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         scale_color_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold'), legend.position='none') +
         labs(title='CAGE footprint', x="3'UTR centers (bp)", y='', subtitle="center")

gg_cage_c <- cage_at_threeUTRs_3ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction, col=direction)) +
         geom_area() + facet_grid(signal ~ sites) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         scale_color_brewer(palette='Set1', direction=-1) +
         cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold')) +
         labs(title='', x="3'UTR 3'ends (bp)", y='', subtitle="3'end")

gg_cage_a + gg_cage_b + gg_cage_c + plot_layout(ncol=3)

# 3b. RNA-seq
rnaseq_at_threeUTRs_5ends   <- tidyMetaProfile(sites=threeUTRs_5ends,   forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
rnaseq_at_threeUTRs_centers <- tidyMetaProfile(sites=threeUTRs_centers, forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
rnaseq_at_threeUTRs_3ends   <- tidyMetaProfile(sites=threeUTRs_3ends,   forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)

gg_rnaseq_a <- rnaseq_at_threeUTRs_5ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites) +
  scale_fill_brewer(palette='Set2') +
  cowplot::theme_cowplot() + theme(legend.position='none', strip.text=element_text(face='bold')) +
  labs(title='', x="3'UTR 5'ends (bp)", y='Average FPM', subtitle="5'end")

gg_rnaseq_b <- rnaseq_at_threeUTRs_centers %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites) +
  scale_fill_brewer(palette='Set2') +
  cowplot::theme_cowplot() + theme(legend.position='none', strip.text=element_text(face='bold')) +
  labs(title='RNA-Seq footprint', x="3'UTR centers (bp)", y='', subtitle="center")

gg_rnaseq_c <- rnaseq_at_threeUTRs_3ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites) +
  scale_fill_brewer(palette='Set2') +
  cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold')) +
  labs(title='', x="3'UTR 3'ends (bp)", y='', subtitle="3'end")

gg_rnaseq_a + gg_rnaseq_b + gg_rnaseq_c + plot_layout(ncol=3)

# 3c. GRO-Seq
groseq_at_threeUTRs_5ends   <- tidyMetaProfile(sites=threeUTRs_5ends,   forward=gro_p, reverse=gro_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
groseq_at_threeUTRs_centers <- tidyMetaProfile(sites=threeUTRs_centers, forward=gro_p, reverse=gro_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)
groseq_at_threeUTRs_3ends   <- tidyMetaProfile(sites=threeUTRs_3ends,   forward=gro_p, reverse=gro_m, upstream=200, downstream=200, trimLower=0.01, trimUpper=0.99)

gg_groseq_a <- groseq_at_threeUTRs_5ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites, scales='free_y') +
  scale_fill_brewer(palette='Set1', direction=-1) +
  cowplot::theme_cowplot() + theme(legend.position='none', strip.text=element_text(face='bold')) +
  labs(title='', x="3'UTR 5'ends (bp)", y='Average signal', subtitle="5'end")

gg_groseq_b <- groseq_at_threeUTRs_centers %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites, scales='free_y') +
  scale_fill_brewer(palette='Set1', direction=-1) +
  cowplot::theme_cowplot() + theme(legend.position='none', strip.text=element_text(face='bold')) +
  labs(title='GRO-Seq footprint', x="3'UTR centers (bp)", y='', subtitle="center")

gg_groseq_c <- groseq_at_threeUTRs_3ends %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
  geom_area() + facet_grid(signal ~ sites, scales='free_y') +
  scale_fill_brewer(palette='Set1', direction=-1) +
  cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold')) +
  labs(title='', x="3'UTR 3'ends (bp)", y='', subtitle="3'end")

gg_groseq_a + gg_groseq_b + gg_groseq_c + plot_layout(ncol=3)

# 3d. MNase & DNase
mnase_at_threeUTRs_5ends   <- tidyMetaProfile(sites=threeUTRs_5ends,   forward=dmnase, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
mnase_at_threeUTRs_centers <- tidyMetaProfile(sites=threeUTRs_centers, forward=dmnase, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
mnase_at_threeUTRs_3ends   <- tidyMetaProfile(sites=threeUTRs_3ends,   forward=dmnase, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)

gg_mnase_a <- ggplot(mnase_at_threeUTRs_5ends, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + ylim(0, NA) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + theme(legend.position='none') +
  labs(title='', x="3'UTR 5'ends (bp)", y='Average signal', subtitle="5'end")

gg_mnase_b <- ggplot(mnase_at_threeUTRs_centers, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + ylim(0, NA) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + theme(legend.position='none') +
  labs(title='MNase footprint', x="3'UTR centers (bp)", y='', subtitle="center")

gg_mnase_c <- ggplot(mnase_at_threeUTRs_3ends, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + ylim(0, NA) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + theme(legend.position='none') +
  labs(title='', x="3'UTR 3'ends (bp)", y='', subtitle="3'end")

gg_mnase_a + gg_mnase_b + gg_mnase_c + plot_layout(ncol=3) & theme(aspect.ratio=1)

# 3e. ChIP-Seq
chipseq_at_threeUTRs_5ends   <- tidyMetaProfile(sites=threeUTRs_5ends,   forward=chipseq, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
chipseq_at_threeUTRs_centers <- tidyMetaProfile(sites=threeUTRs_centers, forward=chipseq, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)
chipseq_at_threeUTRs_3ends   <- tidyMetaProfile(sites=threeUTRs_3ends,   forward=chipseq, upstream=400, downstream=400, trimLower=0.01, trimUpper=0.99)

gg_chipseq_a <- ggplot(chipseq_at_threeUTRs_5ends, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + geom_hline(yintercept=0) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + ylim(0, NA) + theme(strip.text=element_text(face='bold'), legend.position='none') +
  labs(title='', x="3'UTR 5'ends (bp)", y='Average signal', subtitle="5'end")

gg_chipseq_b <- ggplot(chipseq_at_threeUTRs_centers, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + geom_hline(yintercept=0) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + ylim(0, NA) + theme(strip.text=element_text(face='bold'), legend.position='none') +
  labs(title='ChIP-Seq footprint', x="3'UTR centers (bp)", y='', subtitle="center")

gg_chipseq_c <- ggplot(chipseq_at_threeUTRs_3ends, aes(x=pos0, y=sense, col=sites)) +
  geom_line() + geom_hline(yintercept=0) + facet_grid(signal~., scales='free_y') +
  cowplot::theme_cowplot() + ylim(0, NA) + theme(strip.text=element_text(face='bold')) +
  labs(title='', x="3'UTR 3'ends (bp)", y='', subtitle="3'end")

gg_chipseq_a + gg_chipseq_b + gg_chipseq_c + plot_layout(ncol=3)




# 4. FOOTPRINT AT aTSSs themselves ####
# -------------------------------------
# get TCs that are antisense to a 3'UTR (again)
aTSSs_gr <- rowRanges(TCs) %>%
            subset(txType_TAIR10extended=='antisense_threeUTR')
# remove those overlapping anyting else
aTSSs_noOverlap_gr <- aTSSs_gr[countOverlaps(resize(swapRanges(aTSSs_gr), width=500, fix='center'), genes, ignore.strand=T) == 1]
# 312 left

# ChIP-seq footprint
chipseq_at_aTSSs <- tidyMetaProfile(sites=invertStrand(swapRanges(aTSSs_noOverlap_gr)), forward=chipseq, upstream=500, downstream=3000, trimLower=0.01, trimUpper=0.99)

# DNase & MNase footprint
dmnase_at_aTSSs <- tidyMetaProfile(sites=invertStrand(swapRanges(aTSSs_noOverlap_gr)), forward=dmnase, upstream=500, downstream=3000, trimLower=0.01, trimUpper=0.99)

# distance between aTSS and next mRNA TSS (that is antisense)
      # 1. get promoter TSSs
      promoter_TSSs_gr <- rowRanges(TCs) %>% subset(txType_TAIR10=='promoter')
      # 2. get the 354 aTSSs and flip their strand
      aTSSs_1bp_inversted_gr <- aTSSs_gr %>% swapRanges() %>% invertStrand()
      # 3. compute distance from aTSSs to next detected cage TSS that is a promoter
      next_promoter_TSS <- precede(aTSSs_1bp_inversted_gr, swapRanges(promoter_TSSs_gr))
      distances <- distance(aTSSs_1bp_inversted_gr, swapRanges(promoter_TSSs_gr[next_promoter_TSS])) %>% enframe(name=NULL, value='distance_to_next')
      # 4. plot
      gg_dist_next <- ggplot(distances, aes(x=distance_to_next)) + geom_density() +
             scale_x_log10() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
             labs(title='Distance to next promoter TSS', x='distance (bp)')

# distance between aTSS and 3'UTR end
      # 1. get overlaps
      overlaps <- findOverlaps(swapRanges(aTSSs_gr), threeUTRs_with_aTSSs, ignore.strand=T)
      # 2. re-order both
      aTSSs_ordered_gr <- aTSSs_gr[queryHits(overlaps)]
      threeUTRs_with_aTSSs_ordered_gr <- threeUTRs_with_aTSSs[subjectHits(overlaps)]
      # 3. make new 3'UTR GR with only their end
      threeUTRs_with_aTSSs_ordered_endOnly_gr <- threeUTRs_with_aTSSs_ordered_gr
      start(threeUTRs_with_aTSSs_ordered_endOnly_gr) <- ifelse(strand(threeUTRs_with_aTSSs_ordered_endOnly_gr)=='+',
                                                               end(threeUTRs_with_aTSSs_ordered_endOnly_gr),
                                                               start(threeUTRs_with_aTSSs_ordered_endOnly_gr))
      end(threeUTRs_with_aTSSs_ordered_endOnly_gr) <- start(threeUTRs_with_aTSSs_ordered_endOnly_gr)
      # 4. get distances (mind strandness!)
      dist_end <- distance(swapRanges(aTSSs_ordered_gr), threeUTRs_with_aTSSs_ordered_endOnly_gr, ignore.strand=T)
      # 4. plot
      gg_dist_end <- ggplot(enframe(dist_end, name=NULL, value='dist_end'), aes(x=dist_end)) +
             geom_density() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
             labs(x="Distance aTSS-3'UTR end (bp)") + xlim(-500, 3000)

# distance between aTSS and 3'UTR start
# 1. make new 3'UTR GR with only their end
threeUTRs_with_aTSSs_ordered_startOnly_gr <- threeUTRs_with_aTSSs_ordered_gr
start(threeUTRs_with_aTSSs_ordered_startOnly_gr) <- ifelse(strand(threeUTRs_with_aTSSs_ordered_startOnly_gr)=='+',
                                                           start(threeUTRs_with_aTSSs_ordered_startOnly_gr),
                                                           end(threeUTRs_with_aTSSs_ordered_startOnly_gr))
end(threeUTRs_with_aTSSs_ordered_startOnly_gr) <- start(threeUTRs_with_aTSSs_ordered_startOnly_gr)
# 2. get distances (mind strandness!)
dist_start <- distance(swapRanges(aTSSs_ordered_gr), threeUTRs_with_aTSSs_ordered_startOnly_gr, ignore.strand=T)
# 4. plot
gg_dist_start <- ggplot(enframe(dist_start, name=NULL, value='dist_start'), aes(x=dist_start)) +
  geom_density() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  labs(x="Distance aTSS-3'UTR start (bp)") + xlim(-500, 3000)


# plot all
gg_dmnase <- ggplot(dmnase_at_aTSSs, aes(x=pos0, y=sense)) +
  geom_line() + geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
  facet_wrap(~signal, scales='free', ncol=1) + ylim(0, NA) +
  cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold'), aspect.ratio=1, legend.position='none') +
  labs(title="Footprint at 3'UTR aTSSs", x="CAGE TSSs antisense to 3'UTR\n(peak, bp)", y='Normalized signal', subtitle="Direction= 5'-3' (as the 3'UTR)")

gg_dist <- ggplot(distances, aes(x=distance)) + geom_density() +
  cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  labs(x='distance (bp)', title="distance: aTSS to 3'UTR end")

gg_chip <- ggplot(chipseq_at_aTSSs, aes(x=pos0, y=sense)) +
       geom_line() + geom_hline(yintercept=0) + geom_vline(xintercept=0, lty=2) +
       facet_wrap(~signal, scales='free', ncol=2) + ylim(0, NA) +
       cowplot::theme_cowplot() + theme(strip.text=element_text(face='bold'), aspect.ratio=1, legend.position='none') +
       labs(title="Footprint at 3'UTR aTSSs", x="CAGE TSSs antisense to 3'UTR\n(peak, bp)", y='Normalized signal', subtitle="Direction= 5'-3' (as the 3'UTR)")

gg_chip + {gg_dist + gg_dmnase + plot_layout(ncol=1, widths=1, heights=c(1, 2))} + plot_layout(widths=c(2,1))


# DNase and distance to 3'UTR end on same scale
gg_mnase <- dmnase_at_aTSSs %>%
  subset(signal=='DNase_leaf') %>%
  ggplot(aes(x=pos0, y=sense)) + geom_line() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(x='position relative to aTSS (bp)', y='DNase signal')

gg_mnase + gg_dist_end + plot_layout(ncol=1)




# 5. HEATMAP OF DNaseI at aTSSs + line for 3'UTR end ####
# -------------------------------------------------------
# 5a. sort aTSSs by small to wider distance to 3'UTR end
  # get increasing order of distances to end
  increasing_order_to_end <- order(dist_end)
  # sort by this order
  aTSSs_ordered_by_dist_gr <- aTSSs_ordered_gr[increasing_order_to_end]
  dist_end_ordered_by_increasing <- dist_end[increasing_order_to_end]
  dist_end_ordered_by_increasing %>%
    enframe(name='aTSS', value='distance_to_end') %>%
    ggplot(aes(x=distance_to_end, y=-aTSS)) + geom_line() +
           geom_hline(yintercept=c(0, -354)) + geom_vline(xintercept=c(-500, 0, 1000, 2000, 3000)) +
           xlim(-500, 3000)
  # compute DNase heatmap
  dnase_mat <- wideMetaProfile(sites=swapRanges(aTSSs_ordered_by_dist_gr) %>% invertStrand(),
                               forward=dmnase$DNase_leaf, reverse=NULL,
                               upstream=500, downstream=3000)
  pheatmap(dnase_mat, cluster_rows=T, cluster_cols=F, show_rownames=F, show_colnames=F)

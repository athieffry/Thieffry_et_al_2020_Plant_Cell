#### 03c. All sort of footprints around TSSs
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
library(refGenome)
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

setwd('~/masked_path/03 - TSS analysis')




# 1. GET INPUT DATA AND FIX SEQINFO ####
# --------------------------------------
# myseqinfo (remove non-canonical chromosomes)
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
# TSS dataset (comparable ones after +/-100bp extension, quantification with CAGE wt T=0 and >= 1TPM in >= 2 libs)
cage_tss <- readRDS('~/masked_path/SE_CAGE_wt_comparable_TSSs_redone_17June2019.rds')
tair_tss <- readRDS('~/masked_path/SE_TAIR10_comparable_TSSs_redone_17June2019.rds')
araport_tss <- readRDS('~/masked_path/SE_ARAPORT11_comparable_TSSs_redone_17June2019.rds')
peat_tss <- readRDS('~/masked_path/SE_PEAT_comparable_TSSs_redone_17June2019.rds')
nanopare_tss <- readRDS('~/masked_path/SE_nanoPARE_comparable_TSSs_redone_17June2019.rds')



# 2. NORMALIZE INTO TPM AND KEEP 99 PERCENTILE ####
# -------------------------------------------------
if(FALSE){
          # ChIP-Seq PlantDHS.org
          plantDHS_buds_H3K27ac  <- import.bw('~/masked_path/plantdhs.org_Ath_buds_H3K27ac.bw')
          plantDHS_leaf_H3K27ac  <- import.bw('~/masked_path/plantdhs.org_Ath_leaf_H3K27ac.bw')
          plantDHS_buds_H3K27me3 <- import.bw('~/masked_path/plantdhs.org_Ath_buds_H3K27me3.bw')
          plantDHS_leaf_H3K27me3 <- import.bw('~/masked_path/plantdhs.org_Ath_leaf_H3K27me3.bw')
          plantDHS_buds_H3K4me1  <- import.bw('~/masked_path/plantdhs.org_Ath_buds_H3K4me1.bw')
          plantDHS_leaf_H3K4me1  <- import.bw('~/masked_path/plantdhs.org_Ath_leaf_H3K4me1.bw')
          # DNaseI PlantDHS.org
          plantDHS_dnase_flower <- import.bw('~/masked_path/plantdhs.org_Ath.flower.DNase.bw')
          plantDHS_dnase_leaf   <- import.bw('~/masked_path/plantdhs.org_Ath.leaf.DNase.bw')
          # Nucleosomes PlantDHS.org
          plantDHS_mnase_buds <- import.bw('~/masked_path/plantdhs.org_Ath_nucleosomes_buds_NPS.bw')
          plantDHS_mnase_leaf <- import.bw('~/masked_path/plantdhs.org_Ath_nucleosomes_leaf_NPS.bw')
              seqlevelsStyle(plantDHS_mnase_buds) <- seqlevelsStyle(myseqinfo)[1]
              seqlevelsStyle(plantDHS_mnase_leaf) <- seqlevelsStyle(myseqinfo)[1]
          # GRO-seq Nature Plants 2018 (Liu & Jacobsen, Nature 2018)
          GROseq_jac_p <- import.bw('~/masked_path/GROseq_col_mergedAveraged_R123456.plus.bw')
          GROseq_jac_m <- import.bw('~/masked_path/GROseq_col_mergedAveraged_R123456.minus.bw')
          # GRO-seq PNAS 2016 (Hetzel, PNAS 2016)
          GROseq_het_p <- import.bw('~/masked_path/GROseq_mergedAveraged_R12.plus.bw')
          GROseq_het_m <- import.bw('~/masked_path/GROseq_mergedAveraged_R12.minus.bw')
          # 5'GRO-seq PNAS 2016 (Hetzel, PNAS 2016)
          GROcap_het_p <- import.bw('~/masked_path/5GROcap.aln.unique.plus.bw')
          GROcap_het_m <- import.bw('~/masked_path/5GROcap.aln.unique.minus.bw')
          }
if(FALSE){
          # 2a) calculate TPM factor
                tpm_factor_plantDHS_buds_H3K27ac  <- sum(plantDHS_buds_H3K27ac$score)  / 1000000
                tpm_factor_plantDHS_leaf_H3K27ac  <- sum(plantDHS_leaf_H3K27ac$score)  / 1000000
                tpm_factor_plantDHS_buds_H3K27me3 <- sum(plantDHS_buds_H3K27me3$score) / 1000000
                tpm_factor_plantDHS_leaf_H3K27me3 <- sum(plantDHS_leaf_H3K27me3$score) / 1000000
                tpm_factor_plantDHS_buds_H3K4me1  <- sum(plantDHS_buds_H3K4me1$score)  / 1000000
                tpm_factor_plantDHS_leaf_H3K4me1  <- sum(plantDHS_leaf_H3K4me1$score)  / 1000000
                
                tpm_factor_plantDHS_dnase_flower <- sum(plantDHS_dnase_flower$score) / 1000000
                tpm_factor_plantDHS_dnase_leaf   <- sum(plantDHS_dnase_leaf$score) / 1000000
                
                tpm_factor_plantDHS_mnase_buds <- sum(plantDHS_mnase_buds$score) / 1000000
                tpm_factor_plantDHS_mnase_leaf <- sum(plantDHS_mnase_leaf$score) / 1000000
                
                tpm_factor_GROseq_jac <- sum(c(GROseq_jac_p$score, GROseq_jac_m$score)) / 1000000
                tpm_factor_GROseq_het <- sum(c(GROseq_het_p$score, GROseq_het_m$score)) / 1000000
                tpm_factor_GROcap_het <- sum(c(GROcap_het_p$score, GROcap_het_m$score)) / 1000000
          # 2b) normalize to TPM
                plantDHS_buds_H3K27ac$score  <- plantDHS_buds_H3K27ac$score  / tpm_factor_plantDHS_buds_H3K27ac
                plantDHS_leaf_H3K27ac$score  <- plantDHS_leaf_H3K27ac$score  / tpm_factor_plantDHS_leaf_H3K27ac
                plantDHS_buds_H3K27me3$score <- plantDHS_buds_H3K27me3$score / tpm_factor_plantDHS_buds_H3K27me3
                plantDHS_leaf_H3K27me3$score <- plantDHS_leaf_H3K27me3$score / tpm_factor_plantDHS_leaf_H3K27me3
                plantDHS_buds_H3K4me1$score  <- plantDHS_buds_H3K4me1$score  / tpm_factor_plantDHS_buds_H3K4me1
                plantDHS_leaf_H3K4me1$score  <- plantDHS_leaf_H3K4me1$score  / tpm_factor_plantDHS_leaf_H3K4me1
                
                plantDHS_dnase_flower$score <- plantDHS_dnase_flower$score / tpm_factor_plantDHS_dnase_flower
                plantDHS_dnase_leaf$score   <- plantDHS_dnase_leaf$score / tpm_factor_plantDHS_dnase_leaf
                
                plantDHS_mnase_buds$score <- plantDHS_mnase_buds$score / tpm_factor_plantDHS_mnase_buds
                plantDHS_mnase_leaf$score <- plantDHS_mnase_leaf$score / tpm_factor_plantDHS_mnase_leaf
                
                GROseq_jac_p$score <- GROseq_jac_p$score / tpm_factor_GROseq_jac
                GROseq_jac_m$score <- GROseq_jac_m$score / tpm_factor_GROseq_jac
                GROseq_het_p$score <- GROseq_het_p$score / tpm_factor_GROseq_het
                GROseq_het_m$score <- GROseq_het_m$score / tpm_factor_GROseq_het
                GROcap_het_p$score <- GROcap_het_p$score / tpm_factor_GROcap_het
                GROcap_het_m$score <- GROcap_het_m$score / tpm_factor_GROcap_het
          # 2c) calculate upper 99% percentile thresholds
                quantile_plantDHS_buds_H3K27ac  <- plantDHS_buds_H3K27ac$score  %>% quantile(0.99)
                quantile_plantDHS_leaf_H3K27ac  <- plantDHS_leaf_H3K27ac$score  %>% quantile(0.99)
                quantile_plantDHS_buds_H3K27me3 <- plantDHS_buds_H3K27me3$score %>% quantile(0.99)
                quantile_plantDHS_leaf_H3K27me3 <- plantDHS_leaf_H3K27me3$score %>% quantile(0.99)
                quantile_plantDHS_buds_H3K4me1  <- plantDHS_buds_H3K4me1$score  %>% quantile(0.99)
                quantile_plantDHS_leaf_H3K4me1  <- plantDHS_leaf_H3K4me1$score  %>% quantile(0.99)
                
                quantile_plantDHS_dnase_flower <- plantDHS_dnase_flower$score %>% quantile(0.99)
                quantile_plantDHS_dnase_leaf   <- plantDHS_dnase_leaf$score %>% quantile(0.99)
                
                quantile_plantDHS_mnase_buds <- plantDHS_mnase_buds$score %>% quantile(0.99)
                quantile_plantDHS_mnase_leaf <- plantDHS_mnase_leaf$score %>% quantile(0.99)
                
                quantile_GROseq_jac <- c(GROseq_jac_p$score, GROseq_jac_m$score) %>% quantile(0.99)
                quantile_GROseq_het <- c(GROseq_het_p$score, GROseq_het_m$score) %>% quantile(0.99)
                quantile_GROcap_het <- c(GROcap_het_p$score, GROcap_het_m$score) %>% quantile(0.99)
          # 2d) subset below threshold
                plantDHS_buds_H3K27ac  %<>% subset(score < quantile_plantDHS_buds_H3K27ac)
                plantDHS_leaf_H3K27ac  %<>% subset(score < quantile_plantDHS_leaf_H3K27ac)
                plantDHS_buds_H3K27me3 %<>% subset(score < quantile_plantDHS_buds_H3K27me3)
                plantDHS_leaf_H3K27me3 %<>% subset(score < quantile_plantDHS_leaf_H3K27me3)
                plantDHS_buds_H3K4me1  %<>% subset(score < quantile_plantDHS_buds_H3K4me1)
                plantDHS_leaf_H3K4me1  %<>% subset(score < quantile_plantDHS_leaf_H3K4me1)
                
                plantDHS_dnase_flower %<>% subset(score < quantile_plantDHS_dnase_flower)
                plantDHS_dnase_leaf   %<>% subset(score < quantile_plantDHS_dnase_leaf)
                
                plantDHS_mnase_buds %<>% subset(score < quantile_plantDHS_mnase_buds)
                plantDHS_mnase_leaf %<>% subset(score < quantile_plantDHS_mnase_leaf)
                
                GROseq_jac_p %<>% subset(score < quantile_GROseq_jac)
                GROseq_jac_m %<>% subset(score < quantile_GROseq_jac)
                GROseq_het_p %<>% subset(score < quantile_GROseq_het)
                GROseq_het_m %<>% subset(score < quantile_GROseq_het)
                GROcap_het_p %<>% subset(score < quantile_GROcap_het)
                GROcap_het_m %<>% subset(score < quantile_GROcap_het)
          
          # 2e) re-export all processed signals
                export.bw(plantDHS_buds_H3K27ac, '~/masked_path/plantDHS_buds_H3K27ac_TPM99pc.bw')
                export.bw(plantDHS_leaf_H3K27ac, '~/masked_path/plantDHS_leaf_H3K27ac_TPM99pc.bw')
                export.bw(plantDHS_buds_H3K27me3, '~/masked_path/plantDHS_buds_H3K27me3_TPM99pc.bw')
                export.bw(plantDHS_leaf_H3K27me3, '~/masked_path/plantDHS_leaf_H3K27me3_TPM99pc.bw')
                export.bw(plantDHS_buds_H3K4me1, '~/masked_path/plantDHS_buds_H3K4me1_TPM99pc.bw')
                export.bw(plantDHS_leaf_H3K4me1, '~/masked_path/plantDHS_leaf_H3K4me1_TPM99pc.bw')
                
                export.bw(plantDHS_dnase_flower, '~/masked_path/plantDHS_dnase_flower_TPM99pc.bw')
                export.bw(plantDHS_dnase_leaf , '~/masked_path/plantDHS_dnase_leaf_TPM99pc.bw')
                
                export.bw(plantDHS_mnase_buds, '~/masked_path/plantDHS_mnase_buds_TPM99pc.bw')
                export.bw(plantDHS_mnase_leaf, '~/masked_path/plantDHS_mnase_leaf_TPM99pc.bw')
                
                export.bw(GROseq_jac_p, '~/masked_path/GROseq_jacobsen_TPM99pc.plus.bw')
                export.bw(GROseq_jac_m, '~/masked_path/GROseq_jacobsen_TPM99pc.minus.bw')
                
                export.bw(GROseq_het_p, '~/masked_path/GROseq_hetzel_TPM99pc.plus.bw')
                export.bw(GROseq_het_m, '~/masked_path/GROseq_hetzel_TPM99pc.minus.bw')
                
                export.bw(GROcap_het_p, '~/masked_path/GROcap_hetzel_TPM99pc.plus.bw')
                export.bw(GROcap_het_m, '~/masked_path/GROcap_hetzel_TPM99pc.minus.bw')
                }

# ALL Histone marks (as provided by mivanov)
# bigwig files
chipseq <- list.files('~/masked_path/bigwigs', pattern='.bw', full.names=T) %>% BigWigFileList()
names(chipseq) <- list.files('~/masked_path/bigwigs', pattern='.bw') %>% str_remove('.bw')
# bigwig normalization factors (to RPM)
chipseq_rpm_factor <- read.table('~/masked_path/total_bedgraph_reads.txt', h=T, sep='\t') %>%
  as_tibble() %>%
  mutate('Sample'=str_remove(Sample, '.bedgraph'), 
         'Author'=str_split(Sample, '_', simplify=T)[,1],
         'Histone_mark'=str_split(Sample, '_', simplify=T)[,2],
         'RPM_factor'=total_reads_in_bedgraph / 1000000)
# MNaseI PlantDHS.org
mnase <- BigWigFile('~/masked_path/plantDHS_mnase_leaf_TPM99pc.bw')


mnase2 <- BigWigFileList(list.files(pattern='mnase', full.names=T))
names(mnase2) <- c('flower', 'leaf')

# DNaseI PlantDHS.org
dnase <- BigWigFile('~/masked_path/plantDHS_dnase_leaf_TPM99pc.bw')
# GRO-seq
gro_p <- BigWigFileList(list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*plus', full.names=T))
gro_m <- BigWigFileList(list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*minus', full.names=T))
names(gro_p) <- names(gro_m) <- list.files('~/masked_path/02 - Enhancers analyses/', pattern='GRO.*plus') %>% str_remove('_TPM99pc.plus.bw')
# make TSS GRangesList
tss_list <-list("CAGE"=rowRanges(cage_tss),
                "nanoPARE"=rowRanges(nanopare_tss),
                "PEAT"=rowRanges(peat_tss),
                "TAIR10"=rowRanges(tair_tss),
                "ARAPORT11"=rowRanges(araport_tss))



# 3. MNASE & GROSEQ FOOTPRINTS ####
# ---------------------------------
gro_signal <- lapply(tss_list, function(x) tidyMetaProfile(sites=x, forward=gro_p, reverse=gro_m, upstream=400, downstream=400))
mnase_signal <- lapply(tss_list, function(x) tidyMetaProfile(sites=x, forward=mnase, reverse=NULL, upstream=400, downstream=400))

mnase_test_signal <- lapply(tss_list[1:3], function(x) tidyMetaProfile(sites=x, forward=mnase2, reverse=NULL, upstream=400, downstream=400))
mnase_test_signal %<>% plyr::ldply() %>% as_tibble()

ggplot(mnase_test_signal, aes(x=pos1, y=sense, col=.id)) +
  geom_vline(xintercept=0, lty=2) +
  geom_line() + facet_wrap(~signal) +
  cowplot::theme_cowplot() +
  theme(aspect.ratio=1) +
  scale_color_brewer(name='TSS set', palette='Set1') +
  labs(x='Position relative to TSS/TC peak (bp)',
       y='Nucleosome normalized signal',
       title='MNase-I at TSSs')

gro_signal %<>% plyr::ldply() %>% as_tibble()
mnase_signal %<>% plyr::ldply() %>% as_tibble()

    # finding respective maxima
    gro_signal %>%
      select(-anti, -pos1) %>%
      group_by(.id, signal) %>%
      top_n(wt=sense, n=1) %>%
      as_tibble()
    
    mnase_signal %>%
      select(-pos1) %>%
      group_by(.id) %>%
      top_n(wt=sense, n=1)

# plot footprint
gro_signal_df <- gro_signal %>% gather(key="direction", value="score", sense, anti, factor_key=T)

rbind(gro_signal_df,
      mnase_signal %>% mutate('signal'='MNase', 'direction'='sense', 'score'=sense, 'sense'=NULL)) %>%
  ggplot(aes(x=pos0, y=score, lty=direction, col=.id)) +
         geom_vline(xintercept=0, lty=2) + geom_line(lwd=.8) +
         facet_wrap(~signal, scales='free_y', nrow=2) +
         cowplot::theme_cowplot() + theme(aspect.ratio=.75, legend.position='bottom', legend.direction=1) +
         scale_color_brewer(palette='Set1', name='TSS sources') +
         labs(title='Footprints at TSSs', x='Distance to TSSs (bp)', y='Normalized signal')

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
#  Extend TSSs +/- 100bp
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
CTSSs <- readRDS('~/masked_path/SE_CTSSs_1count_min3lib_TSSstory.rds')
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

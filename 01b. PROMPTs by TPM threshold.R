#### Arabidopsis
#### Axel Thieffry - revised May 2018
set.seed(42)
library(tidyverse)
library(magrittr)
library(stringr)
library(reshape2)
library(tidyquant)
library(CAGEfightR)
library(TeMPO)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(BiocParallel)
library(patchwork)
library(Gviz)
library(gt)
library(patchwork)
register(MulticoreParam(workers=4))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(Matrix)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
library(seqLogo)
library(JASPAR2016)
library(TFBSTools)
library(seqPattern)
odb <- org.At.tair.db
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/01b - PROMPTs by TPM threshold')



# 1. PREPARATION ####
# -------------------
### get seqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
### get sample names
sampleNames <- list.files(path='~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_raw_chrWchr/', pattern='_0_.*.plus.chrWchr.bw', full.names=F) %>% str_replace('.raw.plus.chrWchr.bw', '')
### colData
colData <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/colData_TSSstory.rds')
### get CTSSs (before they were support-filtered)
CTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CTSSs_1count_min3lib_TSSstory.rds')
# RNA-Seq BIGWIG FILES
rnaseq_forward <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123_fixed', pattern='Forward', full.names=T))
rnaseq_reverse <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123_fixed', pattern='Reverse', full.names=T))
rnaseq_names <- list.files(path='~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names
# CAGE TCs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')
    # make thick.star a IRanges
    rowRanges(TCs)$thick.start <- IRanges(start=rowRanges(TCs)$thick.start, width=1L)
    rowRanges(TCs)$thick.end <- rowRanges(TCs)$thick.width <- NULL
    rowRanges(TCs)


# 2. PROMPTs AS A FUNCTION OF TPM THRESHOLD #### # TO DEBUG
# ----------------------------------------------
### 4a. get "reverse" (that are "intergenic" on the other strand) and on canonical chromosomes
  reverse <- TCs[rowRanges(TCs)$txType=='reverse' & rowRanges(TCs)$txType_extended=='intergenic'] # 135
  seqlevels(reverse, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC')) # 116 reverse TCs
### 4b. remove reverse TCs that overlap any other genomic feature
  genes <- genes(txdb)
      seqlevelsStyle(genes) <- seqlevelsStyle(reverse)[1]
      seqlevels(genes, pruning.mode='coarse') <- seqlevels(reverse)
      seqinfo(genes) <- seqinfo(reverse)
  reverse <- reverse[countOverlaps(reverse, genes, ignore.strand=T) == 0] # 110
  if(FALSE){export.bed(rowRanges(reverse), '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/01b - PROMPTs by TPM threshold/all_reverse_PROMPT_candidates.bed')}
### 4c. quantify TPM per genotype (average replicates)
  rowRanges(reverse)$wt_TPM   <- subset(reverse, select=genotype=='wt')   %>% assay('TPM') %>% rowMeans()
  rowRanges(reverse)$hen2_TPM <- subset(reverse, select=genotype=='hen2') %>% assay('TPM') %>% rowMeans()
  rowRanges(reverse)$rrp4_TPM <- subset(reverse, select=genotype=='rrp4') %>% assay('TPM') %>% rowMeans()
### 4d. get TCs that were found DE UP in any mutant
  hen2_de_up <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotypehen2.rds') %>% subset(logFC >= 0)
  rrp4_de_up <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotyperrp4.rds') %>% subset(logFC >= 0)
  # keep the reverse ones only
  hen2_reverse_de <- subset(hen2_de_up, txType_TAIR10 == 'reverse') # 63
  rrp4_reverse_de <- subset(rrp4_de_up, txType_TAIR10 == 'reverse') # 94
  # make GR with PROMPTs only
  prompts_ids <- unique(c(hen2_reverse_de$id, rrp4_reverse_de$id))
  prompts_gr <- TCs[names(TCs) %in% prompts_ids] %>% rowRanges()
### 4d. barplot by TPM threshold
  thresholds <- c(0.5, 0.75, 1, 1.25, 1.5, 2, 5)
  # WT: number of PROMPTs remaining with threshold
  wt_df <- data.frame('isDE'=FALSE,
                      'thresholds'=as.factor(thresholds),
                      'genotype'='wt',
                      'nTCs'=sapply(thresholds, function(x) sum(rowRanges(reverse)$wt_TPM >= x)))
  # hen2: number of PROMPTs remaining with threshold
  hen2_df <- lapply(thresholds, function(x) names(subset(rowRanges(reverse), hen2_TPM >= x))) %>%
    set_names(thresholds) %>%
    sapply(function(x) table(x %in% hen2_reverse_de$id)) %>%
    melt() %>%
    set_colnames(c('isDE', 'thresholds', 'nTCs')) %>%
    mutate('genotype'='hen2')
  # rrp4: number of PROMPTs remaining with threshold
  rrp4_df <- lapply(thresholds, function(x) names(subset(rowRanges(reverse), rrp4_TPM >= x))) %>%
    set_names(thresholds) %>%
    sapply(function(x) table(x %in% rrp4_reverse_de$id)) %>%
    melt() %>%
    set_colnames(c('isDE', 'thresholds', 'nTCs')) %>%
    mutate('genotype'='rrp4')
  
  rbind(wt_df, hen2_df, rrp4_df) %>%
  ggplot(aes(x=as.factor(thresholds), y=nTCs, fill=interaction(isDE, genotype))) +
         geom_bar(stat='identity', col='black', position=position_dodge()) +
         geom_text(aes(label=nTCs, y=nTCs-(nTCs/15)), position=position_dodge(width=0.9), angle=90) +
         cowplot::theme_cowplot() + theme(legend.position='bottom', legend.direction=1) +
         labs(x='TPM threshold', y='PROMPT-like TCs') +
         scale_fill_brewer(palette='Paired')


  
# 5. AWTAA / 5'SS signal at PROMPTs ####
# compare PROMPTs vs classical mRNA TSSs AWTAA signal upstream of them
# --------------------------------------------------------------------
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome <- rep('TAIR10', 7)
seqinfo(genome)

# get classical mRNA CAGE TSSs
promoters_300300_seq <- rowRanges(TCs) %>%
  subset(txType_TAIR10=='promoter') %>%
  swapRanges(inputColumn='thick.start') %>%
  promoters(300, 300) %>%
  remove_out_of_bound() %>%
  getSeq(genome, .)

# get PROMPT TSSs
prompts_gr <- rowRanges(TCs) %>%
  subset(txType_TAIR10=='reverse') %>%
  subset(genotypehen2 ==1 | genotyperrp4 ==1)

prompts_300300_seq <- prompts_gr %>%
  swapRanges(inputColumn='thick.start') %>%
  promoters(300, 300) %>%
  remove_out_of_bound() %>%
  getSeq(genome, .)

# plot AWTAAA motifs (.png)
plotPatternDensityMap(regionsSeq=promoters_300300_seq,
                      patterns=c("AWTAAA"),
                      flankUp=300, flankDown=300, 
                      color='cyan', outFile='promoter_TSS_300bp.awtaaa', 
                      plotScale=F, labelCol='black',  addPatternLabel=F,
                      xLabel='', useMulticore=T, nrCores=4)

plotPatternDensityMap(regionsSeq=prompts_300300_seq,
                      patterns=c("AWTAAA"),
                      flankUp=300, flankDown=300,
                      color='cyan', outFile='PROMPTs_TSS_300bp.awtaaa', 
                      plotScale=F, addPatternLabel=F,
                      useMulticore=T, nrCores=4, 
                      xLabel='', bandWidth=c(1, 1))


par(mfrow=c(2, 1), pty='s')
seqPattern::plotPatternOccurrenceAverage(regionsSeq=promoters_300300_seq, patterns=c("AWTAA"), flankUp=300, flankDown=300, xLabel='promoter TSSs (N=17,925)', useMulticore=T, nrCores=4, smoothingWindow=5, color='blue')
seqPattern::plotPatternOccurrenceAverage(regionsSeq=prompts_300300_seq, patterns=c("AWTAA"), flankUp=300, flankDown=300, xLabel='PROMPT TSSs (N=96)', useMulticore=T, nrCores=4, smoothingWindow=5, color='red')
dev.off()



# 5. RNA-seq at PROMPTs ####
# --------------------------
# make PROMPT regions (400bp upstream of TSS and reverse strand)
  prompts_regions_gr <- promoters(txdb, upstream=400, downstream=0) %>% trim() %>% invertStrand()
    export.bed(prompts_regions_gr, '~/Desktop/prompts_regions_gr.bed')
# keep only PROMPTs regions with an up-regulated PROMPT-TC in any mutant
  prompts_TCs_gr <- rowRanges(TCs)[names(rowRanges(TCs)) %in% unique(c(hen2_reverse_de$id, rrp4_reverse_de$id))]
    export.bed(prompts_TCs_gr, '~/Desktop/prompts_TCs_gr.bed')
  prompts_regions_gr <- subsetByOverlaps(prompts_regions_gr, prompts_TCs_gr) %>% unique() # 111
  resize(prompts_regions_gr, width=1, fix='start') %>% export.bed('~/Desktop/prompts_TCs_1bp_gr.bed')
# quantify RNA-seq at those PROMPT regions (signal is SENSE since I inverted strand initially)
  rnaseq_signal <- mapply(function(x, y) wideMetaProfile(sites=resize(prompts_regions_gr, width=1, fix='start'), forward=x, reverse=y, upstream=1, downstream=400), rnaseq_forward, rnaseq_reverse, SIMPLIFY=F)
  rnaseq_at_prompts_df <- tibble('wt'=rowSums(rnaseq_signal$WT$sense),
                                 'lsm8'=rowSums(rnaseq_signal$LSM8$sense),
                                 'rrp4'=rowSums(rnaseq_signal$RRP4$sense),
                                 'lsm8_rrp4'=rowSums(rnaseq_signal$DM$sense))
  rnaseq_at_prompts_df
# plot
  ggplot(rnaseq_at_prompts_df, aes(x=wt+1, y=rrp4+1)) +
         geom_point() + geom_abline(intercept=0, slope=1) +
         scale_x_log10() + scale_y_log10() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(title='RNA-seq support of PROMPTs', x='', y='rrp4-2\naverage FPM') +
  ggplot(rnaseq_at_prompts_df, aes(x=wt+1, y=lsm8+1)) +
         geom_point() + geom_abline(intercept=0, slope=1) +
         scale_x_log10() + scale_y_log10() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(x='', y='lsm8-2\naverage FPM') +
  ggplot(rnaseq_at_prompts_df, aes(x=wt+1, y=lsm8_rrp4+1)) +
         geom_point() + geom_abline(intercept=0, slope=1) +
         scale_x_log10() + scale_y_log10() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(y='lsm8-2 rrp4-2\naverage FPM', x='wt average FPM') +
  plot_layout(ncol=1)
  
  
  
# 6. splice site in PROMPT regions ####
# -------------------------------------
# read Nucleotide composition of Arabidopsis introns at the 5' splice site
# from : https://www.arabidopsis.org/info/splice_site.pdf
ss5_pfm <- readJASPARMatrix(fn='Athal_intron_at_5ss.mat')
ss3_pfm <- readJASPARMatrix(fn='Athal_intron_at_3ss.mat')
ss5_pwm <- toPWM(ss5_pfm)
ss3_pwm <- toPWM(ss3_pfm)

# make seqLogo of 3SS & 5SS
proportion <- function(x){
                          rs <- sum(x);
                          return(x / rs);
}

ss5_pfm_prop <- apply(ss5_pfm@profileMatrix, 2, proportion) %>% makePWM()
ss3_pfm_prop <- apply(ss3_pfm@profileMatrix, 2, proportion) %>% makePWM()

par(pty='s')
seqLogo::seqLogo(ss5_pfm_prop)
seqLogo::seqLogo(ss3_pfm_prop)

# plot Densities
plotMotifDensityMap(regionsSeq=promoters_300300_seq, plotScale=F, addReferenceLine=T,
                    motifPWM=ss5_pwm@profileMatrix, color='blue',
                    minScore='80%', flankUp=300, flankDown=300, outFile='promoters_5SS.png')

plotMotifDensityMap(regionsSeq=prompts_300300_seq, plotScale=F, addReferenceLine=T,
                    motifPWM=ss5_pwm@profileMatrix, color='blue',
                    minScore='80%', flankUp=300, flankDown=300, outFile='prompts_5SS.png')

plotMotifDensityMap(regionsSeq=promoters_300300_seq, plotScale=F, addReferenceLine=T,
                    motifPWM=ss3_pwm@profileMatrix, color='red',
                    minScore='80%', flankUp=300, flankDown=300, outFile='promoters_3SS.png')

plotMotifDensityMap(regionsSeq=prompts_300300_seq, plotScale=F, addReferenceLine=T,
                    motifPWM=ss3_pwm@profileMatrix, color='red',
                    minScore='80%', flankUp=300, flankDown=300, outFile='prompts_3SS.png')



# 7. cumulative as in Ntini et al., NSBM 2013 ####
# percentage of cases having >= 1 match after X bp from TC, moving along & cumulative
# ------------------------------------------------
# do walking promoter & PROMPT sequence from CAGE TC, in PolII direction
promoters_500_walkingDir_seq <- rowRanges(TCs) %>%
  subset(txType_TAIR10=='promoter') %>%
  swapRanges(inputColumn='thick.start') %>%
  promoters(upstream=0, downstream=500) %>%
  remove_out_of_bound() %>%
  getSeq(genome, .)

prompts_500_walkingDir_seq <- prompts_gr %>%
  swapRanges(inputColumn='thick.start') %>%
  promoters(upstream=0, downstream=500) %>%
  remove_out_of_bound() %>%
  getSeq(genome, .)

# scan scores for AWTAA, 5'SS & 3'SS
AWTAA_promoters <- getPatternOccurrenceList(regionsSeq=promoters_500_walkingDir_seq, patterns='AWTAA', useMulticore=T, nrCores=4)$AWTAA %>% as_tibble() %>% select(-value)
AWTAA_prompts <- getPatternOccurrenceList(regionsSeq=prompts_500_walkingDir_seq, patterns='AWTAA', useMulticore=T, nrCores=4)$AWTAA %>% as_tibble() %>% select(-value)

score_5SS_promoters <- motifScanScores(regionsSeq=promoters_500_walkingDir_seq, motifPWM=ss5_pwm@profileMatrix, asPercentage=T)
score_5SS_prompts <- motifScanScores(regionsSeq=prompts_500_walkingDir_seq, motifPWM=ss5_pwm@profileMatrix, asPercentage=T)

score_3SS_promoters <- motifScanScores(regionsSeq=promoters_500_walkingDir_seq, motifPWM=ss3_pwm@profileMatrix, asPercentage=T)
score_3SS_prompts <- motifScanScores(regionsSeq=prompts_500_walkingDir_seq, motifPWM=ss3_pwm@profileMatrix, asPercentage=T)

# binarize matrices depending on score threshold
score_threshold = 80

binarize_by_threshold <- function(x, threshold=80) {
  x[x < threshold] <- 0
  x[x >= threshold] <- 1
  x
}

score_5SS_promoters_bin <- binarize_by_threshold(score_5SS_promoters)
score_5SS_prompts_bin <- binarize_by_threshold(score_5SS_prompts)
score_3SS_promoters_bin <- binarize_by_threshold(score_3SS_promoters)
score_3SS_prompts_bin <- binarize_by_threshold(score_3SS_prompts)

# get position of first match
AWTAA_promoters_fm <- AWTAA_promoters %>% group_by(sequence) %>% filter(position==min(position)) %>% ungroup()
AWTAA_prompts_fm <- AWTAA_prompts %>% group_by(sequence) %>% filter(position==min(position)) %>% ungroup()
score_5SS_promoters_fm <- apply(score_5SS_promoters_bin, 1, function(x) dplyr::first(which(x==1))) %>% as.vector()
score_5SS_prompts_fm <- apply(score_5SS_prompts_bin, 1, function(x) dplyr::first(which(x==1))) %>% as.vector()
score_3SS_promoters_fm <- apply(score_3SS_promoters_bin, 1, function(x) dplyr::first(which(x==1))) %>% as.vector()
score_3SS_prompts_fm <- apply(score_3SS_prompts_bin, 1, function(x) dplyr::first(which(x==1))) %>% as.vector()

# K-S tests
    # AWTAA
        # assess distributions
        ggplot() +
          geom_density(data=AWTAA_prompts_fm, aes(x=position), col='red') +
          geom_density(data=AWTAA_promoters_fm, aes(x=position), col='blue') +
          geom_text(label='TC source: AWTAA at promoters', col='blue', aes(x=100, y=0.003), hjust=0) +
          geom_text(label='TC source: AWTAA at PROMPTs', col='red', aes(x=250, y=0.0025), hjust=0) +
          geom_text(label='K-S test:\nalternative H: CDF of PROMPTs lies above that of promoters):\nD=0.11289\np-value=0.1699', aes(x=30, y=0.0004), hjust=0) +
          labs(x='Position of AWTAA match from TC source (bp)', title='AWTAA distribution') +
          cowplot::theme_cowplot()
        # KS test
        ks.test(x=AWTAA_prompts_fm$position, y=AWTAA_promoters_fm$position, alternative='greater')
    
    # 3SS
        # assess distributions
        data.frame('source'=c(rep('promoters', length(score_3SS_promoters_fm)),
                              rep('prompts', length(score_3SS_prompts_fm))),
                   'position'=c(score_3SS_promoters_fm, score_3SS_prompts_fm)) %>%
          as_tibble() %>%
          ggplot(aes(x=position, col=source, fill=source)) +
                 geom_density(alpha=.5) +
                 labs(x='Position of 3SS match from TC source (bp)',
                      title='3SS distribution') +
                 cowplot::theme_cowplot()
        # KS test
        ks.test(x=score_3SS_prompts_fm, y=score_3SS_promoters_fm, alternative='greater')

# remove NA's
score_5SS_promoters_fm[is.na(score_5SS_promoters_fm)] <- 0
score_5SS_prompts_fm[is.na(score_5SS_prompts_fm)] <- 0
score_3SS_promoters_fm[is.na(score_3SS_promoters_fm)] <- 0
score_3SS_prompts_fm[is.na(score_3SS_prompts_fm)] <- 0

# make dataframe with 0's and 1's
AWTAA_promoters_fm <- plyr::ldply(AWTAA_promoters_fm$position, function(x) c(rep(0, x), rep(1, 496-x)))
AWTAA_promoters_fmdf <- apply(AWTAA_promoters_fm, 2, function(x) sum(x)/17925) # give number manually as some promoters disappeared (no match)

AWTAA_prompts_fm <- plyr::ldply(AWTAA_prompts_fm$position, function(x) c(rep(0, x), rep(1, 496-x)))
AWTAA_prompts_fmdf <- apply(AWTAA_prompts_fm, 2, function(x) sum(x)/96) # give number manually as some prompts disappeared (no match)

score_5SS_promoters_fm <- plyr::ldply(score_5SS_promoters_fm, function(x) c(rep(0, x), rep(1, ncol(score_5SS_promoters_bin)-x)))
score_5SS_promoters_fmdf <- apply(score_5SS_promoters_fm, 2, function(x) sum(x)/nrow(score_5SS_promoters_bin))

score_5SS_prompts_fm <- plyr::ldply(score_5SS_prompts_fm, function(x) c(rep(0, x), rep(1, ncol(score_5SS_prompts_bin)-x)))
score_5SS_prompts_fmdf <- apply(score_5SS_prompts_fm, 2, function(x) sum(x)/nrow(score_5SS_prompts_bin))

score_3SS_promoters_fm <- plyr::ldply(score_3SS_promoters_fm, function(x) c(rep(0, x), rep(1, ncol(score_3SS_promoters_bin)-x)))
score_3SS_promoters_fmdf <- apply(score_3SS_promoters_fm, 2, function(x) sum(x)/nrow(score_3SS_promoters_bin))

score_3SS_prompts_fm <- plyr::ldply(score_3SS_prompts_fm, function(x) c(rep(0, x), rep(1, ncol(score_3SS_prompts_bin)-x)))
score_3SS_prompts_fmdf <- apply(score_3SS_prompts_fm, 2, function(x) sum(x)/nrow(score_3SS_prompts_bin))

# put together and plot
relfreq_5SS <- data.frame('position'=1:length(score_5SS_promoters_fmdf),
           'promoters'=score_5SS_promoters_fmdf,
           'PROMPTs'=score_5SS_prompts_fmdf) %>%
  melt(id.vars='position', variable.name='set') %>%
  as_tibble()

relfreq_3SS <- data.frame('position'=1:length(score_3SS_promoters_fmdf),
                          'promoters'=score_3SS_promoters_fmdf,
                          'PROMPTs'=score_3SS_prompts_fmdf) %>%
  melt(id.vars='position', variable.name='set') %>%
  as_tibble()

relfreq_AWTAA <- data.frame('position'=1:length(AWTAA_promoters_fmdf),
                            'promoters'=AWTAA_promoters_fmdf,
                            'prompts'=AWTAA_prompts_fmdf) %>%
  melt(id.vars='position', variable.name='set') %>%
  as_tibble()

# plot
ggplot(relfreq_5SS, aes(x=position, color=set, y=value)) +
       geom_line() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
       scale_color_brewer(palette='Dark2', name='') +
       labs(title="5' Splice Site", x='Position relative to\nmRNA or PROMPT TSSs (nt)',
            y='Fraction of regions\nwith pattern (cumulative)') +
ggplot(relfreq_3SS, aes(x=position, color=set, y=value)) +
       geom_line() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
       scale_color_brewer(palette='Dark2', name='') +
       labs(title="3' SS", x='Position relative to\nmRNA or PROMPT TSSs (nt)',
            y='') +
ggplot(relfreq_AWTAA, aes(x=position, color=set, y=value)) +
       geom_line() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       scale_color_brewer(palette='Dark2', name='') +
       labs(title='AWTAA', x='Position relative to\nmRNA or PROMPT TSSs (nt)',
            y='') +
plot_layout(ncol=3)


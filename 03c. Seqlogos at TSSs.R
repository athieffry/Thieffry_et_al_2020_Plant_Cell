#### 03d. Seqlogos at TSSs
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

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/03 - TSS analysis')


# 1. GET ALL INPUT DATA ####
# --------------------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
    # remove non-canonical chromosomes
    seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))

# TSS dataset (comparable ones after +/-100bp extension, quantification with CAGE wt T=0 and >= 1TPM in >= 2 libs)
cage_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CAGE_wt_comparable_TSSs_redone_17June2019.rds')
tair_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TAIR10_comparable_TSSs_redone_17June2019.rds')
araport_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_ARAPORT11_comparable_TSSs_redone_17June2019.rds')
peat_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_PEAT_comparable_TSSs_redone_17June2019.rds')
nanopare_tss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_nanoPARE_comparable_TSSs_redone_17June2019.rds')




# 2. SEQLOGO AT COMPARABLE TSSs ####
# ----------------------------------
#  Read TAIR10 genome
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack  

# make promoter sequence regions for (shift everything by -1 to account for the safe-trimming)
upstream=40
downstream=10

# CAGE WT ctrl
cage_seq <- cage_tss %>%
  rowRanges() %>%
  flank(start=T, width=1) %>% # shift -1 bp: 'shift' is not used since it acts independently of strandedness
  promoters(upstream=upstream, downstream=downstream) %>% # make promoters GR regions
  getSeq(genome, .) # fetch DNAStringSet

cage_pfm <- cage_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='CAGE wt')

cage_pfm@background <- colSums(letterFrequency(cage_seq, DNA_BASES)) / sum(colSums(letterFrequency(cage_seq, DNA_BASES)))

# PEAT
peat_seq <- peat_tss %>%
  rowRanges() %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

peat_pfm <- peat_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='PEAT')

peat_pfm@background <- colSums(letterFrequency(peat_seq, DNA_BASES)) / sum(colSums(letterFrequency(peat_seq, DNA_BASES)))

# nanoPARE
nanopare_seq <- nanopare_tss %>%
  rowRanges() %>%
  #flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

nanopare_pfm <- nanopare_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='nanoPARE')

nanopare_pfm@background <- colSums(letterFrequency(nanopare_seq, DNA_BASES)) / sum(colSums(letterFrequency(nanopare_seq, DNA_BASES)))

# TAIR10
tair_seq <- tair_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

tair_pfm <- tair_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='TAIR10')

tair_pfm@background <- colSums(letterFrequency(tair_seq, DNA_BASES)) / sum(colSums(letterFrequency(tair_seq, DNA_BASES)))

# ARAPORT11
araport_seq <- araport_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

araport_pfm <- araport_seq %>%
  consensusMatrix(as.prob=T) %>%
  subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
  as.matrix() %>%
  new('pfm', mat=., name='ARAPORT11')

araport_pfm@background <- colSums(letterFrequency(araport_seq, DNA_BASES)) / sum(colSums(letterFrequency(araport_seq, DNA_BASES)))

# plot logos
mylist <- list('CAGE wt (21,221)'=cage_pfm$mat,
               'PEAT (8,952)'=peat_pfm$mat,
               'nanoPARE (16,310)'=nanopare_pfm$mat,
               'TAIR10 (22,394)'=tair_pfm$mat,
               'ARAPORT11 (24,877)'=araport_pfm$mat)

plusone <- seq(-upstream, downstream-1, 10)
plusone[plusone==0] <- '+1'

ggseqlogo(data=mylist, ncol=1) +
         geom_vline(xintercept=c(1, 11, 21, 31, 41), lty=2) +
         scale_x_continuous(breaks=seq(1, upstream + downstream, 10), labels=plusone) +
         theme(aspect.ratio=.25) +
         labs(x='Distance to TSSs (bp)')

# Special one-shot seqlogo for the 96 exosome-sensitive PROMPT TCs.
# (a bit of back-and-forth because some SE have DE data, some have swapRanges working, whatever.)
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')

prompt_TCs <- rowRanges(TCs) %>%
              subset(txType_TAIR10 == 'reverse') %>%
              subset(genotypehen2 == 1 | genotyperrp4 == 1) %>%
              names()

TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_TPM1_min3lib_TSSstory.rds')

prompt_TCs <- rowRanges(TCs) %>%
              subset(names %in% prompt_TCs) %>%
              swapRanges()

prompt_seq <- prompt_TCs %>%
              promoters(upstream=upstream, downstream=downstream) %>%
              getSeq(genome, .)

prompt_pfm <- prompt_seq %>%
              consensusMatrix(as.prob=T) %>%
              subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
              as.matrix() %>%
              new('pfm', mat=., name='PROMPT TCs')

prompt_pfm@background <- colSums(letterFrequency(prompt_seq, DNA_BASES)) / sum(colSums(letterFrequency(prompt_seq, DNA_BASES)))

ggseqlogo(data=list('PROMPT TCs (96)'=prompt_pfm$mat), ncol=1) +
  geom_vline(xintercept=40, lty=2) +
  scale_x_continuous(breaks=seq(0, upstream + downstream, 10), labels=c(plusone, "+10")) +
  theme(aspect.ratio=.4) +
  labs(x='Distance to TC peaks (bp)')
upstream
downstream

# 4. DINUCLEOTIDE AND TRINUCLEOTIDE FREQUENCIES ####
# --------------------------------------------------
# make nucs, dinucs and trinucs
nucs <- c('A', 'T', 'C', 'G')
dinucs <- expand.grid(nucs, nucs) %>% apply(1, function(x) paste0(x[1], x[2]))
trinucs <- c('AAA', 'TTT', 'CCC', 'GGG')

# calculate frequencies
cage_seq_nucs <- lapply(c(nucs, dinucs, trinucs), function(x) vmatchPattern(pattern=x, subject=cage_seq, fixed=F) %>%
                      as.data.frame() %>%
                      group_by(start) %>%
                      summarise('ratio'=n()/length(cage_seq)) %>%
                      mutate('pattern'=x)) %>%
                      do.call(rbind, .)

tair_seq_nucs <- lapply(c(nucs, dinucs, trinucs), function(x) vmatchPattern(pattern=x, subject=tair_seq, fixed=F) %>%
                      as.data.frame() %>%
                      group_by(start) %>%
                      summarise('ratio'=n()/length(tair_seq)) %>%
                      mutate('pattern'=x)) %>%
                      do.call(rbind, .)

araport_seq_nucs <- lapply(c(nucs, dinucs, trinucs), function(x) vmatchPattern(pattern=x, subject=araport_seq, fixed=F) %>%
                      as.data.frame() %>%
                      group_by(start) %>%
                      summarise('ratio'=n()/length(araport_seq)) %>%
                      mutate('pattern'=x)) %>%
                      do.call(rbind, .)

nanopare_seq_nucs <- lapply(c(nucs, dinucs, trinucs), function(x) vmatchPattern(pattern=x, subject=nanopare_seq, fixed=F) %>%
                      as.data.frame() %>%
                      group_by(start) %>%
                      summarise('ratio'=n()/length(nanopare_seq)) %>%
                      mutate('pattern'=x)) %>%
                      do.call(rbind, .)

peat_seq_nucs <- lapply(c(nucs, dinucs, trinucs), function(x) vmatchPattern(pattern=x, subject=peat_seq, fixed=F) %>%
                      as.data.frame() %>%
                      group_by(start) %>%
                      summarise('ratio'=n()/length(peat_seq)) %>%
                      mutate('pattern'=x)) %>%
                      do.call(rbind, .)

cage_seq_nucs$set <- "CAGE wt (21,221)" %>% as.factor()
tair_seq_nucs$set <- "TAIR10 (22,394)" %>% as.factor()
araport_seq_nucs$set <- "ARAPORT11 (24,877)" %>% as.factor()
nanopare_seq_nucs$set <- "nanoPARE (16,309)" %>% as.factor()
peat_seq_nucs$set <- "PEAT (8,952)" %>% as.factor()

nucs_df <- rbind(cage_seq_nucs, tair_seq_nucs, araport_seq_nucs, nanopare_seq_nucs, peat_seq_nucs)
nucs_df$start <- nucs_df$start-40
nucs_df$pattern %<>% factor()
nucs_df$length <- ifelse(nchar(as.character(nucs_df$pattern))==1, 'nucs', 
                         ifelse(nchar(as.character(nucs_df$pattern))==2, 'dinucs', 'trinucs')) %>%
                                factor()

# nucleotides
gg_nuc <- ggplot(subset(nucs_df, length=='nucs'), aes(x=start, y=ratio, col=pattern)) +
       geom_vline(xintercept=1, lty=2) + geom_line(size=.6) +
       facet_wrap(~set, ncol=1) +
       cowplot::theme_cowplot() + theme(strip.text.x=element_text(face='bold')) +
       scale_color_manual(values=brewer.pal(name='Set1', n=5)[c(3, 2, 5, 1)], name='Nucleotides') +
       labs(titles='Nucleotide frequencies at comparable TSSs', x='Distance to TSSs (bp)', y='Nucleotide ratio (%)')

# dinucleotides
gg_dinuc <- ggplot(subset(nucs_df, length=='dinucs'), aes(x=start, y=ratio, col=set)) +
       geom_vline(xintercept=1, lty=2) + geom_line() +
       facet_wrap(~pattern, scales='free_y') +
       cowplot::theme_cowplot() + theme(strip.text.x=element_text(face='bold'), legend.position='none') +
       labs(title='Dinucleotide frequencies at comparable TSSs',
            x='Distance to TSSs (bp)', y='Dinucleotide ratio (%)') +
       scale_color_manual(values=brewer.pal(name='Paired', n=5), name='TSSs dataset')

# trinucleotids
gg_trinuc <- ggplot(subset(nucs_df, length=='trinucs'), aes(x=start, y=ratio, col=set)) +
       geom_vline(xintercept=1, lty=2) + geom_line(lwd=0.8) +
       cowplot::theme_cowplot() + theme(strip.text.x=element_text(face='bold'), legend.position='bottom', legend.direction=2) +
       facet_wrap(~pattern) +
       scale_color_manual(values=brewer.pal(name='Paired', n=5), name='') +
       labs(title='Trinucleotide frequencies at comparable TSSs',
            x='Distance to TSSs (bp)', y='Trinucleotide ratio (%)')

Rmisc::multiplot(gg_nuc, gg_trinuc, gg_dinuc, layout=matrix(c(1,2,3,3), nrow=2, ncol=2, byrow=T))



# 5. SCAN FOR TATA BOX BY CAGE QUANTILE ####
# ------------------------------------------
# redo promoter sequences with different window
upstream=100
downstream=50

cage_seq <- cage_tss %>%
  rowRanges() %>%
  flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

tair_seq <- tair_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

araport_seq <- araport_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

nanopare_seq <- nanopare_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

peat_seq <- peat_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

# classify TCs TPM expression by rank
    # recalculate score
    rowRanges(cage_tss)$score <- assay(cage_tss, 'TPM') %>% rowMeans()
    # divide 
    quant_cage_tss <- rowRanges(cage_tss)$score %>% quantile(probs=seq(0, 1, 0.2))
    rowRanges(cage_tss)$TPMquantile <- fancycut(x=rowRanges(cage_tss)$score,
                                               '[0-20%)'='[0, 2.2814050)',
                                               '[20-40%)'='[2.2814050, 5.6592173)',
                                               '[40-60%)'='[5.6592173, 14.1196835)',
                                               '[60-80%)'='[14.1196835, 38.6368844)',
                                               '[80-100%]'='[38.6368844, 11935.41541874]',
                                               out.as.factor=FALSE)
    # sanity check
    table(rowRanges(cage_tss)$TPMquantile, useNA='always') # no NA: all ok

# get all PFM POL matrices from JASPAR 2018
matrices <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/03 - TSS analysis/JASPAR2018_POL_matrices/', full.names=T)
matrices_names <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/03 - TSS analysis/JASPAR2018_POL_matrices/') %>% str_replace('.jaspar', '')
matrices <- lapply(matrices, readJASPARMatrix)
names(matrices) <- matrices_names
# make PWM
matrices <- lapply(matrices, toPWM)
# calculate matrices length (need to remove matric_lenght/2 from x axis, later on)
matrices_lengths <- lapply(matrices, function(x) ncol(x@profileMatrix))
# compute scores
scores_cage <- lapply(matrices, function(x) motifScanScores(regionsSeq=cage_seq, asPercentage=T, motifPWM=x@profileMatrix))
scores_tair <- lapply(matrices, function(x) motifScanScores(regionsSeq=tair_seq, asPercentage=T, motifPWM=x@profileMatrix))
scores_araport <- lapply(matrices, function(x) motifScanScores(regionsSeq=araport_seq, asPercentage=T, motifPWM=x@profileMatrix))
scores_nanopare <- lapply(matrices, function(x) motifScanScores(regionsSeq=nanopare_seq, asPercentage=T, motifPWM=x@profileMatrix))
scores_peat <- lapply(matrices, function(x) motifScanScores(regionsSeq=peat_seq, asPercentage=T, motifPWM=x@profileMatrix))
warnings()
# process scores
scores_cage_df <- mapply(function(x, y)
  data.frame('position'= -(upstream-y/2) : (downstream-y/2),
             'CAGE wt'=colMeans(x)) %>%
    melt(id.vars='position', variable.name='TSSs'),
  scores_cage, matrices_lengths, SIMPLIFY=F) %>%
  plyr::ldply()

scores_tair_df <- mapply(function(x, y)
  data.frame('position'= -(upstream-y/2) : (downstream-y/2),
             'TAIR10'=colMeans(x)) %>%
    melt(id.vars='position', variable.name='TSSs'),
  scores_tair, matrices_lengths, SIMPLIFY=F) %>%
  plyr::ldply()

scores_araport_df <- mapply(function(x, y)
  data.frame('position'= -(upstream-y/2) : (downstream-y/2),
             'ARAPORT11'=colMeans(x)) %>%
    melt(id.vars='position', variable.name='TSSs'),
  scores_araport, matrices_lengths, SIMPLIFY=F) %>%
  plyr::ldply()

scores_nanopare_df <- mapply(function(x, y)
  data.frame('position'= -(upstream-y/2) : (downstream-y/2),
             'nanoPARE'=colMeans(x)) %>%
    melt(id.vars='position', variable.name='TSSs'),
  scores_nanopare, matrices_lengths, SIMPLIFY=F) %>%
  plyr::ldply()

scores_peat_df <- mapply(function(x, y)
  data.frame('position'= -(upstream-y/2) : (downstream-y/2),
             'PEAT'=colMeans(x)) %>%
    melt(id.vars='position', variable.name='TSSs'),
  scores_peat, matrices_lengths, SIMPLIFY=F) %>%
  plyr::ldply()

# merge dataframes and plot
all_scores_df <- rbind(scores_cage_df, scores_tair_df, scores_araport_df, scores_nanopare_df, scores_peat_df)

ggplot(all_scores_df, aes(x=position, y=value, col=TSSs)) +
       geom_vline(xintercept=+1, lty=2) +
       geom_line(alpha=.5) +
       facet_wrap(~.id, scales='free_y', ncol=3) +
       cowplot::theme_cowplot() + theme(aspect.ratio=1, strip.text.x=element_text(face='bold'), legend.position=c(0.7, 0.1)) +
       scale_color_brewer(palette='Set2', direction=-1) +
       labs(x='Distance to TSSs (bp)', y='Average score',
            title='POL motifs by TSS datasets',
            caption='matrices from JASPAR2018')

# TATA in CAGE WT ctrl by TPM rank
scores_cage$TATA %>%
  as.data.frame() %>%
  split(rowRanges(cage_tss)$TPMquantile) %>%
  lapply(colMeans) %>%
  plyr::ldply() %>%
  melt(id.vars='.id') %>%
  mutate('position'=as.numeric(variable)) %>%
  dplyr::select(-variable) %>%
  mutate('position'=position - (upstream-matrices_lengths$TATA/2)) %>%
    ggplot(aes(x=position, y=value, col=.id, group=.id)) +
           geom_vline(xintercept=1, lty=2) + geom_line(size=.8) +
           cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position=c(.03, .81)) +
           scale_color_brewer(palette='RdYlBu', direction=-1, name='TPM expression ranks') +
           labs(x='Distance to CAGE TSSs (bp)', y='Average score',
                title='TATA-box at CAGE TSSs by expression rank')

# How many CAGE WT TSSs have a TATA box between -40 to -10 to TSS?
# redo promoter sequences with different window
upstream=200
downstream=200

cage_seq <- cage_tss %>%
  rowRanges() %>%
  flank(start=T, width=1) %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

tair_seq <- tair_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

araport_seq <- araport_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

peat_seq <- peat_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

nanopare_seq <- nanopare_tss %>%
  rowRanges() %>%
  promoters(upstream=upstream, downstream=downstream) %>%
  getSeq(genome, .)

# bla <- motifScanHits(regionsSeq=cage_seq, motifPWM=as.matrix(matrices$TATA), minScore='64%')
# length(unique(bla$sequence)) / length(cage_seq) * 100

plotMotifOccurrenceAverage(regionsSeq=cage_seq , motifPWM=as.matrix(matrices$TATA) , minScore='65%', flankUp=200, flankDown=200,
                          smoothingWindow=3, color=c('darkgreen'), cex.axis=0.9)

plotMotifOccurrenceAverage(regionsSeq=tair_seq , motifPWM=as.matrix(matrices$TATA) , minScore='65%', flankUp=200, flankDown=200,
                           smoothingWindow=3, color=c('blue3'), add=T)

plotMotifOccurrenceAverage(regionsSeq=araport_seq , motifPWM=as.matrix(matrices$TATA) , minScore='65%', flankUp=200, flankDown=200,
                           smoothingWindow=3, color=c('red3'), add=T)

plotMotifOccurrenceAverage(regionsSeq=peat_seq , motifPWM=as.matrix(matrices$TATA) , minScore='65%', flankUp=200, flankDown=200,
                           smoothingWindow=3, color=c('orange'), add=T)

plotMotifOccurrenceAverage(regionsSeq=nanopare_seq , motifPWM=as.matrix(matrices$TATA) , minScore='65%', flankUp=200, flankDown=200,
                           smoothingWindow=3, color=c('violet'), add=T)

legend('topright', legend=c('CAGE wt', 'TAIR10', 'ARAPORT11', 'PEAT', 'nanoPARE'), col=c('darkgreen', 'blue3', 'red3', 'orange', 'violet'), bty='n', lwd=1)
title(main='TATA MotifOccurenceAverage with
      65% minscore')

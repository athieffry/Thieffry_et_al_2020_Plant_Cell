#### 03b. Make comparable set of TSSs
#### Axel Thieffry - July 2018
set.seed(42)
library(tidyverse)
library(tidylog)
library(reshape2)
library(patchwork)
library(magrittr)
library(CAGEfightR)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(BiocParallel)
register(MulticoreParam(workers=4))
library(TxDb.Athaliana.BioMart.plantsmart28)
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
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
    # remove non-canonical chromosomes
    seqlevels(myseqinfo) <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
# CAGE TSSs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_ctrl.rds')
    # remove non-canonical chromosomes
    seqlevels(TCs, pruning.mode='coarse') <- seqlevels(myseqinfo)
# CAGE CTSSs
CTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CTSSs_1count_min3lib_TSSstory.rds')
CTSSs_wt <- subset(CTSSs, select=genotype=='wt')
# TAIR10 TSSs
txdb_tair <- TxDb.Athaliana.BioMart.plantsmart28
    # remove non-canonical chromosomes
    seqlevelsStyle(txdb_tair) <- seqlevelsStyle(myseqinfo)[1]
    seqlevels(txdb_tair, pruning.mode='coarse') <- seqlevels(myseqinfo)
# ARAPORT11 TSSs
araport <- import.gff3('~/Dropbox/Axel_Arabidopsis_Flagellin/ARAPORT11/Araport11_GFF3_genes_transposons.201606.gff3')
    # fix seqinfo
    seqlevels(araport, pruning.mode='coarse') <- seqlevels(myseqinfo)
    seqinfo(araport) <- myseqinfo
    # make TxDb object
    txdb_araport <- suppressWarnings( makeTxDbFromGRanges(araport) )

# PEAT from Morton & al.: as provided by Schon et al.
peat <- read.table('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6a. comparison PEAT Morton/NanoPARE data/supplemental_code_S1/data_tables/PEAT_peaks.bed') %>% as_tibble()
    # clean
    colnames(peat) <- c('chr', 'start', 'end', 'geneID', 'score', 'strand', 'peak')
    peat$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr')
    # make start as peak
    peat %<>% mutate('start'=start+peak, 'end'=start)
    # make GenomicRanges with 1bp as Mode (peak)
    peat_gr <- makeGRangesFromDataFrame(df=peat, keep.extra.columns=F, seqinfo=myseqinfo,
                                        seqnames.field='chr', strand.field='strand',
                                        start.field='start', end.field='start')


# NanoPARE TSSs from Schon et al. (5ng)
nanopare <- read.table('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6a. comparison PEAT Morton/NanoPARE data/supplemental_code_S1/data_tables/fb.W.5P.bed') %>% as_tibble()
    # clean
    colnames(nanopare) <- c('chr', 'start', 'end', 'type', 'peak', 'strand', 'geneID', 'bla', 'blu', 'blo')
    nanopare$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr') %>% str_replace('Chrc', 'ChrC') %>% str_replace('Chrm', 'ChrM')
    nanopare %<>% select(-bla:-blo)
    # look at the different types
    nanopare %<>% separate(type, c('bla', 'type', 'n'), sep='\\.')
    nanopare %<>% select(-bla, -n)
    # make peak realative to start
    nanopare %<>% mutate('start'=start+peak+1, 'end'=start)
    # remove non-canonical chromosomes
    nanopare %<>% subset(chr %in% paste0('Chr', 1:5))
    # make GenomicRanges with 1bp as Mode (peak)
    nanopare_gr <- makeGRangesFromDataFrame(df=nanopare, keep.extra.columns=F, seqinfo=myseqinfo,
                                            seqnames.field='chr', strand.field='strand',
                                            start.field='start', end.field='start', starts.in.df.are.0based=F)
    subset(nanopare_gr, strand=='+')
    subset(nanopare_gr, strand=='-')


# 2. MAKE PROMOTER REGIONS TO BE QUANTIFIED ####
# ----------------------------------------------
# NB: the whole process is done with the CAGE WT as well just for fun/sanity check

# TSSs +/- 100 bp
cage_extended_tss <- rowRanges(TCs) # CAGE is not really extended, rather intrinsic range of TC is considered # 21,221
tair_extended_tss <- txdb_tair %>% promoters(upstream=100, downstream=101) %>% trim() # 41,392
araport_extended_tss <- txdb_araport %>% promoters(upstream=100, downstream=101) %>% trim() # 59,907
peat_extended_tss <- promoters(peat_gr, upstream=100, downstream=101) # 9,326
nanopare_extended_tss <- promoters(nanopare_gr, upstream=100, downstream=101) %>% trim() # 22,364

# TSSs 1bp
cage_tss_1bp <- rowRanges(TCs) %>% swapRanges() # 21,221
tair_tss_1bp <- txdb_tair %>% promoters(upstream=0, downstream=1) # 41,392
araport_tss_1bp <- txdb_araport %>% promoters(upstream=0, downstream=1) # 59,907
peat_tss_1bp <- peat_gr # 9,326
nanopare_tss_1bp <- nanopare_gr # 22,364



# 3. TPM-QUANTIFY ANNOTATION PROMOTERS WITH WT CTRL SAMPLE ####
# -------------------------------------------------------------
# function for quantifying any features (GR) across all samples
quantifyAcrossSample <- function(features, scoring_SE, inputAssay='TPM', sample){
                                 hits <- as(findOverlaps(query=features, subject=rowRanges(scoring_SE)), "List")
                                 newSE <- suppressWarnings( swapScores(object=scoring_SE, outputColumn='score', inputAssay=inputAssay, sample=sample) )
                                 sum(extractList(rowRanges(newSE)$score, hits))
                                  }

# quantify TSSs across samples
sampleNames <- CTSSs_wt$Name %>% as.character()

cage_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=cage_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
tair_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=tair_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
araport_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=araport_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
peat_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=peat_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))
nanopare_extended_tss_TPM <- sapply(sampleNames, function(x) quantifyAcrossSample(features=nanopare_extended_tss, scoring_SE=CTSSs_wt, inputAssay='TPM', sample=x))

# make SE with 1bp rowRanges
cage_1bp_tss_TPM_SE <- SummarizedExperiment(assays=list('TPM'=cage_extended_tss_TPM), rowRanges=cage_tss_1bp, colData=colData(CTSSs_wt)) # 21,221
tair_1bp_tss_TPM_SE <- SummarizedExperiment(assays=list('TPM'=tair_extended_tss_TPM), rowRanges=tair_tss_1bp, colData=colData(CTSSs_wt)) # 41,392
araport_1bp_tss_TPM_SE <- SummarizedExperiment(assays=list('TPM'=araport_extended_tss_TPM), rowRanges=araport_tss_1bp, colData=colData(CTSSs_wt)) # 59,907
peat_1bp_tss_TPM_SE <- SummarizedExperiment(assays=list('TPM'=peat_extended_tss_TPM), rowRanges=peat_tss_1bp, colData=colData(CTSSs_wt)) # 9,326
nanopare_1bp_tss_TPM_SE <- SummarizedExperiment(assays=list('TPM'=nanopare_extended_tss_TPM), rowRanges=nanopare_tss_1bp, colData=colData(CTSSs_wt)) # 22,364

# look at expression distribution of TSSs before any filtering
mean_TPM_all_TSSs_df <- rbind(tibble('nanoPARE'=rowMeans(nanopare_extended_tss_TPM)) %>% melt(variable.name='source', value.name='meanTPM'),
                              tibble('CAGE_wt'=rowMeans(cage_extended_tss_TPM)) %>% melt(variable.name='source', value.name='meanTPM'),
                              tibble('TAIR10'=rowMeans(tair_extended_tss_TPM)) %>% melt(variable.name='source', value.name='meanTPM'),
                              tibble('ARAPORT11'=rowMeans(araport_extended_tss_TPM)) %>% melt(variable.name='source', value.name='meanTPM'),
                              tibble('PEAT_seq'=rowMeans(peat_extended_tss_TPM)) %>% melt(variable.name='source', value.name='meanTPM')) %>%
  mutate('source'=factor(source, levels=c('TAIR10', 'PEAT_seq', 'nanoPARE', 'CAGE_wt', 'ARAPORT11')))


ggplot(mean_TPM_all_TSSs_df, aes(x=meanTPM, col=source, y=..count..)) +
       geom_density(lwd=1) +
       scale_x_log10() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       scale_color_brewer(palette='Set1', direction=-1) +
       labs(x='Expression\n(CAGE wt TPM, average rep.)', title='')


# 4. SUBSET ANNOTATION PROMOTERS ####
# -----------------------------------
# calculate TPM support for all 3 wt t0 replicates
cage_1bp_tss_TPM_SE %<>% calcSupport(inputAssay='TPM', unexpressed=1, outputColumn='support')
tair_1bp_tss_TPM_SE %<>% calcSupport(inputAssay='TPM', unexpressed=1, outputColumn='support')
araport_1bp_tss_TPM_SE %<>% calcSupport(inputAssay='TPM', unexpressed=1, outputColumn='support')
peat_1bp_tss_TPM_SE %<>% calcSupport(inputAssay='TPM', unexpressed=1, outputColumn='support')
nanopare_1bp_tss_TPM_SE %<>% calcSupport(inputAssay='TPM', unexpressed=1, outputColumn='support')

# plot supports
support_df <- rbind(data.frame('source'='CAGE_wt', 'support'=rowData(cage_1bp_tss_TPM_SE)$support),
                    data.frame('source'='TAIR10', 'support'=rowData(tair_1bp_tss_TPM_SE)$support),
                    data.frame('source'='ARAPORT11', 'support'=rowData(araport_1bp_tss_TPM_SE)$support),
                    data.frame('source'='PEAT', 'support'=rowData(peat_1bp_tss_TPM_SE)$support),
                    data.frame('source'='nanoPARE', 'support'=rowData(nanopare_1bp_tss_TPM_SE)$support)) %>% as_tibble()

ggplot(support_df, aes(x=support, fill=source)) +
       geom_bar(stat='count', position=position_dodge()) +
       geom_vline(aes(xintercept=1.5), col='red') +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       scale_fill_brewer(palette='Paired', name='Annotation', direction=-1) +
       labs(title='CAGE wt support of TSS datasets')

# subset by support
cage_tss_new <- cage_1bp_tss_TPM_SE[rowRanges(cage_1bp_tss_TPM_SE)$support >= 2] # 21,221
tair_tss_new <- tair_1bp_tss_TPM_SE[rowRanges(tair_1bp_tss_TPM_SE)$support >= 2] # 22,394
araport_tss_new <- araport_1bp_tss_TPM_SE[rowRanges(araport_1bp_tss_TPM_SE)$support >= 2] # 24,877
peat_tss_new <- peat_1bp_tss_TPM_SE[rowRanges(peat_1bp_tss_TPM_SE)$support >= 2] # 8,952
nanopare_tss_new <- nanopare_1bp_tss_TPM_SE[rowRanges(nanopare_1bp_tss_TPM_SE)$support >= 2] # 16,309

# save
saveRDS(cage_tss_new, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CAGE_wt_comparable_TSSs_redone_17June2019.rds')
saveRDS(tair_tss_new, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TAIR10_comparable_TSSs_redone_17June2019.rds')
saveRDS(araport_tss_new, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_ARAPORT11_comparable_TSSs_redone_17June2019.rds')
saveRDS(peat_tss_new, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_PEAT_comparable_TSSs_redone_17June2019.rds')
saveRDS(nanopare_tss_new, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_nanoPARE_comparable_TSSs_redone_17June2019.rds')

rowRanges(nanopare_tss_new) %>%
  subset(strand=='-')
rowRanges(nanopare_tss_new) %>%
  subset(strand=='+')


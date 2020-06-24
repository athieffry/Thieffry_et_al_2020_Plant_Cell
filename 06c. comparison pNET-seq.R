#### Arabidopsis flg22 : pNET-seq pausing at 3'UTRs (with and without aTSSs)
#### Axel Thieffry - April 2019
set.seed(42)
library(WriteXLS)
library(patchwork)
library(RColorBrewer)
library(edgeR)
library(limma)
library(DESeq2)
library(magrittr)
library(reshape2)
library(tidyquant)
library(GenomicFeatures)
library(GenomicRanges)
library(TeMPO)
library(rtracklayer)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(CAGEfightR)
library(tidyverse)
library(tidylog)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { GR[-idx]}
                                     else {GR}}

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6c. pNET-seq')



# 1. LOAD ALL INPUT FILES ####
# ----------------------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
# 3'UTR aTSSs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_TPM1_min3lib.rds')
aTSSs <- rowRanges(TCs) %>% subset(txType_TAIR10extended == 'antisense_threeUTR')
# 3'UTRs
threeUTRs_with_aTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_with_aTSS_gr.rds')
threeUTRs_wout_aTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_without_aTSS_gr.rds')
# pNET signal as BigWigFileList
pNET_plus <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6c. pNET-seq/averaged_bigwigs_tpm', pattern='_plus_tpm.bw', full.names=T) %>% BigWigFileList()
pNET_minus <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6c. pNET-seq/averaged_bigwigs_tpm', pattern='_minus_tpm.bw', full.names=T) %>% BigWigFileList()
filenames <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6c. pNET-seq/averaged_bigwigs_tpm', pattern='_plus_tpm.bw', full.names=F) %>% str_remove('_r12_plus_tpm.bw') %>% str_remove('_plus_tpm.bw')
names(pNET_plus) <- filenames
names(pNET_minus) <- filenames


# 2. READ AND NORMALIZE pNET-seq signals ####
# -------------------------------------------
if(FALSE){
# make GRanges and compute TPM
makeTPM <- function(x) { tpm_factor <- sum(x$count) / 1000000
                         x$tpm <- x$count / tpm_factor
                         x}

Ser2P_FP_r1   <- read.table('pNET_Ser2P_FP_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser2P_FP_r2   <- read.table('pNET_Ser2P_FP_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser2P_mock_r1 <- read.table('pNET_Ser2P_mock_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser2P_mock_r2 <- read.table('pNET_Ser2P_mock_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser2P_r12     <- read.table('pNET_Ser2P_rep1-2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()

Ser5P_FP_r1   <- read.table('pNET_Ser5P_FP_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser5P_FP_r2   <- read.table('pNET_Ser5P_FP_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser5P_mock_r1 <- read.table('pNET_Ser5P_mock_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
Ser5P_mock_r2 <- read.table('pNET_Ser5P_mock_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()

unph_FP_r1    <- read.table('pNET_unph_FP_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
unph_FP_r2    <- read.table('pNET_unph_FP_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
unph_mock_r1  <- read.table('pNET_unph_mock_rep1.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()
unph_mock_r2  <- read.table('pNET_unph_mock_rep2.txt', h=T) %>% makeGRangesFromDataFrame(seqinfo=myseqinfo, seqnames.field='chromosome', start.field='position', end.field='position', strand.field='strand', keep.extra.columns=T) %>% makeTPM()


# separate by strand and make coverage
Ser2P_FP_r1 %<>% split(strand(Ser2P_FP_r1)) %>% .[1:2]
Ser2P_FP_r2 %<>% split(strand(Ser2P_FP_r2)) %>% .[1:2]
Ser2P_mock_r1 %<>% split(strand(Ser2P_mock_r1)) %>% .[1:2]
Ser2P_mock_r2 %<>% split(strand(Ser2P_mock_r2)) %>% .[1:2]

Ser2P_r12 %<>% split(strand(Ser2P_r12)) %>% .[1:2]

Ser5P_FP_r1 %<>% split(strand(Ser5P_FP_r1)) %>% .[1:2]
Ser5P_FP_r2 %<>% split(strand(Ser5P_FP_r2)) %>% .[1:2]
Ser5P_mock_r1 %<>% split(strand(Ser5P_mock_r1)) %>% .[1:2]
Ser5P_mock_r2 %<>% split(strand(Ser5P_mock_r2)) %>% .[1:2]

unph_FP_r1 %<>% split(strand(unph_FP_r1)) %>% .[1:2]
unph_FP_r2 %<>% split(strand(unph_FP_r2)) %>% .[1:2]
unph_mock_r1 %<>% split(strand(unph_mock_r1)) %>% .[1:2]
unph_mock_r2 %<>% split(strand(unph_mock_r2)) %>% .[1:2]

# compute coverages
Ser2P_FP_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser2P_FP_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser2P_mock_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser2P_mock_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser2P_r12 %<>% lapply(function(x) coverage(x, weight='tpm'))

Ser5P_FP_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser5P_FP_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser5P_mock_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
Ser5P_mock_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))

unph_FP_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
unph_FP_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))
unph_mock_r1 %<>% lapply(function(x) coverage(x, weight='tpm'))
unph_mock_r2 %<>% lapply(function(x) coverage(x, weight='tpm'))

# merge replicates
Ser2P_FP_plus <- (Ser2P_FP_r1$`+` + Ser2P_FP_r2$`+`) / 2
Ser2P_FP_minus <- (Ser2P_FP_r1$`-` + Ser2P_FP_r2$`-`) / 2
Ser2P_mock_plus <- (Ser2P_mock_r1$`+` + Ser2P_mock_r2$`+`) / 2
Ser2P_mock_minus <- (Ser2P_mock_r1$`-` + Ser2P_mock_r2$`-`) / 2

Ser5P_FP_plus <- (Ser5P_FP_r1$`+` + Ser5P_FP_r2$`+`) / 2
Ser5P_FP_minus <- (Ser5P_FP_r1$`-` + Ser5P_FP_r2$`-`) / 2
Ser5P_mock_plus <- (Ser5P_mock_r1$`+` + Ser5P_mock_r2$`+`) / 2
Ser5P_mock_minus <- (Ser5P_mock_r1$`-` + Ser5P_mock_r2$`-`) / 2

unph_FP_plus <- (unph_FP_r1$`+` + unph_FP_r2$`+`) / 2
unph_FP_minus <- (unph_FP_r1$`-` + unph_FP_r2$`-`) / 2
unph_mock_plus <- (unph_mock_r1$`+` + unph_mock_r2$`+`) / 2
unph_mock_minus <- (unph_mock_r1$`-` + unph_mock_r2$`-`) / 2

# ouput as bigwig
export.bw(Ser2P_FP_plus, 'averaged_bigwigs_tpm/Ser2P_FP_r12_plus_tpm.bw')
export.bw(Ser2P_FP_minus, 'averaged_bigwigs_tpm/Ser2P_FP_r12_minus_tpm.bw')
export.bw(Ser2P_mock_plus, 'averaged_bigwigs_tpm/Ser2P_mock_r12_plus_tpm.bw')
export.bw(Ser2P_mock_minus, 'averaged_bigwigs_tpm/Ser2P_mock_r12_minus_tpm.bw')
export.bw(Ser2P_r12$`+`, 'averaged_bigwigs_tpm/Ser2P_r1-2_plus_tpm.bw')
export.bw(Ser2P_r12$`-`, 'averaged_bigwigs_tpm/Ser2P_r1-2_minus_tpm.bw')
export.bw(Ser5P_FP_plus, 'averaged_bigwigs_tpm/Ser5P_FP_r12_plus_tpm.bw')
export.bw(Ser5P_FP_minus, 'averaged_bigwigs_tpm/Ser5P_FP_r12_minus_tpm.bw')
export.bw(Ser5P_mock_plus, 'averaged_bigwigs_tpm/Ser5P_mock_r12_plus_tpm.bw')
export.bw(Ser5P_mock_minus, 'averaged_bigwigs_tpm/Ser5P_mock_r12_minus_tpm.bw')
export.bw(unph_FP_plus, 'averaged_bigwigs_tpm/unph_FP_r12_plus_tpm.bw')
export.bw(unph_FP_minus, 'averaged_bigwigs_tpm/unph_FP_r12_minus_tpm.bw')
export.bw(unph_mock_plus, 'averaged_bigwigs_tpm/unph_mock_r12_plus_tpm.bw')
export.bw(unph_mock_minus, 'averaged_bigwigs_tpm/unph_mock_r12_minus_tpm.bw')

}




# 3. FOOTPRINT AT 3'UTRs aTSSs ####
# ---------------------------------
footprint_aTSSs <- tidyMetaProfile(sites=swapRanges(aTSSs), forward=pNET_plus, reverse=pNET_minus, upstream=300, downstream=300, trimLower=0.01, trimUpper=0.99)

gg_atss <- footprint_aTSSs %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  separate(signal, c('Ser', 'type'), sep='_') %>%
  subset(type != 'r1-2') %>%
  ggplot(aes(x=pos0, y=score, col=type, group=interaction(type, direction))) +
         geom_vline(xintercept=0, lty=2) + geom_hline(yintercept=0) +
         geom_line(alpha=.5) + geom_ma(n=10, lwd=1, lty=1, ma_fun=WMA) +
         facet_wrap(~Ser, nrow=1, ncol=3) +
         theme(strip.text=element_text(face='bold', size=14), strip.background=element_blank()) +
         labs(x="Position relative to 3'UTR aTSSs (bp)",
              y="Normalized signal (15bp moving average)",
              title="pNET-seq at 3'UTR aTSSs",
              caption='Data from Dong et al., 2018')



# 4. FOOTPRINT AT PROMOTER TSSs (sanity-check) ####
# -------------------------------------------------
prom_tss <- rowRanges(TCs) %>% subset(txType_TAIR10=='promoter')
footprint_prom <- tidyMetaProfile(sites=swapRanges(prom_tss), forward=pNET_plus, reverse=pNET_minus, upstream=300, downstream=300, trimLower=0.01, trimUpper=0.99)

gg_tss <- footprint_prom %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=T) %>%
  separate(signal, c('Ser', 'type'), sep='_') %>%
  subset(type != 'r1-2') %>%
  ggplot(aes(x=pos0, y=score, col=type, group=interaction(type, direction))) +
         geom_vline(xintercept=0, lty=2) + geom_hline(yintercept=0) +
         geom_line(lwd=1) +
         facet_wrap(~Ser, nrow=1, ncol=3) +
         theme(strip.text=element_text(face='bold', size=14), strip.background=element_blank()) +
         labs(x="Position relative to promoter TSSs (bp)",
              y="Normalized signal (15bp moving average)",
              title="pNET-seq at promoter TSSs (sanity check)")

gg_tss + gg_atss + plot_layout(nrow=2)

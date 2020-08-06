#### HOW MANY PEAT/NANOPARE TSSs AROUND OUR 96 CAGE PROMPT TCs?
#### Axel Thieffry
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(reshape2)
library(CAGEfightR)
library(patchwork)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(RColorBrewer)
library(BiocParallel)
library(seqPattern)
library(TFBSTools)
library(ggseqlogo)
library(motifStack)
library(JASPAR2016)
library(BSgenome.Athaliana.TAIR.TAIR9)
genome <- BSgenome.Athaliana.TAIR.TAIR9
genome@seqinfo@genome[] <- "TAIR10" # little hack  
register(MulticoreParam(workers=3))
options(scipen=999) # disable scientific notation
'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}
setwd('~/masked_path/11 - PEAT & nanoPARE support of PROMPT TCs')


# 1. GET INPUT DATA ####
# ----------------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')

# CAGE TCs (work with that one because it has IRanges as peaks)
TCs <- readRDS('~/masked_path/SE_TCs_TPM1_min3lib_TSSstory.rds')

# CAGE TCs ids (get those one first to select TC_ids that are PROMPT: )
PROMPT_ids <- readRDS('~/masked_path/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds') %>%
              rowRanges() %>%
              subset(txType_TAIR10=='reverse' & (genotypehen2==1 | genotyperrp4==1)) %>%
              names()

# nanoPARE TSSs
nanopare <- read.table('~/masked_path/fb.W.5P.bed') %>% as_tibble()
    # clean
    colnames(nanopare) <- c('chr', 'start', 'end', 'type', 'peak', 'strand', 'geneID', 'bla', 'blu', 'blo')
    nanopare$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr') %>% str_replace('Chrc', 'ChrC') %>% str_replace('Chrm', 'ChrM')
    nanopare %<>% select(-bla:-blo)
    # look at the different types
    nanopare %<>% separate(type, c('bla', 'type', 'n'), sep='\\.')
    nanopare %<>% select(-bla, -n)
    # make peak realative to start
    nanopare %<>% mutate('start'=start+peak+1, 'end'=start)
    # make GenomicRanges with 1bp as Mode (peak)
    nanopare_gr <- makeGRangesFromDataFrame(df=nanopare, keep.extra.columns=F, seqinfo=myseqinfo,
                                            seqnames.field='chr', strand.field='strand',
                                            start.field='start', end.field='start', starts.in.df.are.0based=F)
    
    # repeat seqlogo at nanoPARE TSSs to see if we're good with positions: all GOOD.
    if(FALSE){
              nanopare_seqs <- nanopare_gr %>%
                promoters(upstream=50, downstream=30) %>%
                remove_out_of_bound() %>%
                getSeq(genome, .)
              
              nanopare_pfm <- nanopare_seqs %>%
                consensusMatrix(as.prob=T) %>%
                subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
                as.matrix() %>%
                new('pfm', mat=., name='nanoPARE')
              
              nanopare_pfm@background <- colSums(letterFrequency(nanopare_seqs, DNA_BASES)) / sum(colSums(letterFrequency(nanopare_seqs, DNA_BASES)))
              
              ggseqlogo(data=list('nanoPARE'=nanopare_pfm@mat)) +
                        scale_x_continuous(breaks=seq(1, 80, 10), labels=seq(-50, 20, 10)) +
                        theme(aspect.ratio=.25) + labs(x='Distance to TSSs (bp)')
            }

# PEAT TSSs
peat <- read.table('~/masked_path/PEAT_peaks.bed') %>% as_tibble()
      # clean
      colnames(peat) <- c('chr', 'start', 'end', 'geneID', 'score', 'strand', 'peak')
      peat$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr')
      # make start as peak
      peat %<>% mutate('start'=start+peak, 'end'=start)
      # make GenomicRanges with 1bp as Mode (peak)
      peat_gr <- makeGRangesFromDataFrame(df=peat, keep.extra.columns=F, seqinfo=myseqinfo,
                                          seqnames.field='chr', strand.field='strand',
                                          start.field='start', end.field='start')
      
      # repeat seqlogo at PEAT TSSs to see if we're good with positions: all GOOD.
      if(FALSE){
        peat_seqs <- peat_gr %>%
                     promoters(upstream=50, downstream=30) %>%
                     remove_out_of_bound() %>%
                     getSeq(genome, .)
        
        peat_pfm <- peat_seqs %>%
                    consensusMatrix(as.prob=T) %>%
                    subset(rownames(.) %in% c('A', 'T', 'C', 'G')) %>%
                    as.matrix() %>%
                    new('pfm', mat=., name='PEAT')
        
        peat_pfm@background <- colSums(letterFrequency(peat_seqs, DNA_BASES)) / sum(colSums(letterFrequency(peat_seqs, DNA_BASES)))
        
        ggseqlogo(data=list('PEAT'=peat_pfm@mat)) +
                  scale_x_continuous(breaks=seq(1, 80, 10), labels=seq(-50, 20, 10)) +
                  theme(aspect.ratio=.25) + labs(x='Distance to TSSs (bp)')
      }

# check that CAGE TCs, nanoPARE TSSs & PEAT TSSs all have same seqinfo
seqinfo(TCs) ; seqinfo(nanopare_gr) ; seqinfo(peat_gr)


# 2. ANALYSIS ####
# ----------------
# 2a. extract PROMPT TCs (n=96)
PROMPT_TCs <- rowRanges(TCs) %>%
              subset(names(.) %in% PROMPT_ids) %>%
              swapRanges()

      # LOOP FOR DIFFERENT DISTANCES AROUND PROMPT TCs
      windows <- seq(0, 150, 10)
      
      support <- data.frame('windows'=windows,
                 'nanoPARE'=sapply(windows, function(x) {
                            PROMPT_TCs_extended <- PROMPT_TCs %>% promoters(upstream=x, downstream=x)
                            # 2c. how many PROMPT TCs have a nanoPARE / PEAT TSSs in their proximity?
                            subsetByOverlaps(PROMPT_TCs_extended, nanopare_gr) %>% length()
                           }),
                 'PEAT'=sapply(windows, function(x) {
                   PROMPT_TCs_extended <- PROMPT_TCs %>% promoters(upstream=x, downstream=x)
                   # 2c. how many PROMPT TCs have a nanoPARE / PEAT TSSs in their proximity?
                   subsetByOverlaps(PROMPT_TCs_extended, peat_gr) %>% length()
                 }))

# plot
support %>%
  melt(id.vars='windows') %>%
  ggplot(aes(x=as.factor(windows), y=value, col=variable)) +
         geom_point(size=3) +
         geom_line(aes(group=variable)) +
         geom_text(aes(label=value, y=value+.75), col='black') +
         cowplot::theme_cowplot() +
         theme(aspect.ratio=1, legend.position=c(.65,.5)) +
         scale_color_brewer(palette='Dark2', name='TSS source') +
         labs(x='Distance from PROMPT CAGE TC peaks (N=96) (bp)',
              y='Number of TSSs (same strandedness)',
              title='Support of PROMPT CAGE TCs\nby PEAT & nanoPARE TSSs')
21/96*100
1/96*100

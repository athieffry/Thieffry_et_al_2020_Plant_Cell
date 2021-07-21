#### comparison with PEAT from Morton & al
#### Same analysis as NanoPARE from Schon & al
#### Axel Thieffry
library(tidyverse)
library(magrittr)
library(reshape2)
library(rtracklayer)
library(CAGEfightR)
library(ggpubr)
library(TeMPO)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(RColorBrewer)
library(BiocParallel)
register(SerialParam())
register(MulticoreParam(workers=3))
library(readxl)
library(WriteXLS)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}

setwd('~/masked_path/6a. comparison PEAT Morton')



# 1. READ DATA ####
# -----------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')

# NanoPARE TSSs from Schon et al. (5ng)
nanopare <- read.table('~/masked_path/fb.W.5P.bed') %>% as_tibble()
    # clean
    colnames(nanopare) <- c('chr', 'start', 'end', 'type', 'peak', 'strand', 'geneID', 'bla', 'blu', 'blo')
    nanopare$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr') %>% str_replace('Chrc', 'ChrC') %>% str_replace('Chrm', 'ChrM')
    nanopare %<>% select(-bla:-blo)
    # look at the different types
    nanopare %<>% separate(type, c('bla', 'type', 'n'), sep='\\.')
    nanopare %<>% select(-bla, -n)
    # look at NanoPARE TSS width distribution
    nanopare %<>% mutate('width'=end-start+1)
    table(nanopare$type)
    gg_width_nanopare <- ggplot(nanopare, aes(x=width)) +
           geom_histogram(binwidth=3) +
           facet_wrap(~type, scales='free_y') +
           xlim(NA, 500) +
           labs(title='NanoPARE TSS width distribution', x='width (bp)', subtitle='Schon et al.', caption='binwidth=3bp') +
           theme(aspect.ratio=1)
    # make peak realative to start
    nanopare %<>% mutate('peak'=start+peak)
    # make GenomicRanges with 1bp as Mode (peak)
    nanopare_gr <- makeGRangesFromDataFrame(df=nanopare, keep.extra.columns=T, seqinfo=myseqinfo,
                                            seqnames.field='chr', strand.field='strand',
                                            start.field='start', end.field='end')
    nanopare_gr

# PEAT from Morton & al.: as provided by Schon et al.
peat <- read.table('~/masked_path/PEAT_peaks.bed') %>% as_tibble()
    # clean
    colnames(peat) <- c('chr', 'start', 'end', 'geneID', 'score', 'strand', 'peak')
    peat$chr %<>% str_remove('Ath_') %>% str_replace('chr', 'Chr')
    # look at PEAT TSS width distribution
    peat %<>% mutate('width'=end-start+1)
    gg_width_peat <- ggplot(peat, aes(x=width)) +
      geom_histogram(binwidth=3) +
      xlim(NA, 500) +
      labs(title='PEAT TSS width distribution', x='width (bp)', subtitle='Morton et al.', caption='binwidth=3bp') +
      theme(aspect.ratio=1)
    # make start as peak
    peat %<>% mutate('peak'=start+peak)
    # make GenomicRanges with 1bp as Mode (peak)
    peat_gr <- makeGRangesFromDataFrame(df=peat, keep.extra.columns=T, seqinfo=myseqinfo,
                                        seqnames.field='chr', strand.field='strand',
                                        start.field='start', end.field='end')

# CAGE wildtype TSSs (t=0)
TCs_ctrl <- readRDS('~/masked_path/SE_TCs_ctrl.rds') # 21,221
    # look at CAGE TSS width distribution
    gg_width_cage <- rowRanges(TCs_ctrl) %>%
                     as.data.frame() %>%
                     as_tibble() %>%
                     ggplot(aes(x=IQR)) +
                            geom_histogram(binwidth=1) +
                            xlim(NA, 500) +
                            labs(title='CAGE wildtype TSS width distribution', x='width (bp)', caption='bindwidth=3bp') +
                            theme(aspect.ratio=1)

Rmisc::multiplot(gg_width_nanopare, gg_width_peat, gg_width_cage, layout=matrix(c(1, 1, 2, 3), byrow=T, nrow=2))
dev.off()



# 2. POSITIONAL OVERLAP BETWEEN TSSs FROM ALL 3 SOURCES ####
# ----------------------------------------------------------
# combine all TSS positions and reduce them
# to avoid un-equal equivalences when overlapping (1 big <=> 2 small)
    # remove mcols
    TCs_ctrl_noMcols <- rowRanges(TCs_ctrl) # 21,221
    nanopare_noMcols_gr <- nanopare_gr # 22,852
    peat_noMcols_gr <- peat_gr # 9,326
    
    mcols(TCs_ctrl_noMcols) <- mcols(nanopare_noMcols_gr) <- mcols(peat_noMcols_gr) <- NULL
    combined_TSSs_gr <- c(TCs_ctrl_noMcols, nanopare_noMcols_gr, peat_noMcols_gr) # 53,399
    combined_TSSs_gr <- reduce(combined_TSSs_gr, ignore.strand=F) # 27,014
    
    # make countOverlap matrix
    overlap_matrix <- data.frame('CAGE_wt'=countOverlaps(combined_TSSs_gr, rowRanges(TCs_ctrl), ignore.strand=F),
                                 'PEAT'=countOverlaps(combined_TSSs_gr, peat_gr, ignore.strand=F),
                                 'NanoPARE'=countOverlaps(combined_TSSs_gr, nanopare_gr, ignore.strand=F))
    par(pty='s')
      limma::vennDiagram(overlap_matrix) # plot with the eulerr.co webapp
    dev.off()



# 3. ANALYSIS: How many PEAT peaks are detected in a 50bp-window around TSSs? ####
# --------------------------------------------------------------------------------
    # TAIR10 models
    txdb <- TxDb.Athaliana.BioMart.plantsmart28
    seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
    
    # ARAPORT11 models
    araport <- import.gff3('~/masked_path/Araport11_GFF3_genes_transposons.201606.gff')
    # clean mcols
    mcols(araport) %<>% as.data.frame() %>% select(type, ID, Name, locus_type)
    seqlevels(araport) <- seqlevels(myseqinfo)
    seqinfo(araport) <- myseqinfo
    
    # make gr versions with PEAKS
    peat_peak_gr <- peat_gr
    start(peat_peak_gr) <- end(peat_peak_gr) <- peat_peak_gr$peak
    
    nanopare_peak_gr <- nanopare_gr
    start(nanopare_peak_gr) <- end(nanopare_peak_gr) <- nanopare_peak_gr$peak
    
# 2a. TAIR10:
    # get 5'ends of genes
    tair_5ends <- promoters(txdb, upstream=50, downstream=51) %>%
                  trim() # 6,961 : 74.64%
    
    # how many PEAT peaks are captured?
    (recall_tair <- subsetByOverlaps(peat_peak_gr, tair_5ends) %>% length())

# 2b. ARAPORT11:
    # get 5'ends of genes
    table(araport$type)
    araport_5ends <- subset(araport, type=='gene') %>%
                     promoters(upstream=50, downstream=51) %>%
                     trim() # 2,252 : 24.14%
    
    # how many PEAT peaks are captured?
    (recall_araport <- subsetByOverlaps(peat_peak_gr, araport_5ends) %>% length())

# 2c. CAGE ctrl TSSs:
    # extend 50bp
    cage_ctrl_5ends <- TCs_ctrl %>%
                       rowRanges() %>%
                       swapRanges() %>%
                       promoters(upstream=50, downstream=51) # 8,322 : 89,2%
    
    # how many PEAT peaks are captured?
    (recall_cage_ctrl <- subsetByOverlaps(peat_peak_gr, cage_ctrl_5ends) %>% length())

# 2d. NanoPARE:
    # extend 50bp
    nanopare_5ends <- nanopare_peak_gr %>%
                      promoters(upstream=50, downstream=51) %>%
                      trim() # 7,620 : 81.7%

    # how many PEAT peaks are captured?
    (recall_nanopare <- subsetByOverlaps(peat_peak_gr, nanopare_5ends) %>% length())

# Summarize recall by TSS source
recall_df <- tibble('TSS_source' = c('TAIR10', 'ARAPORT11', 'CAGE_wt', 'Schon et al.'),
                    'recall' = c(recall_tair, recall_araport, recall_cage_ctrl, recall_nanopare),
                    'type' = c(rep('annotation', 2), 'CAGE', 'NanoPARE'))

recall_df$TSS_source %<>% factor(levels=c('ARAPORT11', 'TAIR10', 'Schon et al.', 'CAGE_wt'))

ggplot(recall_df, aes(x=TSS_source, y=recall, fill=type)) +
       geom_bar(stat='identity', col='black', lwd=.2) +
       cowplot::theme_cowplot() + theme(aspect.ratio=1.8, axis.text.x=element_text(angle=45, hjust=1)) +
       geom_text(aes(label=paste0(round(recall/length(peat_gr)*100, digits=2), '%'), vjust=-1)) +
       scale_fill_brewer(palette='Paired') +
       labs(title='Morton et al. PEAT TSSs: recall fractions',
            subtitle='How many PEAT TSSs are found within\n50bp around different TSSs sources?',
            x='', y='Number of PEAT TSSs within 50bp (100%=9,326)')




# 3. CUMULATIVE FREQUENCY DISTRIBUTION OF POSITIONAL ERROR ####
# -------------------------------------------------------------
# for all features withing 200 bp of a PEAT peak
# 3a. make PEAT peaks +/- 200bp
peat_200_gr <- promoters(peat_peak_gr, upstream=200, downstream=201)

# 3b. keep TSSs peaks of all sources within peat_200_gr
tair_tss_1bp <- txdb %>% promoters(upstream=0, downstream=1)
araport_tss_1bp <- araport %>% subset(type=='five_prime_UTR') %>% promoters(upstream=0, downstream=1)
CAGEctrl_tss_1bp <- TCs_ctrl %>% rowRanges() %>% swapRanges()

tair_tss_1bp_within200 <- subsetByOverlaps(tair_tss_1bp, peat_200_gr)
araport_tss_1bp_within200 <- subsetByOverlaps(araport_tss_1bp, peat_200_gr)
CAGEctrl_tss_1bp_within200 <- subsetByOverlaps(CAGEctrl_tss_1bp, peat_200_gr)
nanopare_tss_1bp_within200 <- subsetByOverlaps(nanopare_peak_gr, peat_200_gr)

# 3c. compute distance to PEAT peak
tair_dist <- distanceToNearest(tair_tss_1bp_within200, peat_peak_gr) %>% mcols() %>% as_tibble()
araport_dist <- distanceToNearest(araport_tss_1bp_within200, peat_peak_gr) %>% mcols() %>% as_tibble()
CAGEctrl_dist <- distanceToNearest(CAGEctrl_tss_1bp_within200, peat_peak_gr) %>% mcols() %>% as_tibble()
nanopare_dist <- distanceToNearest(nanopare_tss_1bp_within200, peat_peak_gr) %>% mcols() %>% as_tibble()

# 3d. plot cumulative frequency of distance
ggplot() + stat_ecdf(data=tair_dist, geom='line', aes(x=distance), col=brewer.pal(n=3, name='Paired')[1], lwd=1.5) +
           stat_ecdf(data=araport_dist, geom='line', aes(x=distance), col=brewer.pal(n=3, name='Paired')[1], lwd=1.5, lty=2) +
           stat_ecdf(data=CAGEctrl_dist, geom='line', aes(x=distance), col=brewer.pal(n=3, name='Paired')[2], lwd=1.5) +
           stat_ecdf(data=nanopare_dist, geom='line', aes(x=distance), col=brewer.pal(n=3, name='Paired')[3], lwd=1.5) +
           theme_bw() + theme(aspect.ratio=1, panel.border=element_rect(fill=NA, size=1), axis.line=element_blank()) +
           scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
           labs(x='Distance to PEAT peak (bp)', y='Cumulative frequency', title='TSS positional error')

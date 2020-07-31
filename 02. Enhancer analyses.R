#### Enhancer analyses at timepoint = 0 only
#### Axel Thieffry - June 2018
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(TeMPO)
library(CAGEfightR)
library(tidyquant)
library(ggpubr)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(BiocParallel)
register(MulticoreParam(workers=4))
library(patchwork)
library(Gviz)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
seqlevels(txdb) <- seqlevels(myseqinfo)
library(Matrix)
library(BSgenome.Athaliana.TAIR.TAIR9)
genome <- BSgenome.Athaliana.TAIR.TAIR9
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]}
                                     else {o <- GR}
                                     o}
setwd('~/masked_path/02 - Enhancers analyses')


# 0. READ INPUT DATA ####
# -----------------------
enhancers <- readRDS('~/masked_path/SE_Enhancers.rds') # 113 enhancers withtout timecourse
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
CTSSs <- readRDS('~/masked_path/SE_CTSSs_1count_min3lib_TSSstory.rds')
dhss <- readRDS('~/masked_path/new_dhssumits_annotated_plantDHSonly.rds')

# CAGE signal
sampleNames <- list.files(path='~/masked_path/bw_files_R123', pattern='_0_R123.minus.tpm.bw', full.names=F) %>% str_remove('_0_R123.minus.tpm.bw')
      # locate files on hard drive & create two named BigWigFileList-objects
      cage_plus <- list.files(path='~/masked_path/bw_files_R123', pattern='_0_R123.plus', full.names=T) %>% BigWigFileList()
      cage_minus <- list.files(path='~/masked_path/bw_files_R123', pattern='_0_R123.minus', full.names=T) %>% BigWigFileList()
      names(cage_plus) <- names(cage_minus) <- sampleNames

# ALL Histone marks (as provided by Maxim Ivanov)
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

# DNaseI & MNase PlantDHS.org
dmnase <- list('~/masked_path/plantdhs.org_Ath.leaf.DNase.bw',
               '~/masked_path/plantdhs.org_Ath_nucleosomes_leaf_NPS.bw') %>% BigWigFileList()
names(dmnase) <- c('DNaseI', 'MNaseI')

# GRO-seq PNAS 2016 (Hetzel)
gro_p <- list.files(path='~/masked_path/03 - TSS analysis', pattern='GRO.*.plus.bw', full.name=T) %>% BigWigFileList()
gro_m <- list.files(path='~/masked_path/03 - TSS analysis', pattern='GRO.*.minus.bw', full.name=T) %>% BigWigFileList()
names(gro_p) <- names(gro_m) <- list.files(path='~/masked_path/03 - TSS analysis', pattern='GRO.*.plus.bw', full.name=F) %>% str_remove('_TPM99pc.plus.bw')

# RNA-Seq BIGWIG FILES
rnaseq_forward <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Forward', full.names=T))
rnaseq_reverse <- BigWigFileList(list.files('~/masked_path/bigwigs_RPM_R123_fixed', pattern='Reverse', full.names=T))
rnaseq_names <- list.files(path='~/masked_path/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names



# 1. DESCRIPTIVE STATISTICS ####
# ------------------------------
### 1a. summary plot of enhancers
enhancers_df <- rowRanges(enhancers) %>% as.data.frame() %>% as_tibble()

txType_labeller <- c(`intron` = paste0('intronic (N=', sum(enhancers_df$txType_TAIR10=='intron'), ')'),
                     `intergenic` = paste0('intergenic (N=', sum(enhancers_df$txType_TAIR10=='intergenic'),')'))

ggplot(enhancers_df, aes(x=balance, y=score, col=bidirectionality)) +
       geom_point() +
       facet_grid(txType_TAIR10~., labeller=as_labeller(txType_labeller)) +
       scale_y_log10() +
       scale_color_continuous(name='Bidirectional\nsupport') +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       labs(title='Enhancer expression and bidirectionality',
            subtitle=paste0('N=', nrow(enhancers_df)),
            x='Transcriptional balance (Bhattacharyya coefficient)', y='Pooled CAGE Expression (TPM)')


### 1b. expression density ridges for all enhancers by sample
enhancers_tpm_df <- assay(enhancers, 'TPM') %>% as.data.frame()
enhancer_annotation <- data.frame('ID'= names(enhancers), 'txType_TAIR10'=rowRanges(enhancers)$txType_TAIR10)

enhancers_tpm_df %>%
  rownames_to_column('ID') %>%
  as_tibble() %>%
  melt(value.name='TPM', variable.name='sample') %>%
  separate(sample, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
  group_by(ID, genotype, timepoint) %>%
  summarise('meanTPM'=mean(TPM)) %>%
  ungroup() %>%
  mutate('bla'=genotype, 'bli'=timepoint) %>%
  unite(sample, c('genotype', 'timepoint'), sep=' ') %>%
  rename('genotype'='bla', 'timepoint'='bli') %>%
  mutate('sample'=factor(sample, levels=c('rrp4 0', 'hen2 0', 'wt 0'))) %>%
  left_join(enhancer_annotation, by='ID') %>%
  ggplot(aes(x=meanTPM+0.25, y=txType_TAIR10, fill=genotype, col=genotype)) +
         geom_density_ridges2(scale=1, alpha=.5) +
         scale_x_log10() +
         theme_ridges(center_axis_labels=T, grid=F) +
         scale_fill_brewer(palette='Set2', name='genotype', direction=-1) +
         scale_color_brewer(palette='Set2', direction=-1) +
         labs(title='Enhancer expression density', subtitle=paste0('N=', nrow(enhancers_df)), y='Genomic location')



# 2. OVERLAP OF CAGEfightR ENHANCERS AND DHSSs ####
# -------------------------------------------------
enhancers_gr <- rowRanges(enhancers)
enhancers_gr$dhss <- countOverlaps(enhancers_gr, dhss) != 0
rowRanges(enhancers)$dhss <- countOverlaps(enhancers_gr, dhss) != 0

enhancers_gr %>%
  as.data.frame() %>%
  ggplot(aes(x=txType_TAIR10, fill=dhss)) +
         geom_bar(stat='count', col='black', lwd=.3) +
         geom_text(stat='count', aes(label=..count..), position=position_stack(vjust=0.5), col=c('white', 'black', 'white', 'black')) +
         scale_fill_manual(values=c('white', 'black')) +
         cowplot::theme_cowplot() + theme(aspect.ratio=.75) +
         labs(title='Enhancer & DHS overlap', x='Enhancer category', y='Nb. enhancers (N=113)')



# 3. EXOSOME SENSITIVITY OF ENHANCERS ####
# ----------------------------------------
# average TPM over replicates
enhancers_meantpm_df <- enhancers_tpm_df %>%
  rownames_to_column('ID') %>%
  melt(id.vars='ID', variable.name='sample', value.name='TPM') %>%
  separate(sample, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
  group_by(ID, genotype) %>%
  summarise('meanTPM'=mean(TPM)) %>%
  as_tibble()

enhancers_meantpm_df2 <- enhancers_meantpm_df %>% spread(genotype, meanTPM)

# log-log plots
ggplot(enhancers_meantpm_df2, aes(x=wt, y=hen2)) +
      geom_abline(intercept=0, slope1) +
      geom_point() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
      scale_x_log10() + scale_y_log10() +
      ggplot(enhancers_meantpm_df2, aes(x=wt, y=rrp4)) +
      geom_abline(intercept=0, slope1) +
      geom_point() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
      scale_x_log10() + scale_y_log10() +
      plot_layout(nrow=2)




# 4. HEATMAP ####
# ---------------
# heatmap annotation dataframe
annot_col_df <- colData(enhancers) %>% as.data.frame() %>% select(genotype)
brewer.pal(name='Set2', n=3)
genotype_col <- c('wt'='#66C2A5', 'hen2'='#8DA0CB', 'rrp4'='#FC8D62')
dhs_col <- c('yes'='grey40', 'no'='grey80')
anno_col_colors <- list('genotype'=genotype_col, 'dhss'=dhs_col)

annot_row_df <- rowRanges(enhancers) %>%
                as.data.frame() %>%
                dplyr::select(balance, bidirectionality, txType_TAIR10, dhss) %>%
                rownames_to_column('id') %>%
                mutate('dhss'=ifelse(dhss, 'yes', 'no')) %>%
                column_to_rownames('id')

pheatmap::pheatmap(log(enhancers_tpm_df+0.1), cluster_cols=T, cluster_rows=T, scale='row', clustering_method='ward.D2',
                   cellwidth=11, show_rownames=F, annotation_col=annot_col_df, annotation_colors=anno_col_colors, border_color=NA,
                   main='Enhancer clustering', annotation_row=annot_row_df, cutree_rows=2, cutree_cols=2)




# 5. GENOME BROWSER SHOTS ####
# ----------------------------
# Genome tracks
axis_track <- GenomeAxisTrack()
seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
tx_track <- GeneRegionTrack(txdb, name="Gene Models", col=NA, fill="bisque4", shape="arrow", showId=T)
# Define plot area (test first enhancer)
plot_region <- enhancers %>% rowRanges() %>% .[1] %>% add(100) %>% unstrand()
# CTSSs track
ctss_track <- CTSSs %>% rowRanges() %>% subsetByOverlaps(plot_region) %>% trackCTSS(name="CTSSs")
# Cluster track
cluster_track <- enhancers %>% subsetByOverlaps(plot_region) %>% trackClusters(name="Enhancers", col=NA, showId=T)
# Plot at tracks together
plotTracks(list(axis_track, ctss_track, cluster_track, tx_track), from=start(plot_region), to=end(plot_region), chromosome=seqnames(plot_region))




# 6. GVIZ LOOP ####
# -----------------
for (i in seq(1:length(enhancers))) {
  print(paste0("Plotting Enhancer ", i, '/', length(enhancers)))
  # Define plot area (test first enhancer)
  plot_region <- enhancers %>% rowRanges() %>% .[i] %>% add(100) %>% unstrand()
  enhancer_name <- names(plot_region) %>% str_replace(':', '-')
  # CTSSs track
  ctss_track <- suppressWarnings(suppressMessages(CTSSs %>% rowRanges() %>% subsetByOverlaps(plot_region) %>% trackCTSS(name="CTSSs")))
  # Cluster track
  cluster_track <- suppressWarnings(suppressMessages(enhancers %>% subsetByOverlaps(plot_region) %>% trackClusters(name="Enhancers", col=NA, showId=T)))
  # Plot at tracks together
  outfile <- paste0('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/02 - Enhancers analyses/enhancer CAGEfightR Gviz screenshots/enhancer_',i, '_', enhancer_name,'.pdf')
  pdf(file=outfile, paper='a4')
  plotTracks(list(axis_track, ctss_track, cluster_track, tx_track), from=start(plot_region), to=end(plot_region), chromosome=seqnames(plot_region))
  dev.off()
}
# retained enhancers from GVIZ
# 34 visual candidates identified after looking manually through 207 Gviz plots
visu_candidates <- c(1, 3, 10:12, 27, 31, 33:35, 43, 45, 49, 51:54, 59:61, 70, 72:74, 79, 83, 92, 94, 102, 103, 106:108, 113)
# look by IGV
rowRanges(enhancers)[visu_candidates] %>% as.data.frame() %>% View()
# 10 nice looking examples to share with Peter & Albin




# 7. FOOTPRINTS ####
# ------------------
# make enhancer midpoints
enh_midpoint <- rowRanges(enhancers) %>% GenomicRanges::resize(width=1, fix='center', use.names=T)
enh_midpoint_byTxType <- split(enh_midpoint, enh_midpoint$txType_TAIR10)

# 7a. CAGE footprint
cage_at_enh <- tidyMetaProfile(sites=enh_midpoint, forward=cage_plus, reverse=cage_minus, upstream=400, downstream=400, trimLower=0.1, trimUpper=0.9)

gg_cage <- cage_at_enh %>%
  gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
  mutate('signal'=factor(signal, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_line(alpha=.4, lwd=1) +
         geom_ma(lwd=.75, n=15, lty=1) +
         facet_grid(.~signal) +
         cowplot::theme_cowplot() +
         theme(aspect.ratio=1, legend.position=c(.02, .8)) + scale_color_brewer(palette='Set1', direction=-1) +
         labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'), y='CAGE TPM\n(15bp moving average)')

# 7b. RNAseq footprint
rnaseq_at_enh <- tidyMetaProfile(sites=enh_midpoint, forward=rnaseq_forward, reverse=rnaseq_reverse, upstream=750, downstream=750, sumFun=colMedians, trimLower=0.01, trimUpper=0.99)

rnaseq_at_enh %>%
  gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
  mutate('signal'=factor(signal, levels=c('WT', 'LSM8', 'RRP4', 'DM'))) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_line(lwd=1, alpha=.4) +
         geom_smooth() +
         facet_grid(~signal) +
         cowplot::theme_cowplot() +
         theme(aspect.ratio=1) + scale_color_brewer(palette='Set1', direction=-1) +
         labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'),
              y='RNAseq RPM (median, smoothened)', title='RNAseq at CAGE enhancers')

# 7b. ALL ChIP-Seq
      # make enhancer_midpoints with seqinfo compatible with chipseq
      enh_midpoint_4chipseq <- enh_midpoint
      seqlevels(enh_midpoint_4chipseq) <- seqlevels(chipseq$Bewick2016_H3K27me3)
      seqlengths(enh_midpoint_4chipseq) <- seqlengths(chipseq$Bewick2016_H3K27me3)
      
      # split by annotation
      enh_midpoint_4chipseq_byTxType <- split(enh_midpoint_4chipseq, enh_midpoint_4chipseq$txType_TAIR10)
      
      # keep only interesting marks
      hmark_albin <- chipseq_rpm_factor$Histone_mark %in% c('H3K4me3', 'H3K4me1', 'H3K27ac', 'PolII')
      
      subset(chipseq_rpm_factor, Histone_mark %in% c('H3K4me3', 'H3K4me1', 'H3K27ac'))
      
      # compute profiles
      histone_at_enh <- tidyMetaProfile(sites=enh_midpoint_4chipseq, forward=chipseq[hmark_albin], reverse=NULL, upstream=1000, downstream=1000, trimLower=0.1, trimUpper=0.9)

      # normalize with RPM factor
      histone_at_enh %<>% left_join(select(chipseq_rpm_factor, Sample, RPM_factor), by=c('signal'='Sample'))
      
      histone_at_enh %<>%
        separate(signal, c('author', 'mark'), sep='_') %>%
        mutate('sense_RPM'=sense/RPM_factor)
      
      # plot
        histone_at_enh %>%
        subset(author %!in% c('Bewick2016', 'Inagaki2017', 'Luo2013', 'Zhu2015')) %>%
        mutate('type'=ifelse(mark == 'PolII', 'PolII', 'Histones marks')) %>%
        subset(type=='Histones marks') %>%
        ggplot(aes(x=pos0, y=sense_RPM, col=mark, group=interaction(author, mark))) +
               geom_line() +
               geom_vline(xintercept=0, lty=2) +
               facet_wrap(~type, scales='free_y') +
               cowplot::theme_cowplot() + theme(aspect.ratio=1, strip.background=element_blank()) +
               scale_color_brewer(palette='Dark2', name='') +
               labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'),
                    y='ChIP-seq RPM')
        
      # recover authors
      histone_at_enh %>%
        subset(author %!in% c('Bewick2016', 'Inagaki2017', 'Luo2013', 'Zhu2015')) %>%
        subset(mark != 'PolII') %$%
        table(author, mark)
        
      histone_at_enh %>%
        subset(author %!in% c('Bewick2016', 'Inagaki2017', 'Luo2013', 'Zhu2015', 'Chen2017')) %>%
        subset(mark != 'PolII') %>%
        subset(author == 'vanDijk2010') %>%
        ggplot(aes(x=pos0, y=sense_RPM, col=mark, group=interaction(author, mark))) +
        geom_line() +
        cowplot::theme_cowplot() +
        scale_color_brewer(palette='Dark2', name='')

# 7c. DNaseI & MNaseI
dnase_at_enh <- tidyMetaProfile(sites=enh_midpoint, forward=dmnase$DNaseI, reverse=NULL, upstream=750, downstream=750, trimLower=0.01, trimUpper=0.99)
mnase_at_enh <- tidyMetaProfile(sites=enh_midpoint, forward=dmnase$MNaseI, reverse=NULL, upstream=750, downstream=750, trimLower=0.01, trimUpper=0.99)

gg_dnase <- dnase_at_enh %>%
  ggplot(aes(x=pos0, y=sense)) +
         geom_line(lwd=1) +
         geom_vline(xintercept=0, lty=2) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'), y='Normalized signal')

gg_mnase <- mnase_at_enh %>%
  ggplot(aes(x=pos0, y=sense)) +
  geom_line(lwd=1) +
  geom_vline(xintercept=0, lty=2) +
  cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'), y='Normalized signal')

gg_dnase + gg_mnase + plot_layout(ncol=1, nrow=2)

# 7d. GRO-seq
groseq_at_enh <- tidyMetaProfile(sites=enh_midpoint, forward=gro_p, reverse=gro_m, upstream=500, downstream=500, trimLower=0.01, trimUpper=0.99)

gg_gro <- groseq_at_enh %>%
  mutate('anti'=-anti) %>%
  gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
  ggplot(aes(x=pos0, y=score, fill=direction)) +
         geom_area(alpha=.7) +
         geom_line(col='black', lwd=.2) + 
         geom_vline(xintercept=0, lty=2) +
         geom_hline(yintercept=0, lwd=.2) +
         facet_wrap(~signal) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         scale_fill_brewer(palette='Set1', direction=-1) +
         labs(x=paste0('Enhancer midpoints (N=', length(enh_midpoint),')'),
              y='Normalized signal')
gg_gro

# 7e. plot all
gg_cage + gg_gro / (gg_dnase + gg_mnase + gg_hist) + plot_layout(nrow=2)



# 8. HOW CLOSE TO EXONS ARE ENHANCERS ? ####
# ------------------------------------------
# 8a. get all exons
txdb <- TxDb.Athaliana.BioMart.plantsmart28
exons <- exons(txdb)
      # fix seqinfo
      seqlevelsStyle(exons) <- seqlevelsStyle(myseqinfo)[1]
      seqlevels(exons, pruning.mode='coarse') <- seqlevels(myseqinfo)
      seqinfo(exons) <- myseqinfo

# 8b. find closest exon on the left side of enhancer midpoint
closest_exon_left <- exons[follow(x=enh_midpoint, subject=exons, ignore.strand=T)]
enhancers_gr$distance_left_exon <- distance(x=enh_midpoint, y=closest_exon_left, ignore.strand=T)

# 8d. find closest exon on the right side of enhancer midpoint
closest_exon_right <- exons[precede(x=enh_midpoint, subject=exons, ignore.strand=T)]
enhancers_gr$distance_right_exon <- distance(x=enh_midpoint, y=closest_exon_right, ignore.strand=T)

# 8e. plot left and right distances
min(enhancers_gr$distance_right_exon)

gg_dist <- enhancers_gr %>%
  as.data.frame() %>%
  ggplot() +
  geom_vline(xintercept=c(-500, 0, 500), lty=2, col=c('grey60', 'black', 'grey60')) +
  geom_density(aes(x=-distance_left_exon, col=txType_TAIR10)) +
  geom_density(aes(x=distance_right_exon, col=txType_TAIR10)) +
  xlim(-1000, 1000) +
  cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  scale_color_brewer(palette='Set1', name='Enhancer location', direction=-1) +
  labs(x='Distance to closest exon (bp)', title='Enhancer surroundings', caption='X-axis limited to +/-1Kb')
  
# 8f. what is the intron length in Arabidopsis?
introns <- intronicParts(txdb) %>%
  as.data.frame() %>%
  subset(width > 1)

gg_intron <- ggplot(introns, aes(x=width)) +
         geom_density() +
         scale_x_log10() +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(title='Intron length in Arabidopsis',
              x='Length (bp)')

gg_dist + gg_intron + plot_layout(nrow=2)

# EOF #

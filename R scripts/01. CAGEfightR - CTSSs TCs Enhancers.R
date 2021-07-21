#### Arabidopsis
#### Axel Thieffry - revised May 2018
set.seed(42)
library(tidyverse)
library(tidylog)
library(magrittr)
library(stringr)
library(reshape2)
library(CAGEfightR)
library(TeMPO)
library(ggpubr)
library(ggridges)
library(ggalluvial)
library(RColorBrewer)
library(BiocParallel)
library(Gviz)
library(gt)
register(MulticoreParam(workers=4))
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(Matrix)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(org.At.tair.db)
odb <- org.At.tair.db
options(scipen=10) # disable scientific notation
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename

setwd('masked_path')

# 1. PREPARATION ####
# -------------------
# 1a. get seqinfo from BSgenome
if(FALSE) {genome <- BSgenome.Athaliana.TAIR.TAIR9
           myseqinfo <- seqinfo(genome)
           genome(myseqinfo) <- rep('TAIR10', 7)
           saveRDS(myseqinfo, file='~/masked_path/myseqinfo.rds')}

myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')

# 1b. get sample names
sampleNames <- list.files(path='~/masked_path/bw_files_raw_chrWchr/', pattern='_0_.*.plus.chrWchr.bw', full.names=F) %>%
               str_replace('.raw.plus.chrWchr.bw', '')

# 1c. make design, and relevel 'wt' and timepoint '0' as references
if(FALSE) {colData <- data.frame(row.names=sampleNames,
                                 'Name'=sampleNames,
                                 'genotype'=c(rep('hen2', 3), rep('rrp4', 3), rep('wt', 3)))
           colData$genotype %<>% relevel(ref='wt')
           saveRDS(colData, file='~/masked_path/colData_TSSstory.rds')
           }

colData <- readRDS('~/masked_path/colData_TSSstory.rds')


# 2. CTSS-LEVEL ANALYSIS ####
# ---------------------------
if(FALSE){
# 2a. Create CTSS object
    # locate files on hard drive
    bw_plus   <- list.files(path='~/masked_path/bw_files_raw_chrWchr', pattern='_0_.*.plus.chrWchr.bw', full.names=T)
    bw_minus  <- list.files(path='~/masked_path/bw_files_raw_chrWchr', pattern='_0_.*.minus.chrWchr.bw', full.names=T)
    # create two named BigWigFileList-objects
    bw_plus %<>% BigWigFileList()
    bw_minus %<>% BigWigFileList()
    names(bw_plus) <- names(bw_minus) <- sampleNames
    # quantify CTSS accross samples
    CTSSs <- quantifyCTSSs(plusStrand=bw_plus, minusStrand=bw_minus, design=colData, genome=myseqinfo, tileWidth=30000000L)

# 2b. remove noise and normalize in TPM
    # calculate support based on readcount
    CTSSs <- calcSupport(CTSSs, inputAssay='counts', outputColumn='support', unexpressed=0)
    # plot CTSS support
    rowRanges(CTSSs)$support %>%
      table() %>%
      as.data.frame() %>%
      set_colnames(c('support', 'Freq')) %>%
      ggplot(aes(x=support, y=Freq/1000000)) +
             geom_bar(stat='identity') +
             geom_vline(xintercept=2.5, col='red', lty=2) + theme(aspect.ratio=1) +
             labs(title='CTSS support', x='support (libraries with readcount > 0)', y='CTSSs (millions)',
                  subtitle='CTSSs kept if supported at least in 3 libraries')
    # remove CTSS supported by less than 3 libraries
    CTSSs <- subset(CTSSs, support >= 3)
    # normalize to TPM: add supported_TPM assay as well as a 'totalTags' column to colData
    CTSSs <- calcTPM(CTSSs, inputAssay='counts', outputAssay='TPM', outputColumn='totalTags')
    # sanitycheck: total TPM is 1 million
    assay(CTSSs, 'TPM') %>% colSums() %>% as.data.frame()
    # plot totalTags (raw and after filtering for support) per library
    colData(CTSSs) %>%
      as.data.frame() %>%
      melt(id.vars=c('Name', 'genotype')) %>%
      ggplot(aes(x=Name, fill=genotype, y=value/1000000)) + geom_bar(stat='identity') +
             theme(aspect.ratio=1, axis.text.x=element_text(angle=90, vjust=0.5)) +
             labs(y='Total tags (millions)', title='Total mapped tags', x='Sample', subtitle='After support filtering') +
             scale_fill_brewer(palette='Dark2')
    # calculate pooled TPM: as score in the rowRanges
    CTSSs <- calcPooled(CTSSs, inputAssay='TPM', outputColumn='score')
    # save CTSSs
    saveRDS(CTSSs, file='~/masked_path/SE_CTSSs_1count_min3lib_TSSstory.rds')
    }

CTSSs <- readRDS('~/masked_path/SE_CTSSs_1count_min3lib_TSSstory.rds')


# 2c. Unidirectional tag clustering (TSSs)
    # TC clustering
    TCs <- clusterUnidirectionally(CTSSs, pooledCutoff=0, mergeDist=20)
    # quantification of TSSs (counts)
    TCs <- quantifyClusters(CTSSs, clusters=TCs, inputAssay='counts')
    # calculate TPM
    TCs <- calcTPM(TCs, totalTags='totalTags', inputAssay='counts', outputAssay='TPM')
    # calculate TPM support
    TCs <- calcSupport(TCs, inputAssay='TPM', outputColumn='TPMsupport', unexpressed=1)
    # plot TC support
    as.data.frame(mcols(TCs)) %>%
      subset(TPMsupport != 0) %>%
      ggplot(aes(x=as.factor(TPMsupport))) + geom_bar(stat='count') +
             geom_vline(xintercept=1.5, col='red', lty=2) +
             theme(aspect.ratio=1) +
             labs(x='Support (libraries with TPM >= 1)', y='TCs',
                  title='TCs support', subtitle='TPM threshold >= 1\nTCs with zero support not shown')
    # save all TCs before support filtering
    export.bed(rowRanges(TCs), 'all_TCs_before_TPMsupport_filtering.bed')
    # remove TCs supported by less than 3 libraries with TPM >= 1
    TCs_all <- TCs
    TCs <- subset(TCs, TPMsupport >= 2)
    
    # PCA on TCs for timepoints = 0
    pca <- assay(TCs, 'TPM')
    pca <- pca[apply(pca, 1, function(x) var(x) != 0), ]
    pca <- prcomp(pca %>% t(), scale.=T, center=T)
      
    pca_importance <- summary(pca)$importance[2,] * 100
    
    blue='#71A0C5'; orange='orange'; red='#BF717F'
    
    pca$x %>%
      as.data.frame() %>%
      rownames_to_column('sample') %>%
      separate('sample', c('genotype', 'time0', 'replicate'), sep='_') %>%
      mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
      ggplot(aes(x=PC1, y=PC2, col=genotype, shape=replicate)) +
             geom_point(size=3) + cowplot::theme_cowplot() + theme(aspect.ratio=1, axis.line=element_blank()) +
             labs(x=paste0('PC1 (', round(pca_importance[1], 1), '%)'), y=paste0('PC2 (', round(pca_importance[2], 1), '%)'),
                  title='PCA of CAGE TCs', subtitle='TPM scaled and centered') +
             scale_color_manual(values=c(blue, orange, red)) + cowplot::panel_border(colour='black', size=1)

# 2d. Bidirectional tag clustering (aka enhancer candidates)
    enhancers <- clusterBidirectionally(CTSSs, window=201)
    # calculate bidirectionality support (must have counts > 0 in both arms)
    enhancers <- calcBidirectionality(enhancers, samples=CTSSs)
    # plot enhancer bidirectionality support
    as.data.frame(mcols(enhancers)) %>%
      ggplot(aes(x=as.factor(bidirectionality))) +
             geom_bar(stat='count') +
             geom_vline(xintercept=3.5, col='red', lty=2) +
             labs(x='Support (libraries with bidirectional signal)', y='Enhancers',
                  title='Enhancer bidirectionality support',
                  subtitle='Sliding window length: 400bp\nPooled bidirectionality threshold: 0.99\nBidirectionality support: count > 0 in both arms') +
             theme(aspect.ratio=1)
    # keep only enhancers that are bidirectional in at least in 3 samples
    enhancers <- subset(enhancers, bidirectionality > 2)
    # quantification of enhancers (counts)
    enhancers <- quantifyClusters(CTSSs, clusters=enhancers, inputAssay='counts')
    # calculate TPM (but don't filter on it)
    enhancers <- calcTPM(enhancers, totalTags='totalTags', inputAssay='counts', outputAssay='TPM')
    # calculate pooled TPM as a new column, as well as TPM for each genotype
    enhancers <- calcPooled(enhancers, inputAssay='TPM', outputColumn='TPMpooled')
    rowRanges(enhancers)$TPM_wt <- calcPooled(subset(enhancers, select=genotype=='wt'), inputAssay='TPM', outputColumn='TPM_wt') %>% rowRanges() %$% .$TPM_wt
    rowRanges(enhancers)$TPM_hen2 <- calcPooled(subset(enhancers, select=genotype=='hen2'), inputAssay='TPM', outputColumn='TPM_hen2') %>% rowRanges() %$% .$TPM_hen2
    rowRanges(enhancers)$TPM_rrp4 <- calcPooled(subset(enhancers, select=genotype=='rrp4'), inputAssay='TPM', outputColumn='TPM_rrp4') %>% rowRanges() %$% .$TPM_rrp4



# 3. ANNOTATION ####
# ------------------
if(FALSE){
# 3a. make annotation from ARAPORT11 GFF3
    araport_txdb <- makeTxDbFromGFF(file='~/masked_path/Araport11_GFF3_genes_transposons.201606.gff3', format='gff3', chrominfo=myseqinfo)

# 3b. Building personal hierarchy of TxTypes (first to appear is higher priority)
# NB: using a GRangesList for the txModels will make assignTxType to not use any argument/option, so everything has to be build here
    # make transcript antisense for TAIR10 and ARAPORT11
    Antis_tair <- transcripts(txdb) %>% invertStrand()
    Antis_araport <- transcripts(araport_txdb) %>% invertStrand()
    # make reverse strand (aka PROMPT regions, antisense to gene and up to 400bp upstream of TSS) for TAIR10 and ARAPORT11
    Reverse_tair <- suppressWarnings( promoters(txdb, upstream=400, downstream=0) %>% trim() %>% invertStrand() )
    Reverse_araport <- suppressWarnings( promoters(araport_txdb, upstream=400, downstream=0) %>% trim() %>% invertStrand() )
    # build custom hierarchy for TAIR10
    custom_hierarchy_tair10 <- GRangesList('promoter' =  suppressWarnings( promoters(txdb, upstream=100, downstream=100) %>% trim() %>% granges() ),
                                           'fiveUTR' =   fiveUTRsByTranscript(txdb) %>% unlist() %>% granges(),
                                           'threeUTR' =  threeUTRsByTranscript(txdb) %>% unlist() %>% granges(),
                                           'CDS' =       cds(txdb) %>% granges(),
                                           'exon' =      exons(txdb) %>% granges(),
                                           'intron' =    intronsByTranscript(txdb) %>% unlist() %>% granges(),
                                           'proximal' =  suppressWarnings( promoters(txdb, upstream=400, downstream=0) %>% trim() %>% granges() ),
                                           'reverse' =   Reverse_tair,
                                           'antisense' = Antis_tair)
      
    # build custom hierarchy for ARAPORT11
    custom_hierarchy_araport <- GRangesList('promoter' =  suppressWarnings( promoters(araport_txdb, upstream=100, downstream=100) %>% trim() %>% granges() ),
                                            'fiveUTR' =   fiveUTRsByTranscript(araport_txdb) %>% unlist() %>% granges(),
                                            'threeUTR' =  threeUTRsByTranscript(araport_txdb) %>% unlist() %>% granges(),
                                            'CDS' =       cds(araport_txdb) %>% granges(),
                                            'exon' =      exons(araport_txdb) %>% granges(),
                                            'intron' =    intronsByTranscript(araport_txdb) %>% unlist() %>% granges(),
                                            'proximal' =  suppressWarnings( promoters(araport_txdb, upstream=400, downstream=0) %>% trim() %>% granges() ),
                                            'reverse' =   Reverse_araport,
                                            'antisense' = Antis_araport)
    # fix seqinfo
    seqlevelsStyle(custom_hierarchy_tair10) <- seqlevelsStyle(custom_hierarchy_araport) <- seqlevelsStyle(myseqinfo)
    seqinfo(custom_hierarchy_tair10) <- seqinfo(custom_hierarchy_araport) <- myseqinfo
    seqinfo(custom_hierarchy_tair10) ; seqinfo(custom_hierarchy_araport)
           
    # save
    saveRDS(custom_hierarchy_tair10, file='~/masked_path/custom_annotation_hierarchy_TAIR10.rds')
    saveRDS(custom_hierarchy_araport, file='~/masked_path/custom_annotation_hierarchy_ARAPORT11.rds')
}
                     
custom_hierarchy_tair10 <- readRDS('~/masked_path/custom_annotation_hierarchy_TAIR10.rds')
custom_hierarchy_araport <- readRDS('~/masked_path/custom_annotation_hierarchy_ARAPORT11.rds')

# build antisense txType custom hierarchy
if(FALSE){
    antisense_tmp <- invertStrand(custom_hierarchy_tair10)
    names(antisense_tmp) <- paste0('antisense_', names(antisense_tmp))
    # merge after sense hierarchy (order is important)
    custom_annotation_hierarchy_extended_TAIR10 <- c(custom_hierarchy_tair10, antisense_tmp)
    # remove unnecessary categories
    custom_annotation_hierarchy_extended_TAIR10$antisense <- NULL # because it will be annotated
    custom_annotation_hierarchy_extended_TAIR10$antisense_antisense <- NULL # because redundant
    custom_annotation_hierarchy_extended_TAIR10$reverse <- NULL # because it is non-genic
    custom_annotation_hierarchy_extended_TAIR10$antisense_promoter <- NULL # because it will either be 5'UTR or reverse
    custom_annotation_hierarchy_extended_TAIR10$antisense_reverse <- NULL # because reverse already captures them
    custom_annotation_hierarchy_extended_TAIR10$antisense_proximal <- NULL # because reverse already captures them
    data.frame('antisense_annot'=names(custom_annotation_hierarchy_extended_TAIR10))
    saveRDS(custom_annotation_hierarchy_extended_TAIR10, '~/masked_path/custom_annotation_hierarchy_TAIR10_extended_antisense.rds')
}
                     
custom_hierarchy_tair10_extended_antisense <- readRDS('~/masked_path/custom_annotation_hierarchy_TAIR10_extended_antisense.rds')

# 3d. annotate TCs
    # sense TCs up to 400bp upstream of the TX start will be assigned to that TX
    seqlevelsStyle(txdb) <- seqlevelsStyle(myseqinfo)
    TCs <- assignTxID(TCs, txModels=txdb, outputColumn='txID', swap='thick', upstream=400, downstream=0)
    TCs_all <- assignTxID(TCs_all, txModels=txdb, outputColumn='txID', swap='thick', upstream=400, downstream=0)
    # assign TxType for TAIR10 and ARAPORT11
    TCs_all <- assignTxType(TCs_all, txModels=custom_hierarchy_tair10, outputColumn='txType_TAIR10', swap='thick', noOverlap='intergenic')
    TCs <- assignTxType(TCs, txModels=custom_hierarchy_tair10, outputColumn='txType_TAIR10', swap='thick', noOverlap='intergenic')
    TCs_all <- assignTxType(TCs_all, txModels=custom_hierarchy_araport, outputColumn='txType_ARAPORT11', swap='thick', noOverlap='intergenic')
    TCs <- assignTxType(TCs, txModels=custom_hierarchy_araport, outputColumn='txType_ARAPORT11', swap='thick', noOverlap='intergenic')
    # assign extended antisense TxType for TAIR10
    TCs <- assignTxType(TCs, txModels=custom_hierarchy_tair10_extended_antisense, outputColumn='txType_TAIR10extended', swap='thick', noOverlap='intergenic')
    TCs_all <- assignTxType(TCs_all, txModels=custom_hierarchy_tair10_extended_antisense, outputColumn='txType_TAIR10extended', swap='thick', noOverlap='intergenic')
    # plot pooled TCs annotation & pooled expression : pan-experiment
    TC_pan_annotation_df <- as.data.frame(rowData(TCs)) %>%
                            select(thick.names, score, TPMsupport, txType_TAIR10, txType_ARAPORT11) %>% as_tibble()
    
    TC_all_pan_annotation_df <- as.data.frame(rowData(TCs_all)) %>%
      select(thick.names, score, TPMsupport, txType_TAIR10, txType_ARAPORT11) %>% as_tibble()
    
    TC_pan_annotation_df %<>%
      melt(id.vars=c('thick.names', 'score', 'TPMsupport'), value.name='txType', variable.name='reference') %>%
      separate(reference, c('bla', 'reference'), sep='_') %>%
      select(-bla) %>% as_tibble()
    
    TC_all_pan_annotation_df %<>%
      melt(id.vars=c('thick.names', 'score', 'TPMsupport'), value.name='txType', variable.name='reference') %>%
      separate(reference, c('bla', 'reference'), sep='_') %>%
      select(-bla) %>% as_tibble()
    
    TC_pan_annotation_df$category <- ifelse(TC_pan_annotation_df$txType %in% c('intergenic', 'reverse', 'proximal'), 'non-genic', 'genic')
    TC_pan_annotation_df$txType %<>% factor(levels=c('intergenic', 'antisense', 'reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))
    
    TC_all_pan_annotation_df$category <- ifelse(TC_all_pan_annotation_df$txType %in% c('intergenic', 'reverse', 'proximal'), 'non-genic', 'genic')
    TC_all_pan_annotation_df$txType %<>% factor(levels=c('intergenic', 'antisense', 'reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))
    
    ggplot(TC_pan_annotation_df, aes(x=txType, fill=reference)) + geom_bar(stat='count', position=position_dodge(), col='black', alpha=.3) +
           theme_ridges(center_axis_labels=T, grid=F) +
           theme(legend.position='top', axis.text.y=element_text(colour=c('red', 'black', rep('red', 2), rep('black', 6)))) +
           scale_fill_brewer(palette='Set1') +
           scale_color_manual(values=c('black', 'red')) +
           labs(title='CAGE t=0 TSSs: annotation and expression', subtitle='TPM >= 1.0 in minimum 3 libraries',
                x='', y=paste0('TSSs (N=', length(TCs), ')') ) +
           coord_flip() +
    
    ggplot(TC_all_pan_annotation_df, aes(x=txType, fill=reference)) + geom_bar(stat='count', position=position_dodge(), col='black', alpha=.3) +
      theme_ridges(center_axis_labels=T, grid=F) +
      theme(legend.position='top', axis.text.y=element_text(colour=c('red', 'black', rep('red', 2), rep('black', 6)))) +
      scale_fill_brewer(palette='Set1') +
      scale_color_manual(values=c('black', 'red')) +
      labs(title='CAGE t=0 TSSs: annotation and expression', subtitle='TPM >= 1.0 in minimum 3 libraries',
           x='', y=paste0('TSSs (N=', length(TCs), ')') ) +
      coord_flip() +
      patchwork::plot_layout(ncol=2)
    
    gg_expression <- ggplot(TC_pan_annotation_df, aes(x=score, y=txType, fill=reference)) + geom_density_ridges(alpha=.3) + scale_x_log10() +
           theme_ridges(center_axis_labels=T, grid=F) +
           theme(legend.position='none', axis.text.y=element_text(colour=c('red', 'black', rep('red', 2), rep('black', 6)))) +
           scale_fill_brewer(palette='Set1') +
           labs(y='', x=paste0('TSS expression (log TPM) (N=', length(TCs), ')') ) +
           scale_y_discrete(expand = c(0.01, 0))
    
    ggarrange(gg_annotation, gg_expression, ncol=2, nrow=1)
    
    # TAIR10 -vs- ARAPORT11: what becomes what? (alluvial plot)
    alluvial_annotation <- TCs %>%
            rowData() %>%
            as.data.frame() %$%
            table(txType_TAIR10, txType_ARAPORT11) %>%
            as.data.frame()

    alluvial_annotation$identical <- with(alluvial_annotation, txType_TAIR10 == txType_ARAPORT11)
    
    alluvial_annotation %>%
      filter(Freq > 100, txType_TAIR10 != txType_ARAPORT11) %>%
      droplevels() %>%
      ggplot(aes(axis1=txType_TAIR10, axis2=txType_ARAPORT11, y=Freq)) +
             geom_alluvium(aes(fill=txType_TAIR10), alpha=.6) + geom_stratum(width=.2) +
             geom_text(stat='stratum', label.strata=T) +
             scale_x_discrete(limits=c('TAIR10', 'ARAPORT11'), expand=c(0, 0)) +
             scale_y_continuous(expand=c(0, 0)) +
             theme(legend.position='none', aspect.ratio=1, panel.grid=element_blank(), panel.background=element_blank()) +
             scale_fill_brewer(palette='Set1', name='TAIR10', direction=-1) +
             labs(x='', y='TSS annotation categories', title='Differences of annotation: TAIR10 -vs- ARAPORT11',
                  subtitle='Only couples with > 100 differences are shown')
    
    # table version of alluvial plot
    TCs %>%
      rowData() %>%
      as.data.frame() %$%
      table(txType_TAIR10, txType_ARAPORT11) %>%
      as.data.frame.matrix() %>%
      rownames_to_column('TAIR10') %>%
      gt() %>% tab_header(title='TAIR10 to ARAPORT11', subtitle='TSS annotation')
    
    # TAIR10 antisense : what becomes what? (alluvial plot)
    alluvial_antisense <- rowData(TCs) %>%
                          as.data.frame() %$%
                          table(txType_TAIR10, txType_TAIR10extended) %>%
                          as.data.frame() %>%
                          subset(Freq != 0)
    
    alluvial_antisense$txType_TAIR10 %<>% as.character()
    alluvial_antisense$txType_TAIR10extended %<>% as.character()
    
    alluvial_antisense %>%
      filter(txType_TAIR10 != txType_TAIR10extended) %>%
      droplevels() %>%
      ggplot(aes(axis1=txType_TAIR10, axis2=str_remove(txType_TAIR10extended, 'antisense_'), y=Freq)) +
             geom_alluvium(aes(fill=txType_TAIR10), alpha=.6) + geom_stratum(width=.2) +
             geom_text(stat='stratum', label.strata=T) +
             scale_x_discrete(limits = c('TAIR10', 'TAIR10 antisense'), expand=c(0, 0)) +
             scale_y_continuous(expand=c(0, 0)) +
             theme(legend.position='none', aspect.ratio=1, panel.grid=element_blank(), panel.background=element_blank()) +
             scale_fill_brewer(palette='Set1', name='Category', direction=-1) +
             labs(x='', y='TSS annotation categories', title='CAGE TSSs: extended TAIR10 antisense annotation')


# 3d. annotate bidirectional clusters (aka enhancer candidates)
    enhancers <- assignTxType(enhancers, txModels=txdb, outputColumn='txType_TAIR10')
    # plot enhancer annotation
    rowRanges(enhancers) %>%
      mcols() %>%
      as.data.frame() %>%
      ggplot(aes(txType_TAIR10)) + geom_bar(stat='count')
    # only keep intronic and intergenic enhancers
    enhancers %<>% subset(txType_TAIR10 %in% c('intron', 'intergenic')) # 125 enhancers
    # only keep enhancers on canonical chromosomes
    seqlevels(enhancers, pruning.mode='coarse') <- paste0('Chr', 1:5) # 113
    # save pooled enhancers as bedfile and Rdata
    saveRDS(enhancers, file='~/masked_path/SE_Enhancers.rds')
    export.bed(rowRanges(enhancers), '~/masked_path/Enhancers_CAGEfightR.bed')
    # export enhancers for supplemental dataset
    rowRanges(enhancers) %>%
      as.data.frame() %>%
      rownames_to_column('TC_id') %>%
      select(-thick.end, -thick.width, -thick.names, -TPMpooled:-TPM_rrp4) %>%
      set_colnames(c('TC_id', 'chr', 'start', 'end', 'width', 'strand', 'pooled_TPM', 'midpoint', 'balance', 'bidirectionality', 'TAIR10_annotation')) %>%
      WriteXLS::WriteXLS('~/masked_path/bidirectional_CAGE_TCs.xlsx', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)


# 4. TSS SHAPE STATISTICS ####
# ----------------------------
TCs <- calcShape(TCs, pooled=CTSSs, outputColumn='IQR', shapeFunction=shapeIQR, lower=0.25, upper=0.75)
TCs <- calcShape(TCs, pooled=CTSSs, outputColumn='Entropy', shapeFunction=shapeEntropy)

rowRanges(TCs) %>%
  as.data.frame() %>%
  ggplot(aes(x=IQR)) +
         geom_histogram(binwidth=1) + xlim(NA, 150) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1) +
         labs(x='IQR (bp)', title='CAGE TCs IQR')
dev.off()




# 5. ANNOTATION AT GENE-LEVEL ####
# --------------------------------
# 5a. Gene ID with TAIR10
  # normal sense genes
  TCs <- assignGeneID(TCs, geneModels=txdb, outputColumn='geneID', upstream=400, downstream=0, swap='thick')
  # antisense gene (must extend to the PROMPT region before inverting strand!)
  extended_genes <- genes(txdb)
  promoters_genes <- promoters(extended_genes, upstream=400, downstream=0) %>% trim()
  start(extended_genes) <- ifelse(strand(extended_genes)=='+', start(promoters_genes) , start(extended_genes))
  end(extended_genes) <- ifelse(strand(extended_genes)=='-', end(promoters_genes), end(extended_genes))
  # IGV export for sanity check
  export.bed(invertStrand(extended_genes), 'extended_genes_antisense.bed')
  TCs <- assignGeneID(TCs, geneModels=invertStrand(extended_genes), outputColumn='geneID_anti', upstream=0, downstream=0, swap='thick')
  # get gene symbols
  symbols <- mapIds(odb, keys=rowRanges(TCs)$geneID, keytype='TAIR', column='SYMBOL')
  symbols_anti <- mapIds(odb, keys=rowRanges(TCs)$geneID_anti, keytype='TAIR', column='SYMBOL')
  rowRanges(TCs)$symbol <- as.character(symbols)
  rowRanges(TCs)$symbol_anti <- as.character(symbols_anti)
  # save pooled TCs as bedfile and Rdata
  saveRDS(TCs, file='~/masked_path/SE_TCs_TPM1_min3lib_TSSstory.rds')
  export.bed(rowRanges(TCs), '~/masked_path/TCs_TPM1_min3lib.bed')

# 5b. Quantify expression at Gene-level
  genelevel <- quantifyGenes(TCs, genes='geneID', inputAssay='counts')
  saveRDS(genelevel, file='~/masked_path/SE_genelevel_TSSstory.rds')

# export TCS for supplementary dataset
TCs_export <- rowRanges(TCs) %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  set_colnames(c('TC_id', 'chr', 'start', 'end', 'width', 'strand', 'pooled_TPM', 'TC_peak', 'TC_peak2', 'thick.width', 'thick.names', 'TPM_support', 'txID_TAIR10', 'txType_TAIR10', 'txType_ARAPORT11', 'txType_TAIR10_anti', 'IQR', 'Entropy', 'geneID_TAIR10', 'geneID_TAIR10_anti', 'gene_TAIR10_symbol', 'gene_TAIR10_symbol_anti')) %>%
  select(-TC_peak2, -thick.width, -thick.names, -Entropy, -gene_TAIR10_symbol, -gene_TAIR10_symbol_anti)

  WriteXLS::WriteXLS(x=TCs_export, ExcelFileName='~/masked_path/CAGE_TCs.xlsx', row.names=F, col.names=T)

# export TPM expression matrix for supplemental dataset
TCs_tpm_export <- assay(TCs, 'TPM') %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  set_colnames(str_remove(colnames(.), '_0'))
  
  write.table(TCs_tpm_export, '~/masked_path/CAGE_TC_expression_matrix_TPM.txt', append=F, quote=F, dec='.', row.names=F, col.names=T, sep='\t')



# 6. FILTERING CLUSTERS BASED ON COMPOSITION - to use for calling DTUs ####
# -------------------------------------------------------------------------
# remove TSSs not belonging to any gene
intragenicTSSs <- subset(TCs, !is.na(geneID))
# calculate composition: the number of samples expression TSSs above 10% of the total gene expression
intragenicTSSs <- calcComposition(intragenicTSSs, inputAssay='counts', outputColumn='composition', unexpressed=0.1, genes='geneID')
# plot
as.data.frame(rowRanges(intragenicTSSs)) %>%
  ggplot(aes(x=composition)) + geom_bar(stat='count') +
         labs(title='TSS composition',
              subtitle='with 10% total gene expression contribution thresholold') +
         theme(aspect.ratio=1)
# subset
intragenicTSSs <- subset(intragenicTSSs, composition >= 3)
# save
saveRDS(intragenicTSSs, file='~/masked_path/SE_intragenicTSSs_for_DTU_TSSstory.rds')

### EOF ###


#### Arabidopsis flg22 : FULL DE at TC level
#### Axel Thieffry - July 2018
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
'rename' <- dplyr::rename
options(scipen=999)
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                    if(length(idx) != 0) { GR[-idx]}
                                    else {GR}}
get_density <- function(x, y, ...) {dens <- MASS::kde2d(x, y, ...)
                                    ix <- findInterval(x, dens$x)
                                    iy <- findInterval(y, dens$y)
                                    ii <- cbind(ix, iy)
                                    return(dens$z[ii])}
setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE')



# 1. LOAD ALL INPUT FILES ####
# ----------------------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
# id mapping
idmapping <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/AGI_ID_mapping.rds')
# CTSSs
CTSSs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_CTSSs_1count_min3lib_TSSstory.rds')
# TCs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_TPM1_min3lib_TSSstory.rds')
# CAGE BIGWIG FILES
cage_p <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', full.names=T, pattern='_0.*plus'))
cage_m <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', full.names=T, pattern='_0.*minus'))
cage_names <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='_0.*plus') %>% str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names
# RNA-Seq BIGWIG FILES
rnaseq_forward <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123_fixed', pattern='Forward', full.names=T))
rnaseq_reverse <- BigWigFileList(list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123_fixed', pattern='Reverse', full.names=T))
rnaseq_names <- list.files(path='~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123', pattern='Forward.RPM.bw', full.names=F) %>% str_remove('_R123.Forward.RPM.bw')
names(rnaseq_forward) <- names(rnaseq_reverse) <- rnaseq_names
# GROseq Jacobsen
grohet_p <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/GRO-seq Nature Plants 2018/GROseq_bigwigs/GROseq_col_mergedAveraged_R123456.plus.bw')
grohet_m <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/GRO-seq Nature Plants 2018/GROseq_bigwigs/GROseq_col_mergedAveraged_R123456.minus.bw')
# GROcap Duttke
grocap_p <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/5GROcap.aln.unique.plus.bw')
grocap_m <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/5GROcap.aln.unique.minus.bw')
# GROseq Duttke
grodut_p <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/GROseq_mergedAveraged_R12.plus.bw')
grodut_m <- BigWigFile('~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/GROseq_mergedAveraged_R12.minus.bw')
# all GROseq at one
gro_p <- BigWigFileList(c('~/Dropbox/Axel_Arabidopsis_Flagellin/GRO-seq Nature Plants 2018/GROseq_bigwigs/GROseq_col_mergedAveraged_R123456.plus.bw',
                          '~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/5GROcap.aln.unique.plus.bw',
                          '~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/GROseq_mergedAveraged_R12.plus.bw'))
gro_m <- BigWigFileList(c('~/Dropbox/Axel_Arabidopsis_Flagellin/GRO-seq Nature Plants 2018/GROseq_bigwigs/GROseq_col_mergedAveraged_R123456.minus.bw',
                          '~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/5GROcap.aln.unique.minus.bw',
                          '~/Dropbox/Axel_Arabidopsis_Flagellin/Duttke/05.STAR_BigWigs/fixed_bws/GROseq_mergedAveraged_R12.minus.bw'))
names(gro_p) <- names(gro_m) <- c('GROseq_Hetzel', 'GROcap_Duttke', 'GROseq_Duttke')
# DHS regions
dhss <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/rrp4_lsm8_novogene/analysis/new_dhssumits_annotated_plantDHSonly.rds')
dhss_promoter <- subset(dhss, txType=='promoter') # 20,198



# 2. MAKE DESIGN AND CONTRAST MATRICES ####
# -----------------------------------------
# design matrix only for wt, hen2 and rrp4 at time 0
if(FALSE) {
          design <- model.matrix(~ genotype, data=colData(TCs)) %>% as.data.frame()
          colnames(design) <- make.names(colnames(design))
          annot_design <- data.frame(row.names  = rownames(design),
                                     'genotype' = c( rep('hen2',3), rep('rrp4',3), rep('wt',3) ))
          annot_design$genotype %<>% as.factor()
          pheatmap(design, cluster_rows=F, cluster_cols=F, cellwidth=30, cellheight=10, color=c('white', 'black'),
                   main='design: ~ genotype', annotation_row=annot_design, legend=F)
          saveRDS(design, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/limma_design_matrix.rds')
          }
design <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/limma_design_matrix.rds')
# contrasts
if(FALSE){
          contrasts <- diag(ncol(design))
          colnames(contrasts) <- colnames(design)
          contrasts <- contrasts[, -1]
          rownames(contrasts) <- colnames(design)
          pheatmap(contrasts, cluster_rows=F, cluster_cols=F, color=c('white', 'darkgreen'),
                   main='Contrast matrix', cellwidth=15, cellheight=15, legend=F, display_numbers=T, number_color='white', number_format='%s')
          saveRDS(contrasts, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/limma_contrast_matrix.rds')
          }
contrast <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/limma_contrast_matrix.rds')



# 3. LIMMA : TSS-level DE ####
# ----------------------------
# 3a. process DE
dge <- DGEList(assay(TCs, 'counts')) %>% calcNormFactors(method='TMM')
v   <- voom(dge, design=design, plot=F)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrast)
eb  <- eBayes(fit, robust=TRUE)
dt  <- decideTests(eb, adjust.method='BH', lfc=1, p.value=0.05)

# diagnostic plots
if(FALSE){
         par(pty='s')
           plotSA(eb, main='SA plot')
         dev.off()
         plotMDS(dge, method = 'bcv')
         }

# get all DEGs
coefs <- colnames(dt)
names(coefs) <- coefs

degs <- map(coefs, topTable, fit=eb, number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>%
  lapply(function(x) rownames_to_column(x, 'TC_id')) %>%
  bind_rows(.id='coef') %>%
  mutate('direction'=ifelse(logFC > 0, 'up', 'down')) %>%
  as_tibble()

# export DEGs for supplementary data
degs_export <- degs
    # rename coefficients
    degs_export$coef <- ifelse(degs_export$coef=='genotypehen2', 'hen2-4 vs. wt', 'rrpp4-2 vs. wt')
    WriteXLS::WriteXLS(x=degs_export, SheetNames='CAGE TCs DE', ExcelFileName='~/Dropbox/Flg22_CAGE/Exosome_TSS paper/Figures and layout/Supplementary Data/CAGE_DE.xlsx', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)

# add to TCs metacolumns
rowRanges(TCs)$genotyperrp4 <- as.data.frame(dt)$genotyperrp4
rowRanges(TCs)$genotypehen2 <- as.data.frame(dt)$genotypehen2

tt_genotypehen2 <- topTable(fit=eb, coef='genotypehen2', number=Inf, adjust.method='BH') %>%
  rownames_to_column('TC_id') %>%
  as_tibble()

tt_genotyperrp4 <- topTable(fit=eb, coef='genotyperrp4', number=Inf, adjust.method='BH') %>%
  rownames_to_column('TC_id') %>%
  as_tibble()

mcols <- mcols(TCs) %>%
  as.data.frame()

mcols %<>% left_join(select(tt_genotypehen2, -AveExpr, -t, -P.Value, -B), by=c('thick.names'='TC_id')) %>%
  rename('logFC_hen2'='logFC', 'adj.P.Val_hen2'='adj.P.Val') %>%
  left_join(select(tt_genotyperrp4, -AveExpr, -t, -P.Value, -B), by=c('thick.names'='TC_id')) %>%
  rename('logFC_rrp4'='logFC', 'adj.P.Val_rrp4'='adj.P.Val')

mcols(TCs) <- mcols

    # see if DOG1 has DE TCs
    rowRanges(TCs) %>%
      subset(geneID=='AT5G45830' | geneID_anti=='AT5G45830') %>%
      as.data.frame() %>%
      View()

# 3b. barplot - amount of DE by coefficient
deByCoef <- summary(dt) %>%
            as.data.frame() %>%
            spread_(key='Var1', value='Freq') %>%
            set_colnames(c('coefficient', 'down', 'nDE', 'up'))

deByCoef$coefficient %<>% factor(levels=c('genotypehen2', 'genotyperrp4'))

melt(deByCoef, id.vars='coefficient', variable.name='direction', value.name='deTSSs') %>%
  subset(direction != 'nDE') %>%
  ggplot(aes(x=coefficient, fill=direction, y=deTSSs)) +
         geom_rect(xmin=0, xmax=3.5, aes(ymin=-Inf, ymax=Inf), fill='blue', alpha=0.006) +
         geom_rect(xmin=3.5, xmax=6.5, aes(ymin=-Inf, ymax=Inf), fill='grey', alpha=0.006) +
         geom_rect(xmin=6.5, xmax=Inf, aes(ymin=-Inf, ymax=Inf), fill='orange', alpha=0.006) +
         geom_bar(stat='identity', position=position_dodge(), col='black') +
         labs(title='Number of DE TSSs',
              x='Coefficients', y='DE TSSs (LogFC > 1, FDR < 0.05)') +
         coord_flip() + theme_bw() + theme(aspect.ratio=1)

# 3c. VennDiagram - DE by coefficient
par(pty='s', mar=rep(0,4))
    vennDiagram(dt[,c('genotypehen2', 'genotyperrp4')], include=c('up', 'down'), counts.col=c('red', 'blue'), circle.col=c('darkviolet', 'darkorange'), names=c("genotype hen2", "genotype rrp4"), mar=rep(0, 4))
dev.off()
    
    # test overlap significant for all UP DE
    de_up_hen2 <- rownames(subset(as.data.frame(dt), genotypehen2 == 1)) # 796   # this is a vector of gene names (AGI)
    de_up_rrp4 <- rownames(subset(as.data.frame(dt), genotyperrp4 == 1)) # 1674  # this too
    go_up <- newGeneOverlap(listA=de_up_hen2, listB=de_up_rrp4, genome.size=n_distinct(names(TCs)))
    go_up <- testGeneOverlap(go_up)
    
    # test overlap significant for all DOWN DE
    de_down_hen2 <- rownames(subset(as.data.frame(dt), genotypehen2 == -1)) # 79
    de_down_rrp4 <- rownames(subset(as.data.frame(dt), genotyperrp4 == -1)) # 1008
    go_down <- newGeneOverlap(listA=de_down_hen2, listB=de_down_rrp4, genome.size=n_distinct(names(TCs)))
    go_down <- testGeneOverlap(go_down)
    
# advanced Venn diagram: logFC vs WT X-Y plot (Figure 6)
    # get all up-regulated TCs
    logFC_xyplot_df <- rowRanges(TCs) %>%
      subset(genotypehen2 == 1 | genotyperrp4 == 1) %>%
      as.data.frame() %>%
      select(genotypehen2, genotyperrp4, logFC_hen2, logFC_rrp4, txType_TAIR10) %>%
      rownames_to_column('TC_id') %>%
      as_tibble()
    # add if up in hen2, rrp4 or both
    logFC_xyplot_df %<>% mutate('DE_in'= case_when(genotypehen2 == 1 & genotyperrp4 != 1 ~ 'hen2_only',
                                                 genotypehen2 != 1 & genotyperrp4 == 1 ~ 'rrp4_only',
                                                 genotypehen2 == 1 & genotyperrp4 == 1 ~ 'both',
                                                 TRUE ~ 'error'))
    # sanity check
    table(logFC_xyplot_df$DE_in, useNA='always') # ok
    table(is.na(logFC_xyplot_df$logFC_hen2)) # no missing FC
    table(is.na(logFC_xyplot_df$logFC_rrp4)) # no missing FC
    
    # logFC XY-plot (all)
    ggplot(logFC_xyplot_df, aes(x=logFC_hen2, y=logFC_rrp4, col=DE_in)) +
      geom_point(alpha=.5) + xlim(-1, NA) + ylim(-1, NA) +
      stat_ellipse(type='norm', level=0.9) +
      geom_vline(xintercept=1, lty=2) +
      geom_hline(yintercept=1, lty=2) +
      geom_abline(intercept=0, slope=1, col='black') +
      cowplot::theme_cowplot() + theme(aspect.ratio=1) +
      scale_color_brewer(palette='Dark2', name='TC up in', direction=-1) +
      labs(title='up-regulated CAGE TCs', 
           subtitle='adjusted P-value <= 0.05',
           x='logFC (hen2-4 / WT)', y='logFC (rrp4-2 / WT )')
    
    # logFC XY-plot (mRNA only: TCs annotated as promoter)
    logFC_xyplot_df %>%
      filter(txType_TAIR10=='promoter') %>%
      ggplot(aes(x=logFC_hen2, y=logFC_rrp4, col=DE_in)) +
        geom_point(alpha=.5) + xlim(-1, NA) + ylim(-1, NA) +
        stat_ellipse(type='norm', level=0.9) +
        geom_vline(xintercept=1, lty=2) +
        geom_hline(yintercept=1, lty=2) +
        geom_abline(intercept=0, slope=1, col='black') +
        cowplot::theme_cowplot() + theme(aspect.ratio=1) +
        scale_color_brewer(palette='Dark2', name='TC up in', direction=-1) +
        labs(title='up-regulated CAGE TCs (mRNAs)', 
             subtitle='adjusted P-value <= 0.05',
             x='logFC (hen2-4 / WT)', y='logFC (rrp4-2 / WT )',
             caption='N(rrp4 only)=367\nN(hen2 only)=19\nN(both)=112')
        # corresponding Venn diagram
        logFC_xyplot_df %>%
          filter(txType_TAIR10=='promoter') %>%
          select(TC_id:genotyperrp4) %>%
          column_to_rownames('TC_id') %>%
          eulerr::euler() %>%
          plot()
        
    
    
    # output hen2-specific up-regulated TCs
    subset(logFC_xyplot_df, DE_in == 'hen2_only') %>%
      select(-genotypehen2, -genotyperrp4, -DE_in) %>%
      left_join(select(tt_genotypehen2, TC_id, adj.P.Val), by='TC_id') %>%
      left_join(select(as.data.frame(rowData(TCs)), thick.names, geneID, txType_TAIR10, symbol, geneID_anti, txType_TAIR10extended, symbol_anti), by=c('TC_id'='thick.names')) %>%
      left_join(select(idmapping, geneID, description), by='geneID') %>%
      left_join(select(idmapping, geneID, 'description_anti'=description), by=c('geneID_anti'='geneID')) %>%
      WriteXLS(ExcelFileName='~/Dropbox/Flg22_CAGE/Exosome_TSS paper/Data/hen2_specific_upregulated_CAGE_TCs.xls', verbose=T, row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)
      
    # output all PROMPTs for GRO-seq support (in another script)
    saveRDS(TCs, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')


# 3d. output bed file for any DE TC in exosome mutant
dt_exosome <- subset(as.data.frame(dt), genotypehen2 != 0 | genotyperrp4 != 0) %>%
  rownames_to_column('TCid') %>%
  as_tibble()
    # add GeneID and other info
    dt_exosome <- left_join(dt_exosome,
                            mcols(TCs) %>% as_tibble() %>% select(thick.names, txType_TAIR10, txType_TAIR10extended, geneID, geneID_anti, symbol, symbol_anti),
                            by=c('TCid'='thick.names'))
    # make BED and output
    dt_exosome.bed <- rowRanges(TCs)[names(rowRanges(TCs)) %in% dt_exosome$TCid]
    export.bed(dt_exosome.bed, 'TCs_DE_in_mutant.bed.txt')

    # how many promoter-DHSs have a TC up-regulated on the reverse strand?
        # get all up-regulated TCs annotated as "reverse" (PROMPTs)
        rrp4_up_prompts <- subset(dt_exosome, genotyperrp4 == 1 & txType_TAIR10=='reverse')
        hen2_up_prompts <- subset(dt_exosome, genotypehen2 == 1 & txType_TAIR10=='reverse')
        # see how many are overlapping with a DHS region
        rowRanges(TCs) %>%
          subset(names %in% unique(c(rrp4_up_prompts$TCid, hen2_up_prompts$TCid))) %>%
          subsetByOverlaps(promoters(dhss_promoter, upstream=400, downstream=400), ignore.strand=TRUE) # 91
        # NO GO ENRICHMENT FOR RRP4 PROMPT genes

    # output all genes (BED) having an antisense TSS up-regulated
        # 549 antisense TCs are up-regulated, they make up 481 genes
        antisense_TC_up_genes <- subset(dt_exosome, genotyperrp4 == 1 | genotypehen2 == 1) %>%
          subset(txType_TAIR10=='antisense') %>%
          .$geneID_anti %>%
          unique()
        # GO enrichment
        universe <- rowRanges(TCs)$geneID %>% na.omit() %>% as.vector()
        genes_antisense_TC_up_GSEA <- gprofiler(query=antisense_TC_up_genes, organism='athaliana', significant=T, correction_method='fdr', custom_bg=universe)
        
        genes_antisense_TC_up_GSEA_top20 <- data.frame('term'=genes_antisense_TC_up_GSEA$term.name,
                                                      'pval'=genes_antisense_TC_up_GSEA$p.value,
                                                      'term.id'=genes_antisense_TC_up_GSEA$term.id) %>%
                                                as_tibble() %>%
                                                mutate('log10pval'=-log(pval, 10)) %>%
                                                arrange(desc(log10pval)) %>%
                                                mutate('term'=factor(term, levels=rev(term)))
        
        genes_antisense_TC_up_GSEA_top20
        # export as BED for deeptools metaplot (RNAseq signal)
        genes <- genes(TxDb.Athaliana.BioMart.plantsmart28)
              seqlevelsStyle(genes) <- seqlevelsStyle(myseqinfo)
              seqlevels(genes) <- seqlevels(myseqinfo)
              seqinfo(genes) <- myseqinfo
              
        export.bed(subset(genes, gene_id %in% antisense_TC_up_genes), '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE/genes_antisense_TC_up.bed.txt')


# 3e. output bed file of genes with a 3'UTR antisense TSS in rrp4 (and/or hen2)
gene_with_up_aTSS_rrp4_3UTR <- subset(dt_exosome, genotyperrp4 == 1 & txType_TAIR10extended=='antisense_threeUTR')$geneID_anti
all_genes_txdb <- genes(TxDb.Athaliana.BioMart.plantsmart28)
  seqlevelsStyle(all_genes_txdb) <- seqlevelsStyle(myseqinfo)
  seqinfo(all_genes_txdb) <- myseqinfo
  subset(all_genes_txdb, gene_id %in% gene_with_up_aTSS_rrp4_3UTR) %>%
    export.bed('genes_with_3UTR_aTSS_up_in_CAGE_rrp4.bed.txt')



# 4. ANNOTATION OF DE TCs ####
# ----------------------------
# 4a. txTypes by DE categories by coefficient
dtByAnot <- select(dt_exosome, TCid:txType_TAIR10) %>%
            melt(id.vars=c('TCid', 'txType_TAIR10'), variable.name='coef') %>%
            filter(value != 0) %>%
            group_by(txType_TAIR10, coef, value) %>%
            summarise('nTCs'=n()) %>%
            ungroup() %>%
            rename('direction'='value')
map <- c('-1'='down', '1'='up')
dtByAnot$direction <- map[as.character(dtByAnot$direction)]
table(dtByAnot$txType_TAIR10)
dtByAnot$txType_TAIR10 %<>% factor(levels=c('intergenic', 'antisense', 'reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))
  
exo_colors2 <- brewer.pal(name='Set3', n=4)[3:4]
dtByAnot$coef %<>% factor(levels=c('genotypehen2', 'genotyperrp4'))
saveRDS(dtByAnot, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/dtByAnot.rds')

ggplot() +
       geom_bar(data=subset(dtByAnot, direction=='up'), aes(x=txType_TAIR10, y=nTCs, fill=coef), stat='identity', position=position_dodge(preserve='single'), col='black', lwd=.2) +
       geom_bar(data=subset(dtByAnot, direction=='down'), aes(x=txType_TAIR10, y=-nTCs, fill=coef), stat='identity', position=position_dodge(preserve='single'), col='black', lwd=.2) +
       geom_hline(yintercept=0) + cowplot::theme_cowplot() + coord_flip() + theme(aspect.ratio=1, legend.position='bottom', legend.direction=1) +
       scale_fill_manual(values=exo_colors2, name='Coefficient') +
       labs(title='Exosome: sense DE TSSs', y='down / up DE TSSs', x='')

    # check some of the promoter down-regulated in genotyperrp4, are they visibly ok?
    subset(dt_exosome, genotyperrp4 == -1 & txType_TAIR10 == 'promoter') %>% head(20)
    # output the gene list for GSEA
    subset(dt_exosome, genotyperrp4 == -1 & txType_TAIR10 == 'promoter', select=geneID) %>%
      write.table('genotyperrp4_all_genes_with_down_promoter_TC.txt', row.names=F, col.names=F, append=F, quote=F)
    subset(dt_exosome, genotyperrp4 == 1 & txType_TAIR10 == 'promoter', select=geneID) %>%
      write.table('genotyperrp4_all_genes_with_up_promoter_TC.txt', row.names=F, col.names=F, append=F, quote=F)
    # also output any gene with a detected promoter TC for the GSEA universe
    subset(rowRanges(TCs), txType_TAIR10 == 'promoter')$geneID %>%
      write.table('GOuniverse_all_detected_genes_with_a_promoter_TC.txt', row.names=F, col.names=F, append=F, quote=F)
    
    # GO on rrp4 DE TSSs promoters genes
    universe <- rowRanges(TCs)$geneID %>% na.omit() %>% as.vector()
    
    genes_prom_down_rrp4 <- subset(rowRanges(TCs), genotyperrp4==-1 & txType_TAIR10=='promoter')$geneID
    genes_prom_up_rrp4 <- subset(rowRanges(TCs), genotyperrp4==1 & txType_TAIR10=='promoter')$geneID
    
    genes_prom_down_rrp4_GSEA <- gprofiler(query=genes_prom_down_rrp4, organism='athaliana', significant=T, correction_method='fdr', custom_bg=universe)
    genes_prom_up_rrp4_GSEA <- gprofiler(query=genes_prom_up_rrp4, organism='athaliana', significant=T, correction_method='fdr', custom_bg=universe)
    
    genes_prom_down_rrp4_GSEA_top20 <- data.frame('term'=genes_prom_down_rrp4_GSEA$term.name,
                                                  'pval'=genes_prom_down_rrp4_GSEA$p.value,
                                                  'term.id'=genes_prom_down_rrp4_GSEA$term.id) %>%
                                       as_tibble() %>%
                                       mutate('log10pval'=-log(pval, 10)) %>%
                                       arrange(desc(log10pval)) %>%
                                       mutate('term'=factor(term, levels=rev(term))) %>%
                                       head(20)
    
    genes_prom_up_rrp4_GSEA_top20 <- data.frame('term'=genes_prom_up_rrp4_GSEA$term.name,
                                                'pval'=genes_prom_up_rrp4_GSEA$p.value,
                                                'term.id'=genes_prom_up_rrp4_GSEA$term.id) %>%
                                      as_tibble() %>%
                                      mutate('log10pval'=-log(pval, 10)) %>%
                                      arrange(desc(log10pval)) %>%
                                      mutate('term'=factor(term, levels=rev(term))) %>%
                                      head(20)
    
    ggplot(genes_prom_down_rrp4_GSEA_top20, aes(x=term, y=log10pval)) +
           geom_bar(stat='identity') +
           geom_text(aes(label=term, y=.1), col='white', hjust=0, size=4) +
           coord_flip() + cowplot::theme_cowplot() + theme(axis.text.y=element_blank(), aspect.ratio=1.5) +
           labs(title='GSEA: rrp4 genes with promoter TC down (Top20 terms)',
                subtitle=paste0('N=', nrow(genes_prom_down_rrp4_GSEA), ' significant terms'),
                caption=paste0('rrp4 genes with TC down (N down=', length(genes_prom_down_rrp4),' genes)'),
                y='-log10(adjusted p.value)')
    
    ggplot(genes_prom_up_rrp4_GSEA_top20, aes(x=term, y=log10pval)) +
           geom_bar(stat='identity') +
           geom_text(aes(label=term, y=.1), col='white', hjust=0, size=6) +
           coord_flip() + cowplot::theme_cowplot() + theme(axis.text.y=element_blank(), aspect.ratio=1.5) +
           labs(title='GSEA: rrp4 genes with promoter TC up (Top20 terms)',
                subtitle=paste0('N=', nrow(genes_prom_up_rrp4_GSEA), ' significant terms'),
                caption=paste0('rrp4 genes with TC up (N up=', length(genes_prom_up_rrp4),' genes)'),
                y='-log10(adjusted p.value)')

# 4c. Expression densities of DE TCs by their annotation, in the mutants only
    # compute TPM of each TC in each genotype
    meanTPM_df <- assay(TCs, 'TPM') %>%
      as.data.frame() %>%
      rownames_to_column('TC_id') %>%
      melt(id.vars='TC_id', variable.name='genotype', value.name='TPM') %>%
      separate(genotype, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
      select(-timepoint, -replicate) %>%
      group_by(TC_id, genotype) %>%
      summarise('meanTPM'=mean(TPM)) %>%
      ungroup() %>%
      spread(genotype, meanTPM) %>%
      set_colnames(c('TC_id', 'meanTPM_hen2', 'meanTPM_rrp4', 'meanTPM_wt')) %>%
      as_tibble()
    # arrange in same order and add to TCs SE
    meanTPM_df <- meanTPM_df[match(names(TCs), meanTPM_df$TC_id), ]
    rowRanges(TCs)$meanTPM_wt <- meanTPM_df$meanTPM_wt
    rowRanges(TCs)$meanTPM_hen2 <- meanTPM_df$meanTPM_hen2
    rowRanges(TCs)$meanTPM_rrp4 <- meanTPM_df$meanTPM_rrp4
    
    # expression boxplots for sense DE TCs
    rbind(rowRanges(TCs) %>%
            as_tibble() %>%
            subset(genotypehen2 != 0) %>%
            select(thick.names, txType_TAIR10, genotypehen2, meanTPM_hen2) %>%
            mutate(genotypehen2 = ifelse(genotypehen2=='-1', 'down-regulated', 'up-regulated'), set='hen2') %>%
            rename('direction'='genotypehen2', 'meanTPM'='meanTPM_hen2'),
          rowRanges(TCs) %>%
            as_tibble() %>%
            subset(genotyperrp4 != 0) %>%
            select(thick.names, txType_TAIR10, genotyperrp4, meanTPM_rrp4) %>%
            mutate(genotyperrp4 = ifelse(genotyperrp4=='-1', 'down-regulated', 'up-regulated'), set='rrp4') %>%
            rename('direction'='genotyperrp4', 'meanTPM'='meanTPM_rrp4')) %>%
      mutate('txType_TAIR10'=factor(txType_TAIR10, levels=c('intergenic', 'reverse', 'antisense', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))) %>%
      ggplot(aes(x=txType_TAIR10, y=meanTPM, fill=set)) +
             geom_boxplot() + scale_y_log10() + facet_wrap(~direction, scales='free') + theme(aspect.ratio=2) + coord_flip() +
             labs(title='Expression of sense DE TCs', y='Expression in respective mutant (logTPM)') +
             scale_fill_manual(values=exo_colors2)
    
    # expression boxplots for antisense DE TCs
    rbind(rowRanges(TCs) %>%
            as_tibble() %>%
            subset(genotypehen2 != 0 & txType_TAIR10 == 'antisense') %>%
            select(thick.names, txType_TAIR10extended, genotypehen2, meanTPM_hen2) %>%
            mutate(genotypehen2 = ifelse(genotypehen2=='-1', 'down-regulated', 'up-regulated'), set='hen2') %>%
            rename('direction'='genotypehen2', 'meanTPM'='meanTPM_hen2'),
          rowRanges(TCs) %>%
            as_tibble() %>%
            subset(genotyperrp4 != 0 & txType_TAIR10 == 'antisense') %>%
            select(thick.names, txType_TAIR10extended, genotyperrp4, meanTPM_rrp4) %>%
            mutate(genotyperrp4 = ifelse(genotyperrp4=='-1', 'down-regulated', 'up-regulated'), set='rrp4') %>%
            rename('direction'='genotyperrp4', 'meanTPM'='meanTPM_rrp4')) %>%
      ggplot(aes(x=txType_TAIR10extended, y=meanTPM, fill=set)) +
             geom_boxplot() + scale_y_log10() + facet_wrap(~direction, scales='free') + theme(aspect.ratio=1.3) + coord_flip() +
             labs(title='Expression of antisense DE TCs', y='Expression in respective mutant (logTPM)') +
             scale_fill_manual(values=exo_colors2)
    
    # save TCs for seqlogo analysis only
    saveRDS(TCs, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_for_seqlogos.rds')




# 4d. 3'UTR LENGTH AND DOWNSTREAM-INTERGENIC LENGTH AT 3'UTR-aTSSs
# ----------------------------------------------------------------
  # there are 354 antisense 3'UTR DE TSSs
  sum(rowRanges(TCs)$txType_TAIR10extended == 'antisense_threeUTR')
  # get those as a GR
  threeUTR_aTSSs_gr <- subset(rowRanges(TCs), txType_TAIR10extended == 'antisense_threeUTR')
  width(threeUTR_aTSSs_gr) %>% summary()
  # who are those genes with 3'UTR aTSS?
  threeUTR_aTSS_genes_df <- as.data.frame(threeUTR_aTSSs_gr) %>% as_tibble() %>% select(txType_TAIR10:meanTPM_rrp4)
  # add their description
  idmapping <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/AGI_ID_mapping.rds')
  threeUTR_aTSS_genes_df %<>% left_join(select(idmapping, -name), by=c('geneID_anti'='geneID'))
  # clean
  threeUTR_aTSS_genes_df %<>% select(-Entropy, -symbol, -geneID)
  threeUTR_aTSS_genes_df$genotyperrp4 <- ifelse(threeUTR_aTSS_genes_df$genotyperrp4==1, 'up', ifelse(threeUTR_aTSS_genes_df$genotyperrp4==-1, 'down', 'not DE'))
  threeUTR_aTSS_genes_df$genotypehen2 <- ifelse(threeUTR_aTSS_genes_df$genotypehen2==1, 'up', ifelse(threeUTR_aTSS_genes_df$genotypehen2==-1, 'down', 'not DE'))
  # add if known TF
  TFs <- read.csv('~/Dropbox/Axel_Arabidopsis_Flagellin/AGRIS/AtTFDB/families_data.txt', sep='\t', h=F, col.names=c('family', 'locus', 'name', 'description', 'bla', 'blaa', 'blaaa', 'blaaaa', 'No')) %>% select(family, locus) %>% as_tibble()
  TFs$locus %<>% toupper()
  threeUTR_aTSS_genes_df %<>% left_join(TFs, by=c('geneID_anti'='locus'))
  # add if DE up in flg22 treatment
  flg22_de <- rbind(readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/DE_geneLevel_timepoint10_df.rds') %>% as_tibble() %>% mutate('coeff'='timepoint10'),
                    readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/DE_geneLevel_timepoint30_df.rds') %>% as_tibble() %>% mutate('coeff'='timepoint30'),
                    readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/DE_geneLevel_wt3010_df.rds') %>% as_tibble() %>% mutate('coeff'='wt3010'))
  flg22_de_up <- subset(flg22_de, logFC > 0)
  flg22_de_genes <- unique(flg22_de_up$id)
  threeUTR_aTSS_genes_df$flg22_up <- ifelse(threeUTR_aTSS_genes_df$geneID_anti %in% flg22_de_genes, 'yes', 'no')
  # make venn
  threeUTR_aTSS_genes_df %>%
    select(genotyperrp4, genotypehen2, family, flg22_up) %>%
    mutate('genotypehen2'=ifelse(genotypehen2=='up', 1, 0)) %>%
    mutate('genotyperrp4'=ifelse(genotyperrp4=='up', 1, 0)) %>%
    mutate('family'=ifelse(is.na(as.character(family)), 0, 1)) %>%
    mutate('flg22_up'=ifelse(flg22_up=='yes', 1, 0)) %>%
    rename('is TF'='family') %>%
    eulerr::euler() %>% plot(quantities=T, main="Genes with 3'UTR aTSSs")
  # for the genes with 3'UTR aTSS that are up in flg22, are their 3'UTR aTSS down-regulated in treatment?
      # get genes with 3'UTR aTSS that are up in flg22
      tmp_genes <- subset(threeUTR_aTSS_genes_df, flg22_up=='yes')$geneID_anti %>% unique()
      # get genes with down-regulated TCs in flg22 treatement
      tmp_down_TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/DE_TC.rds') %>% subset(timepoint10 == -1 | timepoint30 == -1 | wt3010 == -1) %>% select(TCid, timepoint10, timepoint30, wt3010, txType)
      tmp_TCs_flg22 <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/SE_TCs_TPM1_min3lib.rds') %>% rowRanges() %>% as.data.frame() %>% as_tibble()
      tmp_down_TCs %<>% left_join(tmp_TCs_flg22, by=c('TCid'='thick.names'))
  
  # export to XLSX
  WriteXLS(threeUTR_aTSS_genes_df, ExcelFileName='~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE/3UTR aTSS genes.xlsx', row.names=F, col.names=T, BoldHeaderRow=T)
  # make sure none of them are on ChrC or ChrM
  table(seqnames(threeUTR_aTSSs_gr))
  # output 3'UTR aTSSs for IGV
  export.bed(threeUTR_aTSSs_gr, 'threeUTR_aTSSs.bed')
  # 197 of them are up-regulated in a mutant
  subset(threeUTR_aTSSs_gr, genotypehen2 == 1 | genotyperrp4 == 1)
  # make venn diagram
  data.frame('all'=names(subset(threeUTR_aTSSs_gr, genotypehen2 == 1 | genotyperrp4 == 1))) %>%
    mutate('hen2'=ifelse(all %in% names(subset(threeUTR_aTSSs_gr, genotypehen2 == 1)), 1, 0)) %>%
    mutate('rrp4'=ifelse(all %in% names(subset(threeUTR_aTSSs_gr, genotyperrp4 == 1)), 1, 0)) %>%
    column_to_rownames('all') %>%
    eulerr::euler() %>% plot(quantities=T, main="up-regulated 3'UTR aTSSs", labels=c('hen2-4', 'rrp4'))

  # 4d.1. get all intergenic regions as a GR
      intergenic_gr <- genes(TxDb.Athaliana.BioMart.plantsmart28) %>%
        GenomicRanges::reduce(ignore.strand=T) %>%
        gaps(start=NA, end=NA)
            # fix seqinfo
            seqlevelsStyle(intergenic_gr) <- seqlevelsStyle(myseqinfo)
            seqinfo(intergenic_gr) <- myseqinfo
            # remove ChrM and ChrC
            seqlevels(intergenic_gr, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
            # output for IGV
            export.bed(intergenic_gr, 'intergenic_gr.bed')
            # plot overall width
            data.frame(intergenic_gr) %>%
              ggplot(aes(x=width, col=seqnames)) + geom_density(lwd=1) +
                     scale_x_log10() + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
                     scale_color_brewer(palette='Paired', name='Chr') +
                     labs(x='width (log bp)', title='intergenic width')
            summary(width(unique(intergenic_gr)))
  # 4d.2. get 3'UTR regions with antisense TSS
      # all 3'UTR regions
      threeUTRs_gr <- threeUTRsByTranscript(TxDb.Athaliana.BioMart.plantsmart28) %>%
                      unlist() %>%
                      GenomicRanges::reduce()
                      # fix seqinfo
                      seqlevelsStyle(threeUTRs_gr) <- seqlevelsStyle(myseqinfo)
                      seqinfo(threeUTRs_gr) <- myseqinfo
                      # remove ChrC or ChrM
                      seqlevels(threeUTRs_gr, pruning.mode='coarse') <- setdiff(seqlevels(myseqinfo), c('ChrM', 'ChrC'))
                      # plot overall width
                      data.frame(threeUTRs_gr) %>%
                        ggplot(aes(x=width, col=seqnames)) + geom_density(lwd=1) +
                               cowplot::theme_cowplot() + theme(aspect.ratio=1) +
                               scale_color_brewer(palette='Paired', name='Chr') +
                               labs(x='width (bp)', title="3'UTR width", caption='limited to 1000 bp') + xlim(NA, 1000)
                      summary(width(threeUTRs_gr))
      # 3'UTR regions with a aTSS
            tmp_ranges <- IRanges(start=threeUTR_aTSSs_gr$thick.start, width=1)
            threeUTR_aTSSs_gr$thick.start <- IRanges(start=threeUTR_aTSSs_gr$thick.start, width=1)
      threeUTRs_with_aTSS_gr <- subsetByOverlaps(threeUTRs_gr, swapRanges(threeUTR_aTSSs_gr, inputColumn='thick.start'), ignore.strand=T)
      unique(threeUTRs_with_aTSS_gr) # 323
      width(unique(threeUTRs_with_aTSS_gr)) %>% summary()
      # 3'UTR regions with a DE-UP aTSS
      threeUTRs_with_aTSS_DE_up_gr <- subsetByOverlaps(threeUTRs_gr, swapRanges(subset(threeUTR_aTSSs_gr, genotypehen2 == 1 | genotyperrp4 == 1), inputColumn='thick.start'), ignore.strand=T)
      unique(threeUTRs_with_aTSS_DE_up_gr) # 181
      width(unique(threeUTRs_with_aTSS_DE_up_gr)) %>% summary()
          # output list of genes with 3UTR aTSS up for Peter
          threeUTR_aTSS_up_genes <- threeUTR_aTSSs_gr %>%
            subset(genotypehen2 == 1 | genotyperrp4 == 1) %$%
            unique(geneID_anti)
          
          subset(idmapping, geneID %in% threeUTR_aTSS_up_genes) %>%
            WriteXLS('~/Dropbox/Flg22_CAGE/Exosome_TSS paper/Data/gene_list_wit_3UTR_antisense_CAGE_TC_upregulated.xls', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)
            
      # 3'UTR regions without any aTSS
      threeUTRs_no_aTSS_gr <- subsetByOverlaps(threeUTRs_gr, threeUTRs_with_aTSS_gr, invert=T)
      unique(threeUTRs_no_aTSS_gr) # 22,366
      width(unique(threeUTRs_no_aTSS_gr)) %>% summary()
      # plot all widths
      ggplot() +
             geom_density(data=as_tibble(threeUTRs_with_aTSS_gr), aes(x=width, col='aTSS (N=323)'), lwd=1) +
             geom_density(data=as_tibble(threeUTRs_with_aTSS_DE_up_gr), aes(x=width, col='aTSS DE up (N=181)'), lwd=1) +
             geom_density(data=as_tibble(threeUTRs_no_aTSS_gr), aes(x=width, col='no aTSS (N=22,366)'), lwd=1) +
             scale_color_brewer(palette='Dark2', name="3'UTR regions with") + cowplot::theme_cowplot() + theme(aspect.ratio=1) +
             labs(title="3'UTRs lengths and antisense TSSs", x='width (bp)', caption='limited to 1200 bp') + xlim(NA, 1200)

  # 4d.3. use "follow()" to select the upstream intergenic regions to a 3'UTR-antisense TSS
      intergenic_gr # 31,068
      # intergenic regions following 3'UTR with antisense TSSs
      intergenic_with_aTSS_gr <- intergenic_gr[follow(threeUTR_aTSSs_gr, intergenic_gr)]
      unique(intergenic_with_aTSS_gr) # 322
      width(unique(intergenic_with_aTSS_gr)) %>% summary()
      # intergenic regions following 3'UTR with antisense TSSs DE UP
      intergenic_with_aTSS_DE_up_gr <- intergenic_gr[follow(subset(threeUTR_aTSSs_gr, genotypehen2 == 1 | genotyperrp4 == 1), intergenic_gr)]
      unique(intergenic_with_aTSS_DE_up_gr) # 181
      width(unique(intergenic_with_aTSS_DE_up_gr)) %>% summary()
      # all other intergenic regions
      intergenic_wo_aTSS_gr <- subsetByOverlaps(intergenic_gr, intergenic_with_aTSS_gr, invert=T)
      unique(intergenic_wo_aTSS_gr) # 30,746
      width(unique(intergenic_wo_aTSS_gr)) %>% summary()
      # plot all widths
      ggplot() +
        geom_density(data=as_tibble(intergenic_with_aTSS_gr), aes(x=width, col='aTSS (N=322)'), lwd=1) +
        geom_density(data=as_tibble(intergenic_with_aTSS_DE_up_gr), aes(x=width, col='aTSS DE up (N=181)'), lwd=1) +
        geom_density(data=as_tibble(intergenic_wo_aTSS_gr), aes(x=width, col='no aTSS (N=30,746)'), lwd=1) +
        scale_color_brewer(palette='Dark2', name="intergenic regions preceding") + cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position=c(.25, .7)) +
        labs(title="Lengths of intergenic regions\npreceding an antisense TSSs", x='width (bp)', caption='limited to 15000') + xlim(NA, 15000)

    
  # 4d.4 Testing statistically (log-transform for normality) or wilcox (non-transformed) + alternative=greater, non-paired)
      # intergenic -vs- intergenic aTSS
      t.test(log(width(intergenic_with_aTSS_gr)), log(width(intergenic_wo_aTSS_gr)), alternative='greater', paired=F)
      wilcox.test(width(intergenic_with_aTSS_gr), width(intergenic_wo_aTSS_gr), altenative='greater', paired=F)
      # intergenic -vs- intergenic aTSS DE UP
      t.test(log(width(intergenic_with_aTSS_DE_up_gr)), log(width(intergenic_wo_aTSS_gr)), alternative='greater', paired=F)
      wilcox.test(width(intergenic_with_aTSS_DE_up_gr), width(intergenic_wo_aTSS_gr), altenative='greater', paired=F)
      # 3'UTR -vs 3'UTR aTSS
      t.test(log(width(threeUTRs_with_aTSS_gr)), log(width(threeUTRs_no_aTSS_gr)), alternative='greater', paired=F)
      wilcox.test(width(threeUTRs_with_aTSS_gr), width(threeUTRs_no_aTSS_gr), altenative='greater', paired=F)
      # 3'UTR -vs 3'UTR aTSS DE UP
      t.test(log(width(threeUTRs_with_aTSS_DE_up_gr)), log(width(threeUTRs_no_aTSS_gr)), alternative='greater', paired=F)
      wilcox.test(width(threeUTRs_with_aTSS_DE_up_gr), width(threeUTRs_no_aTSS_gr), altenative='greater', paired=F)

  # 4d.5 any GO enrichment for the genes with a 3'UTR aTSS?
      genes_threeUTR_aTSS_GSEA <- gprofiler(query=unique(threeUTR_aTSSs_gr$geneID_anti), organism='athaliana', significant=T, correction_method='fdr', custom_bg=universe)
      tibble('term'=genes_threeUTR_aTSS_GSEA$term.name, 'pval'=genes_threeUTR_aTSS_GSEA$p.value, 'term.id'=genes_threeUTR_aTSS_GSEA$term.id)
      # gene list of TFs out for PETER
      genes_threeUTR_aTSS_GSEA_TFs <- genes_threeUTR_aTSS_GSEA$intersection %>% str_split(',', simplify=T) %>% as.vector()
      subset(idmapping, geneID %in% genes_threeUTR_aTSS_GSEA_TFs) %>% WriteXLS('~/Dropbox/Flg22_CAGE/Exosome_TSS paper/Data/gene_list_3UTR_antisense_upregulated_RNAseq_TFs.xls', row.names=F, col.names=T, BoldHeaderRow=T)
      

# 4e. HOW MANY 3'UTR-aTSSs HAVE A NORMAL TSS (ON THE OPPOSITE STRAND) IN A REASONNABLE DISTANCE
      # FORMING A TSS CONSTELLATION AKIN OF PROMPTs-mRNA-TSS?
      # 4e.1 get promoter TSSs
      promoter_TSSs_gr <- rowRanges(TCs) %>% subset(txType_TAIR10=='promoter')
      # 4e.2 get 3'UTR aTSSs
      threeUTR_aTSSs_gr
      # 4e.3 expend 3'UTR aTSSs 300bp upstream
      threeUTR_aTSSs_300upstream_gr <- promoters(swapRanges(threeUTR_aTSSs_gr), upstream=300, downstream=1)
      # 4e.4 flip strand of 3'UTR aTSSs
      threeUTR_aTSSs_1bp_inverted_gr <- threeUTR_aTSSs_gr %>%
        swapRanges() %>%
        invertStrand()
      # 4e.5 compute distance following 3'UTR aTSSs inverted
      next_promoter_TSS <- precede(threeUTR_aTSSs_1bp_inverted_gr, swapRanges(promoter_TSSs_gr))
      distances_df <- distance(threeUTR_aTSSs_1bp_inverted_gr, swapRanges(promoter_TSSs_gr[next_promoter_TSS])) %>%
        enframe(name=NULL, value='distance_bp')
      table(distances_df$distance_bp < 500)
      
      threeUTR_aTSSs_1bp_inverted_gr[2]
      swapRanges(promoter_TSSs_gr[next_promoter_TSS[2]])

      gg_a <- ggplot(distances_df, aes(x=distance_bp)) + geom_histogram() + labs(title="distance from 3'UTR aTSS to mRNA TSS", subtitle='on opposite strand')
      gg_b <- ggplot(distances_df, aes(x=distance_bp)) + geom_density()
      gg_a + gg_b + plot_layout(ncol=1)
      
      ggplot(distances_df, aes(x=distance_bp)) +
             geom_histogram(binwidth=2500) +
             facet_zoom(x=distance_bp < 1000, binwidth=100) +
             cowplot::theme_cowplot() +
             labs(title="Distance from 3'UTR aTSS to mRNA TSS")

      ggplot() + 
        geom_histogram(aes(x=distance_bp), mutate(distances_df, z=F), binwidth=1000, col='black', lwd=.3, fill='white') +
        geom_histogram(aes(x=distance_bp), mutate(distances_df, z=T), binwidth=50, col='black', lwd=.3, fill='white') +
        facet_zoom(xlim=c(0, 1025), ylim=c(0, 7.5), zoom.data=z, horizontal=F) +
        cowplot::theme_cowplot() + theme(zoom.y=element_blank(), validate=F) +
        labs(title="Distance from 3'UTR aTSS to mRNA TSS", x='Distance (bp)')

# 4f. WHERE DO THE aTSSs FALL WITHIN THE 3'UTR REGION
      threeUTR_aTSSs_gr # 354 aTSSs
      threeUTRs_gr # 22,689
      # match the 3'UTRs with the 3'UTR aTSSs so that they are in the same order
      overlap_hits <- findOverlaps(swapRanges(threeUTR_aTSSs_gr), threeUTRs_gr, ignore.strand=T)
      threeUTR_aTSSs_ordered_gr <- threeUTR_aTSSs_gr[queryHits(overlap_hits)] # 354
      threeUTRs_ordered_gr <- threeUTRs_gr[subjectHits(overlap_hits)] # 354
      # get width
      threeUTRs_widths <- width(threeUTRs_ordered_gr)
      # get distance from start of 3'UTR
      dist_from_threeUTR_start <- start(swapRanges(threeUTR_aTSSs_ordered_gr)) - start(threeUTRs_ordered_gr)
      # plot distance as percent
      dist_from_threeUTR_start <- enframe(dist_from_threeUTR_start / threeUTRs_widths, name=NULL, value='relative_pos')
      dist_from_threeUTR_start$isDE <- (threeUTR_aTSSs_ordered_gr$genotypehen2 == 1 | threeUTR_aTSSs_ordered_gr$genotyperrp4 ==1)
      # plot
      ggplot(dist_from_threeUTR_start, aes(x=relative_pos, col=isDE)) + geom_density() +
             labs(x='Relative position (%)', title="Position of aTSSs in 3'UTRs") +
      ggplot(dist_from_threeUTR_start, aes(x=relative_pos, fill=isDE)) + geom_histogram(binwidth=.1, position=position_dodge()) +
             labs(x='Relative position (%)') +
      plot_layout(ncol=1)
      # absolute distance to 3'UTR end
      absolute_dist_from_threeUTR_end <- (start(swapRanges(threeUTR_aTSSs_ordered_gr)) - end(threeUTRs_ordered_gr)) %>%
                                         enframe(name=NULL, value='dist_3end') %>%
                                         as_tibble()
      absolute_dist_from_threeUTR_end$isDE <- (threeUTR_aTSSs_ordered_gr$genotypehen2 == 1 | threeUTR_aTSSs_ordered_gr$genotyperrp4 ==1)
      
      ggplot(absolute_dist_from_threeUTR_end, aes(x=abs(dist_3end), col=isDE)) + geom_density() +
        labs(x="Distance from 3'end (abs bp)", title="aTSSs position in 3'UTRs") +
        cowplot::theme_cowplot() + theme(aspect.ratio=1) + scale_color_brewer(palette='Dark2', direction=-1) +
      ggplot() +
        geom_density(data=enframe(threeUTRs_widths, name=NULL, value='width'), aes(x=width, col="Yes")) +
        geom_density(data=enframe(width(threeUTRs_no_aTSS_gr), name=NULL, value='width'), aes(x=width, col="No")) + xlim(NA, 1000) +
        labs(x="3'UTR width (bp)", title="3'UTRs widths") + cowplot::theme_cowplot() + theme(aspect.ratio=1) + scale_color_manual(values=c('black', 'red'), name='has aTSS') +
        plot_layout(ncol=1)
      
      cbind(absolute_dist_from_threeUTR_end, threeUTRs_widths) %>%
        as_tibble() %>%
        ggplot(aes(x=threeUTRs_widths, y=abs(dist_3end), col=isDE)) +
               geom_point() + geom_abline(slope=1, intercept=0, lty=2) +
               cowplot::theme_cowplot() + theme(legend.position='bottom', legend.direction=2) +
               labs(x="3'UTR width (bp)", y="aTSS distance to 3'end (bp)") + coord_equal()
      
      # export dataset for 3'UTR aTSSs
      saveRDS(threeUTRs_with_aTSS_gr, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_with_aTSS_gr.rds')
      saveRDS(threeUTRs_no_aTSS_gr, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/threeUTRstory_3utrs_without_aTSS_gr.rds')
      

### Do the "promoter TSSs" that are down in rrp4 are actually overlapping the gene?
### to see if it could be why RNAseq DE doesn't capture them so much
      # get those TSSs
      prom_down_rrp4_gr <-  rowRanges(TCs) %>% subset(genotyperrp4 == -1 & txType_TAIR10 == 'promoter')
      # get gene bodies
      txdb <- TxDb.Athaliana.BioMart.plantsmart28
      genes <- genes(txdb)
          # fix seqinfo
          seqlevelsStyle(genes) <- seqlevelsStyle(cage_p$hen2)
          seqlevels(genes) <- seqlevels(cage_p$hen2)
          seqinfo(genes) <- seqinfo(cage_p$hen2)
      # look at overlap
      countOverlaps(swapRanges(prom_down_rrp4_gr), genes) %>% table() # yes=476, no=205



    # make PROMPT regions
    txdb <- TxDb.Athaliana.BioMart.plantsmart28
    genes <- genes(txdb)
          seqlevelsStyle(genes) <- seqlevelsStyle(cage_p$hen2)
          seqlevels(genes) <- seqlevels(cage_p$hen2)
          seqinfo(genes) <- seqinfo(cage_p$hen2)
    # remove OoB regions
    prompts <- promoters(genes, upstream=600, downstream=0)
    (OoB_prompts <- GenomicRanges:::get_out_of_bound_index(prompts))
    genes_noOoB <- genes[-OoB_prompts]
    prompts_noOoB <- prompts[-OoB_prompts]
    # remove PROMPTs overlapping another PROMPT
    prompts_noOoB_noSelfOverlap <- prompts_noOoB[countOverlaps(prompts_noOoB, ignore.strand=TRUE) == 1]
    genes_noOoB_noSelfOverlap <- genes_noOoB[countOverlaps(prompts_noOoB, ignore.strand=TRUE) == 1]
    # remove PROMPTs overlapping another gene
    prompts_noOoB_noOverlap <- prompts_noOoB_noSelfOverlap[countOverlaps(prompts_noOoB_noSelfOverlap, genes, ignore.strand=T) == 0]
    genes_noOoB_noOverlap <- genes_noOoB_noSelfOverlap[countOverlaps(prompts_noOoB_noSelfOverlap, genes, ignore.strand=T) == 0]
    # compute signals at those PROMPTs
    upstream=600 ; downstream=1
    
    test <- data.frame('geneID'    =genes_noOoB_noOverlap$gene_id,
                       'cage_wt'   =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=cage_p$wt,   reverse=cage_m$wt,   upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'cage_hen2' =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=cage_p$hen2, reverse=cage_m$hen2, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'cage_rrp4' =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=cage_p$rrp4, reverse=cage_m$rrp4, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'rna_wt'    =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=rnaseq_forward$WT,   reverse=rnaseq_reverse$WT, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'rna_lsm8'  =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=rnaseq_forward$LSM8,   reverse=rnaseq_reverse$LSM8, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'rna_dm'    =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=rnaseq_forward$DM,   reverse=rnaseq_reverse$DM, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'rna_rrp4'  =wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=rnaseq_forward$RRP4, reverse=rnaseq_reverse$RRP4, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'GROseq_het'=wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=grohet_p,    reverse=grohet_m, upstream=upstream, downstream=downstream)$anti %>% rowSums(),
                       'GROcap_dut'=wideMetaProfile(sites=resize(genes_noOoB_noOverlap, width=1, fix='start'), forward=grocap_p,    reverse=grocap_m, upstream=upstream, downstream=downstream)$anti %>% rowSums()) %>%
            as_tibble()

    # get PROMPTs: UP TCs annotated as "reverse"
    hen2_prompts <- rowRanges(TCs) %>% subset(genotypehen2 == 1) %>% subset(txType_TAIR10 == 'reverse') # 63
    rrp4_prompts <- rowRanges(TCs) %>% subset(genotyperrp4 ==1) %>% subset(txType_TAIR10 == 'reverse') # 94
    all_prompts  <- rowRanges(TCs) %>% subset(genotypehen2 == 1 | genotyperrp4 ==1) %>% subset(txType_TAIR10 == 'reverse') # 96
    # export all PROMPTs
    all_prompts_df <- all_prompts %>%
      as.data.frame() %>%
      as_tibble() %>%
      select(-score, -thick.end, -thick.width, -txType_ARAPORT11, -Entropy) %>%
      rename('peak'='thick.start', 'TC_id'='thick.names', 'CAGE_DE_hen2'='genotypehen2', 'CAGE_DE_rrp4'='genotyperrp4')
    View(all_prompts_df)
    
    export.bed(all_prompts, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE/CAGE_PROMPTs.bed')
    
    # how many genes with CAGE PROMPT are found in Chekanova 2007?: NONE
    chekanova_2007_unts <- readxl::read_xlsx('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/06 - Comparisons/6d. Chekanova UNTs/Chekanova 2007 UNTs in rrp4est.xlsx') %>% as_tibble()
    table(all_prompts_df$geneID_anti %in% chekanova_2007_unts$AGI)

    # HERE HERE HERE
    # HERE HERE HERE
    # HERE HERE HERE
    # HERE HERE HERE
    
    # add if gene has PROMPT
    test$PROMPT <- ifelse(test$geneID %in% all_prompts$geneID_anti, TRUE, FALSE)
    
    gg_a <- ggplot(test, aes(x=GROseq_het+1, y=cage_rrp4+1, col=PROMPT)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() + scale_y_log10() + #geom_abline(intercept=0, slope=1, col='red') +
      labs(x='GRO-seq Hetzel', y='CAGE rrp4', title='Antisense signal at all PROMPT regions', subtitle='-600bp upstream of TSSs, not overlapping anything else') +
      scale_color_brewer(palette='Dark2', direction=-1) + theme(aspect.ratio=1)
    
    gg_b <- ggplot(subset(test, PROMPT==TRUE), aes(x=GROseq_het+1, y=cage_rrp4+1)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() + scale_y_log10() + #geom_abline(intercept=0, slope=1, col='red') +
      labs(x='GRO-seq Hetzel', y='CAGE rrp4', title='Antisense signal at validated PROMPTs', subtitle='-600bp upstream of TSSs, not overlapping anything else') +
      theme(aspect.ratio=1)
    
    test %>%
      melt(id.vars=c('geneID', 'PROMPT')) %>%
      ggplot(aes(x=value, col=PROMPT, fill=PROMPT)) + geom_density(alpha=.6) + facet_wrap(~variable, scales='free') + scale_x_log10() +
      cowplot::theme_cowplot() + theme(aspect.ratio=1, strip.text=element_text(face='bold')) +
      labs(title='Antisense signal at PROMPT regions', x='Normalized signal (log)')
    
    test %>%
      subset(PROMPT==T) %>%
      select(geneID:rna_rrp4) %>%
      melt(id.vars='geneID') %>%
      separate(variable, c('source', 'sample'), sep='_') %>%
      ggplot(aes(x=value + .1, col=sample)) +
             geom_density(alpha=.6) + facet_wrap(~source, scales='free') + scale_x_log10() +
             cowplot::theme_cowplot() + theme(aspect.ratio=1, strip.text=element_text(face='bold')) +
             labs(title='Antisense signal at PROMPT regions', x='Normalized signal (log)')
    
    ggarrange(gg_a, gg_b, ncol=1, nrow=2, align='hv')
    
    # make venn diagram: hen2 prompts / rrp4 prompts / DHS regions
    tmp_prompts_in_dhs <- names(subsetByOverlaps(all_prompts, promoters(dhss, upstream=400, downstream=400), ignore.strand=T))
    tmp_prompts_in_hen2 <- names(hen2_prompts)
    tmp_prompts_in_rrp4 <- names(rrp4_prompts)
    data.frame(row.names=names(all_prompts),
               'hen2_prompts' = names(all_prompts) %in% tmp_prompts_in_hen2,
               'rrp4_prompts' = names(all_prompts) %in% tmp_prompts_in_rrp4,
               'dhs_regions' = names(all_prompts) %in% tmp_prompts_in_dhs) %>%
      limma::vennDiagram()
    rm(tmp_prompts_in_dhs, tmp_prompts_in_hen2, tmp_prompts_in_rrp4)
    # get TSS of genes having PROMPTs (1bp)
    all_genes_TSSs <- genes %>% resize(width=1, fix='start') # 33602
    # remove genes that will have PROMPT out of bound
        oob <- suppressWarnings(promoters(all_genes_TSSs, upstream=401, downstream=2) %>% GenomicRanges:::get_out_of_bound_index())
        all_genes_TSSs <- all_genes_TSSs[-oob]
    # compute total signal on PROMPT regions (400 bp upstream of TSSs)
    df_anti <- data.frame('geneID'=all_genes_TSSs$gene_id,
                          'cage_wt' = wideMetaProfile(sites=all_genes_TSSs, forward=cage_p$wt, reverse=cage_m$wt, upstream=400, downstream=1)$anti %>% rowSums(),
                          'cage_rrp4' = wideMetaProfile(sites=all_genes_TSSs, forward=cage_p$rrp4, reverse=cage_m$rrp4, upstream=400, downstream=1)$anti %>% rowSums(),
                          'gro_hetzel'= wideMetaProfile(sites=all_genes_TSSs, forward=grohet_p, reverse=grohet_m, upstream=400, downstream=1)$anti %>% rowSums(),
                          'gro_duttke'= wideMetaProfile(sites=all_genes_TSSs, forward=grodut_p, reverse=grodut_m, upstream=400, downstream=1)$anti %>% rowSums(),
                          'grocap_duttke'= wideMetaProfile(sites=all_genes_TSSs, forward=grocap_p, reverse=grocap_m, upstream=400, downstream=1)$anti %>% rowSums()) %>%
               as_tibble()
    
    # add if gene has PROMPT or not
    df_anti$PROMPT <- ifelse(df_anti$geneID %in% all_prompts$geneID_anti, TRUE, FALSE)
    
    df_anti <- data.frame('geneID'=names(TCs),
                          'cage_wt' = wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=cage_p$wt, reverse=cage_m$wt, upstream=20, downstream=20)$anti %>% rowSums(),
                          'cage_rrp4' = wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=cage_p$rrp4, reverse=cage_m$rrp4, upstream=20, downstream=20)$anti %>% rowSums(),
                          'gro_hetzel'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grohet_p, reverse=grohet_m, upstream=20, downstream=20)$anti %>% rowSums(),
                          'gro_duttke'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grodut_p, reverse=grodut_m, upstream=20, downstream=20)$anti %>% rowSums(),
                          'grocap_duttke'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grocap_p, reverse=grocap_m, upstream=20, downstream=20)$anti %>% rowSums()) %>%
      as_tibble()
    df_sense <- data.frame('geneID'=names(TCs),
                           'cage_wt' = wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=cage_p$wt, reverse=cage_m$wt, upstream=20, downstream=20)$sense %>% rowSums(),
                           'cage_rrp4' = wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=cage_p$rrp4, reverse=cage_m$rrp4, upstream=20, downstream=20)$sense %>% rowSums(),
                           'gro_hetzel'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grohet_p, reverse=grohet_m, upstream=20, downstream=20)$sense %>% rowSums(),
                           'gro_duttke'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grodut_p, reverse=grodut_m, upstream=20, downstream=20)$sense %>% rowSums(),
                           'grocap_duttke'= wideMetaProfile(sites=swapRanges(rowRanges(TCs)), forward=grocap_p, reverse=grocap_m, upstream=20, downstream=20)$sense %>% rowSums()) %>%
      as_tibble()
    
    df_anti$PROMPT <- ifelse(df_anti$geneID %in% names(subset(rowRanges(TCs), txType_TAIR10=='reverse' & (genotypehen2 == 1 | genotyperrp4==1))), TRUE, FALSE)
    df_sense$PROMPT <- ifelse(df_anti$geneID %in% names(subset(rowRanges(TCs), txType_TAIR10=='reverse' & (genotypehen2 == 1 | genotyperrp4==1))), TRUE, FALSE)
    Rmisc::multiplot(ggplot(df_anti, aes(x=cage_wt, y=cage_rrp4, col=PROMPT)) + geom_point() + scale_x_log10() + scale_y_log10() + labs(title='Antisense'),
                     ggplot(df_sense, aes(x=cage_wt, y=cage_rrp4, col=PROMPT)) + geom_point() + scale_x_log10() + scale_y_log10() + labs(title='Sense'),
                     cols=1)

    Rmisc::multiplot(
      ggplot(subset(df_anti, PROMPT==TRUE), aes(x=cage_rrp4, y=gro_hetzel)) +
             geom_point() + scale_x_log10() + scale_y_log10() +
             scale_color_manual(values=c('black', 'red')) + theme(aspect.ratio=1),
      
      ggplot(subset(df_anti, PROMPT==TRUE), aes(x=cage_rrp4, y=gro_duttke)) +
             geom_point() + scale_x_log10() + scale_y_log10() +
             scale_color_manual(values=c('black', 'red')) + theme(aspect.ratio=1),
      
      ggplot(subset(df_anti, PROMPT==TRUE), aes(x=cage_rrp4, y=grocap_duttke)) +
             geom_point() + scale_x_log10() + scale_y_log10() +
             scale_color_manual(values=c('black', 'red')) + theme(aspect.ratio=1),
      cols=1)

    # correlation between WT TSSs and GRO-Seq at CAGE WT TSSs +/- 150bp
    TCs_ctrl <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_ctrl.rds')
    
    sites <- rowRanges(TCs_ctrl) %>% swapRanges()
    upstream=10 ; downstream=200
    grodf_sense <- data.frame('geneID' = names(sites),
                              'GROseq_duttke'=wideMetaProfile(sites=sites, forward=grodut_p,  reverse=grodut_m,  upstream=upstream, downstream=downstream)$sense %>% rowSums(),
                              'GROseq_hetzel'=wideMetaProfile(sites=sites, forward=grohet_p,  reverse=grohet_m,  upstream=upstream, downstream=downstream)$sense %>% rowSums(),
                              'GROcap_duttke'=wideMetaProfile(sites=sites, forward=grocap_p,  reverse=grocap_m,  upstream=upstream, downstream=downstream)$sense %>% rowSums(),
                              'CAGE_wt'      =wideMetaProfile(sites=sites, forward=cage_p$wt, reverse=cage_m$wt, upstream=upstream, downstream=downstream)$sense %>% rowSums()) %>%
                              as_tibble()
    
    grodf_anti <- data.frame('geneID' = names(sites),
                             'GROseq_duttke'=wideMetaProfile(sites=sites, forward=grodut_p,  reverse=grodut_m,  upstream=downstream, downstream=upstream)$anti %>% rowSums(),
                             'GROseq_hetzel'=wideMetaProfile(sites=sites, forward=grohet_p,  reverse=grohet_m,  upstream=downstream, downstream=upstream)$anti %>% rowSums(),
                             'GROcap_duttke'=wideMetaProfile(sites=sites, forward=grocap_p,  reverse=grocap_m,  upstream=downstream, downstream=upstream)$anti %>% rowSums(),
                             'CAGE_wt'      =wideMetaProfile(sites=sites, forward=cage_p$wt, reverse=cage_m$wt, upstream=downstream, downstream=upstream)$anti %>% rowSums()) %>%
                             as_tibble()
      
    grodf_sense %>%
      melt(id.vars='geneID') %>%
      ggplot(aes(x=value, col=variable)) +
             geom_density(lwd=1) + scale_x_log10() +
             theme(aspect.ratio=1) + labs(title='Signal density at CAGE WT TSSs -10/+200bp', subtitle='sense only', x='TPM (log)')
      
    grodf_anti %>%
      melt(id.vars='geneID') %>%
      ggplot(aes(x=value, col=variable)) +
             geom_density(lwd=1) + scale_x_log10() +
             theme(aspect.ratio=1) + labs(title='Signal density at CAGE WT TSSs -200/+10bp', subtitle='antisense only', x='TPM (log)')
    
    ggpairs(na.omit(log(select(grodf_sense, -geneID))), showStrips=T, xlab='Normalized signal (log)', ylab='Normalized signal (log)', title='Signal at CAGE WT TSSs -10/+200bp (sense)')
    ggpairs(na.omit(log(select(grodf_anti, -geneID))),  showStrips=T, xlab='Normalized signal (log)', ylab='Normalized signal (log)', title='Signal at CAGE WT TSSs -200/+10bp (antisense)')

# 4f. log-log plot DE TSSs
      # get TPM matrix and average replicates
      TPMaveraged <- assay(TCs, 'TPM') %>%
        as.data.frame() %>%
        rownames_to_column('id') %>%
        melt(id.vars='id', value.name='TPM') %>%
        separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
        select(-timepoint) %>%
        group_by(id, genotype) %>%
        summarise('meanTPM'=mean(TPM)) %>%
        dcast(id ~ genotype, value.var='meanTPM') %>%
        as_tibble()
      
      # add txType
      TPMaveraged <- left_join(TPMaveraged, select(as.data.frame(rowData(TCs)), thick.names, txType_TAIR10), by=c('id'='thick.names'))

      # log-log hen2 / wt
      tt_genotypehen2 <- topTable(fit=eb, coef='genotypehen2', number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>% rownames_to_column('id') %>% as_tibble()
      tt_genotypehen2 %<>% left_join(as.data.frame(rowRanges(TCs)) %>% select(thick.names, txType_TAIR10, geneID, symbol), by=c('id'='thick.names'))
      saveRDS(tt_genotypehen2, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotypehen2.rds')
      
      volc_hen2 <- topTable(fit=eb, coef='genotypehen2', number=Inf, adjust.method='BH') %>%
                   rownames_to_column('id') %>%
                   ggplot(aes(x=logFC, y=-log10(adj.P.Val), col=(adj.P.Val < 0.05 & abs(logFC) > 1))) +
                          geom_point() + geom_vline(xintercept = c(-1, 1), lty=2) +
                          cowplot::theme_cowplot() + theme(aspect.ratio=1) + labs(title='Volcano: hen2-4') +
                          scale_color_brewer(palette='Set1', name='ajd.P.val < 0.5\nabs(logFC) > 1', direction=-1)
      
      gg_scatter_hen2 <- subset(TPMaveraged, id %in% tt_genotypehen2$id) %>%
        ggplot(aes(x=wt, y=hen2, col=txType_TAIR10)) +
               geom_point(size=.7) +
               geom_abline(slope=1, intercept=0, col='black', lty=2) +
               facet_wrap(~txType_TAIR10, nrow=1) +
               cowplot::theme_cowplot() +
               theme(aspect.ratio=1, legend.position='none') +
               scale_x_log10(labels=scales::comma_format()) + scale_y_log10(labels=scales::comma_format()) +
               labs(x="WT TPM (log)", y="hen2-4 TPM (log)", title="DE TSSs: genotype hen2")
      
      # log-log rrp4 / wt
      tt_genotyperrp4 <- topTable(fit=eb, coef='genotyperrp4', number=Inf, adjust.method='BH', p.value=0.05, lfc=1) %>% rownames_to_column('id') %>% as_tibble()
      tt_genotyperrp4 %<>% left_join(as.data.frame(rowRanges(TCs)) %>% select(thick.names, txType_TAIR10, geneID, symbol), by=c('id'='thick.names'))
      saveRDS(tt_genotyperrp4, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotyperrp4.rds')
      
      volc_rrp4 <- topTable(fit=eb, coef='genotyperrp4', number=Inf, adjust.method='BH') %>%
                   rownames_to_column('id') %>%
                   ggplot(aes(x=logFC, y=-log10(adj.P.Val), col=(adj.P.Val < 0.05 & abs(logFC) > 1))) +
                          geom_point() + geom_vline(xintercept = c(-1, 1), lty=2) +
                          cowplot::theme_cowplot() + theme(aspect.ratio=1) + labs(title='Volcano: rrp4-1') +
                          scale_color_brewer(palette='Set1', name='ajd.P.val < 0.5\nabs(logFC) > 1', direction=-1)
      
      gg_scatter_rrp4 <- subset(TPMaveraged, id %in% tt_genotyperrp4$id) %>%
        ggplot(aes(x=wt, y=rrp4, col=txType_TAIR10)) +
               geom_point(size=.7) +
               geom_abline(slope=1, intercept=0, col='black', lty=2) +
               facet_wrap(~txType_TAIR10, nrow=1) +
               cowplot::theme_cowplot() +
               theme(aspect.ratio=1, legend.position='none') +
               scale_x_log10(labels=scales::comma_format()) + scale_y_log10(labels=scales::comma_format()) +
               labs(x="WT TPM (log)", y="rrp4-1 TPM (log)", title="DE TSSs: genotype rrp4")
      
      # combined log-log plots
      ggarrange(gg_scatter_hen2, gg_scatter_rrp4, nrow=2, ncol=1)
      # combined volcano-plots
      ggarrange(volc_hen2, volc_rrp4, ncol=1, nrow=2, common.legend=T, legend='bottom')
      
      # FoldChange by TxType for SENSE DE
      FC_hen2_df <- rowRanges(TCs) %>%
        subset(genotypehen2 != 0) %>%
        as_tibble() %>%
        left_join(select(tt_genotypehen2, id, logFC) %>% rename('logFC_hen2'='logFC'), by=c('thick.names'='id')) %>%
        mutate('txType_TAIR10'=factor(txType_TAIR10, levels=c('intergenic', 'reverse', 'antisense', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR')))
        
      FC_rrp4_df <- rowRanges(TCs) %>%
        subset(genotyperrp4 != 0) %>%
        as_tibble() %>%
        left_join(select(tt_genotyperrp4, id, logFC) %>% rename('logFC_rrp4'='logFC'), by=c('thick.names'='id')) %>%
        mutate('txType_TAIR10'=factor(txType_TAIR10, levels=c('intergenic', 'reverse', 'antisense', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR')))
      
      ggplot() +
             geom_boxplot(data=FC_hen2_df, aes(x=txType_TAIR10, y=logFC_hen2, group=interaction(genotypehen2, txType_TAIR10), fill='hen2'), position=position_nudge(x=-.2), width=.4) +
             geom_boxplot(data=FC_rrp4_df, aes(x=txType_TAIR10, y=logFC_rrp4, group=interaction(genotyperrp4, txType_TAIR10), fill='rrp4'), position=position_nudge(x=.2), width=.4) +
             geom_hline(yintercept=c(-1, 1), lty=2) + coord_flip() +
             labs(title='logFC of DE TCs (sense)', y='logFC') + theme(aspect.ratio=1.3) +
             scale_fill_manual(values=exo_colors2)
      
      FC_anti_hen2_df <- rowRanges(TCs) %>%
        subset(genotypehen2 != 0 & txType_TAIR10 == 'antisense') %>%
        as_tibble() %>%
        left_join(select(tt_genotypehen2, id, logFC) %>% rename('logFC_hen2'='logFC'), by=c('thick.names'='id'))
      
      FC_anti_rrp4_df <- rowRanges(TCs) %>%
        subset(genotyperrp4 != 0 & txType_TAIR10 == 'antisense') %>%
        as_tibble() %>%
        left_join(select(tt_genotyperrp4, id, logFC) %>% rename('logFC_rrp4'='logFC'), by=c('thick.names'='id'))

      ggplot() +
        geom_boxplot(data=FC_anti_hen2_df, aes(x=txType_TAIR10extended, y=logFC_hen2, group=interaction(genotypehen2, txType_TAIR10extended), fill='hen2'), position=position_nudge(x=-.2), width=.4) +
        geom_boxplot(data=FC_anti_rrp4_df, aes(x=txType_TAIR10extended, y=logFC_rrp4, group=interaction(genotyperrp4, txType_TAIR10extended), fill='rrp4'), position=position_nudge(x=.2), width=.4) +
        geom_hline(yintercept=c(-1, 1), lty=2) + coord_flip() +
        labs(title='logFC of DE TCs (antisense)', y='logFC') + theme(aspect.ratio=1) +
        scale_fill_manual(values=exo_colors2)

# 4g. antisense DE TSSs annotation
      # annotation categories of antisense DE UP in hen2 and rrp4
            # make DF with TCs id and both txTypes
            TSSs_annotated_df <- TCs %>%
                           rowRanges() %>%
                           as.data.frame() %>%
                           rownames_to_column('id') %>%
                           select(id, txType_TAIR10extended, geneID_anti) %>% as_tibble()
            # add both annotation to tt's and also which coeff
            tt_genotypehen2 <- left_join(tt_genotypehen2, TSSs_annotated_df, by='id') %>% mutate('coef'='genotypehen2')
            tt_genotyperrp4 <- left_join(tt_genotyperrp4, TSSs_annotated_df, by='id') %>% mutate('coef'='genotyperrp4')
            
            # put all tts together
            tt_all <- rbind(tt_genotypehen2, tt_genotyperrp4)
            # add coefficient information and effect
            tt_all$coef %<>% factor(levels=c('genotypehen2', 'genotyperrp4'))
            # add DE direction
            tt_all$direction <- ifelse(tt_all$logFC > 0, 'up', 'down')
            # plot txTypeExtended of antisense TSSs that are up-regulated
            dtByAnot_antisense <- tt_all %>%
              subset(txType_TAIR10 == 'antisense') %>%
              select(txType_TAIR10extended, coef, direction) %>%
              group_by_all() %>%
              summarise('count'=n())
            
            saveRDS(dtByAnot_antisense, '~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/dtByAnot_antisense.rds')
            
            ggplot(dtByAnot_antisense, aes(x=txType_TAIR10extended, y=ifelse(direction=='down', -count, count), fill=coef)) +
                   geom_bar(stat='identity', position=position_dodge(), col='black', lwd=.2) +
                   geom_label(aes(label=count, group=coef), position=position_dodge(width=.9), fill='white', fontface='bold') +
                   geom_hline(yintercept=0, lty=2) +
                   cowplot::theme_cowplot() + theme(aspect.ratio=1) + coord_flip() +
                   scale_fill_manual(values=exo_colors2, name='Coefficient') +
                   labs(title='Extended annotation of antisense TSSs',
                        x='', y='down / up DE TSSs')
            
      # are the up-regulated 'reverse' TSSs shared between hen2 and rrp4?
      reverse_TSSs_up_in_hen2 <- subset(tt_genotypehen2, logFC > 0 & txType_TAIR10=='reverse')$id
      reverse_TSSs_up_in_rrp4 <- subset(tt_genotyperrp4, logFC > 0 & txType_TAIR10=='reverse')$id
      dev.off() ; VennDiagram::draw.pairwise.venn(area1=length(reverse_TSSs_up_in_hen2),
                                                  area2=length(reverse_TSSs_up_in_rrp4),
                                                  cross.area=length(intersect(reverse_TSSs_up_in_hen2, reverse_TSSs_up_in_rrp4)),
                                                  category=c('genotype hen2', 'genotype rrp4'), fill=c('#9698B7', '#303467'), alpha=.9, cex=3)
      
      # are the up-regulated 'antisense' TSSs shared between hen2 and rrp4 ?
      antisense_TSSs_up_in_hen2 <- subset(tt_genotypehen2, logFC > 0 & txType_TAIR10=='antisense')$id
      antisense_TSSs_up_in_rrp4 <- subset(tt_genotyperrp4, logFC > 0 & txType_TAIR10=='antisense')$id
      dev.off() ; VennDiagram::draw.pairwise.venn(area1=length(antisense_TSSs_up_in_hen2),
                                                  area2=length(antisense_TSSs_up_in_rrp4),
                                                  cross.area=length(intersect(antisense_TSSs_up_in_hen2, antisense_TSSs_up_in_rrp4)),
                                                  category=c('genotype hen2', 'genotype rrp4'), fill=c('#9698B7', '#303467'), alpha=.9, cex=3)
      
      # are the up-regulated 'intergenic' TSSs shared between hen2 and rrp4 ?
      intergenic_TSSs_up_in_hen2 <- subset(tt_genotypehen2, logFC > 0 & txType_TAIR10=='intergenic')$id
      intergenic_TSSs_up_in_rrp4 <- subset(tt_genotyperrp4, logFC > 0 & txType_TAIR10=='intergenic')$id
      dev.off() ; VennDiagram::draw.pairwise.venn(area1=length(intergenic_TSSs_up_in_hen2),
                                                  area2=length(intergenic_TSSs_up_in_rrp4),
                                                  cross.area=length(intersect(intergenic_TSSs_up_in_hen2, intergenic_TSSs_up_in_rrp4)),
                                                  category=c('genotype hen2', 'genotype rrp4'), fill=c('#9698B7', '#303467'), alpha=.9, cex=3)
      dev.off()

      # get all genes with any sort of antisense DE TC in an exosome mutant
      table(tt_all$txType_TAIR10extended) %>% as.tibble()
      exo_antisense_up <- subset(TCs,
                                 names(TCs) %in% unique(subset(tt_all, effect=='exosome' & direction=='up-regulated' & txType_TAIR10 %in% c('reverse', 'antisense'))$id)) %>%
                          rowRanges()
      
      as.tibble(mcols(exo_antisense_up)) %>%
        ggplot(aes(x=txType_TAIR10extended)) +
               geom_bar(stat='count') + coord_flip() +
               labs(title='DE TSSs: exosome up-regulated', subtitle='being reverse or antisense')
          
          # For T=0
          # CAGE footprint around those TCs (need to invert strand as the initial TC are antisense and we want to see them as antisense)
          cage_tmP_exo_anti_up <- tidyMetaProfile(sites=split(swapRanges(invertStrand(exo_antisense_up)), exo_antisense_up$txType_TAIR10extended, drop=T),
                                                  forward=cage_p[c(1,4,6)], reverse=cage_m[c(1,4,6)],
                                                  upstream=150, downstream=150,
                                                  trimLower=0.1, trimUpper=0.9)
          
          cage_tmP_exo_anti_up_count <- as.tibble(mcols(exo_antisense_up)) %>%
                          subset(txType_TAIR10extended != 'intergenic') %>%
                          group_by(txType_TAIR10extended) %>%
                          summarise('n'=n()) %>%
                          mutate('txType_TAIR10extended'=str_remove(txType_TAIR10extended, 'antisense_')) %>%
                          arrange(txType_TAIR10extended) %$%
                          rep(.$n, 4)
          
          cage_tmP_exo_anti_up %>%
            mutate('anti'=-anti, 'signal'=factor(str_remove(signal, '_0'), levels=c('rrp4', 'hen2', 'wt'))) %>%
            gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
            separate(sites, c('antisense', 'txType'), sep='_') %>%
            ggplot(aes(x=pos0*-1, y=score, col=direction)) + # * -1 to invert direction as I inverted strand before
                   geom_vline(xintercept=0, lwd=.2) + geom_line(lwd=.8) +
                   annotate(geom='text', x=80, y=0.4, label=cage_tmP_exo_anti_up_count) +
                   facet_grid(signal~txType) + cowplot::theme_cowplot() +
                   theme(legend.position='bottom', axis.text=element_text(size=8), strip.text=element_text(size=10), aspect.ratio=1) +
                   scale_color_brewer(palette='Set1', direction=-1) +
                   labs(title='up-regulated antisense TSSs in exosome effect',
                        x='position relative to antisense CAGE TSS (bp)', y='CAGE TPM (0.1-0.9%tile)')
          
          rnaseq_tmP_exo_anti_up <- tidyMetaProfile(sites=split(swapRanges(invertStrand(exo_antisense_up)), exo_antisense_up$txType_TAIR10extended, drop=T),
                                                  forward=rnaseq_forward, reverse=rnaseq_reverse,
                                                  upstream=500, downstream=100,
                                                  trimLower=0.1, trimUpper=0.9)
          
          rnaseq_tmP_exo_anti_up %>%
            mutate('anti'=-anti) %>%
            gather(key="direction", value="score", sense, anti, factor_key=TRUE) %>%
            separate(sites, c('antisense', 'txType'), sep='_') %>%
            ggplot(aes(x=pos0, y=score, fill=direction)) +
                   geom_vline(xintercept=0, col='grey60', lwd=.1) +
                   geom_area(alpha=.75) + annotate(geom='text', x=-300, y=3.5, label=cage_tmP_exo_anti_up_count) +
                   facet_grid(signal~txType) + cowplot::theme_cowplot() + theme(legend.position='bottom', axis.text=element_text(size=8), strip.text=element_text(size=10), aspect.ratio=1) +
                   scale_fill_brewer(palette='Set1', direction=-1) +
                   labs(title='up-regulated antisense TSSs in exosome effect (N=849)',
                        x='position relative to antisense CAGE TSS (bp)', y='RNAseq FPM (0.1-0.9%tile)')

        # get PROMPTs
        subset(exo_antisense_up, txType_TAIR10extended=='antisense_proximal')$geneID_anti













# 5. ALTERNATIVE TSS ANALYSIS DURING TIMECOURSE #### USE THE OTHER TC DATASET MADE IN PART ONE!!
# --------------------------------------------------
# 4a) geneLevel TPM expression for treatment
load('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES/RDATA/SE_geneLevel.Rdata', verbose=T)
# let's try when considering gene expression only in WT for the treatment
gene_expression <- data.frame('wt_0_tpm_gene'    = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_0')) %>% rowMeans(),
                              'wt_10_tpm_gene'   = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_10')) %>% rowMeans(),
                              'wt_30_tpm_gene'   = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_30')) %>% rowMeans(),
                              'hen2_0_tpm_gene'  = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_0')) %>% rowMeans(),
                              'hen2_10_tpm_gene' = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_10')) %>% rowMeans(),
                              'hen2_30_tpm_gene' = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_30')) %>% rowMeans(),
                              'rrp4_0_tpm_gene'  = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('rrp4_0')) %>% rowMeans(),
                              'rrp4_30_tpm_gene' = assay(geneLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('rrp4_30')) %>% rowMeans())
gene_expression$geneID <- rownames(gene_expression)
rownames(gene_expression) <- NULL
# add TC expression ratio to gene
tc_expression <- data.frame('wt_0_tpm_tc'    = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_0')) %>% rowMeans(),
                            'wt_10_tpm_tc'   = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_10')) %>% rowMeans(),
                            'wt_30_tpm_tc'   = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('wt_30')) %>% rowMeans(),
                            'hen2_0_tpm_tc'  = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_0')) %>% rowMeans(),
                            'hen2_10_tpm_tc' = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_10')) %>% rowMeans(),
                            'hen2_30_tpm_tc' = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('hen2_30')) %>% rowMeans(),
                            'rrp4_0_tpm_tc'  = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('rrp4_0')) %>% rowMeans(),
                            'rrp4_30_tpm_tc' = assay(txLevel) %>% cpm() %>% as.data.frame() %>% select(., matches('rrp4_30')) %>% rowMeans())
tc_expression$TCid <- rownames(tc_expression)
tc_expression$geneID <- rowRanges(txLevel)$geneID
rownames(tc_expression) <- NULL
# merge both
tpm_expression <- left_join(tc_expression, gene_expression, by='geneID')
rm(tc_expression, gene_expression)
# compute average TPM expression by timepoint, independent of genotype
tpm_expression$aveGeneExp_0 <- select(tpm_expression, matches('_0_tpm_gene')) %>% rowMeans()
tpm_expression$aveGeneExp_10 <- select(tpm_expression, matches('_10_tpm_gene')) %>% rowMeans()
tpm_expression$aveGeneExp_30 <- select(tpm_expression, matches('_30_tpm_gene')) %>% rowMeans()    
# compute average TC expression by timepoint, independent of genotype
tpm_expression$aveTCexp_0 <- select(tpm_expression, matches('_0_tpm_tc')) %>% rowMeans()
tpm_expression$aveTCexp_10 <- select(tpm_expression, matches('_10_tpm_tc')) %>% rowMeans()
tpm_expression$aveTCexp_30 <- select(tpm_expression, matches('_30_tpm_tc')) %>% rowMeans()
# reorder columns & save
tpm_expression <- select(tpm_expression, TCid, geneID,
                                         wt_0_tpm_tc, wt_10_tpm_tc, wt_30_tpm_tc, hen2_0_tpm_tc, hen2_10_tpm_tc, hen2_30_tpm_tc, rrp4_0_tpm_tc, rrp4_30_tpm_tc,
                                         wt_0_tpm_gene, wt_10_tpm_gene, wt_30_tpm_gene, hen2_0_tpm_gene, hen2_10_tpm_gene, hen2_30_tpm_gene, rrp4_0_tpm_gene, rrp4_30_tpm_gene,
                                         aveGeneExp_0, aveGeneExp_10, aveGeneExp_30,
                                         aveTCexp_0, aveTCexp_10, aveTCexp_30)
save(tpm_expression, file='~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES/RDATA/panExperiment_tc_and_gene_TPM_expression_mean_R123.Rdata')


# 5b) number of TSS per gene and their DE status
    # keep only DE info for timecourse
    TSS_per_gene <- tcdf %>% as.data.frame() %>% select(timepoint10, timepoint30, wt3010) %>% mutate('TCid'=rownames(.))
    # add gene for each TSS
    TSS_per_gene %<>% mutate('geneID'=rowRanges(txLevel)$geneID)
    # remove TSS not attributed to a gene
    TSS_per_gene %<>% subset(!is.na(geneID)) # 24,291 gene-attributed TSSs
    # add DE info: is the TSS DE at any treatment timepoint, and in any direction?
    TSS_per_gene$isTimecourseDE <- select(TSS_per_gene, -TCid, -geneID) %>% apply(1, function(x) ifelse(any(x != 0), 'DE', 'not DE'))
    TSS_per_gene$isTimecourseDE %<>% factor(levels=c('not DE', 'DE'))
    # add total number of known TSSs in the gene
    TSS_per_gene %<>% left_join(., table(TSS_per_gene$geneID) %>% as.data.frame() %>% set_colnames(c('geneID', 'knownTSSnb')), by='geneID')
    # add total number of DE TSSs in the gene
    tmp <- TSS_per_gene %>% group_by(geneID) %>% summarise('deTSSnb'=sum(isTimecourseDE=='DE'))
    TSS_per_gene %<>% left_join(., tmp, by='geneID')
    rm(tmp)
    # add structure (single/multi-TSS as a factor)
    TSS_per_gene %<>% mutate('structure'=ifelse(.$knownTSSnb > 1, 'multi-TSS', 'single-TSS') %>% as.factor())
    # plot
    head(TSS_per_gene)
    ggplot(TSS_per_gene, aes(x=as.factor(knownTSSnb), fill=interaction(structure, isTimecourseDE, lex.order=T, sep=' '))) +
           geom_bar() + labs(x='Gene structure', y='Nb. of TSSs', title='Flagellin-induced TSSs by gene structure') +
           scale_fill_brewer(palette='Paired', name='') + theme(aspect.ratio=1.5, legend.position='bottom') + guides(fill=guide_legend(ncol=2))

# 5c) Limma DTU detection : topSplice FDR setting is buggy, better to subset on FDR manually
ex <- diffSplice(fit=fit, geneid=rowRanges(txLevel)$geneID, exonid=names(rowRanges(txLevel)), robust=TRUE, verbose=TRUE)
# retrieve exon-level DTUs, using the t-statistic
dtu_timepoint10  <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='timepoint10'), FDR <= 0.05) ; dim(dtu_timepoint10) # 178
dtu_wt3010       <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='wt3010'),      FDR <= 0.05) ; dim(dtu_wt3010) # 367
dtu_timepoint30  <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='timepoint30'), FDR <= 0.05) ; dim(dtu_timepoint30) # 848
dtu_genotypehen2 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='genotypehen2'), FDR <= 0.05) ; dim(dtu_genotypehen2) # 237
dtu_genotyperrp4 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='genotyperrp4'), FDR <= 0.05) ; dim(dtu_genotyperrp4) # 498
dtu_genotypehen2.timepoint10 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='genotypehen2.timepoint10'), FDR <= 0.05) ; dim(dtu_genotypehen2.timepoint10) # 0
dtu_genotypehen2.timepoint30 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='genotypehen2.timepoint30'), FDR <= 0.05) ; dim(dtu_genotypehen2.timepoint30) # 0
dtu_genotyperrp4.timepoint30 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='genotyperrp4.timepoint30'), FDR <= 0.05) ; dim(dtu_genotyperrp4.timepoint30) # 2
dtu_hen3010 <- subset(topSplice(ex, test='t', FDR=1, number=Inf, coef='hen3010'), FDR <= 0.05) ; dim(dtu_hen3010) # 0
    # plot number of DTU genes per coefficient
    tot_dtu_per_coef <- data.frame('coefficient'=c('t=10', 't=30/10', 't=30', 'hen2', 'rrp4'),
                                   'DTU genes'=c( length(unique(dtu_timepoint10$GeneID)), length(unique(dtu_wt3010$GeneID)), length(unique(dtu_timepoint30$GeneID)), 
                                                  length(unique(dtu_genotypehen2$GeneID)), length(unique(dtu_genotyperrp4$GeneID))),
                                   'Base effect'=c(rep('Flg22', 3), rep('Exosome', 2)))
    tot_dtu_per_coef$coefficient %<>% factor(levels=c('t=10', 't=30/10', 't=30', 'hen2', 'rrp4'))
    ggplot(tot_dtu_per_coef, aes(x=coefficient, y=DTU.genes, fill=Base.effect)) + geom_bar(stat='identity') +
           labs(x='', y='Nb. DTU genes (1165)') + scale_fill_brewer(palette='Set1', name='Base effect', direction=-1) +
           geom_text(aes(label=DTU.genes), vjust=-1) + theme(aspect.ratio=1.5)

# merge all treatment DTUs
dtus <- rbind(dtu_timepoint10, dtu_timepoint30, dtu_wt3010) ; dim(dtus) # 1393
# add coefficient
dtus$coefficient <- c(rep('timepoint10', nrow(dtu_timepoint10)),
                      rep('timepoint30', nrow(dtu_timepoint30)),
                      rep('wt3010', nrow(dtu_wt3010)))
# add UP/DOWN factor for DE
dtus$direction <- ifelse(dtus$logFC > 0, 'up', 'down') %>% as.factor()
# add WT TPM expression of TCs and gene
dtus <- left_join(dtus, dplyr::select(tpm_expression, TCid, wt_0_tpm_tc, wt_0_tpm_gene, wt_10_tpm_tc, wt_10_tpm_gene, wt_30_tpm_tc, wt_30_tpm_gene), by=c('ExonID'='TCid'))
head(dtus)
subset(dtus, coefficient=='timepoint10' & direction=='up')[with(subset(dtus, coefficient=='timepoint10' & direction=='up'), aveTCexp_10 <= aveTCexp_0), ]




# # # # # # # # # #
# Actual Analysis #
# # # # # # # # # #
t10   <- subset(dtus, coefficient=='timepoint10') ; dim(t10) # 178
t3010 <- subset(dtus, coefficient=='wt3010') ; dim(t3010) # 367
t30   <- subset(dtus, coefficient=='timepoint30') ; dim(t30) # 848
# a) gene must have at least 2 TSSs: OK
all(t10$GeneID %in% multiTSSgenes$geneID) #ok
all(t3010$GeneID %in% multiTSSgenes$geneID) #ok
all(t30$GeneID %in% multiTSSgenes$geneID) #ok
# b) at least on TSS must be up-regulated
  # list of genes with one up-regulated TC in timepoint10
  lst_10 <- subset(df, timepoint10 == 1)$geneID %>% unique() # 552
  lst_3010 <- subset(df, wt3010 == 1)$geneID %>% unique() # 1456
  lst_30 <- subset(df, timepoint30 == 1)$geneID %>% unique() # 1847
  t10 <- subset(t10, GeneID %in% lst_10) # 90
  t3010 <- subset(t3010, GeneID %in% lst_3010) # 231
  t30 <- subset(t30, GeneID %in% lst_30) # 506
# c) the up-regulated TSS must contribute at least 10% of gene expression
  # calculate gene expression ratio
  t10 <- mutate(t10, 'ratio_t0'=aveTCexp_0/aveGeneExp_0*100,
                     'ratio_t10'=aveTCexp_10/aveGeneExp_10*100)
  t3010 <- mutate(t3010, 'ratio_t10'=aveTCexp_10/aveGeneExp_10*100,
                         'ratio_t30'=aveTCexp_30/aveGeneExp_30*100)
  t30 <- mutate(t30, 'ratio_t0'=aveTCexp_0/aveGeneExp_0*100,
                     'ratio_t30'=aveTCexp_30/aveGeneExp_30*100)
  # make list of genes with TSS up >= 10%
  lst_10 <- subset(t10, direction == 'up' & ratio_t10 >= 10)$GeneID %>% unique() # 32
  lst_3010 <- subset(t3010, direction == 'up' & ratio_t30 >= 10)$GeneID %>% unique() # 90
  lst_30 <- subset(t30, direction == 'up' & ratio_t30 >= 10)$GeneID %>% unique() # 200
  t10 <- subset(t10, GeneID %in% lst_10) # 61
  t3010 <- subset(t3010, GeneID %in% lst_3010) # 186
  t30 <- subset(t30, GeneID %in% lst_30) # 396

t10$GeneID %>% unique() # 32
t3010$GeneID %>% unique() # 90
t30$GeneID %>% unique() # 200

# Venn diagram of DTUs
venn(x=list(t10$GeneID, t3010$GeneID, t30$GeneID),
     snames=c('DTUgenes_10', 'DTUgenes_30/10','DTUgenes_30'),
     cexil=2, cexsn=1,
     zcolor=c('green', 'blue', 'orange'))


# d) categorization of alternative TSSs:
  # one TSS up, another is down = proper TSS switch
  dim(dtus)
  unique(dtus$GeneID) %>% length()
  test <- dtus %>% group_by(GeneID, coefficient) %>% summarise('any_up'=any(direction=='up'), 'any_down'=any(direction=='down')) %>% as.data.frame()
  subset(dtus, GeneID %in% subset(test, any_up==T & any_down==T)$GeneID)
  subset(test, any_up==T & any_down==F) %>% dim()
  categories <- table('timepoint'=test$coefficient, 'UP'=test$any_up, 'DOWN'=test$any_down) %>% melt() %>% as.data.frame()
  sum(categories$value)
  # one TSS up, no other not down = alternative alternative TSS (by default those that don't qualify in the previous criteria)


  
# EXPORT DATA
# ALL PROMPTS: 96 TCs that are annotated as 'reverse' and up-regualted in hen2 or rrp4
rowRanges(TCs) %>%
  subset(txType_TAIR10=='reverse' & (genotypehen2==1 | genotyperrp4==1)) %>%
  as.data.frame() %>%
  left_join(idmapping, by=c('geneID_anti'='geneID')) %>%
  WriteXLS::WriteXLS(ExcelFileName='~/Desktop/prompts.xlsx', row.names=F, col.names=T, AdjWidth=T, BoldHeaderRow=T)

readxl::read_xlsx('~/Dropbox/Flg22_CAGE/Exosome_TSS paper/PLANT CELL/Supplementary Data S1 - CAGE TCs.xlsx', sheet='Unidirectional TCs', col_names=T) %>%
  left_join(idmapping, by=c('geneID_TAIR10'='geneID')) %>%
  rename('gene_symbol'='name', 'gene_description'='description') %>%
  left_join(idmapping, by=c('geneID_TAIR10_anti'='geneID')) %>%
  rename('gene_anti_symbol'='name', 'gene_anti_descrition'='description') %>%
  WriteXLS::WriteXLS('~/Desktop/Uni_TCs.xlsx', row.names=F, col.names=T)

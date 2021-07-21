#### Arabidopsis flg22 : Alternative TSSs
#### Axel Thieffry
set.seed(42)
library(WriteXLS)
library(ggpubr)
library(RColorBrewer)
library(edgeR)
library(limma)
library(DESeq2)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(org.At.tair.db)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(CAGEfightR)
library(tidyverse)
library(venn)
'select' <- dplyr::select
'%!in%' <- function(x,y)!('%in%'(x,y))
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                    if(length(idx) != 0) { GR[-idx]}
                                    else {GR}}
setwd('~/masked_path/04 - TSS_Level DE/')


# 1. LOAD ALL INPUT FILES ####
# ----------------------------
# myseqinfo
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')

# INTRAGENIC TSSs for DTU ANALYSIS
# Those are TSSs associated to a gene and contributing at least 10% to the total gene expression
intraTSSs <- readRDS('~/masked_path/SE_intragenicTSSs_for_DTU.rds')
    # drop empty annotation levels
    rowRanges(intraTSSs)$txType_TAIR10 %<>% droplevels()
    # sanity check
    rowRanges(intraTSSs) %$% table(.$txType_TAIR10) %>% as.data.frame()

# AGI ID mapping
idmapping <- readRDS('~/masked_path/AGI_ID_mapping.rds')

# design and contrast
design <- readRDS('~/masked_pathA/limma_design_matrix.rds')
contrasts <- readRDS('~/masked_path/limma_contrast_matrix.rds')

# DE TSSs
dt <- readRDS('~/masked_path/DE_TCs_summary.rds')
tt_timepoint10 <- readRDS('~/masked_path/DE_TSSs_topTable_timepoint10.rds')
tt_wt3010 <- readRDS('~/masked_path/DE_TSSs_topTable_wt3010.rds')
tt_timepoint30 <- readRDS('~/masked_path/DE_TSSs_topTable_timepoint30.rds')



# 3. ALTERNATIVE TSSs ####
# ------------------------
# 3a. TSSs per gene with 10% expression contribution
rowRanges(intraTSSs) %>%
  as.data.frame() %$%
  table(.$geneID) %>%
  as.data.frame() %$%
  table(.$Freq) %>%
  as.data.frame() %>%
  mutate('Var1'=as.integer(Var1)) %>%
  ggplot(aes(x=Var1, y=Freq, fill=ifelse(Var1==1, 'Single-TSS gene', '>1 TSS gene'))) + geom_bar(stat='identity') +
         labs(title='Gene structure', x='Number of TSSs per gene', y=paste0('Number of Genes / ', length(intraTSSs)), subtitle='Intragenic TSSs with 10%\nexpression contribution') +
         geom_text(aes(label=Freq), vjust=-.5, size=3) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1.5, legend.position='bottom') +
         scale_fill_manual(values=c('orange', 'grey40'), name='')


# 3b. annotation of TSSs per :
# single/multiple-TSS gene
# dominant/non-dominant TSS

  # 3b.1 simple annotation of single-TSS genes
      # get single-TSS genes
      singleTSS_genes <- rowRanges(intraTSSs) %>%
              as.tibble() %$%
              table(.$geneID) %>%
              as.tibble() %>%
              subset(n==1) %>%
              set_colnames(c('geneID', 'nTSSs')) %$%
              .$geneID
      # annotation barplot
      rowRanges(intraTSSs) %>%
        as.data.frame() %>%
        subset(geneID %in% singleTSS_genes) %>%
        ggplot(aes(x=txType_TAIR10)) +
               geom_bar(stat='count') + cowplot::theme_cowplot() + theme(legend.position='bottom', aspect.ratio=10) +
               facet_grid(~txType_TAIR10, scales='free_x') + ylim(NA, 15500)
    
  # 3b.2 annotation of multi-TSS genes by dominant/non-dominant TSS
      # get multi-TSS genes
      multiTSS_genes <- rowRanges(intraTSSs) %>%
        as.tibble() %$%
        table(.$geneID) %>%
        as.tibble() %>%
        subset(n!=1) %>%
        set_colnames(c('geneID', 'nTSSs')) %$%
        .$geneID
      # get dominant TSSs
      domTSS <- subset(intraTSSs, geneID %in% multiTSS_genes) %>%
        rowRanges() %>%
        as.tibble() %>%
        group_by(geneID) %>%
        filter(score==max(score)) %$%
        .$thick.names
      # annotation barplot
      rowRanges(intraTSSs) %>%
        as.data.frame() %>%
        subset(geneID %in% multiTSS_genes) %>%
        ggplot(aes(x=txType_TAIR10, fill=ifelse(thick.names %in% domTSS, 'Dominant TSS', 'Non-Dominant TSS'))) +
               geom_bar(stat='count') + cowplot::theme_cowplot() + theme(legend.position='bottom', aspect.ratio=10) +
               scale_fill_manual(values=c('darkorange', 'gold'), name='') +
               facet_grid(~txType_TAIR10, scales='free_x') + ylim(NA, 15000)


# 2c. annotation of alternative TSSs with number of DE (direction does not matter)
de_TSSs_any_timecourse <- dt %>%
                          as.data.frame() %>%
                          select(timepoint10, wt3010, timepoint30) %>%
                          subset(apply(., 1, function(x) sum(x!=0) != 0)) %>%
                          rownames_to_column('ID')

de_TSSs_any_exosome <- dt %>%
                       as.data.frame() %>%
                       select(genotypehen2, genotyperrp4) %>%
                       subset(apply(., 1, function(x) sum(x!=0) != 0)) %>%
                       rownames_to_column('ID')

rowRanges(intraTSSs) %>%
  as.data.frame() %>%
  rownames_to_column('ID') %>%
  as.tibble() %>%
  subset(txType_TAIR10 != 'promoter') %>%
  mutate('txType_TAIR10' = factor(txType_TAIR10, levels=c('proximal', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))) %>%
  mutate('isTimecourseDE' = ifelse(ID %in% de_TSSs_any_timecourse$ID, TRUE, FALSE)) %>%
  ggplot(aes(x=txType_TAIR10, fill=isTimecourseDE)) + geom_bar(stat='count') +
         geom_text(stat='count', aes(label=..count..), position=position_stack(.5), size=3) +
         labs(title='Alternative TSSs in timecourse', x='', y='count',
              caption='') +
         scale_fill_brewer(palette='Set2', direction=-1) +
         theme(aspect.ratio=1.5, legend.position=c(.8, .9), axis.text.x=element_text(angle=45, hjust=1))

gg_d <- rowRanges(intraTSSs) %>%
  as.data.frame() %>%
  rownames_to_column('ID') %>%
  subset(txType_TAIR10 != 'promoter') %>%
  mutate('txType_TAIR10' = factor(txType_TAIR10, levels=c('proximal', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR'))) %>%
  mutate('isExosomeDE' = ifelse(ID %in% de_TSSs_any_exosome$ID, TRUE, FALSE)) %>%
  ggplot(aes(x=txType_TAIR10, fill=isExosomeDE)) + geom_bar(stat='count') +
         geom_text(stat='count', aes(label=..count..), position=position_stack(.5), size=3) +
         labs(title='Alternative TSSs in exosome', x='', y='count',
              caption='all TSSs here are >= 10% contribution to total gene expression') +
         scale_fill_brewer(palette='Accent', direction=-1) +
         theme(aspect.ratio=1.5, legend.position=c(.8, .9), axis.text.x=element_text(angle=45, hjust=1))

ggpubr::ggarrange(plotlist=list(gg_c, gg_d), align='hv')



# 2. LIMMA DIFFSPLICE ####
# ------------------------
dge <- DGEList(assay(intraTSSs, 'counts')) %>% calcNormFactors(method='TMM')
v   <- voom(dge, design=design, plot=FALSE)
fit <- lmFit(v, design=design)
fit <- contrasts.fit(fit, contrasts)
eb  <- eBayes(fit, robust=TRUE)

ex <- diffSplice(fit=eb, geneid=rowRanges(intraTSSs)$geneID, exonid=names(rowRanges(intraTSSs)), robust=TRUE, verbose=TRUE)

# retrieve exon-level DTUs, using the t-statistic (FDR is the BH correction)
dtu_timepoint10  <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='timepoint10') %>% mutate('coefficient'=factor('timepoint10')) %>% as.tibble() #97
dtu_wt3010       <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='wt3010') %>% mutate('coefficient'=factor('wt3010')) %>% as.tibble() # 172
dtu_timepoint30  <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='timepoint30') %>% mutate('coefficient'=factor('timepoint30')) %>% as.tibble() # 331
dtu_genotypehen2 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='genotypehen2') %>% mutate('coefficient'=factor('genotypehen2')) %>% as.tibble() # 65
dtu_genotyperrp4 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='genotyperrp4') %>% mutate('coefficient'=factor('genotyperrp4')) %>% as.tibble() # 145
dtu_genotypehen2.timepoint10 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='genotypehen2.timepoint10') %>% mutate('coefficient'=factor('genotypehen2.timepoint10')) %>% as.tibble() # 0
dtu_hen3010 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='hen3010') %>% mutate('coefficient'=factor('hen3010')) %>% as.tibble() # 1
dtu_genotypehen2.timepoint30 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='genotypehen2.timepoint30') %>% mutate('coefficient'=factor('genotypehen2.timepoint30')) %>% as.tibble() # 0
dtu_genotyperrp4.timepoint30 <- topSplice(ex, test='simes', FDR=0.05, number=Inf, coef='genotyperrp4.timepoint30') %>% mutate('coefficient'=factor('genotyperrp4.timepoint30')) %>% as.tibble() # 6

# plot number of DTU genes per coefficient
coef_levels <- factor(c('timepoint10', 'wt3010', 'timepoint30', 'genotypehen2', 'genotyperrp4', 'genotypehen2.timepoint10', 'hen3010', 'genotypehen2.timepoint30', 'genotyperrp4.timepoint30'))

dtu_per_coef <- rbind(dtu_timepoint10, dtu_wt3010, dtu_timepoint30,
                      dtu_genotypehen2, dtu_genotyperrp4,
                      dtu_genotypehen2.timepoint10, dtu_hen3010, dtu_genotypehen2.timepoint30, dtu_genotyperrp4.timepoint30) %>%
                  mutate('coefficient'=factor(.$coefficient, levels=coef_levels)) %>%
                  group_by(coefficient) %>%
                  summarise('count'=n()) %>%
                  complete(coefficient, fill=list(count=0)) %>%
                  mutate('effect' = c(rep('timecourse', 3), rep('exosome', 2), rep('interaction', 4)) %>%
                           factor(levels=c('timecourse', 'exosome', 'interaction'))) %>%
                  ungroup()

ggplot(dtu_per_coef, aes(x=coefficient, y=count, fill=effect)) +
       geom_bar(stat='identity') +
       geom_text(data=subset(dtu_per_coef, count!=0), aes(label=count), vjust=-1) + 
       theme(axis.text.x=element_text(angle=45, hjust=1)) +
       scale_fill_brewer(palette='Set1') +
       labs(title='DTU genes per coefficient')

# merge results
dtus <- rbind(dtu_timepoint10, dtu_wt3010, dtu_timepoint30,
              dtu_genotypehen2, dtu_genotyperrp4,
              dtu_genotypehen2.timepoint10, dtu_hen3010, dtu_genotypehen2.timepoint30, dtu_genotyperrp4.timepoint30)

dtus %>%
  group_by(coefficient) %>%
  summarise('n_DTU'=n()) %>%
  mutate('effect'=c(rep('timecourse', 3), rep('exosome mutant', 2), rep('interaction', 2))) %>%
  ggplot(aes(x=coefficient, y=n_DTU, fill=effect)) + geom_bar(stat='identity') +
         theme(axis.text.x=element_text(angle=45, hjust=1), aspect.ratio=1) +
         labs(title='Differential Transcript Usage per coefficient', x='', y='Nb. Genes') +
         scale_fill_brewer(palette='Set2', direction=-1)

# Venn diagram of DTUs - Gene level
venn(x=list(unique(dtu_timepoint10$GeneID), unique(dtu_wt3010$GeneID), unique(dtu_timepoint30$GeneID)),
     snames=c('timepoint10', 'contrast 30/10','timepoint30'),
     cexil=2, cexsn=1,
     zcolor=c('green', 'blue', 'orange'))

# Venn diagram of DTUs - Gene level
venn(x=list(unique(dtu_genotypehen2$GeneID), unique(dtu_genotyperrp4$GeneID)),
     snames=c('genotypehen2', 'genotyperrp4'),
     cexil=2, cexsn=1,
     zcolor=c('green', 'blue', 'orange'))

# how many DTU per gene?
table(dtus$GeneID) %>%
  as.data.frame() %$%
  table(.$Freq) %>%
  as.data.frame() %>%
  set_colnames(c('DTU_TSSs_per_gene', 'count')) %>%
  ggplot(aes(x=DTU_TSSs_per_gene, y=count)) + geom_bar(stat='identity') +
         theme(aspect.ratio=1) + 
         labs(title='DTU TSSs frequency per gene', x='DTU TSSs', y='N(genes)')

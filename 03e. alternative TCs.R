#### 03d. Alternative TCS 
#### Axel Thieffry - July 2019
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
register(MulticoreParam(workers=3))
library(TxDb.Athaliana.BioMart.plantsmart28)
options(scipen=999) # disable scientific notation
'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))
setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/03 - TSS analysis')




# 1. GET INPUT DATA ####
# ----------------------
# CAGE TCs
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_TPM1_min3lib_TSSstory.rds')
# CAGE genes (timepoint 0 only)
geneLevel <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_genelevel_TSSstory.rds')
# CAGE DE (timepoint 0 only)
de_hen2 <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotypehen2.rds')
de_rrp4 <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/DE_TSSs_topTable_genotyperrp4.rds')


# 2. COMPUTE TPM EXPRESSION ####
# ------------------------------
# TCs: compute average TPM expression per genotype
rowRanges(TCs)$wt_TC_TPM <- TCs %>%
  subset(select=genotype=='wt') %>%
  assay('TPM') %>%
  rowMeans()

rowRanges(TCs)$hen2_TC_TPM <- TCs %>%
  subset(select=genotype=='hen2') %>%
  assay('TPM') %>%
  rowMeans()

rowRanges(TCs)$rrp4_TC_TPM <- TCs %>%
  subset(select=genotype=='rrp4') %>%
  assay('TPM') %>%
  rowMeans()

# geneLevel: compute average TPM expression per genotype
geneLevel %<>% calcTPM(inputAssay='counts', outputAssay='TPM', outputColumn='score', totalTags='totalTags')

gene_df <- data.frame('geneID'=names(geneLevel),
                      'wt_gene_TPM'=geneLevel %>% subset(select=genotype=='wt') %>% assay('TPM') %>% rowMeans(),
                      'hen2_gene_TPM'=geneLevel %>% subset(select=genotype=='hen2') %>% assay('TPM') %>% rowMeans(),
                      'rrp4_gene_TPM'=geneLevel %>% subset(select=genotype=='rrp4') %>% assay('TPM') %>% rowMeans()) %>%
  as_tibble()

# merge all together: TC, TC TPMs, gene, gene TPMs, TC annotation
data <- rowRanges(TCs) %>%
  as.data.frame() %>%
  select(thick.names, txType_TAIR10, geneID, symbol, wt_TC_TPM:rrp4_TC_TPM) %>%
  subset(!is.na(geneID)) %>%
  as_tibble()

data %<>% left_join(gene_df, by='geneID')
data %<>% select(thick.names, txType_TAIR10, geneID, symbol, wt_TC_TPM, wt_gene_TPM, hen2_TC_TPM, hen2_gene_TPM, rrp4_TC_TPM, rrp4_gene_TPM)

# calculate gene expression contribution (as percentage)
data %<>% mutate('wt_contrib'=wt_TC_TPM / wt_gene_TPM * 100,
                 'hen2_contrib'=hen2_TC_TPM / hen2_gene_TPM * 100,
                 'rrp4_contrib'=rrp4_TC_TPM / rrp4_gene_TPM * 100)



# 3. NUMBER OF TCs PER GENE ####
# ------------------------------
# filter for 10% contribution to total gene expression (WT only)
data_wt10 <- subset(data, wt_contrib >= 10.0)

# in wt only, filter for 10% contribution
data_wt10 %>%
  select(thick.names:wt_gene_TPM, wt_contrib) %>%
  group_by(geneID) %>%
  summarise('n_TCs'=n()) %>%
  ungroup() %>%
  ggplot(aes(x=as.factor(n_TCs), fill=n_TCs==1)) +
         geom_bar(stat='count', col='black', lwd=.3) +
         geom_text(stat='count', aes(label=..count..), angle=90, hjust=-.1) +
         cowplot::theme_cowplot() +
         scale_fill_manual(values=c('orange', 'grey40'), name='') +
         labs(title='Gene structure in wt (t=0)',
              subtitle='Intragenic TSSs',
              x='Nb. CAGE wt TCs per gene',
              y=paste0('Nb. Genes (N=', n_distinct(data_wt10$geneID), ')'))

# same graph but comparing with filtering for: 10%, 5%, 1%
data_wt5 <- subset(data, wt_contrib >= 5.0)
data_wt1 <- subset(data, wt_contrib >= 1.0)

rbind(data_wt10 %>%
        select(thick.names:wt_gene_TPM, wt_contrib) %>%
        group_by(geneID) %>%
        summarise('n_TCs'=n()) %>%
        ungroup() %>%
        mutate('set'='10%'),
      data_wt5 %>%
        select(thick.names:wt_gene_TPM, wt_contrib) %>%
        group_by(geneID) %>%
        summarise('n_TCs'=n()) %>%
        ungroup() %>%
        mutate('set'='5%'),
  data_wt1 %>%
        select(thick.names:wt_gene_TPM, wt_contrib) %>%
        group_by(geneID) %>%
        summarise('n_TCs'=n()) %>%
        ungroup() %>%
        mutate('set'='1%')) %>%
  mutate('set'=factor(set, levels=c('1%', '5%', '10%'))) %>%
  ggplot(aes(x=as.factor(n_TCs), fill=n_TCs==1)) +
         geom_bar(stat='count', col='black', lwd=.3) +
         geom_text(stat='count', aes(label=..count..), angle=90, hjust=-.1) +
         cowplot::theme_cowplot() +
         facet_grid(set~.) +
         scale_fill_manual(values=c('orange', 'grey40'), name='') +
         labs(title='Gene structure in wt (t=0)',
              subtitle='Intragenic TSSs',
              x='Nb. CAGE wt TCs per gene',
              y=paste0('Nb. Genes (N=', n_distinct(data_wt10$geneID), ')'))


# total detected genes in wt: 18,330
n_distinct(data_wt10$geneID)
# 1,839 have more than one TC with 10% min: exactly 10 %
table(data_wt10$geneID) %>%
  enframe(name='geneID', value='nb_TCs') %>%
  subset(nb_TCs > 1)



# 4. ANNOTATION OF SINGLE AND MULTI TSS GENES ####
# ------------------------------------------------
# add if TC is dominant or marginal (in expression)
wt10_dominant_TCs <- data_wt10 %>%
  group_by(geneID) %>%
  filter(wt_contrib==max(wt_contrib)) %>%
  ungroup() %>%
  .$thick.names

wt5_dominant_TCs <- data_wt5 %>%
  group_by(geneID) %>%
  filter(wt_contrib==max(wt_contrib)) %>%
  ungroup() %>%
  .$thick.names

wt1_dominant_TCs <- data_wt1 %>%
  group_by(geneID) %>%
  filter(wt_contrib==max(wt_contrib)) %>%
  ungroup() %>%
  .$thick.names

data_wt10 %<>% mutate('wt_dominance'=ifelse(thick.names %in% wt10_dominant_TCs, 'major', 'minor'))
data_wt10 %<>% mutate('wt_dominance'=factor(wt_dominance, levels=c('major', 'minor')))

data_wt5 %<>% mutate('wt_dominance'=ifelse(thick.names %in% wt5_dominant_TCs, 'major', 'minor'))
data_wt5 %<>% mutate('wt_dominance'=factor(wt_dominance, levels=c('major', 'minor')))

data_wt1 %<>% mutate('wt_dominance'=ifelse(thick.names %in% wt1_dominant_TCs, 'major', 'minor'))
data_wt1 %<>% mutate('wt_dominance'=factor(wt_dominance, levels=c('major', 'minor')))

# add if gene is single or multi tss
wt10_singleTC_genes <- table(data_wt10$geneID) %>%
  enframe(value='nb_TCs', name='geneID') %>%
  subset(nb_TCs==1) %>%
  .$geneID

wt5_singleTC_genes <- table(data_wt5$geneID) %>%
  enframe(value='nb_TCs', name='geneID') %>%
  subset(nb_TCs==1) %>%
  .$geneID

wt1_singleTC_genes <- table(data_wt1$geneID) %>%
  enframe(value='nb_TCs', name='geneID') %>%
  subset(nb_TCs==1) %>%
  .$geneID

data_wt10 %<>% mutate('gene_structure'=ifelse(geneID %in% wt10_singleTC_genes, 'single_TSS', 'multi_TSS'))
data_wt10 %<>% mutate('gene_structure'=factor(gene_structure, levels=c('single_TSS', 'multi_TSS')))

data_wt5 %<>% mutate('gene_structure'=ifelse(geneID %in% wt5_singleTC_genes, 'single_TSS', 'multi_TSS'))
data_wt5 %<>% mutate('gene_structure'=factor(gene_structure, levels=c('single_TSS', 'multi_TSS')))

data_wt1 %<>% mutate('gene_structure'=ifelse(geneID %in% wt1_singleTC_genes, 'single_TSS', 'multi_TSS'))
data_wt1 %<>% mutate('gene_structure'=factor(gene_structure, levels=c('single_TSS', 'multi_TSS')))

# plot
tricolors <- brewer.pal(n=3, name='Paired')[c(1, 2, 3)]

ggplot(data_wt10, aes(x=txType_TAIR10, fill=interaction(wt_dominance, gene_structure))) +
       geom_bar(stat='count', col='black', lwd=.3, position=position_dodge()) +
       geom_text(stat='count', aes(label=..count..), angle=90, hjust=-.1, vjust=0.5, position=position_dodge(width=0.9)) +
       cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='none') +
       scale_fill_manual(values=tricolors, name='') +
       labs(title='10%', y='Nb. CAGE wt TCs', x='')

rbind(data_wt10 %>% mutate('set'='10%'),
      data_wt5 %>% mutate('set'='5%'),
      data_wt1 %>% mutate('set'='1%')) %>%
  mutate('set'=factor(set, levels=c('1%', '5%', '10%'))) %>%
  ggplot(aes(x=txType_TAIR10, fill=interaction(wt_dominance, gene_structure))) +
         geom_bar(stat='count', col='black', lwd=.3, position=position_dodge()) +
         geom_text(stat='count', aes(label=..count..), angle=90, hjust=-.1, vjust=0.5, position=position_dodge(width=0.9)) +
         cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position='bottom') +
         scale_fill_manual(values=tricolors, name='') +
         labs(title='10%', y='Nb. CAGE wt TCs', x='') +
         facet_grid(~set)
  
  

# get examples of genes with two TCs: one marginal in promoter and one secondary in CDS
bi_TSS_genes <- data_wt10 %>%
  select(thick.names:wt_gene_TPM, wt_contrib) %>%
  group_by(geneID) %>%
  summarise('n_TCs_10pc'=n()) %>%
  ungroup() %>%
  subset(n_TCs_10pc == 2) %$%
  .$geneID

subset(data_wt10, geneID %in% bi_TSS_genes) %$%
  table(.$geneID, .$txType_TAIR10) %>%
  as.data.frame() %>%
  as_tibble() %>%
  spread(Var2, Freq) %>%
  rename('geneID'='Var1') %>%
  subset(promoter==1 & CDS == 1) %>% View()
  
  



# 5. EXOSOME SENSITIVITY OF ALTERNATIVE TCs ####
# ----------------------------------------------
de_hen2_up <- subset(de_hen2, logFC > 0)
de_rrp4_up <- subset(de_rrp4, logFC > 0)
de_any_up <- unique(c(de_hen2_up$id, de_rrp4_up$id))

# add if up-regulated in (any) exosome
data_wt10$exosome_up <- ifelse(data_wt10$thick.names %in% de_any_up, 'sensitive', 'insensitive') %>% factor(levels=c('sensitive', 'insensitive'))

# add number of TCs in gene
wt10_nbTCs_inGene <- data_wt10 %>%
  group_by(geneID) %>%
  summarise('n_TCs'=n())

data_wt10 %<>% left_join(wt10_nbTCs_inGene, by='geneID')

# barplot number of genes
data_wt10 %>%
  ggplot(aes(x=as.factor(n_TCs), fill=exosome_up)) +
         geom_bar(stat='count', col='black', lwd=.3) +
         geom_text(stat='count', aes(label=..count.., col=exosome_up), position=position_stack(vjust=1.1)) +
         cowplot::theme_cowplot() +
         scale_fill_brewer(palette='Set2', direction=-1) +
         scale_color_brewer(palette='Set2', direction=-1) +
         labs(x='Nb. CAGE wt TCs per gene',
              y=paste0('Nb. Genes (N=', n_distinct(data_wt10$geneID), ')'),
              title='Exosome-sensitive TCs per gene')

# number of TCs in multi-TC genes that are not exosome sensitive
data_wt10 %>%
  subset(n_TCs > 1)

# plot density ridges
View(data_wt10)

data_wt10 %>%
  subset(exosome_up == 'sensitive') %>%
  select(n_TCs, wt_TC_TPM, hen2_TC_TPM, rrp4_TC_TPM) %>%
  melt(id.vars='n_TCs', variable.name='genotype') %>%
  separate(genotype, 'genotype', sep='_', extra='drop') %>%
  mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=value, y=as.factor(n_TCs), col=genotype, fill=genotype)) +
         geom_density_ridges(alpha=.3) +
         cowplot::theme_cowplot() +
         scale_x_log10() +
         scale_fill_brewer(palette='Dark2') +
         scale_color_brewer(palette='Dark2') +
         labs(title='Exosome sensitive CAGE TCs',
              x='CAGE expression (TPM)',
              y='Nb. TCs in gene')

data_wt10 %>%
  subset(exosome_up == 'sensitive') %$%
  table(.$n_TCs)

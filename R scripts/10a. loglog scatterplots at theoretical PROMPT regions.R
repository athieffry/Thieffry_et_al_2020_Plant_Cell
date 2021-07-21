#### Arabidopsis TSS/exosome: scatter plot
#### Axel Thieffry
set.seed(42)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(TeMPO)
library(rtracklayer)
library(viridis)
library(stringr)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(BSgenome.Athaliana.TAIR.TAIR9)
library(CAGEfightR)
library(tidyverse)
library(tidylog)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
remove_out_of_bound <- function(GR) {idx=GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { o <- GR[-idx]
                                                            message(paste0('*** removed ', length(idx), ' out-of-bound regions ***'))
                                     }
else {o <- GR}
o}

setwd('~/masked_path/together')

# 1. LOAD ALL INPUT FILES ####
# ----------------------------
myseqinfo <- readRDS('~/masked_path/myseqinfo.rds')
idmapping <- readRDS('~/masked_path/AGI_ID_mapping.rds')
TCs <- readRDS('~/masked_path/SE_TCs_TPM1_min3lib_TSSstory.rds')
seqlevels(TCs) <- c(paste0('Chr', 1:5), 'ChrC', 'ChrM')

# small RNA-seq BIGWIGs (Dec 2019) : all lengths - RAW COUNTS
    # leaf (rosette)
    sRNA_plus_counts_leaf <- list.files('~/masked_path/3.bigwigs_alltogether', pattern='*Forward', full.names=T) %>% BigWigFileList()
    sRNA_minus_counts_leaf <- list.files('~/masked_path/3.bigwigs_alltogether', pattern='*Reverse', full.names=T) %>% BigWigFileList()
    sRNA_names_counts_leaf <- list.files('~/masked_path/3.bigwigs_alltogether', pattern='*Forward') %>%
                              str_remove('.Forward.raw.bw')
    names(sRNA_plus_counts_leaf) <- names(sRNA_minus_counts_leaf) <- sRNA_names_counts_leaf
    # flowers
    sRNA_plus_counts_flower <- list.files('~/masked_path/2.bigWigs_rawcounts_all_lengths', pattern='*Forward', full.names=T) %>% BigWigFileList()
    sRNA_minus_counts_flower <- list.files('~/masked_path/2.bigWigs_rawcounts_all_lengths', pattern='*Reverse', full.names=T) %>% BigWigFileList()
    sRNA_names_counts_flower <- list.files('~/masked_path/2.bigWigs_rawcounts_all_lengths', pattern='*Forward') %>%
                                str_remove('.Forward.bw')
    names(sRNA_plus_counts_flower) <- names(sRNA_minus_counts_flower) <- sRNA_names_counts_flower

# small RNA-seq stats
    # leaf (rosette)
    stats_leaf <- readxl::read_xlsx('~/masked_path/Sample_info.xlsx', col_names=T, trim_ws=T) %>%
                  as_tibble() %>%
                  select(-Raw_read_length)
    # flower
    stats_flower <- readxl::read_xlsx('~/masked_path/Sample_info.xlsx', col_names=T, trim_ws=T) %>%
                    as_tibble() %>%
                    select(-Raw_read_length)

# get genes hosting a PROMPT CAGE TC
# get PROMPT TC GR (N=96)
PROMPT_TCs <- readRDS('~/masked_path/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds') %>%
  rowRanges() %>%
  subset(genotypehen2==1 | genotyperrp4==1) %>%
  subset(txType_TAIR10=='reverse')
# get PROMPT TC ids (N=96)
PROMPT_TC_ids <- names(PROMPT_TCs)
# get ids of genes hosting a PROMPT TC (N=94)
PROMPT_TC_genes <- PROMPT_TCs$geneID_anti %>% unique()



# 2. PREPARE NORMALIZATION FACTORS ####
# -------------------------------------
normFact_leaf <- with(stats_leaf, data.frame('Sample'=Sample,
                                             'total_raw_reads_normFact' = Total_raw_reads / 1000000,
                                             'total_clean_reads_normFact' = Passing_filters / 1000000,
                                             'unique_mappers_normFact' = Unique_mappers / 1000000,
                                             'all_mappers_normFact' = (Unique_mappers + Multi_mappers) / 1000000))

normFact_flower <- with(stats_flower, data.frame('Sample'=str_remove(Sample, '-5'),
                                                 'total_raw_reads_normFact' = Total_raw_reads / 1000000,
                                                 'total_clean_reads_normFact' = Passing_filters / 1000000,
                                                 'unique_mappers_normFact' = Unique_mappers / 1000000,
                                                 'all_mappers_normFact' = (Unique_mappers + Multi_mappers) / 1000000))



# 3. MAKE THEORETICAL PROMPT REGIONS: -500 bp from annotated TSSs ###
# -------------------------------------------------------------------
# 3a) get PROMPT regions (-500 bp from TSS) (N=33542)
theo_prompts <- suppressWarnings( genes(txdb) %>%
                                    promoters(upstream=500, downstream=0) %>%
                                    remove_out_of_bound() %>% # remove out of bound regions
                                    unique() ) # keep unique regions only
      # fix seqinfo & keep only canonical chrs
      seqlevelsStyle(theo_prompts) <- seqlevelsStyle(myseqinfo)
      seqinfo(theo_prompts) <- myseqinfo
      seqlevels(theo_prompts, pruning.mode='coarse') <- paste0('Chr', 1:5) # (N=33623)

# 3b) remove PROMPT regions that overlap any other PROMT
theo_prompts_not_overlapping_any_other_theo_prompt <- countOverlaps(theo_prompts, ignore.strand=T) %>%
  enframe() %>%
  subset(value == 1) %>%
  pull(name)

theo_prompts %<>% subset(gene_id %in% theo_prompts_not_overlapping_any_other_theo_prompt) # (N=26033)

# 3c) remove PROMPT regions that overlap any other annotated feature (gene)
genes <- genes(txdb)
      # fix seqinfo
      seqlevelsStyle(genes) <- seqlevelsStyle(myseqinfo)
      seqinfo(genes) <- myseqinfo

theo_prompts_not_overlapping_any_other_annotated_feature <- countOverlaps(theo_prompts, genes, ignore.strand=T) %>%
                                                            enframe() %>%
                                                            subset(value == 0) %>% # not sure if should be zero or <= 1
                                                            pull(name)

theo_prompts %<>% subset(gene_id %in% theo_prompts_not_overlapping_any_other_annotated_feature) # (N=20659)

# 3d) make 1-bp version
theo_prompts_1bp <- resize(theo_prompts, width=1, fix='start')



# 4. COMPUTE smRNA COVERAGE AT THEORETICAL PROMPT REGIONS ####
# ------------------------------------------------------------
if(FALSE) {
      # 4a) compute sum of signal (must be downstream, see above in 4)
      wideMeta_signal_theo_prompts_leaf_altogether <- mapply(function(x, y) wideMetaProfile(sites=theo_prompts_1bp,
                                                                                            forward=x, reverse=y,
                                                                                            downstream=500, upstream=1),
                                                             sRNA_plus_counts_leaf, sRNA_minus_counts_leaf, SIMPLIFY=F)
      
      wideMeta_signal_theo_prompts_flower_altogether <- mapply(function(x, y) wideMetaProfile(sites=theo_prompts_1bp,
                                                                                              forward=x, reverse=y,
                                                                                              downstream=500, upstream=1),
                                                               sRNA_plus_counts_flower, sRNA_minus_counts_flower, SIMPLIFY=F)
      
      
      
      
      # 4b) sum up sense and antisense signals
      signal_theo_prompts_leaf_altogether <- lapply(wideMeta_signal_theo_prompts_leaf_altogether, function(x) rowSums(x$sense) + rowSums(x$anti))
      signal_theo_prompts_flower_altogether <- lapply(wideMeta_signal_theo_prompts_flower_altogether, function(x) rowSums(x$sense) + rowSums(x$anti))
      
      # 4c) save
      saveRDS(signal_theo_prompts_leaf_altogether, '~/masked_path/scatter_signal_theo_prompts_leaf_altogether.rds')
      saveRDS(signal_theo_prompts_flower_altogether, '~/masked_path/scatter_signal_theo_prompts_flower_altogether.rds')
}

signal_theo_prompts_leaf_altogether <- readRDS('~/masked_path/scatter_signal_theo_prompts_leaf_altogether.rds')
signal_theo_prompts_flower_altogether <- readRDS('~/masked_path/scatter_signal_theo_prompts_flower_altogether.rds')



# 5. LOG-LOG SCATTER PLOT ####
# ----------------------------
loglog_flower_df <- signal_theo_prompts_flower_altogether %>%
                    plyr::ldply(.id='sample') %>%
                    melt(id.vars='sample', variable.name='GeneID', value.name='counts') %>%
                    left_join(select(normFact_flower, Sample, total_clean_reads_normFact), by=c('sample'='Sample')) %>%
                    mutate('RPM'=counts/total_clean_reads_normFact) %>% # normalize readcounts into RPM
                    select(-counts, -total_clean_reads_normFact) %>%
                    separate(sample, c('genotype', 'replicate'), sep='_') %>%
                    group_by(genotype, GeneID) %>%
                    summarise('meanRPM'=mean(RPM)) %>% # average RPM accross replicates
                    ungroup() %>%
                    spread(genotype, meanRPM)

loglog_leaf_df <- signal_theo_prompts_leaf_altogether %>%
                  plyr::ldply(.id='sample') %>%
                  melt(id.vars='sample', variable.name='GeneID', value.name='counts') %>%
                  left_join(select(normFact_leaf, Sample, total_clean_reads_normFact), by=c('sample'='Sample')) %>%
                  mutate('RPM'=counts/total_clean_reads_normFact) %>% # normalize readcounts into RPM
                  select(-counts, -total_clean_reads_normFact) %>%
                  separate(sample, c('genotype', 'replicate'), sep='_') %>%
                  group_by(genotype, GeneID) %>%
                  summarise('meanRPM'=mean(RPM)) %>% # average RPM accross replicates
                  ungroup() %>%
                  spread(genotype, meanRPM)

loglog_df <- rbind(loglog_flower_df %>% mutate('tissue'='flower') %>% select(-hen2),
                   loglog_leaf_df %>% mutate('tissue'='leaf')) %>%
             mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))

# plot rrp4's (flower + leaves)
ggplot() +
  geom_point(data=subset(loglog_df, myPROMPTs==FALSE), aes(x=log(wt + 1), y=log(rrp4 + 1)), size=.2, col='black') +
  geom_point(data=subset(loglog_df, myPROMPTs==TRUE), aes(x=log(wt + 1), y=log(rrp4 + 1)), size=1.2, col='red') +
  geom_abline(slope=1, intercept=0, col='red') +
  facet_grid(tissue~.) +
  cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
  labs(title='Normalized small RNAs at theoretical PROMPT regions',
       x='wt log(RPM+1)', y='rrp4 log(RPM+1)')

# plot hen2 (flower only)
loglog_flower_df %<>% mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))
 
ggplot() +
  geom_point(data=subset(loglog_flower_df, myPROMPTs==FALSE), aes(x=wt + 1, y=hen2 + 1), size=.2, col='black') +
  geom_point(data=subset(loglog_flower_df, myPROMPTs==TRUE), aes(x=wt + 1, y=hen2 + 1), size=1.2, col='red') +
  scale_x_log10() + scale_y_log10() +
  geom_abline(slope=1, intercept=0, col='red') +
  cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
  labs(title='Normalized small RNAs at theoretical PROMPT regions',
       x='wt log(RPM+1)', y='rrp4 log(RPM+1)')



# 6. BOXPLOT VERSION ###
# ----------------------
# boxplot with all values
loglog_full_df <- rbind(loglog_flower_df %>%
                          mutate('tissue'='flower') %>%
                          melt(id.vars=c('GeneID', 'myPROMPTs', 'tissue'), variable.name='genotype', value.name='coverage'),
                        loglog_leaf_df %>%
                          mutate('tissue'='leaf', 'myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE)) %>%
                          melt(id.vars=c('GeneID', 'myPROMPTs', 'tissue'), variable.name='genotype', value.name='coverage'))

loglog_full_df %>%
  mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=myPROMPTs, y=coverage+1, fill=genotype)) +
         geom_boxplot() +
         scale_y_log10() +
         facet_grid(~tissue) +
         scale_fill_manual(values=c('aquamarine3', 'darkorchid3', 'brown3')) +
         cowplot::theme_cowplot() +
         labs(y='small RNA-seq (log coverage+1)',
              x='Region has a CAGE PROMPT TC',
              title='small RNA-seq at theoretical PROMPT regions',
              subtitle='-500 bp from annotated TAIR10 TSSs')


# bowplot with trimmed top 1%
loglog_flower_trimmed_df <- loglog_flower_df %>% mutate('hen2'=ifelse(hen2 > quantile(hen2, 0.99), NA, hen2),
                                                        'rrp4'=ifelse(rrp4 > quantile(rrp4, 0.99), NA, rrp4),
                                                        'wt'=ifelse(wt > quantile(wt, 0.99), NA, wt))

loglog_leaf_trimmed_df <- loglog_leaf_df %>% mutate('rrp4'=ifelse(rrp4 > quantile(rrp4, 0.99), NA, rrp4),
                                                    'wt'=ifelse(wt > quantile(wt, 0.99), NA, wt))

loglog_full_trimmed_df <- rbind(loglog_flower_trimmed_df %>%
                                      mutate('tissue'='flower') %>%
                                      melt(id.vars=c('GeneID', 'myPROMPTs', 'tissue'), variable.name='genotype', value.name='coverage'),
                                loglog_leaf_trimmed_df %>%
                                      mutate('tissue'='leaf', 'myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE)) %>%
                                      melt(id.vars=c('GeneID', 'myPROMPTs', 'tissue'), variable.name='genotype', value.name='coverage'))

loglog_full_trimmed_df %>%
  subset(!is.na(coverage)) %>%
  mutate('genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
  ggplot(aes(x=myPROMPTs, y=coverage+1, fill=genotype)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(~tissue) +
  scale_fill_manual(values=c('aquamarine3', 'darkorchid3', 'brown3')) +
  cowplot::theme_cowplot() +
  labs(y='small RNA-seq (log coverage+1, trimmed 0.99)',
       x='Region has a CAGE PROMPT TC',
       title='small RNA-seq at theoretical PROMPT regions',
       subtitle='-500 bp from annotated TAIR10 TSSs')

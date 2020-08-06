#### Arabidopsis TSS/exosome: boxplot of smRNA at CAGE TC PROMPT regions -vs- mRNA region by quantiles
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
library(ggrastr)
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

# cage data
cage_p <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.plus', full.names=T) %>% BigWigFileList()
cage_m <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.minus', full.names=T) %>% BigWigFileList()
cage_names <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.plus') %>%
              str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names

# small RNA-seq BIGWIGs (Dec 2019) : individual lengths - RAW COUNTS
  # leaf (rosette)
  sRNA_plus_counts_leaf_byLength <- list.files('~/masked_path/4.bigwigs_by_length', pattern='*Forward', full.names=T) %>% BigWigFileList()
  sRNA_minus_counts_leaf_byLength <- list.files('~/masked_path/4.bigwigs_by_length', pattern='*Reverse', full.names=T) %>% BigWigFileList()
  sRNA_names_counts_leaf_byLength <- list.files('~/masked_path/4.bigwigs_by_length', pattern='*Forward') %>%
                                     str_remove('.Forward.bw') %>% str_replace('\\.', '_')
  names(sRNA_plus_counts_leaf_byLength) <- names(sRNA_minus_counts_leaf_byLength) <- sRNA_names_counts_leaf_byLength
  # flower
  sRNA_plus_counts_flower_byLength <- list.files('~/masked_path/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Forward', full.names=T) %>% BigWigFileList()
  sRNA_minus_counts_flower_byLength <- list.files('~/masked_path/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Reverse', full.names=T) %>% BigWigFileList()
  sRNA_names_counts_flower_byLength <- list.files('~/masked_path/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Forward') %>%
                                     str_remove('.raw.Forward.bw') %>% str_replace('\\.', '_')
  names(sRNA_plus_counts_flower_byLength) <- names(sRNA_minus_counts_flower_byLength) <- sRNA_names_counts_flower_byLength
  
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



# 3. MAKE THEORETICAL PROMPT REGIONS HAVING A PROMPT CAGE TC: -500 bp from annotated TSSs ###
# -------------------------------------------------------------------------------------------
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

# 3d) keep only those with a CAGE PROMPT TC (N=60)
theo_prompts_with_PROMPT_TC <- subset(theo_prompts, gene_id %in% PROMPT_TC_genes)

# 3d) make 1-bp version
theo_prompts_with_PROMPT_TC_1bp <- resize(theo_prompts_with_PROMPT_TC, width=1, fix='start')



# 4. COMPUTE smRNA COVERAGE AT THEORETICAL PROMPT REGIONS ####
# ------------------------------------------------------------
# 4a) compute sum of signal (must be downstream, see above in 4)
wideMeta_signal_theo_prompts_leaf_byLength <- mapply(function(x, y) wideMetaProfile(sites=theo_prompts_with_PROMPT_TC_1bp,
                                                                                    forward=x, reverse=y,
                                                                                    downstream=500, upstream=1),
                                                     sRNA_plus_counts_leaf_byLength, sRNA_minus_counts_leaf_byLength, SIMPLIFY=F)

wideMeta_signal_theo_prompts_flower_byLength <- mapply(function(x, y) wideMetaProfile(sites=theo_prompts_with_PROMPT_TC_1bp,
                                                                                      forward=x, reverse=y,
                                                                                      downstream=500, upstream=1),
                                                       sRNA_plus_counts_flower_byLength, sRNA_minus_counts_flower_byLength, SIMPLIFY=F)

# 4b) sum up signal on sense and antisense strands
signal_theo_prompts_leaf_byLength <- lapply(wideMeta_signal_theo_prompts_leaf_byLength, function(x) rowSums(x$sense) + rowSums(x$anti))
signal_theo_prompts_flower_byLength <- lapply(wideMeta_signal_theo_prompts_flower_byLength, function(x) rowSums(x$sense) + rowSums(x$anti))

# 4c) normalize raw read count to number of inputted clean reads to STAR
    # match normalization factors with signal names
    match_vec_prompts_leaf <- match(str_remove(names(signal_theo_prompts_leaf_byLength), '_[0-9][0-9]$'), normFact_leaf$Sample)
    data.frame('sample'=names(signal_theo_prompts_leaf_byLength), 'norm'=normFact_leaf$Sample[match_vec_prompts_leaf])
    
    match_vec_prompts_flower <- match(str_remove(names(signal_theo_prompts_flower_byLength), '_[0-9][0-9]$'), normFact_flower$Sample)
    data.frame('sample'=names(signal_theo_prompts_flower_byLength), 'norm'=normFact_flower$Sample[match_vec_prompts_flower])

signal_theo_prompts_leaf_byLength <- mapply(function(x, y) x/y, signal_theo_prompts_leaf_byLength, normFact_leaf$total_clean_reads_normFact[match_vec_prompts_leaf])
signal_theo_prompts_flower_byLength <- mapply(function(x, y) x/y, signal_theo_prompts_flower_byLength, normFact_flower$total_clean_reads_normFact[match_vec_prompts_flower])

# 4d) compute average between replicates, by length
signal_theo_prompts_leaf_byLength %<>% as.data.frame() %>%
                                    rownames_to_column('GeneID') %>%
                                    melt(id.vars='GeneID') %>%
                                    separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                                    unite(sample, c(genotype, length), sep='_') %>%
                                    group_by(GeneID, sample) %>%
                                    summarise('meanValue'=mean(value)) %>%
                                    ungroup() %>%
                                    spread(key=sample, value=meanValue)

signal_theo_prompts_flower_byLength %<>% as.data.frame() %>%
                                      rownames_to_column('GeneID') %>%
                                      melt(id.vars='GeneID') %>%
                                      separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                                      unite(sample, c(genotype, length), sep='_') %>%
                                      group_by(GeneID, sample) %>%
                                      summarise('meanValue'=mean(value)) %>%
                                      ungroup() %>%
                                      spread(key=sample, value=meanValue)

# 4e) add which gene has a PROMPT TC
signal_theo_prompts_leaf_byLength %<>% mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))
signal_theo_prompts_flower_byLength %<>% mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))



# 5. MAKE mRNA REGION FOR GENES WITHOUT A CAGE PROMPT TC ####
# -----------------------------------------------------------
# 5a. get all annotated TSSs and extend them 500 bp downstream (N=33323)
mRNAs <- genes(txdb) %>%
         promoters(upstream=1, downstream=500) %>%
         remove_out_of_bound()
      # fix seqinfo
      seqlevelsStyle(mRNAs) <- seqlevelsStyle(myseqinfo)
      seqlevels(mRNAs) <- seqlevels(myseqinfo)
      seqinfo(mRNAs) <- myseqinfo
      seqlevels(mRNAs, pruning.mode='coarse') <- paste0('Chr', 1:5)

# 5b. remove those that have a PROMPT TC (N=33230)
mRNAs <- subset(mRNAs, gene_id %!in% PROMPT_TC_genes)

# 5c. remove those overlapping any other known gene (except themselves) (N=31289)
mRNAs_not_overlapping_anything_else <- countOverlaps(mRNAs, genes, ignore.strand=T) %>%
                                       enframe() %>%
                                       subset(value==1) %>%
                                       pull(name)

mRNAs <- subset(mRNAs, gene_id %in% mRNAs_not_overlapping_anything_else)

# 5d. export for IGV check
export.bed(mRNAs, 'mRNA_500bp_regions.bed')

# 5e. make 1-bp version
mRNAs_1bp <- resize(mRNAs, width=1, fix='start')



# 6. COMPUTE smRNA COVERAGE AT mRNA REGIONS ####
# ----------------------------------------------
# 6a) compute sum of signal (must be downstream)
if(FALSE){
      wideMeta_signal_mRNAs_leaf_byLength <- mapply(function(x, y) wideMetaProfile(sites=mRNAs_1bp,
                                                                                   forward=x, reverse=y,
                                                                                   downstream=500, upstream=1),
                                                    sRNA_plus_counts_leaf_byLength, sRNA_minus_counts_leaf_byLength, SIMPLIFY=F)
      
      wideMeta_signal_mRNAs_flower_byLength <- mapply(function(x, y) wideMetaProfile(sites=mRNAs_1bp,
                                                                                     forward=x, reverse=y,
                                                                                     downstream=500, upstream=1),
                                                      sRNA_plus_counts_flower_byLength, sRNA_minus_counts_flower_byLength, SIMPLIFY=F)
      
      saveRDS(wideMeta_signal_mRNAs_leaf_byLength, 'wideMeta_signal_mRNAs_leaf_byLength.rds')
      saveRDS(wideMeta_signal_mRNAs_flower_byLength, 'wideMeta_signal_mRNAs_flower_byLength.rds')
}

wideMeta_signal_mRNAs_leaf_byLength <- readRDS('wideMeta_signal_mRNAs_leaf_byLength.rds')
wideMeta_signal_mRNAs_flower_byLength <- readRDS('wideMeta_signal_mRNAs_flower_byLength.rds')

# 6b) sum up signal on sense and antisense strands
signal_mRNAs_leaf_byLength <- lapply(wideMeta_signal_mRNAs_leaf_byLength, function(x) rowSums(x$sense) + rowSums(x$anti))
signal_mRNAs_flower_byLength <- lapply(wideMeta_signal_mRNAs_flower_byLength, function(x) rowSums(x$sense) + rowSums(x$anti))

# 6c) normalize raw read count to number of inputted clean reads to STAR
      # match normalization factors with signal names
      match_vec_mRNA_leaf <- match(str_remove(names(signal_mRNAs_leaf_byLength), '_[0-9][0-9]$'), normFact_leaf$Sample)
      data.frame('sample'=names(signal_mRNAs_leaf_byLength), 'norm'=normFact_leaf$Sample[match_vec_mRNA_leaf])
      
      match_vec_mRNA_flower <- match(str_remove(names(signal_mRNAs_flower_byLength), '_[0-9][0-9]$'), normFact_flower$Sample)
      data.frame('sample'=names(signal_mRNAs_flower_byLength), 'norm'=normFact_flower$Sample[match_vec_mRNA_flower])
      
signal_mRNAs_leaf_byLength <- mapply(function(x, y) x/y, signal_mRNAs_leaf_byLength, normFact_leaf$total_clean_reads_normFact[match_vec_mRNA_leaf])
signal_mRNAs_flower_byLength <- mapply(function(x, y) x/y, signal_mRNAs_flower_byLength, normFact_flower$total_clean_reads_normFact[match_vec_mRNA_flower])

# 6d) compute average between replicates, by length
signal_mRNAs_leaf_byLength %<>% as.data.frame() %>%
                               rownames_to_column('GeneID') %>%
                               melt(id.vars='GeneID') %>%
                               separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                               unite(sample, c(genotype, length), sep='_') %>%
                               group_by(GeneID, sample) %>%
                               summarise('meanValue'=mean(value)) %>%
                               ungroup() %>%
                               spread(key=sample, value=meanValue)

signal_mRNAs_flower_byLength %<>% as.data.frame() %>%
                                  rownames_to_column('GeneID') %>%
                                  melt(id.vars='GeneID') %>%
                                  separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                                  unite(sample, c(genotype, length), sep='_') %>%
                                  group_by(GeneID, sample) %>%
                                  summarise('meanValue'=mean(value)) %>%
                                  ungroup() %>%
                                  spread(key=sample, value=meanValue)

# 6e) add which gene has a PROMPT TC (should be ZERO)
signal_mRNAs_leaf_byLength %<>% mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))
signal_mRNAs_leaf_byLength$myPROMPTs %>% table(useNA='always')

signal_mRNAs_flower_byLength %<>% mutate('myPROMPTs'=ifelse(GeneID %in% PROMPT_TC_genes, TRUE, FALSE))
signal_mRNAs_flower_byLength$myPROMPTs %>% table(useNA='always')



# 7. SPLIT mRNAs by THEIR EXPRESSION RANKS ####
# ---------------------------------------------
# 7a. get gene expression: from promoter TCs
# get TCs that are promoters
promoter_TCs_df <- rowRanges(TCs) %>%
                   subset(txType_TAIR10=='promoter') %>%
                   as.data.frame() %>%
                   select('TC_id'=thick.names, 'GeneID'=geneID) %>%
                   as_tibble()
# 7b. get TPM of TCs that are promoters
promoter_TC_meanTPM <- assay(TCs, 'TPM') %>%
                       as.data.frame() %>%
                       rownames_to_column('TC_id') %>%
                       subset(TC_id %in% promoter_TCs_df$TC_id) %>%
                       melt(id.vars='TC_id') %>%
                       separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
                       select(-timepoint) %>%
                       group_by(TC_id, genotype) %>%
                       summarise('meanTPM'=sum(value)) %>%
                       ungroup() %>%
                       spread(key=genotype, value=meanTPM)
# 7c. add gene ID to TC meanTPM
promoter_TC_meanTPM %<>% left_join(promoter_TCs_df, by='TC_id')

# 7d. compute sum of promoter TC TPM per gene (for those who have more than 1 promoter TC)
promoter_TC_meanTPM %<>% group_by(GeneID) %>%
                         summarise('hen2_promoter'=sum(hen2),
                                   'rrp4_promoter'=sum(rrp4),
                                   'wt_promoter'=sum(wt)) %>%
                         ungroup()

# 7e. keep only the genes that are in my mRNA region set
promoter_TC_meanTPM <- subset(promoter_TC_meanTPM, GeneID %in% mRNAs$gene_id)

# 7f. stratify genes by their expression: make 20% percentiles
# quantiles are calculated from the mean expressions across genotypes (wt + hen2 + rrp4)/3
qtiles <- seq(from=0, to=1, by=.2)
average_expression_across_genotypes <- (promoter_TC_meanTPM$hen2_promoter + promoter_TC_meanTPM$rrp4_promoter + promoter_TC_meanTPM$wt_promoter) / 3
unified_breaks_promoter <- quantile(average_expression_across_genotypes, qtiles)

promoter_TC_rank_unified <- data.frame('GeneID'=promoter_TC_meanTPM$GeneID,
                                       'rank'=cut(average_expression_across_genotypes, breaks=unified_breaks_promoter, labels=1:5, include.lowest=T)) %>%
                            as_tibble()

# 7g. put together everything
# - mRNA region coverage
# - PROMPT region coverage
# - ranking (override for PROMPT regions)

# flower
flower_data <- rbind(signal_mRNAs_flower_byLength,                    # mRNAs regions
                     signal_theo_prompts_flower_byLength) %>%         # PROMPT regions
               left_join(promoter_TC_rank_unified, by='GeneID') %>%   # mRNA ranking
               mutate('rank'=ifelse(myPROMPTs, 0, rank))              # override rank to zero if myPROMPTs


# leaf
leaf_data <- rbind(signal_mRNAs_leaf_byLength,
                   signal_theo_prompts_leaf_byLength) %>%
             left_join(promoter_TC_rank_unified, by='GeneID') %>%
             mutate('rank'=ifelse(myPROMPTs, 0, rank))

# 7h. how many of myPROMPTs are left when I remove the genes with no rank? (N=60)
flower_data %>% subset(!is.na(rank)) %>% pull(myPROMPTs) %>% table(useNA='always')
leaf_data %>% subset(!is.na(rank)) %>% pull(myPROMPTs) %>% table(useNA='always')




# 8. BOXPLOT OF smRNA COVERAGE AT THEORETICAL PROMPT REGIONS BY GENE EXPRESSION STRATA ####
# -----------------------------------------------------------------------------------------
# 8a. melt dataframes
leaf_data_melted <- leaf_data %>%
                    melt(id.vars=c('GeneID', 'myPROMPTs', 'rank'), value.name='coverage') %>%
                    separate(variable, c('genotype', 'length'), sep='_') %>%
                    mutate('myPROMPTs'=ifelse(myPROMPTs, 'PROMPT\nregion', 'mRNA region (TSS + 500 bp)'),
                           'genotype'=factor(genotype, levels=c('wt', 'rrp4'))) %>%
                    subset(!is.na(rank)) %>%
                    as_tibble()

flower_data_melted <- flower_data %>%
                      melt(id.vars=c('GeneID', 'myPROMPTs', 'rank'), value.name='coverage') %>%
                      separate(variable, c('genotype', 'length'), sep='_') %>%
                      mutate('myPROMPTs'=ifelse(myPROMPTs, 'PROMPT\nregion', 'mRNA region (TSS + 500 bp)'),
                             'genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
                      subset(!is.na(rank)) %>%
                      as_tibble()

# 8b. plot
    pseudo <- 1
    # percentile values of non-PROMPT-host genes
    unified_breaks_promoter
    # average expression of PROMPT TCs
    subset(TCs, names %in% PROMPT_TC_ids) %>%
      assay('TPM') %>%
      as.data.frame() %>%
      rownames_to_column('TC_id') %>%
      melt(id.vars='TC_id') %>%
      separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
      select(-timepoint) %>%
      group_by(TC_id, genotype) %>%
      summarise('meanTPM'=mean(value)) %>%
      ungroup() %>%
      group_by(genotype) %>%
      summarise('PROMPT_TPM_mean'=mean(meanTPM))


leaf_data_melted %>%
  ggplot(aes(x=as.factor(rank), y=coverage + pseudo, fill=genotype)) +
         geom_boxplot(alpha=.75, outlier.size=.5) +
         scale_y_log10() +
         cowplot::theme_cowplot() +
         facet_grid(length ~ myPROMPTs, scales='free_x', space='free_x') +
         scale_fill_manual(values=c('aquamarine3', 'darkorchid4', 'brown3')) +
         labs(title='small RNA coverage - LEAF',
              y='Normalized small RNA read count (log RPM+1)',
              x='mRNA expression percentiles (by 20%)\nfrom promoter CAGE TCs')

flower_data_melted %>%
  ggplot(aes(x=as.factor(rank), y=coverage + pseudo, fill=genotype)) +
  geom_boxplot(alpha=.75, outlier.size=.5) +
  scale_y_log10() +
  cowplot::theme_cowplot() +
  facet_grid(length ~ myPROMPTs, scales='free_x', space='free_x') +
  scale_fill_manual(values=c('aquamarine3', 'darkorchid4', 'brown3')) +
  labs(title='small RNA coverage - FLOWER',
       y='Normalized small RNA read count (log RPM+1)',
       x='mRNA expression percentiles (by 20%)\nfrom promoter CAGE TCs')




# 7. BOXPLOT OR RATIO: between PROMPT TCs and non-PROMPT TC regions ####
# ----------------------------------------------------------------------
# 7a) compute ratios to wt
ratio_df <- total_signal_theo_prompts %>%
            melt(id.Vars=c('GeneID', 'myPROMPTs')) %>%
            separate(variable, c('genotype', 'length'), sep='_') %>%
            spread(key=genotype, value=value) %>%
            mutate('hen2_ratio' = (hen2 + 1) / (wt + 1),
                   'rrp4_ratio' = (rrp4 + 1) / (wt + 1)) %>%
            select(-hen2, -rrp4, -wt) %>%
            melt(id.vars=c('GeneID', 'myPROMPTs', 'length')) %>%
            as_tibble()

# 7b) plot boxplot
ratio_df

ggplot(ratio_df, aes(x=variable, y=value, fill=myPROMPTs)) +
       geom_boxplot() +
       scale_y_log10() +
       facet_grid(.~length) +
       cowplot::theme_cowplot() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
       labs(y='Normalized read count ratio log(mut+1/wt+1)',
            subtitle=paste0('N(regions with PROMPT TC) = ', sum(ratio_df$myPROMPTs==TRUE)/5/2, '\nN(regions without PROMPT TC) = ', sum(ratio_df$myPROMPTs==FALSE)/5/2), 
            title='Ratio of small RNAs at theoretical PROMPT regions',
            x='Compared mutant')

# 7c) find best cases with highest ratio
subset(ratio_df, length==21) %>%
        arrange(desc(value)) %>%
        filter(!duplicated(GeneID))



# 8. +/- 750bp PROFILE AROUND PROMPT TCs ####
# -------------------------------------------
# 8a) get PROMPT TC peaks (1bp)
PROMPT_TCs_1bp <- rowRanges(TCs) %>%
                  subset(names %in% names(PROMPT_TCs)) %>%
                  swapRanges()

# export PROMPT TC extended for IGV sanity check
export.bed(promoters(PROMPT_TCs_1bp, upstream=750, downstream=750), 'IGV sanity check - extended PROMPT TCs 750bp.bed')

# 8b) compute signal: /!\ invert PROMPT TC strand so that we look in the same direction as mRNA TSS
PROMPT_TC_profile <- tidyMetaProfile(sites=invertStrand(PROMPT_TCs_1bp),
                                     forward=sRNA_plus_counts_byLength, reverse=sRNA_minus_counts_byLength,
                                     upstream=750, downstream=750, trimLower=0.05, trimUpper=0.95)

# 8c) plot
PROMPT_TC_profile %>%
  select(-pos1) %>%
  separate(signal, c('genotype', 'replicate', 'length'), sep='_') %>%
  unite('sample', c('genotype', 'replicate'), sep='_') %>%
  left_join(select(normFact, Sample, total_clean_reads_normFact), by=c('sample'='Sample')) %>%
  mutate('sense'= sense / total_clean_reads_normFact,
          'anti'= anti / total_clean_reads_normFact) %>%
  select(-total_clean_reads_normFact) %>%
  separate(sample, c('genotype', 'replicate'), by='_') %>%
  group_by(genotype, length, pos0) %>%
  summarise('sense'=mean(sense),
            'anti'=mean(anti)) %>%
  ungroup() %>%
  mutate(anti=-anti) %>%
  gather(key='direction', value='score', sense, anti, factor_key=T) %>%
  ggplot(aes(x=pos0, y=score, col=direction)) +
         geom_line() +
         geom_hline(yintercept=0) +
         geom_vline(xintercept=0, lty=2) +
         facet_grid(genotype~length) +
         cowplot::theme_cowplot() +
         labs(title='Profile plot at PROMPT TCs',
              x='Position relative to PROMPT TC peak (N=96) (bp)',
              y='Normalized small RNA read counts (0.05-0.95% trimming)') +
         scale_color_brewer(palette='Set1', direction=-1, name='Direction')



# 9. PROFILE AROUND TSS OF GENES HOSTING A PROMPT TC ####
# -------------------------------------------------------
# 9a) get all annotated TSSs (1bp)
Gene_TSSs_1bp <- genes(txdb) %>%
                 resize(width=1, fix='start')
      # fix seqinfo
      seqlevels(Gene_TSSs_1bp) <- seqlevels(myseqinfo)
      seqinfo(Gene_TSSs_1bp) <- myseqinfo

# 9b) keep only those with PROMPT TC (N=94)
Gene_TSSs_1bp_with_PROMPT_TC <- subset(Gene_TSSs_1bp, gene_id %in% PROMPT_TC_genes)

# 9c) compute smRNA signal
profile_Gene_TSSs_1bp_with_PROMPT_TC <- tidyMetaProfile(sites=Gene_TSSs_1bp_with_PROMPT_TC,
                                                        forward=sRNA_plus_counts_byLength, reverse=sRNA_minus_counts_byLength,
                                                        upstream=1000, downstream=1000,
                                                        trimLower=0.05, trimUpper=0.95)

# 9d) plot
profile_Gene_TSSs_1bp_with_PROMPT_TC %>%
    select(-pos1) %>%
    separate(signal, c('genotype', 'replicate', 'length'), sep='_') %>%
    unite('sample', c('genotype', 'replicate'), sep='_') %>%
    left_join(select(normFact, Sample, total_clean_reads_normFact), by=c('sample'='Sample')) %>%
    mutate('sense'= sense / total_clean_reads_normFact,
           'anti'= anti / total_clean_reads_normFact) %>%
    select(-total_clean_reads_normFact) %>%
    separate(sample, c('genotype', 'replicate'), by='_') %>%
    group_by(genotype, length, pos0) %>%
    summarise('sense'=mean(sense),
              'anti'=mean(anti)) %>%
    ungroup() %>%
    mutate('anti'=-anti) %>%
    gather(key='direction', value='score', sense, anti, factor_key=T) %>%
    ggplot(aes(x=pos0, y=score, col=direction)) +
           geom_rect(aes(xmin=-400, xmax=0, ymin=-Inf, ymax=Inf), fill='grey80', col=NA) +
           geom_hline(yintercept=0) +
           geom_vline(xintercept=0, lty=2) +
           geom_line() +
           facet_grid(genotype~length) +
           cowplot::theme_cowplot() +
           labs(title='Profile at TSSs of genes hosting a PROMPT TC (N=94)',
                y='Normalized small RNA read count (0.05-0.95%)',
                x='TSSs of gene hosting a PROMPT TC (bp)') +
           scale_color_brewer(palette='Set1', direction=-1, name='Direction')



# 10. PROMPT -vs- mRNA TSS ####
# ----------------------------
# ----------------------------

#            PROMPT TC       mRNA TSS
#         < - - - - - |      | - - - - >
#                     |      |=> => => => => => => =>
#    [  500 bp max ]      [ 500 bp ]
#      
#      sum(smRNAseq)/sum(CAGE)      sum(smRNAseq)/sum(CAGE)

# PART 1: PROMPT SIDE
# -------------------
# 10.1a) get PROMPT TCs extended 500 bp downstream
PROMPTs_down500 <- PROMPT_TCs_1bp %>%
                   promoters(upstream=1, downstream=500)
    # IGV check
    export.bed(PROMPTs_down500, 'prompts_downstream_500.bed')

# 10.1b) keep those NOT overlapping a known gene on the same strand as the PROMPT (N=94)
PROMPTs_down500_noOverlap_SameStrand <- subsetByOverlaps(PROMPTs_down500, genes, invert=T)

# 10.1c) keep those NOT overlapping themselves
tmp <- countOverlaps(PROMPTs_down500_noOverlap_SameStrand, ignore.strand=T) %>%
  enframe() %>%
  subset(value==1) %>%
  pull(name)

PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand <- PROMPTs_down500_noOverlap_SameStrand[names(PROMPTs_down500_noOverlap_SameStrand) %in% tmp]
rm(tmp)

# 10.1d) Get average TPM CAGE signal for those
PROMPT_TC_CAGE_meanTPM <- assay(TCs, 'TPM') %>%
  subset(rownames(.) %in% names(PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand)) %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  melt(id.vars='TC_id') %>%
  separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
  group_by(TC_id, genotype, timepoint) %>%
  summarise('meanTPM'=mean(value)) %>%
  ungroup() %>%
  unite('sample', c('genotype', 'timepoint'), sep='_') %>%
  spread(key=sample, value=meanTPM)
  
# 10.1e) Compute total small RNA signal
PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand_1bp <- resize(PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand, width=1, fix='start')

smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand <- mapply(function(x, y) wideMetaProfile(sites=PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand_1bp,
                                                                                                        forward=x, reverse=y,
                                                                                                        downstream=500, upstream=1),
                                                                         sRNA_plus_counts_byLength, sRNA_minus_counts_byLength,
                                                                         SIMPLIFY=F)

smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand %<>% lapply(function(sample) rowSums(sample$sense) + rowSums(sample$anti))

# 10.1f) normalize raw read count to number of inputted clean reads to STAR
smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand <- mapply(function(x, y) x/y, smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand, rep(normFact$total_clean_reads_normFact, each=5))

# 10.1g) compute average between replicates, by length
smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand %<>% as.data.frame() %>%
                                                                    rownames_to_column('TC_id') %>%
                                                                    melt(id.vars='TC_id') %>%
                                                                    separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                                                                    unite(sample, c(genotype, length), sep='_') %>%
                                                                    group_by(TC_id, sample) %>%
                                                                    summarise('meanValue'=mean(value)) %>%
                                                                    ungroup() %>%
                                                                    spread(key=sample, value=meanValue)

# 10.1h). normalize smRNA to total CAGE expression (rrp4, t=0 only)
prompt_500bp_all_signal <- left_join(smRNA_PROMPTs_down500_noOverlap_geneSameStrand_PROMPTanyStrand,
                                     PROMPT_TC_CAGE_meanTPM %>% rename('cage_hen2'='hen2_0', 'cage_rrp4'='rrp4_0', 'cage_wt'='wt_0'),
                                     by='TC_id')
pseudocount <- 1

rrp4_smRNA_expNorm <- with(prompt_500bp_all_signal,
                           data.frame(TC_id,
                                      'rrp4_21_PROMPT_smRNA_normExp'=(rrp4_21 + pseudocount)/(cage_rrp4 + pseudocount),
                                      'rrp4_22_PROMPT_smRNA_normExp'=(rrp4_22 + pseudocount)/(cage_rrp4 + pseudocount),
                                      'rrp4_23_PROMPT_smRNA_normExp'=(rrp4_23 + pseudocount)/(cage_rrp4 + pseudocount),
                                      'rrp4_24_PROMPT_smRNA_normExp'=(rrp4_24 + pseudocount)/(cage_rrp4 + pseudocount)))

                                                                         

# PART 2: GENE SIDE
# -----------------
# 10.2a) get 500 bp dowstream from annotated TSS (into gene body)
genes_500 <- genes(txdb) %>%
             promoters(upstream=1, downstream=500) %>%
             remove_out_of_bound()
      # fix seqinfo
      seqlevelsStyle(genes_500) <- seqlevelsStyle(myseqinfo)
      seqlevels(genes_500) <- seqlevels(myseqinfo)

# 10.2b) keep only those that have a PROMPT TC
genes_500 %<>% subset(gene_id %in% PROMPT_TC_genes)

# 10.2c) keep those not overlapping each other (there are none)
tmp <- countOverlaps(genes_500, ignore.strand=T) %>%
       enframe() %>%
       subset(value==1) %>%
       pull(name)

genes_500 %<>% subset(gene_id %in% tmp)
rm(tmp)

# 10.3c) keep those not overlapping any other feature
tmp <- countOverlaps(genes_500, genes, ignore.strand=T) %>%
  enframe() %>%
  subset(value==1) %>%
  pull(name)

genes_500 %<>% subset(gene_id %in% tmp)
rm(tmp)

# 10.4c) export for sanity check
export.bed(genes_500, 'genes_500bp.bed')

# 10.4d) compute smRNA signal
genes_500_1bp <- resize(genes_500, width=1, fix='start')

smRNA_genes_500 <- mapply(function(x, y) wideMetaProfile(sites=genes_500_1bp,
                                                         forward=x, reverse=y,
                                                         downstream=500, upstream=1),
                          sRNA_plus_counts_byLength, sRNA_minus_counts_byLength,
                          SIMPLIFY=F)

smRNA_genes_500 %<>% lapply(function(sample) rowSums(sample$sense) + rowSums(sample$anti))

# 10.4e) normalize raw read count to number of inputted clean reads to STAR
smRNA_genes_500 <- mapply(function(x, y) x/y, smRNA_genes_500, rep(normFact$total_clean_reads_normFact, each=5))

# 10.4f) compute average between replicates, by length
smRNA_genes_500 %<>% as.data.frame() %>%
                     rownames_to_column('TC_id') %>%
                     melt(id.vars='TC_id') %>%
                     separate(variable, c('genotype', 'replicate', 'length'), sep='_') %>%
                     unite(sample, c(genotype, length), sep='_') %>%
                     group_by(TC_id, sample) %>%
                     summarise('meanValue'=mean(value)) %>%
                     ungroup() %>%
                     spread(key=sample, value=meanValue)

# 10.4g) get CAGE expression on the 500 bps
cage_genes_500 <- mapply(function(x, y) wideMetaProfile(sites=genes_500_1bp,
                                                        forward=x, reverse=y,
                                                        downstream=500, upstream=1),
                         cage_p, cage_m,
                         SIMPLIFY=F)

cage_genes_500 %<>% sapply(function(sample) rowSums(sample$sense) + rowSums(sample$anti)) %>% as.data.frame() %>% rownames_to_column('gene_ID')
cage_genes_500 %<>% rename('cage_hen2'='hen2', 'cage_rrp4'='rrp4', 'cage_wt'='wt')

# 10.4h) merge smRNA and CAGE at genes 500
smRNA_genes_500 <- left_join(smRNA_genes_500, cage_genes_500, by=c('TC_id'='gene_ID'))

# 10.4i) normalize smRNA to CAGE
pseudocount <- 1

rrp4_smRNA_expNorm_geneside <- with(smRNA_genes_500,
                                    data.frame(TC_id,
                                               'rrp4_21_gene_smRNA_normExp'=(rrp4_21 + pseudocount)/(cage_rrp4 + pseudocount),
                                               'rrp4_22_gene_smRNA_normExp'=(rrp4_22 + pseudocount)/(cage_rrp4 + pseudocount),
                                               'rrp4_23_gene_smRNA_normExp'=(rrp4_23 + pseudocount)/(cage_rrp4 + pseudocount),
                                               'rrp4_24_gene_smRNA_normExp'=(rrp4_24 + pseudocount)/(cage_rrp4 + pseudocount)))


# PART 3: PUT BOTH SIDES TOGETHER AND COMPUTE RATIO
# -------------------------------------------------
# PROMPT-side
rrp4_smRNA_expNorm %<>% as_tibble()
# GENE-side
rrp4_smRNA_expNorm_geneside %<>% as_tibble()

# 10.3a) make PROMPT TC / GeneID mapping table
prompt_TC_gene_ID_mapping <- PROMPT_TCs %>%
                             as.data.frame() %>%
                             select(thick.names, geneID_anti) %>%
                             set_rownames(NULL) %>%
                             as_tibble()
prompt_TC_gene_ID_mapping

# 10.3b) merge both using the GENE/TC mapping table
rrp4_smRNA_expNorm %<>% left_join(prompt_TC_gene_ID_mapping, by=c('TC_id'='thick.names'))

final <- left_join(rrp4_smRNA_expNorm, rrp4_smRNA_expNorm_geneside, by=c('geneID_anti'='TC_id'))

# 10.3c) plot scatter plot
final %>%
  select(-geneID_anti) %>%
  melt(id.vars='TC_id') %>%
  separate(variable, c('rrp4', 'length', 'side', 'bla', 'blu'), sep='_') %>%
  select(-rrp4, -bla, -blu) %>%
  spread(key=side, value=value) %>%
  ggplot(aes(x=PROMPT, y=gene), col=length) +
         geom_abline(intercept=0, slope=1, lty=2, lwd=.3, col='black') +
         geom_point() +
         scale_x_log10() + scale_y_log10() +
         facet_wrap(~length, ncol=2, nrow=2) +
         cowplot::theme_cowplot() + cowplot::panel_border(colour='black') + cowplot::background_grid(colour.major='grey80') +
         scale_color_brewer(palette='Set1', name='read\nlength (bp)') +
         labs(x='PROMPT small RNAs normalized to PROMPT expression (log)',
              y='GENE small RNAs normalized to GENE expression (log)',
              title='Data for RRP4')

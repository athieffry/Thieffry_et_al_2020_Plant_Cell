#### Arabidopsis TSS/exosome: boxplot of smRNA at CAGE TC PROMPT regions -vs- mRNA region by quantiles (RATIO)
#### Axel Thieffry - November 2019
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

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/10 - smRNA at PROMPTs/together')



# 1. LOAD ALL INPUT FILES ####
# ----------------------------
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
idmapping <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSES_v2/00 - RDATA/AGI_ID_mapping.rds')
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_TPM1_min3lib_TSSstory.rds')
seqlevels(TCs) <- c(paste0('Chr', 1:5), 'ChrC', 'ChrM')

# cage data
cage_p <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.plus', full.names=T) %>% BigWigFileList()
cage_m <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.minus', full.names=T) %>% BigWigFileList()
cage_names <- list.files('~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123', pattern='*_0_R123.plus') %>%
              str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names

# small RNA-seq BIGWIGs (Dec 2019) : individual lengths - RAW COUNTS
  # leaf (rosette)
  sRNA_plus_counts_leaf_byLength <- list.files('~/Dropbox/Maria_smallRNA_re-analysis_Dec2019/4.bigwigs_by_length', pattern='*Forward', full.names=T) %>% BigWigFileList()
  sRNA_minus_counts_leaf_byLength <- list.files('~/Dropbox/Maria_smallRNA_re-analysis_Dec2019/4.bigwigs_by_length', pattern='*Reverse', full.names=T) %>% BigWigFileList()
  sRNA_names_counts_leaf_byLength <- list.files('~/Dropbox/Maria_smallRNA_re-analysis_Dec2019/4.bigwigs_by_length', pattern='*Forward') %>%
                                     str_remove('.Forward.bw') %>% str_replace('\\.', '_')
  names(sRNA_plus_counts_leaf_byLength) <- names(sRNA_minus_counts_leaf_byLength) <- sRNA_names_counts_leaf_byLength
  # flower
  sRNA_plus_counts_flower_byLength <- list.files('~/Dropbox/Peter_exosome_re-analysis/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Forward', full.names=T) %>% BigWigFileList()
  sRNA_minus_counts_flower_byLength <- list.files('~/Dropbox/Peter_exosome_re-analysis/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Reverse', full.names=T) %>% BigWigFileList()
  sRNA_names_counts_flower_byLength <- list.files('~/Dropbox/Peter_exosome_re-analysis/3.bigWigs_rawcounts_splitted_by_lengths', pattern='*Forward') %>%
                                     str_remove('.raw.Forward.bw') %>% str_replace('\\.', '_')
  names(sRNA_plus_counts_flower_byLength) <- names(sRNA_minus_counts_flower_byLength) <- sRNA_names_counts_flower_byLength
  
# small RNA-seq stats
  # leaf (rosette)
  stats_leaf <- readxl::read_xlsx('~/Dropbox/Maria_smallRNA_re-analysis_Dec2019/Sample_info.xlsx', col_names=T, trim_ws=T) %>%
                as_tibble() %>%
                select(-Raw_read_length)
  # flower
  stats_flower <- readxl::read_xlsx('~/Dropbox/Peter_exosome_re-analysis/Sample_info.xlsx', col_names=T, trim_ws=T) %>%
                  as_tibble() %>%
                  select(-Raw_read_length)

# get genes hosting a PROMPT CAGE TC
    # get PROMPT TC GR (N=96)
    PROMPT_TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds') %>%
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

#                              TSS
#         PROMPT TC           |
#             |                => => => => => => => GENE BODY => => => => => =>
#     -> -> -> -> -> -> -> ->
#    |     PROMPT region (same orientation as gene)
#    1bp position

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



# 7. GET EXPRESSION FROM PROMOTER TCs ####
# ----------------------------------------
# 7a. get TCs that are promoters
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

# 7f. put together smRNA and expression to compute ratio
# - mRNA region coverage & mRNA promoter expression
pseudo <- 1

ratio_flower_mRNA <- left_join(signal_mRNAs_flower_byLength, promoter_TC_meanTPM, by='GeneID') %>%
                     mutate('hen2_21_ratio' = (hen2_21 + pseudo) / (hen2_promoter + pseudo),
                            'hen2_22_ratio' = (hen2_22 + pseudo) / (hen2_promoter + pseudo),
                            'hen2_23_ratio' = (hen2_23 + pseudo) / (hen2_promoter + pseudo),
                            'hen2_24_ratio' = (hen2_24 + pseudo) / (hen2_promoter + pseudo),
                            'rrp4_21_ratio' = (rrp4_21 + pseudo) / (rrp4_promoter + pseudo),
                            'rrp4_22_ratio' = (rrp4_22 + pseudo) / (rrp4_promoter + pseudo),
                            'rrp4_23_ratio' = (rrp4_23 + pseudo) / (rrp4_promoter + pseudo),
                            'rrp4_24_ratio' = (rrp4_24 + pseudo) / (rrp4_promoter + pseudo),
                            'wt_21_ratio' = (wt_21 + pseudo) / (wt_promoter + pseudo),
                            'wt_22_ratio' = (wt_22 + pseudo) / (wt_promoter + pseudo),
                            'wt_23_ratio' = (wt_23 + pseudo) / (wt_promoter + pseudo),
                            'wt_24_ratio' = (wt_24 + pseudo) / (wt_promoter + pseudo)) %>%
                      select(GeneID, myPROMPTs, matches('ratio'))

ratio_leaf_mRNA <- left_join(signal_mRNAs_leaf_byLength, promoter_TC_meanTPM, by='GeneID') %>%
                   mutate('rrp4_21_ratio' = (rrp4_21 + pseudo) / (rrp4_promoter + pseudo),
                          'rrp4_22_ratio' = (rrp4_22 + pseudo) / (rrp4_promoter + pseudo),
                          'rrp4_23_ratio' = (rrp4_23 + pseudo) / (rrp4_promoter + pseudo),
                          'rrp4_24_ratio' = (rrp4_24 + pseudo) / (rrp4_promoter + pseudo),
                          'wt_21_ratio' = (wt_21 + pseudo) / (wt_promoter + pseudo),
                          'wt_22_ratio' = (wt_22 + pseudo) / (wt_promoter + pseudo),
                          'wt_23_ratio' = (wt_23 + pseudo) / (wt_promoter + pseudo),
                          'wt_24_ratio' = (wt_24 + pseudo) / (wt_promoter + pseudo)) %>%
                  select(GeneID, myPROMPTs, matches('ratio'))

# - PROMPT region coverage & PROMPT TC expression
    # get PROMPT TC expression first
    TC_id_gene_id <- rowRanges(TCs) %>%
                     as.data.frame() %>%
                     select('TC_id'=thick.names, 'GeneID'=geneID_anti) %>%
                     as_tibble()

    PROMPT_expression <- subset(TCs, names %in% PROMPT_TC_ids) %>%
      assay('TPM') %>%
      as.data.frame() %>%
      rownames_to_column('TC_id') %>%
      melt(id.vars='TC_id') %>%
      separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
      select(-timepoint) %>%
      group_by(TC_id, genotype) %>%
      summarise('meanTPM'=mean(value)) %>%
      ungroup() %>%
      spread(genotype, meanTPM) %>%
      left_join(TC_id_gene_id, by='TC_id')

ratio_flower_prompt <- left_join(signal_theo_prompts_flower_byLength, PROMPT_expression, by='GeneID') %>%
                       mutate('hen2_21_ratio' = (hen2_21 + pseudo) / (hen2 + pseudo),
                              'hen2_22_ratio' = (hen2_22 + pseudo) / (hen2 + pseudo),
                              'hen2_23_ratio' = (hen2_23 + pseudo) / (hen2 + pseudo),
                              'hen2_24_ratio' = (hen2_24 + pseudo) / (hen2 + pseudo),
                              'rrp4_21_ratio' = (rrp4_21 + pseudo) / (hen2 + pseudo),
                              'rrp4_22_ratio' = (rrp4_22 + pseudo) / (hen2 + pseudo),
                              'rrp4_23_ratio' = (rrp4_23 + pseudo) / (hen2 + pseudo),
                              'rrp4_24_ratio' = (rrp4_24 + pseudo) / (hen2 + pseudo),
                              'wt_21_ratio' = (wt_21 + pseudo) / (hen2 + pseudo),
                              'wt_22_ratio' = (wt_22 + pseudo) / (hen2 + pseudo),
                              'wt_23_ratio' = (wt_23 + pseudo) / (hen2 + pseudo),
                              'wt_24_ratio' = (wt_24 + pseudo) / (hen2 + pseudo)) %>%
                       select(GeneID, myPROMPTs, matches('ratio'))

ratio_leaf_prompt <- left_join(signal_theo_prompts_leaf_byLength, PROMPT_expression, by='GeneID') %>%
                     mutate('rrp4_21_ratio' = (rrp4_21 + pseudo) / (hen2 + pseudo),
                            'rrp4_22_ratio' = (rrp4_22 + pseudo) / (hen2 + pseudo),
                            'rrp4_23_ratio' = (rrp4_23 + pseudo) / (hen2 + pseudo),
                            'rrp4_24_ratio' = (rrp4_24 + pseudo) / (hen2 + pseudo),
                            'wt_21_ratio' = (wt_21 + pseudo) / (hen2 + pseudo),
                            'wt_22_ratio' = (wt_22 + pseudo) / (hen2 + pseudo),
                            'wt_23_ratio' = (wt_23 + pseudo) / (hen2 + pseudo),
                            'wt_24_ratio' = (wt_24 + pseudo) / (hen2 + pseudo)) %>%
                     select(GeneID, myPROMPTs, matches('ratio'))


# 7g. put mRNA and PROMPT together
ratio_flower <- rbind(ratio_flower_mRNA, ratio_flower_prompt)
ratio_leaf <- rbind(ratio_leaf_mRNA, ratio_leaf_prompt)



# 8. BOXPLOT OF smRNA COVERAGE AT THEORETICAL PROMPT REGIONS BY GENE EXPRESSION STRATA ####
# -----------------------------------------------------------------------------------------
# 8a. trim ratio to 0.95 %
    # trimming threshold
    tt <- 0.99

ratio_flower %<>% mutate('hen2_21_ratio'=ifelse(hen2_21_ratio > quantile(hen2_21_ratio, na.rm=T, tt), NA, hen2_21_ratio),
                         'hen2_22_ratio'=ifelse(hen2_22_ratio > quantile(hen2_22_ratio, na.rm=T, tt), NA, hen2_22_ratio),
                         'hen2_23_ratio'=ifelse(hen2_23_ratio > quantile(hen2_23_ratio, na.rm=T, tt), NA, hen2_23_ratio),
                         'hen2_24_ratio'=ifelse(hen2_24_ratio > quantile(hen2_24_ratio, na.rm=T, tt), NA, hen2_24_ratio),
                         'rrp4_21_ratio'=ifelse(rrp4_21_ratio > quantile(rrp4_21_ratio, na.rm=T, tt), NA, rrp4_21_ratio),
                         'rrp4_22_ratio'=ifelse(rrp4_22_ratio > quantile(rrp4_22_ratio, na.rm=T, tt), NA, rrp4_22_ratio),
                         'rrp4_23_ratio'=ifelse(rrp4_23_ratio > quantile(rrp4_23_ratio, na.rm=T, tt), NA, rrp4_23_ratio),
                         'rrp4_24_ratio'=ifelse(rrp4_24_ratio > quantile(rrp4_24_ratio, na.rm=T, tt), NA, rrp4_24_ratio),
                         'wt_21_ratio'=ifelse(wt_21_ratio > quantile(wt_21_ratio, na.rm=T, tt), NA, wt_21_ratio),
                         'wt_22_ratio'=ifelse(wt_22_ratio > quantile(wt_22_ratio, na.rm=T, tt), NA, wt_22_ratio),
                         'wt_23_ratio'=ifelse(wt_23_ratio > quantile(wt_23_ratio, na.rm=T, tt), NA, wt_23_ratio),
                         'wt_24_ratio'=ifelse(wt_24_ratio > quantile(wt_24_ratio, na.rm=T, tt), NA, wt_24_ratio))

ratio_leaf %<>% mutate('rrp4_21_ratio'=ifelse(rrp4_21_ratio > quantile(rrp4_21_ratio, na.rm=T, tt), NA, rrp4_21_ratio),
                       'rrp4_22_ratio'=ifelse(rrp4_22_ratio > quantile(rrp4_22_ratio, na.rm=T, tt), NA, rrp4_22_ratio),
                       'rrp4_23_ratio'=ifelse(rrp4_23_ratio > quantile(rrp4_23_ratio, na.rm=T, tt), NA, rrp4_23_ratio),
                       'rrp4_24_ratio'=ifelse(rrp4_24_ratio > quantile(rrp4_24_ratio, na.rm=T, tt), NA, rrp4_24_ratio),
                       'wt_21_ratio'=ifelse(wt_21_ratio > quantile(wt_21_ratio, na.rm=T, tt), NA, wt_21_ratio),
                       'wt_22_ratio'=ifelse(wt_22_ratio > quantile(wt_22_ratio, na.rm=T, tt), NA, wt_22_ratio),
                       'wt_23_ratio'=ifelse(wt_23_ratio > quantile(wt_23_ratio, na.rm=T, tt), NA, wt_23_ratio),
                       'wt_24_ratio'=ifelse(wt_24_ratio > quantile(wt_24_ratio, na.rm=T, tt), NA, wt_24_ratio))


# 8b. melt dataframes
ratio_flower_melted <- ratio_flower %>%
                       melt(id.vars=c('GeneID', 'myPROMPTs'), value.name='ratio') %>%
                       separate(variable, c('genotype', 'length'), sep='_', extra='drop') %>%
                       mutate('myPROMPTs'=ifelse(myPROMPTs, 'PROMPT\nregion', 'mRNA region\n(TSS + 500 bp)'),
                              'genotype'=factor(genotype, levels=c('wt', 'hen2', 'rrp4'))) %>%
                       as_tibble()


ratio_leaf_melted <- ratio_leaf %>%
                     melt(id.vars=c('GeneID', 'myPROMPTs'), value.name='ratio') %>%
                     separate(variable, c('genotype', 'length'), sep='_', extra='drop') %>%
                     mutate('myPROMPTs'=ifelse(myPROMPTs, 'PROMPT\nregion', 'mRNA region\n(TSS + 500 bp)'),
                            'genotype'=factor(genotype, levels=c('wt', 'rrp4'))) %>%
                     as_tibble()

# 8b. plot
ratio_flower_melted %>%
  subset(!is.na(ratio)) %>%
  ggplot(aes(x=myPROMPTs, y=ratio+1, fill=genotype)) +
         geom_boxplot(alpha=.75, outlier.size=.5) +
         scale_y_log10() +
         cowplot::theme_cowplot() +
         facet_grid(length~., scales='free_x', space='free_x') +
         scale_fill_manual(values=c('aquamarine3', 'darkorchid4', 'brown3')) +
         labs(title='Ratio: small RNA / expression - FLOWER',
              y='Ratio of small RNA versus PROMPT TC expression or mRNA expression\nlog-scale, pseudocount +1',
              x='',
              subtitle='pseudocount: +1')

ratio_leaf_melted %>%
  subset(!is.na(ratio)) %>%
  ggplot(aes(x=myPROMPTs, y=ratio+1, fill=genotype)) +
         geom_boxplot(alpha=.75, outlier.size=.5) +
         scale_y_log10() +
         cowplot::theme_cowplot() +
         facet_grid(length~., scales='free_x', space='free_x') +
         scale_fill_manual(values=c('aquamarine3', 'brown3')) +
         labs(title='Ratio: small RNA / expression - LEAF',
              y='Ratio of small RNA versus PROMPT TC expression or mRNA expression\nlog-scale, pseudocount +1',
              x='',
              subtitle='pseudocount: +1')




# 9 . FIND EXAMPLES ####
# ----------------------
ratio_flower %>%
  mutate('ratio_mut_wt'=rrp4_22_ratio / wt_22_ratio) %>%
  select(GeneID, myPROMPTs, ratio_mut_wt) %>%
  arrange(desc(ratio_mut_wt)) %>%
  subset(myPROMPTs)

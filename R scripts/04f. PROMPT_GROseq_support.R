#### Arabidopsis TSS story: GRO-seq support of PROMPTs
#### Axel Thieffry
set.seed(42)
library(ggridges)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(CAGEfightR)
library(TeMPO)
library(tidyverse)
library(tidylog)
options(scipen=999)

'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename

remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { GR[-idx]}
                                     else {GR}}

setwd('~/masked_path/04 - TSS_Level DE')



# 1. GET DATA ####
# ----------------
# CAGE TCs
TCs <- readRDS('~/masked_path/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')
# TxDB
txdb <- TxDb.Athaliana.BioMart.plantsmart28
# GROseq: 5'GRO-cap (Hetzel), GRO-seq (Hetzel), GRO-seq (Jacobsen)
gro_p <- list.files('~/masked_path/03 - TSS analysis', pattern='GRO.*plus.bw', full.names=T) %>% BigWigFileList()
gro_m <- list.files('~/masked_path/03 - TSS analysis', pattern='GRO.*minus.bw', full.names=T) %>% BigWigFileList()
names(gro_p) <- names(gro_m) <- c('GROcap_Hetzel', 'GROseq_Hetzel', 'GROseq_Jacobsen')
# CAGE signal
cage_p <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*plus'))
cage_m <- BigWigFileList(list.files('~/masked_path/bw_files_R123', full.names=T, pattern='_0.*minus'))
cage_names <- list.files('~/masked_path/bw_files_R123', pattern='_0.*plus') %>% str_remove('_0_R123.plus.tpm.bw')
names(cage_p) <- names(cage_m) <- cage_names
# make seqinfo from the CAGE signal
myseqinfo <- seqinfo(cage_p$wt)



# 2. DEFINE ALL PROMPT REGIONS FROM TXDB ####
# -------------------------------------------
# 2a. get genes
genes <- genes(txdb)
  seqlevelsStyle(genes) <- seqlevelsStyle(myseqinfo)
  seqlevels(genes) <- seqlevels(myseqinfo)
  seqinfo(genes) <- myseqinfo
  
# 2b. get PROMPT regions
athal_prompts <- promoters(genes, upstream=400, downstream=0, use.names=T) %>% remove_out_of_bound()
  table(width(athal_prompts))
  
# 2c. invert strandness of PROMPT regions
athal_prompts %<>% invertStrand()

# 2d. sanity-check IGV (all OK)
if(FALSE) {export.bed(genes, '~/masked_path/genes.bed')
           export.bed(athal_prompts, '~/masked_path/Athal_PROMPTs.bed')}

# 2e. indicate if PROMPT overlap a gene (strand-independent)
df_prompts_gene_countOverlap <- countOverlaps(athal_prompts, genes, ignore.strand=T) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_colnames(c('geneID', 'counts')) %>%
  mutate('overlap_gene'=ifelse(counts==0, 'no', 'yes')) %>%
  as_tibble()

athal_prompts$overlap_gene <- ifelse(athal_prompts$gene_id %in% subset(df_prompts_gene_countOverlap, overlap_gene=='no')$geneID, 'no', 'yes')
  stopifnot(table(df_prompts_gene_countOverlap$overlap_gene) == table(athal_prompts$overlap_gene))

# 2f. indicate if PROMPTs overlaps another PROMPT (strand-independent)
df_prompts_prompts_countOverlap <- countOverlaps(athal_prompts, athal_prompts, ignore.strand=T) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_colnames(c('geneID', 'counts')) %>%
  mutate('overlap_prompt'=ifelse(counts == 1, 'no', 'yes')) %>%
  as_tibble()

athal_prompts$overlap_prompt <- ifelse(athal_prompts$gene_id %in% subset(df_prompts_prompts_countOverlap, overlap_prompt=='no')$geneID, 'no', 'yes')
  stopifnot(table(df_prompts_prompts_countOverlap$overlap_prompt) == table(athal_prompts$overlap_prompt))

# 2g. indicate if among my 96 exosome-sensitive PROMPTs
my_PROMPTs <- rowRanges(TCs) %>%
  subset(genotypehen2 == 1 | genotyperrp4 == 1) %>%
  subset(txType_TAIR10 == 'reverse')

  # check that all my exosome-sensitive PROMPTs are in the Atlas
  stopifnot(all(my_PROMPTs$geneID_anti %in% athal_prompts$gene_id))
  # indicate in the Atlas if belong to my exosome-sensitive PROMPTs
  athal_prompts$my_PROMPT <- ifelse(athal_prompts$gene_id %in% my_PROMPTs$geneID_anti, 'yes', 'no')
  stopifnot(all(my_PROMPTs$geneID_anti %in% subset(athal_prompts, my_PROMPT=='yes')$gene_id))
  # NB: there are 94 gene for 96 PROMPTs because some gene have the same PROMPT region

# 2h. add PROMPT start as IRanges in the mcols (for TeMPO compatibility)
athal_prompts$PROMPT_start <- ifelse(strand(athal_prompts)=='+', start(athal_prompts), end(athal_prompts))
athal_prompts$PROMPT_start <- IRanges(start=athal_prompts$PROMPT_start, width=1)



# 3. COMPUTE CAGE SIGNAL AT PROMPTS ####
# --------------------------------------
# compute sense and antisense CAGE signal for PROMPTs
# then keep only sense signal because PROMPT strand is antisense to the gene
# then compute rowSums to get total signal per genotype per PROMPT
# do not SIMPLIFY so that it gets out as a matrix, then converted into a tibble
cage_signal <- mapply(function(x, y)wideMetaProfile(sites=swapRanges(athal_prompts, inputColumn='PROMPT_start'), forward=x, reverse=y, upstream=1, downstream=400)$sense %>% rowSums(), cage_p, cage_m)
cage_signal %<>% as.data.frame() %>% rownames_to_column('gene_id') %>% as_tibble()
# add CAGE signal to master GR (need to re-order first then add, left_join will not work)
cage_signal <- cage_signal[match(athal_prompts$gene_id, cage_signal$gene_id), ]
  stopifnot(athal_prompts$gene_id == cage_signal$gene_id)
mcols(athal_prompts) <- cbind(mcols(athal_prompts), select(cage_signal, -gene_id))
  # check all PROMPTs have signal
  stopifnot(any(!is.na(athal_prompts$wt)), any(!is.na(athal_prompts$hen2)), any(!is.na(athal_prompts$rrp4)))



# 4. COMPUTE GRO-SEQ SIGNAL AT PROMPTS ####
# # ---------------------------------------
# 4a. GRO-seq Hetzel
groseq_signal <- mapply(function(x, y) wideMetaProfile(sites=swapRanges(athal_prompts, inputColumn='PROMPT_start'), forward=x, reverse=y, upstream=1, downstream=400)$sense %>% rowSums(), gro_p, gro_m)
groseq_signal %<>% as.data.frame() %>% rownames_to_column('gene_id') %>% as_tibble()
# add CAGE signal to master GR (need to re-order first then add, left_join will not work)
groseq_signal <- groseq_signal[match(athal_prompts$gene_id, groseq_signal$gene_id), ]
  stopifnot(athal_prompts$gene_id == groseq_signal$gene_id)
mcols(athal_prompts) <- cbind(mcols(athal_prompts), select(groseq_signal, -gene_id))
  # check all PROMPTs have signal
  stopifnot(any(!is.na(athal_prompts$GROscap_Hetzel)), any(!is.na(athal_prompts$GROseq_Hetzel)), any(!is.na(athal_prompts$GROseq_Jacobsen)))



# 5. INVESTIGATE SIGNALS ####
# ---------------------------
# 5a. calculate mean GROX signal accross PROMPTs, genome-wide, after filtering for zero overlapping
mean_GRO_df <- athal_prompts_df %>%
    subset(overlap_gene=='no' & overlap_prompt=='no') %>%
    select(GROcap_Hetzel, GROseq_Hetzel, GROseq_Jacobsen) %>%
    colMeans() %>%
    t() %>%
    as.data.frame()

  
# 5b. number of PROMPT regions above 2 * mean genome-wide signal
table(subset(athal_prompts_df, GROcap_Hetzel >= 2 * mean_GRO_df$GROcap_Hetzel)$my_PROMPT)[2] / 94 * 100
table(subset(athal_prompts_df, GROseq_Hetzel >= 2 * mean_GRO_df$GROseq_Hetzel)$my_PROMPT)[2] / 94 * 100
table(subset(athal_prompts_df, GROseq_Jacobsen >= 2 * mean_GRO_df$GROseq_Jacobsen)$my_PROMPT)[2] / 94 * 100

data.frame('GRO'=c('GROcap_Hetzel', 'GROseq_Hetzel', 'GROseq_Jacobsen'),
           'yes'=c(63, 50, 36),
           'no'=c(94-63, 94-50, 94-36)) %>%
  melt(id.vars='GRO') %>%
  mutate('variable'=factor(variable, levels=c('no', 'yes'))) %>%
  mutate('GRO'=str_replace(GRO, '_', ' ')) %>%
  mutate('GRO'=str_replace(GRO, 'GRO', 'GRO-')) %>%
  ggplot(aes(x=GRO, fill=variable, y=value)) +
         geom_bar(stat='identity', lwd=.2, col='black') +
         geom_text(aes(label=paste0(round(value/94*100, 1), '%')), position=position_stack(), vjust=1.5) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, axis.text.x=element_text(angle=45, hjust=1)) +
         scale_fill_brewer(palette='Set2', direction=-1, name='supported') +
         labs(title='GRO-X support of\nexosome-sensitive PROMPTs',
              x='', y='Nb. exosome-sensitive\nPROMPTs regions',
              caption='Support= higher GRO-X signal than\ntwice the genome-wide average')

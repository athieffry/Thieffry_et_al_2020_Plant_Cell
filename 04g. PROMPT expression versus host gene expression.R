#### Arabidopsis: PROMPT expression versus host gene expression
#### Axel Thieffry - December 2019
set.seed(42)
library(patchwork)
library(RColorBrewer)
library(magrittr)
library(reshape2)
library(GenomicFeatures)
library(GenomicRanges)
library(ggpubr)
library(TeMPO)
library(rtracklayer)
library(BiocParallel)
register(MulticoreParam(4), default=TRUE)
library(CAGEfightR)
library(tidyverse)
library(tidylog)
'%!in%' <- function(x,y)!('%in%'(x,y))
'select' <- dplyr::select
'rename' <- dplyr::rename
setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE')



# 1. LOAD ALL INPUT FILES ####
# ----------------------------
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
TCs <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_TCs_with_all_data_for_PROMPT_GROseq_support.rds')
geneLevel <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/SE_genelevel_TSSstory.rds')

# 2. GET PROMPTs only ####
# ------------------------
# prompt GR
prompts <- TCs %>%
  rowRanges() %>%
  subset(genotypehen2 == 1 | genotyperrp4 == 1) %>%
  subset(txType_TAIR10 == 'reverse')

# prompt ids
prompt_ids <- names(prompts) %>% unique()

# prompt genes
host_gene_ids <- prompts %>%
                 as.data.frame() %>%
                 select(thick.names, geneID_anti)


# 3. GET PROMPTs expression in average TPM for each sample ####
# -------------------------------------------------------------
prompt_tpm <- assay(TCs, 'TPM') %>%
  as.data.frame() %>%
  rownames_to_column('TC_id') %>%
  melt(id.vars='TC_id') %>%
  separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
  select(-timepoint) %>%
  group_by(TC_id, genotype) %>%
  summarise('meanTPM'=mean(value)) %>%
  ungroup() %>%
  spread(key=genotype, value=meanTPM) %>%
  subset(TC_id %in% prompt_ids)

# rename to know it's prompt TPM
prompt_tpm %<>% rename('prompt_hen2'='hen2', 'prompt_rrp4'='rrp4', 'prompt_wt'='wt')

# add host gene
prompt_tpm %<>% left_join(host_gene_ids, by=c('TC_id'='thick.names')) %>% head(30)



# 4. GET Gene TPM expression ####
# -------------------------------
# calculate TPM
geneLevel <- calcTPM(geneLevel, inputAssay='counts', outputAssay='TPM', totalTags='totalTags')

# get average TPM for host genes
gene_tpm <- assay(geneLevel, 'TPM') %>%
  as.data.frame() %>%
  rownames_to_column('gene_id') %>%
  melt(id.vars='gene_id') %>%
  separate(variable, c('genotype', 'timepoint', 'replicate'), sep='_') %>%
  select(-timepoint) %>%
  group_by(gene_id, genotype) %>%
  summarise('meanTPM'=mean(value)) %>%
  ungroup() %>%
  spread(key=genotype, value=meanTPM) %>%
  subset(gene_id %in% host_gene_ids$geneID_anti)

# rename to know it's gene TPM
gene_tpm %<>% rename('gene_hen2'='hen2', 'gene_rrp4'='rrp4', 'gene_wt'='wt')

# merge PROMPT & gene expression
data <- left_join(prompt_tpm, gene_tpm, by=c('geneID_anti'='gene_id'))

# re-order
data %<>% select('PROMPT'=TC_id, 'host_gene'=geneID_anti, prompt_hen2, gene_hen2, prompt_rrp4, gene_rrp4, prompt_wt, gene_wt)


# 5. LOg-LOg plots ####
# ---------------------
# wt
gg_wt <- ggplot(data, aes(x=gene_wt +1, y=prompt_wt +1)) +
       geom_point() +
       geom_smooth(method='lm', se=F) +
       stat_cor(output.type='text', method='spearman') +
       scale_x_log10() + scale_y_log10() +
       cowplot::theme_cowplot() + theme(aspect.ratio=1) +
       labs(x='wt host gene\nlog(TPM+1)',
            y='wt PROMPT TC\nlog(TPM+1)')

# hen2
gg_hen2 <- ggplot(data, aes(x=gene_hen2 +1, y=prompt_hen2 +1)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_cor(output.type='text', method='spearman') +
  scale_x_log10() + scale_y_log10() +
  cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  labs(x='hen2 host gene\nlog(TPM+1)',
       y='hen2 PROMPT TC\nlog(TPM+1)')

# rrp4
gg_rrp4 <- ggplot(data, aes(x=gene_rrp4 +1, y=prompt_rrp4 +1)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_cor(output.type='text', method='spearman') +
  scale_x_log10() + scale_y_log10() +
  cowplot::theme_cowplot() + theme(aspect.ratio=1) +
  labs(x='rrp4 host gene\nlog(TPM+1)',
       y='rrp4 PROMPT TC\nlog(TPM+1)')

gg_wt + gg_hen2 + gg_rrp4 + plot_layout(ncol=1, nrow=3)


# 6. CORRELATIONS ####
# --------------------
with(na.omit(data), cor(log(prompt_wt+1), log(gene_wt+1), method='spearman'))
with(na.omit(data), cor(prompt_hen2, gene_hen2))
with(na.omit(data), cor(prompt_rrp4, gene_rrp4))

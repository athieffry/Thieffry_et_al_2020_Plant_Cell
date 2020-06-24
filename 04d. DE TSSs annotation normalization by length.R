#### Arabidopsis flg22 : DE TSSs by annotation, normalized by annotation length
#### Axel Thieffry - Mars 2019
set.seed(42)
library(WriteXLS)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggalluvial)
library(ggridges)
library(pheatmap)
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
library(patchwork)
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
remove_out_of_bound <- function(GR) {idx = GenomicRanges:::get_out_of_bound_index(GR)
                                     if(length(idx) != 0) { GR[-idx]}
                                     else {GR}}

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/04 - TSS_Level DE')


# 1. LOAD ALL INPUT FILES ####
# ----------------------------
# myseqinfo
myseqinfo <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/myseqinfo.rds')
# colors
exo_colors2 <- brewer.pal(name='Set3', n=4)[3:4]
# DE TSSs: annotation DF
sense_annot <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/dtByAnot.rds')
anti_annot <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/dtByAnot_antisense.rds') %>% ungroup()
# annotation hierarchies
sense_hierarchy <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/custom_annotation_hierarchy_TAIR10.rds')
anti_hierarchy  <- readRDS('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - RDATA/custom_annotation_hierarchy_TAIR10_extended_antisense.rds')



# 2. GET EXACT LENGTH OF SENSE ANNOTATION CATEGORIES ####
# -------------------------------------------------------
# promoter lengths
promoter_gr <- sense_hierarchy$promoter %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(promoter_gr, '04d. bed sanity checks/promoters.bed')
promoter_length <- promoter_gr %>% width() %>% sum()

# 5'UTR lengths (5'UTRs - overlap with promoters)
fiveUTR_gr <- sense_hierarchy$fiveUTR %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(fiveUTR_gr, promoter_gr)
tmp_grl  <- extractList(promoter_gr, as(tmp_hits, "List"))
fiveUTR_gr <- psetdiff(fiveUTR_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(fiveUTR_gr, '04d. bed sanity checks/5UTRs.bed')
fiveUTR_length <- fiveUTR_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs
    prom_five_gr <- GenomicRanges::reduce(c(promoter_gr, fiveUTR_gr))

# 3'UTR lengths (3'UTRs - overlap with promoters and 5'UTRs)
threeUTR_gr <- sense_hierarchy$threeUTR %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(threeUTR_gr, prom_five_gr)
tmp_grl  <- extractList(prom_five_gr, as(tmp_hits, "List"))
threeUTR_gr <- psetdiff(threeUTR_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(threeUTR_gr, '04d. bed sanity checks/3UTRs.bed')
threeUTR_length <-  threeUTR_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs and 3'UTRs
    prom_five_three_gr <- GenomicRanges::reduce(c(prom_five_gr, threeUTR_gr), ignore.strand=T)

# CDS lengths (CDS - overlap with promoters and 5'UTRs and 3'UTRs)
CDS_gr <- sense_hierarchy$CDS %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(CDS_gr, prom_five_three_gr)
tmp_grl  <- extractList(prom_five_three_gr, as(tmp_hits, "List"))
CDS_gr <- psetdiff(CDS_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(CDS_gr, '04d. bed sanity checks/CDS.bed')
CDS_length <- CDS_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs and 3'UTRs and CDS
    prom_five_three_CDS_gr <- GenomicRanges::reduce(c(prom_five_three_gr, CDS_gr), ignore.strand=T)

# exons lengths (exons - overlap with promoters and 5'UTRs and 3'UTRs and CDS)
exon_gr <- sense_hierarchy$exon %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(exon_gr, prom_five_three_CDS_gr)
tmp_grl  <- extractList(prom_five_three_CDS_gr, as(tmp_hits, "List"))
exon_gr <- psetdiff(exon_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(exon_gr, '04d. bed sanity checks/exon.bed')
exon_length <- exon_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs and 3'UTRs and CDS and exons
    prom_five_three_CDS_exon_gr <- GenomicRanges::reduce(c(prom_five_three_CDS_gr, exon_gr), ignore.strand=T)

# intron lengths (intron - overlap with promoters and 5'UTRs and 3'UTRs and CDS and exon)
intron_gr <- sense_hierarchy$intron %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(intron_gr, prom_five_three_CDS_exon_gr)
tmp_grl  <- extractList(prom_five_three_CDS_exon_gr, as(tmp_hits, "List"))
intron_gr <- psetdiff(intron_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(intron_gr, '04d. bed sanity checks/intron.bed')
intron_length <- intron_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs and 3'UTRs and CDS and exons and introns
    prom_five_three_CDS_exon_introns_gr <- GenomicRanges::reduce(c(prom_five_three_CDS_exon_gr, intron_gr), ignore.strand=T)

# proximal lengths (proximal - overlap with promoters and 5'UTRs and 3'UTRs and CDS and exon and introns)
proximal_gr <- sense_hierarchy$proximal %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(proximal_gr, prom_five_three_CDS_exon_introns_gr)
tmp_grl  <- extractList(prom_five_three_CDS_exon_introns_gr, as(tmp_hits, "List"))
proximal_gr <- psetdiff(proximal_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(proximal_gr, '04d. bed sanity checks/proximal.bed')
proximal_length <- proximal_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge promoters and 5'UTRs and 3'UTRs and CDS and exons and introns and proximal
    prom_five_three_CDS_exon_introns_proximal_gr <- GenomicRanges::reduce(c(prom_five_three_CDS_exon_introns_gr, proximal_gr), ignore.strand=T)

# intergenic lengths (intergenic - overlap with promoters and 5'UTRs and 3'UTRs and CDS and exon and introns and intergenic and proximal)
intergenic_gr <- sense_hierarchy[names(sense_hierarchy) != 'antisense'] %>% # get sense regions only
                 unlist() %>%
                 GenomicRanges::reduce(ignore.strand=T) %>%
                 gaps(start=NA, end=NA)
export.bed(intergenic_gr, '04d. bed sanity checks/intergenic.bed')
intergenic_length <- intergenic_gr %>% width() %>% sum()

# antisense lengths
antisense_gr <- sense_hierarchy$antisense %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(antisense_gr, '04d. bed sanity checks/antisense.bed')
antisense_lengths <- antisense_gr %>% width() %>% sum()

# reverse lengths
reverse_gr <- sense_hierarchy$reverse %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(reverse_gr, '04d. bed sanity checks/reverse.bed')
reverse_lengths <- reverse_gr %>% width() %>% sum()

# combine all in one df
df_length_sense <- data.frame('annotation'=c('promoter', 'fiveUTR', 'threeUTR', 'CDS', 'exon', 'intron', 'proximal', 'intergenic', 'antisense', 'reverse'),
                             'length_bp'=c(promoter_length, fiveUTR_length, threeUTR_length, CDS_length, exon_length, intron_length, proximal_length, intergenic_length, antisense_lengths, reverse_lengths)) %>%
                   as_tibble()

levels_sense=c('intergenic', 'antisense', 'reverse', 'proximal', 'promoter', 'fiveUTR', 'CDS', 'exon', 'intron', 'threeUTR')

df_length_sense %<>% mutate('length_Mbp'=length_bp/1000000)
df_length_sense %<>% mutate('annotation'=factor(annotation, levels=levels_sense))
df_length_sense




# 3. GET EXACT LENGTH OF ANTISENSE ANNOTATION CATEGORIES ####
# -----------------------------------------------------------
# anti-5'UTR lengths
anti_fiveUTR_gr <- anti_hierarchy$antisense_fiveUTR %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(anti_fiveUTR_gr, '04d. bed sanity checks/anti_5UTR.bed')
anti_fiveUTR_length <- anti_fiveUTR_gr %>% width() %>% sum()

# anti-3'UTR lengths (3'UTRs - overlap with 5'UTRs)
anti_threeUTR_gr <- anti_hierarchy$antisense_threeUTR %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(anti_threeUTR_gr, anti_fiveUTR_gr)
tmp_grl  <- extractList(anti_fiveUTR_gr, as(tmp_hits, "List"))
anti_threeUTR_gr <- psetdiff(anti_threeUTR_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(anti_threeUTR_gr, '04d. bed sanity checks/anti_3UTR.bed')
anti_threeUTR_length <- anti_threeUTR_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge 5'UTRs and 3'UTRs (anti)
    anti_five_three_gr <- GenomicRanges::reduce(c(anti_fiveUTR_gr, anti_threeUTR_gr))

# anti-CDS lengths (CDS - overlap with 5'UTRs and 3'UTRs anti)
anti_CDS_gr <- anti_hierarchy$antisense_CDS %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(anti_CDS_gr, anti_five_three_gr)
tmp_grl  <- extractList(anti_five_three_gr, as(tmp_hits, "List"))
anti_CDS_gr <- psetdiff(anti_CDS_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(anti_CDS_gr, '04d. bed sanity checks/anti_CDS.bed')
anti_CDS_length <- anti_CDS_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge 5'UTRs and 3'UTRs and CDS (anti)
    anti_five_three_CDS_gr <- GenomicRanges::reduce(c(anti_five_three_gr, anti_CDS_gr))

# anti-exon lengths (exon - overlap with 5'UTRs and 3'UTRs and CDS anti)
anti_exon_gr <- anti_hierarchy$antisense_exon %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(anti_exon_gr, anti_five_three_CDS_gr)
tmp_grl  <- extractList(anti_five_three_CDS_gr, as(tmp_hits, "List"))
anti_exon_gr <- psetdiff(anti_exon_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(anti_exon_gr, '04d. bed sanity checks/anti_exon.bed')
anti_exon_length <- anti_exon_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

    # merge 5'UTRs and 3'UTRs and CDS and exon (anti)
    anti_five_three_CDS_exon_gr <- GenomicRanges::reduce(c(anti_five_three_CDS_gr, anti_exon_gr))

# anti-intron lengths (intron - overlap with 5'UTRs and 3'UTRs and CDS and intron anti)
anti_intron_gr <- anti_hierarchy$antisense_intron %>% GenomicRanges::reduce(ignore.strand=T)
tmp_hits <- findOverlaps(anti_intron_gr, anti_five_three_CDS_exon_gr)
tmp_grl  <- extractList(anti_five_three_CDS_exon_gr, as(tmp_hits, "List"))
anti_intron_gr <- psetdiff(anti_intron_gr, tmp_grl) %>% unlist() %>% GenomicRanges::reduce(ignore.strand=T)
export.bed(anti_intron_gr, '04d. bed sanity checks/anti_intron.bed')
anti_intron_length <- anti_intron_gr %>% width() %>% sum()
rm(tmp_hits, tmp_grl)

# combine all in one df
anti_annot$txType_TAIR10extended %<>% droplevels()

anti_annot$txType_TAIR10extended %>% unique() %>% as.data.frame()
df_length_anti <- data.frame('annotation'=c(paste0('antisense_', c('fiveUTR', 'threeUTR', 'CDS', 'exon', 'intron'))),
                             'length_bp'=c(anti_fiveUTR_length, anti_threeUTR_length, anti_CDS_length, anti_exon_length, anti_intron_length)) %>%
                  as_tibble()

levels_anti=c('antisense_threeUTR', 'antisense_intron', 'antisense_fiveUTR', 'antisense_exon', 'antisense_CDS')
df_length_anti %<>% mutate('length_Mbp'=length_bp/1000000)
df_length_anti %<>% mutate('annotation'=factor(annotation, levels=levels_anti))
df_length_anti


  
# 4. PLOT SENSE ####
# ------------------
# sense annotation
gg_sense_a <- ggplot(sense_annot, aes(x=txType_TAIR10, y=ifelse(direction=='down', -nTCs, nTCs), fill=coef)) +
       geom_bar(stat='identity', position=position_dodge(), col='black', lwd=.2) +
       geom_text(aes(y=ifelse(direction=='down', -nTCs-60, nTCs+60), label=nTCs, group=coef), position=position_dodge(width=.9), size=3) +
       geom_hline(yintercept=0, lty=2) +
       cowplot::theme_cowplot() + coord_flip() +
       theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10)) +
       scale_fill_manual(values=exo_colors2, name='Coefficient') +
       labs(title='Sense DE TSSs', y='down / up DE TSSs', x='')

# total lengths
gg_sense_b <- ggplot(df_length_sense, aes(x=annotation, y=length_Mbp)) +
       geom_bar(stat='identity') +
       labs(title='Lengths of sense annotation', x='', y='Total region length (Mbp)') +
       coord_flip() + cowplot::theme_cowplot() +
       theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10))

# combine and normalize
sense_annot %<>% left_join(df_length_sense, by=c('txType_TAIR10'='annotation')) %>%
  mutate('fraction'=nTCs/length_Mbp) %>%
  mutate('txType_TAIR10'=factor(txType_TAIR10, levels=levels_sense))

# sense annotation normalized
gg_sense_c <- ggplot(sense_annot, aes(x=txType_TAIR10, y=ifelse(direction=='down', -fraction, fraction), fill=coef)) +
       geom_bar(stat='identity', position=position_dodge(), col='black', lwd=.2) +
       geom_text(aes(y=ifelse(direction=='down', -fraction-10, fraction+10), label=round(fraction, 1), group=coef), position=position_dodge(width=.9), size=3) +
       geom_hline(yintercept=0, lty=2) +
       cowplot::theme_cowplot() + coord_flip() +
       theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10)) +
       scale_fill_manual(values=exo_colors2, name='Coefficient') +
       labs(title='Sense DE TSSs - normalized', y='down / up DE TSSs\nNormalized to total length of annotation', x='')

gg_sense_a + gg_sense_b + gg_sense_c + plot_layout(ncol=1)




# 5. PLOT ANTISENSE ####
# ----------------------
# antisense annotation
gg_anti_a <- ggplot(anti_annot, aes(x=txType_TAIR10extended, y=ifelse(direction=='down', -count, count), fill=coef)) +
  geom_bar(stat='identity', position=position_dodge(), col='black', lwd=.2) +
  geom_text(aes(y=ifelse(direction=='down', -count-15, count+20), label=count, group=coef), position=position_dodge(width=.9), size=3) +
  geom_hline(yintercept=0, lty=2) +
  cowplot::theme_cowplot() + coord_flip() +
  theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10)) +
  scale_fill_manual(values=exo_colors2, name='Coefficient') +
  labs(title='Antisense DE TSSs', subtitle='extended annotation',
       x='', y='down / up DE TSSs')

# total lengths
gg_anti_b <- ggplot(df_length_anti, aes(x=annotation, y=length_Mbp)) +
  geom_bar(stat='identity') +
  labs(title='Lengths of antisense annotation', x='', y='Total region length (Mbp)') +
  coord_flip() + cowplot::theme_cowplot() +
  theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10))

# combine and normalize
anti_annot %<>% left_join(df_length_anti, by=c('txType_TAIR10extended'='annotation')) %>%
  mutate('fraction'=count/length_Mbp) %>%
  mutate('txType_TAIR10extended'=factor(txType_TAIR10extended, levels=levels_anti))

# sense annotation normalized
gg_anti_c <- ggplot(anti_annot, aes(x=txType_TAIR10extended, y=ifelse(direction=='down', -fraction, fraction), fill=coef)) +
       geom_bar(stat='identity', position=position_dodge(), col='black', lwd=.2) +
       geom_text(aes(y=ifelse(direction=='down', -fraction-3, fraction+3), label=round(fraction, 1), group=coef), position=position_dodge(width=.9), size=3) +
       geom_hline(yintercept=0, lty=2) +
       cowplot::theme_cowplot() + coord_flip() +
       theme(axis.text=element_text(size=10), text=element_text(size=10), plot.title=element_text(size=10)) +
       scale_fill_manual(values=exo_colors2, name='Coefficient') +
       labs(title='Antisense DE TSSs - normalized', y='down / up DE TSSs\nNormalized to total length of annotation', x='')


# PLOT ALL OF THEM TOGETHER
plot_spacer() / (gg_sense_a | gg_sense_b | gg_sense_c) / (gg_anti_a | gg_anti_b | gg_anti_c) / plot_spacer() + plot_layout(ncol=1, nrow=4, heights=c(1,3,3,1)) & theme(legend.position='none')

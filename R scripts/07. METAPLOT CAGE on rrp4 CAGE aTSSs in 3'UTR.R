#### DeepTools metaplot: CAGE signal on CAGE rrp4 3'UTR aTSSs
#### Axel Thieffry - June 2019
library(tidyverse)
library(magrittr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(robustbase)

'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))

setwd('~/masked_path/07 - metagene plots antisense')

# CAGE SCALED METAPLOT ON RNASEQ RRP4 ADEGS (FROM DEEPTOOLS) ####
# ---------------------------------------------------------------
pp <- read.delim('~/masked_path/matrix_CAGE_plus_BED_plus.mat', skip=2)[, -1]
mm <- read.delim('~/masked_path/matrix_CAGE_minus_BED_minus.mat', skip=2)[, -1]
pm <- read.delim('~/masked_path/matrix_CAGE_plus_BED_minus.mat', skip=2)[, -1]
mp <- read.delim('~/masked_path/matrix_CAGE_minus_BED_plus.mat', skip=2)[, -1]

(samples <- colnames(pp) %>% str_remove('\\..*') %>% unique())

matrix_sense <- list(rbind(pp[, 1:240], mm[, 1:240]),
                     rbind(pp[, 241:480], mm[, 241:480]),
                     rbind(pp[, 481:720], mm[, 481:720]))

matrix_anti <- list(rbind(pm[, 1:240], mp[, 1:240]),
                    rbind(pm[, 241:480], mp[, 241:480]),
                    rbind(pm[, 481:720], mp[, 481:720]))

names(matrix_sense) <- samples
names(matrix_anti) <- samples

# compute average
mat_sense_ave <- sapply(matrix_sense, function(x) colMeans(as.matrix(x))) %>% as_tibble()
mat_anti_ave <- sapply(matrix_anti, function(x) colMeans(as.matrix(x))) %>% as_tibble()

# add position
mat_sense_ave %<>% mutate('position'=1:240, 'direction'='sense')
mat_anti_ave %<>% mutate('position'=1:240, 'direction'='anti')

# plot average
gg_ave_sense <- rbind(mat_sense_ave, mat_anti_ave) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='averageRPM') %>%
  subset(direction=='sense') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -averageRPM, averageRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(20, 220), lty=2) +
         geom_line(lwd=.65) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position=c(.5, .5)) +
         labs(title='CAGE at rrp4 CAGE 3UTR aTSSs (N=176)',
              x='', y='Average CAGE TPM') +
         scale_color_brewer(palette='Set2', name='', direction=-1)

gg_ave_anti <- rbind(mat_sense_ave, mat_anti_ave) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='averageRPM') %>%
  subset(direction=='anti') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -averageRPM, averageRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(20, 220), lty=2) +
         geom_line(lwd=.65) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
         labs(title='CAGE at rrp4 CAGE 3UTR aTSSs (N=176)',
              x='', y='Average CAGE TPM') +
         scale_color_brewer(palette='Set2', name='', direction=-1)

gg_ave_sense + gg_ave_anti + plot_layout(ncol=1)

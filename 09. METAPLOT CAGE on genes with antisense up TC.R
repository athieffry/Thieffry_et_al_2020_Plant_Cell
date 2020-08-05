#### DeepTools metaplot: CAGE signal on any gene with aTSS up
#### Axel Thieffry - June 2019
library(tidyverse)
library(magrittr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(tidyquant)
library(RColorBrewer)
library(patchwork)
library(robustbase)

'select' <- dplyr::select
'rename' <- dplyr::rename
'%!in%' <- function(x,y)!('%in%'(x,y))

setwd('~/masked_path/09 - metagene plots - CAGE on genes with any antisense up CAGE TC/')

# CAGE SCALED METAPLOT ON RNASEQ RRP4 ADEGS (FROM DEEPTOOLS) ####
# ---------------------------------------------------------------
pp <- read.delim('~/masked_path/matrix_CAGE_plus_BED_plus.mat', skip=2)[, -1]
mm <- read.delim('~/masked_path/matrix_CAGE_minus_BED_minus.mat', skip=2)[, -1]
pm <- read.delim('~/masked_path/matrix_CAGE_plus_BED_minus.mat', skip=2)[, -1]
mp <- read.delim('~/masked_path/matrix_CAGE_minus_BED_plus.mat', skip=2)[, -1]

(samples <- colnames(pp) %>% str_remove('\\..*') %>% unique())
matdim <- ncol(pp)
nsamples=3
matdim/nsamples

matrix_sense <- list(rbind(pp[, 1:280], mm[, 1:280]),
                     rbind(pp[, 281:560], mm[, 281:560]),
                     rbind(pp[, 561:840], mm[, 561:840]))

matrix_anti <- list(rbind(pm[, 1:280], mp[, 1:280]),
                    rbind(pm[, 281:560], mp[, 281:560]),
                    rbind(pm[, 561:840], mp[, 561:840]))

names(matrix_sense) <- samples
names(matrix_anti) <- samples

# compute average
mat_sense_ave <- sapply(matrix_sense, function(x) colMeans(as.matrix(x))) %>% as_tibble()
mat_anti_ave <- sapply(matrix_anti, function(x) colMeans(as.matrix(x))) %>% as_tibble()

# compute median
mat_sense_med <- sapply(matrix_sense, function(x) colMedians(as.matrix(x))) %>% as_tibble()
mat_anti_med <- sapply(matrix_anti, function(x) colMedians(as.matrix(x))) %>% as_tibble()

# add position
mat_sense_ave %<>% mutate('position'=1:280, 'direction'='sense')
mat_anti_ave %<>% mutate('position'=1:280, 'direction'='anti')

mat_sense_med %<>% mutate('position'=1:280, 'direction'='sense')
mat_anti_med %<>% mutate('position'=1:280, 'direction'='anti')

# make colors
plotcol <- c(brewer.pal(n=4, name='Set1')[4], brewer.pal(n=3, name='Set1')[1], brewer.pal(n=4, name='Paired')[4])

# plot average
gg_ave_sense <- rbind(mat_sense_ave, mat_anti_ave) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='averageRPM') %>%
  subset(direction=='sense') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -averageRPM, averageRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(40, 240), lty=2) +
         geom_line(lwd=.65) +
         cowplot::theme_cowplot() + theme(legend.position=c(.5, .5)) +
         labs(title='CAGE at any gene with aTSS up',
              x='', y='Average CAGE TPM') +
         scale_color_manual(values=plotcol)

gg_ave_anti <- rbind(mat_sense_ave, mat_anti_ave) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='averageRPM') %>%
  subset(direction=='anti') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -averageRPM, averageRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(40, 240), lty=2) +
         geom_line(lwd=.65) +
         cowplot::theme_cowplot() + theme(legend.position='none') +
         labs(title='CAGE at any gene with aTSS up',
              x='', y='Average CAGE TPM') +
         scale_color_manual(values=plotcol)

gg_ave_sense + gg_ave_anti + plot_layout(ncol=1, heights=c(1,1), widths=c(1,1))

# plot median
gg_med_sense <- rbind(mat_sense_med, mat_anti_med) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='medianRPM') %>%
  subset(direction=='sense') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -medianRPM, medianRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(40, 240), lty=2) +
         #geom_line(lwd=.65, alpha=.3) +
         geom_ma(n=10, lty=1) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position=c(.5, .5)) +
         labs(title='CAGE at any gene with aTSS up',
              x='', y='CAGE TPM (median, ma)') +
         scale_color_manual(values=plotcol)

gg_med_anti <- rbind(mat_sense_med, mat_anti_med) %>%
  melt(id.vars=c('position', 'direction'), variable.name='genotype', value.name='medianRPM') %>%
  subset(direction=='anti') %>%
  ggplot(aes(x=position, y=ifelse(direction=='anti', -medianRPM, medianRPM), group=interaction(direction, genotype), col=genotype)) +
         geom_vline(xintercept=c(40, 240), lty=2) +
         #geom_line(lwd=.65, alpha=.3) +
         geom_ma(n=10, lty=1) +
         cowplot::theme_cowplot() + theme(aspect.ratio=1, legend.position='none') +
         labs(title='', x='', y='CAGE TPM (median, ma)') +
         scale_color_manual(values=plotcol)

gg_med_sense + gg_med_anti + plot_layout(ncol=1)

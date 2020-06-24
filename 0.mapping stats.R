# Analysis of mapping stats
# Axel Thieffry - February 2019
library(tidyverse)
library(magrittr)

setwd('~/Dropbox/Axel_Arabidopsis_Flagellin/ANALYSIS_TSSstory/00 - mapping stats/')


stats <- read.table('mapping.stats', h=T) %>% as_tibble()
stats %<>% separate(sample, c('genotype', 'timepoint', 'replicate'), sep='_')

stats_zero <- subset(stats, timepoint==0)

# average unique mappers per library (in millions reads)
mean(stats_zero$unique) / 1000000
# same number in percent
mean(stats_zero$unique_pc)
  

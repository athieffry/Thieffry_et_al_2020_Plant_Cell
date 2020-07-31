# Analysis of stats from mapping CAGE reads to Arabidopsis thaliana reference genome
# Axel Thieffry - February 2019
library(tidyverse)
library(tidylog)
library(magrittr)

setwd('masked_path')

stats <- read.table('mapping.stats', h=T) %>% as_tibble()
stats %<>% separate(sample, c('genotype', 'timepoint', 'replicate'), sep='_')

stats_zero <- subset(stats, timepoint==0)

# average unique mappers per library (in millions reads and percentage)
mean(stats_zero$unique) / 1000000
mean(stats_zero$unique_pc)
  

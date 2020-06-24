#!/bin/bash

# computeMatrix of CAGE - on BED +
computeMatrix scale-regions \
        --binSize 5 -p 4 -a 100 -b 100 -m 1000 --skipZeros --metagene --missingDataAsZero \
        --samplesLabel 'hen2' 'rrp4' 'wt' \
        -R genes_with_3UTR_aTSS_up_in_CAGE_rrp4.bed.plus.txt \
        -S ~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123/*_0_R123.minus.tpm.bw \
        -o matrix_CAGE_minus_BED_plus \
        --outFileNameMatrix matrix_CAGE_minus_BED_plus.mat


# plot average profile
plotProfile -m matrix_CAGE_minus_BED_plus \
        -o profiles_CAGE_minus_BED_plus.pdf \
        --averageType mean \
        --perGroup \
        --plotTitle 'CAGE - on BED +' \
        -y 'Normalized Read Count (mean)' \
        --numPlotsPerRow 4 \
        --colors red violet blue

# plot heatmap
#plotHeatmap -m matrix_CAGE_minus_BED_plus \
#        -o heatmap_CAGE_minus_BED_plus.pdf \
#        --sortRegions descend \
#        --sortUsing mean \
#        --averageTypeSummaryPlot mean \
#        --missingDataColor white \
#        --startLabel TSS \
#        --endLabel PAS \
#        --regionsLabel GeneBody \
#        --yAxisLabel NormReadCount \
#        --plotFileFormat pdf \
#        --colorList white,blue,red \
#        --samplesLabel 'lsm8-2 rrp4' 'lsm8-2' 'rrp4' 'wt' \
#        --heatmapHeight 14 --heatmapWidth 8

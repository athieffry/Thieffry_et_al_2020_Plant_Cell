#!/bin/bash

# computeMatrix of CAGE - on BED -
computeMatrix scale-regions \
        --binSize 5 -p 4 -a 200 -b 200 -m 1000 --skipZeros --metagene --missingDataAsZero \
        --samplesLabel 'hen2' 'rrp4' 'wt' \
        -R genes_antisense_TC_up.bed.minus.txt \
        -S ~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123/*_0_R123.minus.tpm.bw \
        -o matrix_CAGE_minus_BED_minus \
        --outFileNameMatrix matrix_CAGE_minus_BED_minus.mat


# plot average profile
plotProfile -m matrix_CAGE_minus_BED_minus \
        -o profiles_CAGE_minus_BED_minus.pdf \
        --averageType mean \
        --perGroup \
        --plotTitle 'CAGE - on BED -' \
        -y 'Normalized Read Count (mean)' \
        --numPlotsPerRow 4 \
        --colors violet red darkgreen

# plot heatmap
#plotHeatmap -m matrix_CAGE_minus_BED_minus \
#        -o heatmap_CAGE_minus_BED_minus.pdf \
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

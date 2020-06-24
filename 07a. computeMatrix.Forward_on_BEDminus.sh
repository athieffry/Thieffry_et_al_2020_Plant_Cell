#!/bin/bash

# computeMatrix of CAGE + on BED -
computeMatrix scale-regions \
        --binSize 5 -p 4 -a 100 -b 100 -m 1000 --skipZeros --metagene --missingDataAsZero \
        --samplesLabel 'hen2' 'rrp4' 'wt' \
        -R genes_with_3UTR_aTSS_up_in_CAGE_rrp4.bed.minus.txt \
        -S ~/Dropbox/Axel_Arabidopsis_Flagellin/CAGE/bw_files_R123/*_0_R123.plus.tpm.bw \
        -o matrix_CAGE_plus_BED_minus \
        --outFileNameMatrix matrix_CAGE_plus_BED_minus.mat


# plot average profile
plotProfile -m matrix_CAGE_plus_BED_minus \
        -o profiles_CAGE_plus_BED_minus.pdf \
        --averageType mean \
	    --perGroup \
        --plotTitle 'CAGE + on BED -' \
	    -y 'Normalized Read Count (mean)' \
	    --numPlotsPerRow 4 \
        --colors red violet blue

# plot heatmap
#plotHeatmap -m matrix_CAGE_plus_BED_minus \
#        -o heatmap_CAGE_plus_BED_minus.pdf \
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
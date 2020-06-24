#!/bin/bash

# computeMatrix of RNAseq + on BED -
computeMatrix scale-regions \
        --binSize 5 -p 4 -a 200 -b 200 -m 1000 --skipZeros --metagene --missingDataAsZero \
        --samplesLabel 'DM' 'LSM8' 'RRP4' 'WT' \
        -R genes_antisense_TC_up.bed.minus.txt \
        -S ../../ANALYSES_v2/rrp4_lsm8_novogene/bigwigs_RPM_R123_fixed/*_Forward_R123_fixed.bw \
        -o matrix_RNAseq_plus_BED_minus \
        --outFileNameMatrix matrix_RNAseq_plus_BED_minus.mat


# plot average profile
plotProfile -m matrix_RNAseq_plus_BED_minus \
        -o profiles_RNAseq_plus_BED_minus.pdf \
        --averageType median \
	    --perGroup \
        --plotTitle 'RNAseq + on BED -' \
        -y 'Normalized Read Count (median)' \
	    --numPlotsPerRow 4 \
        --colors grey violet red green

# plot heatmap
#plotHeatmap -m matrix_RNAseq_plus_BED_minus \
#        -o heatmap_RNAseq_plus_BED_minus.pdf \
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
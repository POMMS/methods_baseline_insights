#!/bin/bash
#SBATCH --job-name=qiime2_C13_core_picrust
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem=10000
#SBATCH --output=log/C13_core_picrust.%j.out
#SBATCH --nodes=1
#SBATCH --error=log/C13_core_picrust.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: C13_core_picrust-metrics.slurm
# alpha, beta diversity and PCoa analysis

echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# core analysis, alpha, beta diversity. 
echo -e "step C13: core-picrust..." 
 
 qiime diversity core-metrics \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --i-table $WKPATH/q2-picrust2_output/ec_metagenome.qza \
 --p-sampling-depth 5098896 \
 --output-dir $WKPATH/${PREFIX}core-picrust-ec 
 
  qiime diversity core-metrics \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --i-table $WKPATH/q2-picrust2_output/ko_metagenome.qza \
 --p-sampling-depth 7681730 \
 --output-dir $WKPATH/${PREFIX}core-picrust-ko
 
  qiime diversity core-metrics \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --i-table $WKPATH/q2-picrust2_output/pathway_abundance.qza \
 --p-sampling-depth 5098896 \
 --output-dir $WKPATH/${PREFIX}core-picrust-path-ab
 
echo -e "finished core..."
echo -e $(date)

# alpha diversity analysis, also can use the other alpha diversity analysis results for visualization
echo -e "step C13B: alpha diversity..." 
 
 qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/faith_pd_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-picrust-ec/faith-pd-group-significance.qzv

 qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/evenness_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-picrust-ec/evenness-group-significance.qzv

  qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/shannon_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-picrust-ec/shannon-group-significance.qzv

 qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/faith_pd_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
 --o-visualization $WKPATH/${PREFIX}core-picrust-ec/faith_pd_vector_correlation.qzv
 
  qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/faith_pd_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
--o-visualization $WKPATH/${PREFIX}core-picrust-ec/evenness_vector_correlation.qzv
 
  qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-picrust-ec/shannon_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
  --o-visualization $WKPATH/${PREFIX}core-picrust-ec/shannon_vector_correlation.qzv
 
echo -e "finished alpha diversity..."
echo -e $(date)

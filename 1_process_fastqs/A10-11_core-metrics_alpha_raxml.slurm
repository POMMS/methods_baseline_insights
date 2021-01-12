#!/bin/bash
#SBATCH --job-name=A10-11_core
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem=10000
#SBATCH --output=log/A10-11_core.%j.out
#SBATCH --nodes=1
#SBATCH --error=log/A10-11_core.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: 10A_core-metrics.slurm
# alpha, beta diversity and PCoa analysis

echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# core analysis, alpha, beta diversity. 
echo -e "step 10A: core-metrics-phylogenetic..." 
 
 qiime diversity core-metrics-phylogenetic \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --i-phylogeny $WKPATH/qza/${PREFIX}09A_seq-aligned-mask-rooted-raxml.nwk.qza  \
 --i-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
 --p-sampling-depth 12374 \
 --output-dir $WKPATH/${PREFIX}core-metrics-results 
 
echo -e "finished core-metrics..."
echo -e $(date)

# alpha diversity analysis, also can use the other alpha diversity analysis results for visualization
echo -e "step 11A: alpha diversity..." 
 
 qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/faith_pd_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-metrics-results/faith-pd-group-significance.qzv

 qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/evenness_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-metrics-results/evenness-group-significance.qzv

  qiime diversity alpha-group-significance \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/shannon_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/${PREFIX}core-metrics-results/shannon-group-significance.qzv

 qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/faith_pd_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
 --o-visualization $WKPATH/${PREFIX}core-metrics-results/faith_pd_vector_correlation.qzv
 
  qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/evenness_vector.qza \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
--o-visualization $WKPATH/${PREFIX}core-metrics-results/evenness_vector_correlation.qzv
 
  qiime diversity alpha-correlation \
 --i-alpha-diversity $WKPATH/${PREFIX}core-metrics-results/shannon_vector.qza  \
 --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --p-method pearson \
  --o-visualization $WKPATH/${PREFIX}core-metrics-results/shannon_vector_correlation.qzv
 
echo -e "finished alpha diversity..."
echo -e $(date)

echo -e "start beta diversity bioenv..."
qiime diversity bioenv \
--i-distance-matrix $WKPATH/${PREFIX}core-metrics-results/bray_curtis_distance_matrix.qza \
--m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
--o-visualization $WKPATH/${PREFIX}core-metrics-results/bray_curtis_distance_matrix_bioenv.qzv

qiime diversity bioenv \
--i-distance-matrix $WKPATH/${PREFIX}core-metrics-results/jaccard_distance_matrix.qza \
--m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
--o-visualization $WKPATH/${PREFIX}core-metrics-results/jaccard_distance_matrix_bioenv.qzv

qiime diversity bioenv \
--i-distance-matrix $WKPATH/${PREFIX}core-metrics-results/unweighted_unifrac_distance_matrix.qza \
--m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
--o-visualization $WKPATH/${PREFIX}core-metrics-results/unweighted_unifrac_distance_matrix_bioenv.qzv

qiime diversity bioenv \
--i-distance-matrix $WKPATH/${PREFIX}core-metrics-results/weighted_unifrac_distance_matrix.qza \
--m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
--o-visualization $WKPATH/${PREFIX}core-metrics-results/weighted_unifrac_distance_matrix_bioenv.qzv

echo -e "finished beta diversity bioenv..."
echo -e $(date)

# rarefaction, check sample depth & seq read depth sufficient or not.
echo -e "step 11A: alpha diversity, rarefaction..." 
qiime diversity alpha-rarefaction \
  --i-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
  --i-phylogeny $WKPATH/qza/${PREFIX}09A_seq-aligned-mask-rooted-raxml.nwk.qza \
  --p-max-depth 100000 \
  --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
  --o-visualization $WKPATH/${PREFIX}core-metrics-results/alpha-rarefaction.qzv

  echo -e "finished rarefaction..."
echo -e $(date)


echo -e "step 11B: exporting core analysis results..."
# beta diversity distance matrix
qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/bray_curtis_distance_matrix.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/distance-matrix.tsv $WKPATH/export/export_txt/${PREFIX}10A_bray_curtis_distance_matrix.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/jaccard_distance_matrix.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/distance-matrix.tsv $WKPATH/export/export_txt/${PREFIX}10A_jaccard_distance_matrix.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/distance-matrix.tsv $WKPATH/export/export_txt/${PREFIX}10A_unweighted_unifrac_distance_matrix.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/distance-matrix.tsv $WKPATH/export/export_txt/${PREFIX}10A_weighted_unifrac_distance_matrix.tsv 

#core analysis results
#beta diversity pcoa table
qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/bray_curtis_pcoa_results.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/ordination.txt $WKPATH/export/export_txt/${PREFIX}10A_bray_curtis_pcoa.txt 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/jaccard_pcoa_results.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/ordination.txt $WKPATH/export/export_txt/${PREFIX}10A_jaccard_pcoa.txt 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/ordination.txt $WKPATH/export/export_txt/${PREFIX}10A_unweighted_unifrac_pcoa.txt 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/ordination.txt $WKPATH/export/export_txt/${PREFIX}10A_weighted_unifrac_pcoa.txt 

#alpha diversity matrix
qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/evenness_vector.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/alpha-diversity.tsv $WKPATH/export/export_txt/${PREFIX}10A_evenness_vector.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/faith_pd_vector.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/alpha-diversity.tsv $WKPATH/export/export_txt/${PREFIX}10A_faith_pd_vector.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/observed_otus_vector.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/alpha-diversity.tsv $WKPATH/export/export_txt/${PREFIX}10A_observed_otus_vector.tsv 

qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/shannon_vector.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/alpha-diversity.tsv $WKPATH/export/export_txt/${PREFIX}10A_shannon_vector.tsv 


#feature table biom
qiime tools export \
  --input-path $WKPATH/${PREFIX}core-metrics-results/rarefied_table.qza \
  --output-path $WKPATH/export/export_biom

mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}10A_rarefied_table.biom   
biom convert -i $WKPATH/export/export_biom/${PREFIX}10A_rarefied_table.biom -o $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table.tsv --to-tsv
  
echo -e "finished exporting..."
echo -e $(date)

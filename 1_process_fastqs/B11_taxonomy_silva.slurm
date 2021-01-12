#!/bin/bash
#SBATCH --job-name=B11_taxonomy_silva
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/B11_taxonomy_silva.%j.out
#SBATCH -c 4
#SBATCH --error=log/B11_taxonomy_silva.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: B11_taxonomy_silva.slurm
# taxonomy & plot + table

echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# clasification, taxonomy
echo -e "step 6B: classify-sklearn ..." 
 qiime feature-classifier classify-sklearn \
  --i-classifier $REFPATH/silva-132-99-515-806-nb-classifier.qza \
  --i-reads $WKPATH/qza/${PREFIX}04_seq.fasta.qza \
  --o-classification $WKPATH/qza/${PREFIX}06B_seq-taxonomy_silva.tsv.qza 

echo -e "finished classify-sklearn..."
echo -e $(date)

echo -e "step 7B: tabulate..."   

# taxonomy table
qiime metadata tabulate \
  --m-input-file $WKPATH/qza/${PREFIX}06B_seq-taxonomy_silva.tsv.qza \
  --o-visualization $WKPATH/qzv/${PREFIX}07B_seq-taxonomy_silva.qzv

echo -e "finished tabulate..."
echo -e $(date)

echo -e "step 8B: taxa barplot..."    

# taxonomy barplot
qiime taxa barplot \
  --i-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
  --i-taxonomy $WKPATH/qza/${PREFIX}06B_seq-taxonomy_silva.tsv.qza \
  --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
  --o-visualization $WKPATH/qzv/${PREFIX}08B_seq_taxa-barplot_silva.qzv 
  
qiime taxa barplot \
  --i-table $WKPATH/${PREFIX}core-metrics-results/rarefied_table.qza \
  --i-taxonomy $WKPATH/qza/${PREFIX}06B_seq-taxonomy_silva.tsv.qza \
  --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
  --o-visualization $WKPATH/qzv/${PREFIX}11B_seq-rarefied_taxa-barplot_silva.qzv   

echo -e "finished taxa barplot..."
echo -e $(date)

# export files
echo -e "exporting taxonomy files..."
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}06B_seq-taxonomy_silva.tsv.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/taxonomy.tsv $WKPATH/export/export_txt/${PREFIX}06B_seq-taxonomy_silva.tsv  

# merge feature/otu & rarefied table with taxonomy info
#get the header from feature table
head -2 $WKPATH/export/export_txt/${PREFIX}05_filterSample-table.tsv > $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-taxon.tsv

#join the feature table and taxonomy, sort first, then join, append to the output file
join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}05_filterSample-table.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}06B_seq-taxonomy_silva.tsv) >> $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-taxon.tsv

#get the header from feature table
head -2 $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table.tsv > $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon.tsv

#join the feature table and taxonomy, sort first, then join, append to the output file
join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}06B_seq-taxonomy_silva.tsv) >> $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon.tsv

echo -e "finished exporting 06B_seq-taxonomy_silva.tsv 05_filterSample-table-taxon.tsv 10A_rarefied_table-taxon.tsv ..."
echo -e $(date)

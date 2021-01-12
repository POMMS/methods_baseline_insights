#!/bin/bash
#SBATCH --job-name=05_featureTabl
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem=10000
#SBATCH --output=log/05_featureTabl.%j.out
#SBATCH --nodes=1
#SBATCH --error=log/05_featureTabl.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: 05_featureTable.slurm
# filter samples and summarize statistics

# record time and method name, env vars
echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# filter samples
echo -e "step 5: filter samples..."

qiime feature-table filter-samples \
  --i-table $WKPATH/qza/${PREFIX}04_table.biom.qza \
  --m-metadata-file $WKPATH/meta/${PREFIX}map.txt \
  --p-where "samplefilter='keep'" \
  --p-min-frequency 5000 \
  --o-filtered-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
  --verbose
 
echo -e "finished filter samples."  
echo -e $(date)

# feature is the unique sequence. summarlize unique sequences
echo -e "step 5b: starting feature-table summarizing...\n"

qiime feature-table summarize \
 --i-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
 --m-sample-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/qzv/${PREFIX}05_filterSample-table.qzv

echo -e "finished feature-table summarize." 
echo -e $(date)

echo -e "step 5c: starting feature-table tabulate-seqs...\n"
qiime feature-table tabulate-seqs \
 --i-data $WKPATH/qza/${PREFIX}04_seq.fasta.qza \
 --o-visualization $WKPATH/qzv/${PREFIX}04_seq.qzv
 
echo -e "finished feature-table tabulate-seqs." 
echo -e $(date)

echo -e "step 5d: starting export 05_filterSample-table.biom to tsv format\n"
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}05_filterSample-table.biom 
biom convert -i $WKPATH/export/export_biom/${PREFIX}05_filterSample-table.biom -o $WKPATH/export/export_txt/${PREFIX}05_filterSample-table.tsv --to-tsv  

echo -e "finished export 05_filterSample-table.biom to tsv format." 
echo -e $(date)


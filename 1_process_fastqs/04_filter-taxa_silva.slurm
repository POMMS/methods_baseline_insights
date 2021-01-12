#!/bin/bash
#SBATCH --job-name=04_filter-taxa_silva
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/04_filter-taxa_silva.%j.out
#SBATCH -c 4
#SBATCH --error=log/04_filter-taxa_silva.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: 04_filter-taxa_silva.slurm
# taxonomy & plot + table

echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# clasification, taxonomy
echo -e "step 04A: classify-sklearn ..."
 qiime feature-classifier classify-sklearn \
  --i-classifier $REFPATH/silva-132-99-515-806-nb-classifier.qza \
  --i-reads $WKPATH/qza/${PREFIX}03_seq.fasta.qza \
  --o-classification $WKPATH/qza/${PREFIX}04_seq-taxonomy_silva.tsv.qza

echo -e "finished classify-sklearn..."
echo -e $(date)

echo -e "step 04B2: filter based on taxanomy, remove mitochondria & chloroplast seqs ..."   
# filter based on taxanomy, remove mitochondria & chloroplast DNAs
qiime taxa filter-table \
  --i-table $WKPATH/qza/${PREFIX}03_table.biom.qza \
  --i-taxonomy $WKPATH/qza/${PREFIX}04_seq-taxonomy_silva.tsv.qza \
  --p-exclude mitochondria,chloroplast,unassigned,eukaryota \
  --o-filtered-table $WKPATH/qza/${PREFIX}04_table.biom.qza

qiime taxa filter-seqs \
  --i-sequences $WKPATH/qza/${PREFIX}03_seq.fasta.qza \
  --i-taxonomy $WKPATH/qza/${PREFIX}04_seq-taxonomy_silva.tsv.qza \
  --p-exclude mitochondria,chloroplast,unassigned,eukaryota \
  --o-filtered-sequences $WKPATH/qza/${PREFIX}04_seq.fasta.qza 
echo -e "finished filter based on taxanomy..."

# exporting
echo -e "step 4C: starting export 04_pair_dada2_tbl.tsv, 04_pair_dada2_seq.fasta ..."
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}04_table.biom.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}04_pair_dada2_tbl.biom
biom convert -i $WKPATH/export/export_biom/${PREFIX}04_pair_dada2_tbl.biom -o $WKPATH/export/export_txt/${PREFIX}04_pair_dada2_tbl.tsv --to-tsv  

qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}04_seq.fasta.qza \
  --output-path $WKPATH/export/export_fast
mv $WKPATH/export/export_fast/dna-sequences.fasta $WKPATH/export/export_fast/${PREFIX}04_pair_dada2_seq.fasta 

echo -e "finished export 04_dada2_stats.tsv, 04_pair_dada2_tbl.tsv, 04_pair_dada2_seq.fasta ..."
echo -e $(date)

echo -e "exporting taxonomy files..."
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}04_seq-taxonomy_silva.tsv.qza \
  --output-path $WKPATH/export/export_txt

mv $WKPATH/export/export_txt/taxonomy.tsv $WKPATH/export/export_txt/${PREFIX}04_seq-taxonomy_silva.tsv  

# merge feature/otu & rarefied table with taxonomy info
#get the header from feature table
head -2 $WKPATH/export/export_txt/${PREFIX}04_pair_dada2_tbl.tsv > $WKPATH/export/export_txt/${PREFIX}04_pair_dada2_tbl-taxon.tsv

# join the feature table and taxonomy, sort first, then join, append to the output file
join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}04_pair_dada2_tbl.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}04_seq-taxonomy_silva.tsv) >> $WKPATH/export/export_txt/${PREFIX}04_pair_dada2_tbl-taxon.tsv


echo -e "finished exporting 04_seq-taxonomy_silva.tsv 04_pair_dada2_tbl-taxon.tsv ..."
echo -e $(date)

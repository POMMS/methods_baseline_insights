#!/bin/bash
#SBATCH --job-name=03_dada2
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/03_dada2.%j.out
#SBATCH -c 4
#SBATCH --error=log/03_dada2.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: 03_pair_dada2.slurm

# record time and method name, env vars
echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX
 
# QC filtering, trim sequencing data. 
echo -e "step 4: starting dada2 denoising..."

qiime dada2 denoise-paired \
 --i-demultiplexed-seqs $WKPATH/qza/${PREFIX}02_pair-demux.fastq.qza \
 --p-trunc-len-f 230 \
 --p-trunc-len-r 210 \
 --p-trim-left-f 25 \
 --p-trim-left-r 25 \
 --o-denoising-stats $WKPATH/qza/${PREFIX}03_dada2_stats.tsv.qza \
 --o-table $WKPATH/qza/${PREFIX}03_table.biom.qza \
 --o-representative-sequences $WKPATH/qza/${PREFIX}03_seq.fasta.qza \
 --p-n-threads 4 \
 --verbose
echo -e "finished dada2 denoise."
echo -e $(date)

echo -e "step 4b: starting export 03_dada2_stats.tsv, 03_pair_dada2_tbl.tsv, 03_pair_dada2_seq.fasta ..."
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}03_dada2_stats.tsv.qza \
  --output-path $WKPATH/export/export_txt
mv $WKPATH/export/export_txt/stats.tsv $WKPATH/export/export_txt/${PREFIX}03_dada2_stats.tsv    

qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}03_table.biom.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}03_pair_dada2_tbl.biom
biom convert -i $WKPATH/export/export_biom/${PREFIX}03_pair_dada2_tbl.biom -o $WKPATH/export/export_txt/${PREFIX}03_pair_dada2_tbl.tsv --to-tsv  

qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}03_seq.fasta.qza \
  --output-path $WKPATH/export/export_fast
mv $WKPATH/export/export_fast/dna-sequences.fasta $WKPATH/export/export_fast/${PREFIX}03_pair_dada2_seq.fasta 

echo -e "finished export 03_dada2_stats.tsv, 03_pair_dada2_tbl.tsv, 03_pair_dada2_seq.fasta ..."
echo -e $(date)

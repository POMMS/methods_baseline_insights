#!/bin/bash
#SBATCH --job-name=01-02_pair-import-demux
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/01-02_pair-import-demux.%j.out
#SBATCH -c 4
#SBATCH --error=log/01-02_pair-import-demux.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: 01-02_pair_import-demux.slurm
# import and demultiplex EMP pair end sequences.

# record time and method name, env vars
echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# importing data
# [fastq.gz]; [EMPPairedEndSequences]
echo -e "step 1: importing paired end data..."

qiime tools import \
 --type EMPPairedEndSequences \
 --input-path $WKPATH/data \
 --output-path $WKPATH/qza/${PREFIX}01_EMP-pair-seqs.qza
 
echo -e "finished importing."
echo -e $(date)

# demultiplex
# [EMPPairedEndSequences];[PairedEndSequenceWithQuality]
echo -e "step 2: starting demux..."

qiime demux emp-paired \
 --i-seqs $WKPATH/qza/${PREFIX}01_EMP-pair-seqs.qza \
 --m-barcodes-file $WKPATH/meta/${PREFIX}map.txt \
 --m-barcodes-column BarcodeSequence \
 --o-per-sample-sequences $WKPATH/qza/${PREFIX}02_pair-demux.fastq.qza \
 --p-rev-comp-mapping-barcodes \
 --verbose
 
echo -e "finished demux"
echo -e $(date)

# summarizing & visualization based on demultiplexing results.
echo -e "step 2b: starting demux summarizing..."

qiime demux summarize \
 --i-data $WKPATH/qza/${PREFIX}02_pair-demux.fastq.qza \
 --o-visualization $WKPATH/qzv/${PREFIX}02_pair-demux.qzv
 
echo -e "finished demux summarizing."
echo -e $(date)

echo -e "step 2c: exporting 02_pair-demux.fastq.qza to fastq format..."
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}02_pair-demux.fastq.qza \
  --output-path $WKPATH/export/export_fast   
echo -e "finished exporting 02_pair-demux.fastq.qza to fastq format."
echo -e $(date)

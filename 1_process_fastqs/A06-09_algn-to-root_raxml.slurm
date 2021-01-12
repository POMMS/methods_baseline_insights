#!/bin/bash
#SBATCH --job-name=A06-A09_align-root
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/A06-A09_align-root.%j.out
#SBATCH --nodes=1
#SBATCH -c 4
#SBATCH --error=log/A06-A09_align-root.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: A06-A09_algn-to-root.slurm
# align and build phylogeny tree

echo -e $(date) 
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX

# align seqs
echo -e "step 6A: mafft..."

qiime alignment mafft \
 --i-sequences $WKPATH/qza/${PREFIX}04_seq.fasta.qza \
 --o-alignment $WKPATH/qza/${PREFIX}06A_seq-aligned.fasta.qza
 
echo -e "finished mafft..."
echo -e $(date)
echo -e "step 7A: mask..."
 
 qiime alignment mask \
 --i-alignment $WKPATH/qza/${PREFIX}06A_seq-aligned.fasta.qza \
 --o-masked-alignment $WKPATH/qza/${PREFIX}07A_seq-aligned-mask.fasta.qza 

echo -e "finished mask..."
echo -e $(date)

#build trees
echo -e "step 8A: raxml..."
 
 qiime phylogeny raxml \
 --i-alignment $WKPATH/qza/${PREFIX}07A_seq-aligned-mask.fasta.qza  \
 --p-seed 58853 \
 --p-n-searches 10 \
 --p-n-threads 4 \
 --o-tree $WKPATH/qza/${PREFIX}08A_seq-aligned-mask-unrooted-raxml.nwk.qza 

echo -e "finished raxml..."
echo -e $(date)
echo -e "step 9A: midpoint-root..."  
 
 qiime phylogeny midpoint-root \
 --i-tree $WKPATH/qza/${PREFIX}08A_seq-aligned-mask-unrooted-raxml.nwk.qza \
 --o-rooted-tree $WKPATH/qza/${PREFIX}09A_seq-aligned-mask-rooted-raxml.nwk.qza 

echo -e "finished midpoint-root..."
echo -e $(date)

echo -e "step 9B:exporting align rooted fasta and nwk files..." 
qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}06A_seq-aligned.fasta.qza \
  --output-path $WKPATH/export/export_fast

mv $WKPATH/export/export_fast/aligned-dna-sequences.fasta $WKPATH/export/export_fast/${PREFIX}06A_seq-aligned.fasta  

qiime tools export \
  --input-path $WKPATH/qza/${PREFIX}09A_seq-aligned-mask-rooted-raxml.nwk.qza \
  --output-path $WKPATH/export/export_nwk

mv $WKPATH/export/export_nwk/tree.nwk $WKPATH/export/export_nwk/${PREFIX}09A_seq-aligned-mask-rooted-raxml.nwk  

echo -e "finished exporting fasta and nwk files..."
echo -e $(date)

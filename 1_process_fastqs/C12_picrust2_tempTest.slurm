#!/bin/bash
#SBATCH --job-name=C12_picrust2
#SBATCH --partition=all
#SBATCH --time=UNLIMITED
#SBATCH --mem-per-cpu=25000
#SBATCH --output=log/C12_picrust2.%j.out
#SBATCH -c 4
#SBATCH --error=log/C12_picrust2.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zwei@duke.edu

# file: C12_picrust2e.slurm
# picrust

# record time and method name, env vars
echo -e $(date)
echo -e "working path: " $WKPATH
echo -e "prefix: " $PREFIX
 
echo -e "step C12B: picrust..."

# picrust2
qiime picrust2 custom-tree-pipeline \
  --i-table $WKPATH/qza/${PREFIX}05_filterSample-table.biom.qza \
  --i-tree $WKPATH/placed_out_for_picrust/tree.qza \
  --output-dir $WKPATH/q2-picrust2_output \
  --p-threads 4 \
  --p-hsp-method mp \
  --p-max-nsti 2
  
echo -e "finished picrust."  
echo -e $(date)

# summarizing visualizations
echo -e "step C12C: starting feature-table summarizing...\n"

qiime feature-table summarize \
 --i-table $WKPATH/q2-picrust2_output/ec_metagenome.qza \
 --m-sample-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/qzv/${PREFIX}C12_picrust_ec.qzv

qiime feature-table summarize \
 --i-table $WKPATH/q2-picrust2_output/ko_metagenome.qza \
 --m-sample-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/qzv/${PREFIX}C12_picrust_ko.qzv
  
qiime feature-table summarize \
 --i-table $WKPATH/q2-picrust2_output/pathway_abundance.qza \
 --m-sample-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/qzv/${PREFIX}C12_picrust_pathway-ab.qzv 
 
qiime feature-table summarize \
 --i-table $WKPATH/q2-picrust2_output/pathway_coverage.qza \
 --m-sample-metadata-file $WKPATH/meta/${PREFIX}map.txt \
 --o-visualization $WKPATH/qzv/${PREFIX}C12_picrust_pathway-co.qzv  
 
echo -e "finished feature-table summarize." 
echo -e $(date)


echo -e "step C12D: starting export picrust results to biom & tsv format\n"
qiime tools export \
  --input-path $WKPATH/q2-picrust2_output/ec_metagenome.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}C12_ec_metagenome.biom 
biom convert -i $WKPATH/export/export_biom/${PREFIX}C12_ec_metagenome.biom -o $WKPATH/export/export_txt/${PREFIX}C12_ec_metagenome.tsv --to-tsv  

qiime tools export \
  --input-path $WKPATH/q2-picrust2_output/ko_metagenome.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}C12_ko_metagenome.biom 
biom convert -i $WKPATH/export/export_biom/${PREFIX}C12_ko_metagenomee.biom -o $WKPATH/export/export_txt/${PREFIX}C12_ko_metagenome.tsv --to-tsv  

qiime tools export \
  --input-path $WKPATH/q2-picrust2_output/pathway_abundance.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}C12_pathway_abundance.biom 
biom convert -i $WKPATH/export/export_biom/${PREFIX}C12_pathway_abundance.biom -o $WKPATH/export/export_txt/${PREFIX}C12_pathway_abundance.tsv --to-tsv

qiime tools export \
  --input-path $WKPATH/q2-picrust2_output/pathway_coveragee.qza \
  --output-path $WKPATH/export/export_biom
mv $WKPATH/export/export_biom/feature-table.biom $WKPATH/export/export_biom/${PREFIX}C12_pathway_coverage.biom 
biom convert -i $WKPATH/export/export_biom/${PREFIX}C12_pathway_coverage.biom -o $WKPATH/export/export_txt/${PREFIX}C12_pathway_coverage.tsv --to-tsv  

echo -e "finished export picrust files to tsv format." 
echo -e $(date)


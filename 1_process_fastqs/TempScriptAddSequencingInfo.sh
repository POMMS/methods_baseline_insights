cat $WKPATH/export/export_fast/${PREFIX}03_pair_dada2_seq.fasta | tr '\n' '\t' | tr '>' '\n' > $WKPATH/export/export_txt/${PREFIX}03_pair_dada2_seq.tsv

head -2 $WKPATH/export/export_txt/${PREFIX}05_filterSample-table.tsv > $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq.tsv

join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}05_filterSample-table.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}03_pair_dada2_seq.tsv) >> $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq.tsv

head -2 $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq.tsv > $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq-taxon.tsv

join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}06B_seq-taxonomy_silva.tsv) >> $WKPATH/export/export_txt/${PREFIX}05_filterSample-table-seq-taxon.tsv

head -2 $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon.tsv > $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon-seq.tsv


join -t $'\t' <(sort  $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon.tsv) <(sort  $WKPATH/export/export_txt/${PREFIX}03_pair_dada2_seq.tsv) >> $WKPATH/export/export_txt/${PREFIX}10A_rarefied_table-taxon-seq.tsv

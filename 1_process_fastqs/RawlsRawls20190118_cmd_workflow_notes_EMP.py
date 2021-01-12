# file:  cmd_workflow_notes_EMP.py
# this is a fake py file, the reason to use py is to take advantage of the color coding from editor. This is just a file to take notes & tracking everything.

# Part 1: Run Qiime2 on Hardec
# Info to collect 
# a.	Library protocol, EMP or not? (515-806, multiple gg reference seq to choose for taxonomy) 
# b.	Reverse complement barcodes or not? (for demux) 
# c.	Read length (for trimming, dada2) 
# d.	Reverse complemented mapping barcodes? (for demux) 
# e.	Mapping files (prefix+map.txt)
#           1) run LeftJoinTablesOnFirstCol.py to merge mapping txt file with experiment design tables.
#                a) save files needed to be merged into one folder, first column contains primary IDs for merging.
#                b) run python script to merge file, when prompt enter folder path.
#                c) check root folder for merged file.
#   34, treatment = NA, Holly lookedup in RedCap, changed NA to lifestyle.
#           2) check experiment design balanceness with TableStats.py, when prompt enter input file path, output csv file path and export pdf file path
#           3) add a group + sample name column, will be useful for graphing purpose.
#           4) https://docs.google.com/spreadsheets -> newSheet -> copy and paste content -> Addon-> Keemi -> validate Qiime2 meta : check validaty of the meta file

# Login to Hardec [hardac-login.genome.duke.edu]

## Part 2: Make folders with 00_initialize.sh: 
newgrp omicscore
umask 007
# to make folder structures for analysis
# download from gitlab and sftp 00_initialization.sh to /data/omicscore/projectFolder
cd /data/omicscore/Rawls-Rawls-20190118
sh 00_initialization.sh

# a.	script, script/log : all slurms store under script folder, output from script will be stored in log folder. 
# b.	rawData/ stores original fastq.gz files; data/"barcodes.fastq.gz", "forward.fastq.gz" and "reverse.fastq.gz" : data folder store raw data, exact filenames. folder should contain only those 3 files.
# c.	meta/${PREFIX}map.txt : store mapping file, sample filter column can be used to filter samples downstream (e.g. separate samples from different projects)
# 		i.	#SampleID, BarcodeSequence, LinkerPrimerSequence, condition/groups, platform (MiSeq150PE), libProtocol(EMP16s_515_806), seqOrder(#), samplefilter (to be added later to filter out samples) 
# 		ii.	#q2:types, numeric, categorical
#       iii. Tab separated 
# d.	qza, qzv : qza will store all data file from Qiime2 methods; qzv folder will store all visualization file.
# e. 	export/: export-txt, export-fastq, export-nwk, export-biom


##Part 3A:  Sftp sequencing data from dnaseqcore server
# go to data folder where the fastq files need to be saved.
# copy the raw data into folder from duke box
cd /data/omicscore/Rawls-Rawls-20190118/data/

# link R1, R2 and I1 to "../data/barcodes.fastq.gz", "forward.fastq.gz" and "reverse.fastq.gz"
ln -s /data/omicscore/Rawls-Rawls-20190118/rawData/4851-P1_S1_L001_R1_001.fastq.gz forward.fastq.gz
ln -s /data/omicscore/Rawls-Rawls-20190118/rawData/4851-P1_S1_L001_R2_001.fastq.gz reverse.fastq.gz
ln -s /data/omicscore/Rawls-Rawls-20190118/rawData/4851-P1_S1_L001_I1_001.fastq.gz barcodes.fastq.gz

##Part 3B: copy mapping file to meta file folder
sftp *map.txt /data/omicscore/projectFolder/meta/

##Part 4: lines to run at the beginning of each run
# Get files ready, data file, map file
newgrp omicscore
umask 007
# go to script folder
cd /data/omicscore/Rawls-Rawls-20190118/script
# load current version of qiime2
source /data/common/qiime2/miniconda/bin/activate qiime2-2018.11
# start sbatch runs
srun -p interactive --pty --mem 4096 /bin/bash 

# define file paths
export WKPATH="/data/omicscore/Rawls-Rawls-20190118"
export PREFIX="Rawls-Rawls-20190118-R24_"
export REFPATH="/data/omicscore/Qiime/reference"

# run sbatch scripts
sbatch 01-02_pair_import-demux.slurm
# check 02_pair-demux.qzv
# 1) sample total # of reads
# 2) parameters for dada2, trim left and right for forward and revesre read
	
sbatch 03_pair_dada2.slurmA

sbatch 04_filter-taxa_silva.slurm

sbatch 05_featureTable.slurm

sbatch A06-09_algn-to-root_raxml.slurm 

sbatch A10-11_core-metrics_alpha_raxml.slurm
# check txt files
# check visualization files

sbatch B11_taxonomy_silva.slurm

(qiime2-2018.11) [zwei@x2-01-1 reference]$ wget http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.fna.qza
--2019-03-15 11:11:15--  http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.fna.qza
Resolving kronos.pharmacology.dal.ca... 129.173.80.34
Connecting to kronos.pharmacology.dal.ca|129.173.80.34|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 4408600 (4.2M)
Saving to: âreference.fna.qzaâ

reference.fna.qza                            100%[============================================================================================>]   4.20M  4.16MB/s    in 1.0s

2019-03-15 11:11:16 (4.16 MB/s) - âreference.fna.qzaâ saved [4408600/4408600]

(qiime2-2018.11) [zwei@x2-01-1 reference]$ wget http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.tre.qza
--2019-03-15 11:11:36--  http://kronos.pharmacology.dal.ca/public_files/tutorial_datasets/picrust2_tutorial_files/reference.tre.qza
Resolving kronos.pharmacology.dal.ca... 129.173.80.34
Connecting to kronos.pharmacology.dal.ca|129.173.80.34|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 191187 (187K)
Saving to: âreference.tre.qzaâ

reference.tre.qza                            100%[============================================================================================>] 186.71K   773KB/s    in 0.2s

2019-03-15 11:11:37 (773 KB/s) - âreference.tre.qzaâ saved [191187/191187]

(qiime2-2018.11) [zwei@x2-01-1 reference]$ qiime fragment-insertion sepp   --i-representative-sequences $WKPATH/qza/${PREFIX}04_seq.fasta.qza   --p-threads 1   --i-reference-alignment $REFPATH/reference.fna.qza   --i-reference-phylogeny $REFPATH/reference.tre.qza   --output-dir $WKPATH/placed_out_for_picrust
Saved Phylogeny[Rooted] to: /data/omicscore/Rawls-Rawls-20190118/placed_out_for_picrust/tree.qza
Saved Placements to: /data/omicscore/Rawls-Rawls-20190118/placed_out_for_picrust/placements.qza

(qiime2-2018.11) [zwei@x2-01-1 script]$ sbatch C12_picrust2_tempTest.slurm

https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial
https://github.com/picrust/picrust2/wiki
https://picrust.github.io/picrust/tutorials/algorithm_description.html

[zwei@hardac-login omicscore]$ srun -p interactive --mem=64G --pty bash
[zwei@x2-01-1 omicscore]$ module load ddsclient
[zwei@x2-01-1 omicscore]$ ddsclient upload -p R24-16S-Results /data/omicscore/R24-16S-Results/

[zwei@x2-01-1 omicscore]$ ddsclient upload -p R24-16S /data/omicscore/R24-16S/
Uploading 0 projects, 0 folders, 1 file.
Progress: 100% - sending R24-16S-20190118__Microbiome_Data_Analysis_Report-ResulDone: 100%


Upload Report for Project: 'R24-16S-Results' 2019-03-26 20:26:16.309123

SENT FILENAME                                                                                              ID                                      SIZE     HASH
/data/omicscore/R24-16S-Results/20190118/R24-16S-20190118__Microbiome_Data_Analysis_Report-Results.docx    4ffb51ab-ae5d-4c98-9aca-27ba8dd97079    65086    0e916dd69ef8e8ccf27a8da1277f50e3


URL to view project: https://dataservice.duke.edu/#/project/dbd4269d-97da-44cf-875e-47f6e44a56b4

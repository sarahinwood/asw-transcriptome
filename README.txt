##to set up directory for assembly, make sure the below has been set up/run in directory:

##make sure trinotate dependencies linked to directory
mkdir bin
ln -s /Volumes/archive/deardenlab/tomharrop/data/trinotate_db/20200203 bin/trinotate_db

##link RNAseq files into data/reads as follows:
mkdir data
mkdir reads
data/reads/FLOWCELLID/FLOWCELLID_myrnaseqsample.fastq.gz

so for CA7HAANXX-2125-04-01-1_S22_L003_R1_001.fastq.gz
data/reads/CA7HAANXX/CA7HAANXX-2125-04-01-1_S22_L003_R1_001.fastq.gz

##make sure a sample key file is present here:
data/sample_key.csv

##sample_key.csv should have the following columns with rows for every sample that should be included in assembly:
OGBF_sample_id (for above example file = 2125_04)
Sample_name (ASW_Abdomen)
Flow_cell (CA7HAANXX)

##download transrate if desired into bin/transrate (not a necessary step)



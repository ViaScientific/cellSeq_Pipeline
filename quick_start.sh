#!/bin/sh

# If you're running for the first time, you need to download genome_data as follows:

# Please replace /path/to/save with your directory and execute following commands:
genomeDir="/path/to/save"
wget https://web.dolphinnext.com/umw_biocore/dnext_data/genome_data/c_elegans/ws279/ -P ${genomeDir}/genome_data/c_elegans/ws279/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 

# download indrop pipeline repository
git -C ${genomeDir} clone https://github.com/ViaScientific/cellSeq_Pipeline.git && export NXF_VER=22.10.7

## Edit nextflow.config file: ##
# Check nextflow config parameters for your run environment. (https://www.nextflow.io/docs/latest/config.html). For example, you can add following parameters for LSF scheduler.
# process.executor = 'lsf'
# process.time = '240m'
# process.cpus = 1
# process.queue = 'short'
# process.memory = '32 GB'


# example run for paired reads
nextflow ${genomeDir}/cellSeq_Pipeline/main.nf  -profile singularity \
--DOWNDIR ${genomeDir} \
--reads '*_{R1,R2}.fastq.gz' \
--cellBarcodeFile  's3://viafoundry/run_data/genome_data_other/CelSeq/bcSet_full.txt'
--mate 'pair' \
--genome_build 'c_elegans_ws279_refseq_PRJNA13758' \
--run_STAR = 'yes' \
--run_IGV_TDF_Conversion  'no' \
--run_RSeQC  'yes' \
--run_BigWig_Conversion  'no' \
--run_Picard_CollectMultipleMetrics  'yes'
--run_Trimmer  'no'
--run_Quality_Filtering  'no'
--mate_split  'single'
--pdfbox_path  '/usr/local/bin/dolphin-tools/pdfbox-app-2.0.0-RC2.jar'
--gtf2bed_path  '/usr/local/bin/dolphin-tools/gtf2bed'
--run_Adapter_Removal  'yes'
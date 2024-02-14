# ViaScientific/cellSeq_Pipeline: Usage

## Introduction
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow main.nf -profile docker --DOWNDIR /path/to/save/genome_data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build c_elegans_ws279_refseq_PRJNA13758
```

If you're running for the first time, you need to download genome_data as follows:

```bash
# Please replace /path/to/save with your directory and execute following commands:

genomeDir="/path/to/save"
wget https://web.dolphinnext.com/umw_biocore/dnext_data/genome_data/c_elegans/ws279/ -P ${genomeDir}/genome_data/c_elegans/ws279/ -l inf -nc -nH --cut-dirs=4 -r --no-parent -R "index.html*" 

nextflow main.nf -profile docker --DOWNDIR ${genomeDir}/genome_data --reads '*_R{1,2}.fastq.gz' --mate 'pair' --genome_build c_elegans_ws279_refseq_PRJNA13758 
```


This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull ViaScientific/cellSeq_Pipeline
```

## Main arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker, test` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls containers from Quay.io
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from Quay.io


### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{R1,R2}.fastq' --mate 'pair'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired read data, the path must use `{R1,R2}` notation to specify reads.


### `--mate`
Two options (pair) available for `--mate` parameter. For example:

```bash
--reads 'path/to/data/sample_*_{R1,R2}.fastq' --mate 'pair'
```


## Reference genomes

### `--genome_build` 
Currently only C. Elegans genome is supported for CellSeq pipeline. To run the pipeline, you must specify which to use with the `--genome_build` flag.

* C. Elegans
  * `--genome_build c_elegans_ws279_refseq_PRJNA13758`

Note: For new genome requests, please send e-mail to support@viascientific.com.


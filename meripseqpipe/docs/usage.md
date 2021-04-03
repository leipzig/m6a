# nf-core/meripseqpipe: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
- [nf-core/meripseqpipe: Usage](#nf-coremeripseqpipe-usage)
  - [Table of contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Running the pipeline](#running-the-pipeline)
    - [Updating the pipeline](#updating-the-pipeline)
    - [Reproducibility](#reproducibility)
  - [Main arguments](#main-arguments)
    - [`-profile`](#profile)
    - [`--reads`](#reads)
    - [`--singleEnd`](#singleend)
    - [`--gzip`](#gzip)
    - [`--stranded`](#stranded)
    - [`--mapq_cutoff [int]`](#mapqcutoff-int)
    - [`--designfile`](#designfile)
    - [`--comparefile`](#comparefile)
    - [`--aligners`, `--peakCalling_mode`, `--peakMerged_mode`, `--expression_analysis_mode`, `--methylation_analysis_mode`](#aligners---peakcallingmode---peakmergedmode---expressionanalysismode---methylationanalysismode)
  - [Reference genomes](#reference-genomes)
    - [`--genome` (using iGenomes)](#genome-using-igenomes)
    - [`--fasta`](#fasta)
    - [`--gtf`](#gtf)
    - [`--tophat2_index`, `hisat2_index`, `--bwa_index`, `--star_index`](#tophat2index-hisat2index---bwaindex---starindex)
    - [`--rRNA_fasta`](#rrnafasta)
    - [`--igenomesIgnore`](#igenomesignore)
  - [Skipping steps](#skipping-steps)
    - [`--skip_metpeak`, `--skip_macs2`, `--skip_matk`, `--skip_meyer`](#skipmetpeak---skipmacs2---skipmatk---skipmeyer)
    - [`--skip_fastp`, `--skip_fastqc`, `--skip_rseqc`, `--skip_createbedgraph`](#skipfastp---skipfastqc---skiprseqc---skipcreatebedgraph)
  - [Job resources](#job-resources)
    - [Automatic resubmission](#automatic-resubmission)
    - [Custom resource requests](#custom-resource-requests)
  - [AWS Batch specific parameters](#aws-batch-specific-parameters)
    - [`--awsqueue`](#awsqueue)
    - [`--awsregion`](#awsregion)
  - [Other command line parameters](#other-command-line-parameters)
    - [`--outdir`](#outdir)
    - [`--email`](#email)
    - [`--email_on_fail`](#emailonfail)
    - [`-name`](#name)
    - [`-resume`](#resume)
    - [`-c`](#c)
    - [`--custom_config_version`](#customconfigversion)
    - [`--custom_config_base`](#customconfigbase)
    - [`--max_memory`](#maxmemory)
    - [`--max_time`](#maxtime)
    - [`--max_cpus`](#maxcpus)
    - [`--plaintext_email`](#plaintextemail)
    - [`--monochrome_logs`](#monochromelogs)
    - [`--multiqc_config`](#multiqcconfig)
<!-- TOC END -->

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO nf-core: Document required command line parameters to run the pipeline-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/meripseqpipe --reads '*_R{1,2}.fastq.gz' --designfile 'path/to/designfile/design.csv' --comparefile 'path/to/designfile/compare.txt' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/meripseqpipe
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/meripseqpipe releases page](https://github.com/nf-core/meripseqpipe/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

- `conda`
  - A generic configuration profile to be used with [conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`nfcore/meripseqpipe`](http://hub.docker.com/r/nfcore/meripseqpipe/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`nfcore/meripseqpipe`](http://hub.docker.com/r/nfcore/meripseqpipe/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

<!-- TODO nf-core: Document required command line parameters -->

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--readPaths`. For example:

```bash
--singleEnd --readPaths 'path/to/data/'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--gzip`

By default, the pipeline expects uncompressed sequencing data( .fastq ). If you have compressed data( .fastq.gz ), you need to specify `--gzip` on the command line when you launch the pipeline. For example:

```bash
--gzip --readPaths 'path/to/data/'
```

### `--stranded`

By default, the pipeline non-strand-specific sequencing data. If you have strand-specific sequencing data data, you need to specify `--stranded [yes/no/reverse]` on the command line when you launch the pipeline. "yes" means foward strand-specific sequencing data while "reverse" means reverse strand-specific sequencing data. "no" means non-strand-specific sequencing data. For example:

```bash
--stranded yes --readPaths 'path/to/data/'
```

### `--mapq_cutoff [int]`

By default, the pipeline will keep multi mapping reads for PeakCalling. If you want to filter reads of mapping quality, you can specify `--mapq_cutoff [0-255]`. "255" means the pipeline will only keeps the unique mapping reads. For example:

```bash
--mapq_cutoff 255
```

### `--designfile`

Edit the nextflow.config and define "readPaths", "aligners", "designfile", "comparefile" and correspondiente alignment index for recommend.
Designfile is just like the following table with a comma (,) separated, which is .csv suffix file. It's recommended edited by Excel and save as .csv suffix file.

```bash
--singleEnd --readPaths 'path/to/data/' --designfile 'path/to/designfile/design.csv'
```

Example:
| Sample_ID | input_FileName | ip_FileName | Group |
| --- | --- | --- | --- |
| H1A_Endo | A | B | group_Endo |
| H1A_ES | C | D | group_ES |
| H1B_Endo | E | F | group_Endo |
| H1B_ES | G | H | group_ES |

>Tips
>
>1. A, B, C... mean the filenames of data, just like A.fastq.gz.
>2. If your data is .fastq.gz suffix file, please add the parameter of gzip, just like "--gzip true".
>3. If your filename of data is "Hela_cell_input.fastq.gz", please write its filename as "Hela_cell_input".

### `--comparefile`

Comparefile is just like the following text which is a "\_vs\_" between two groups.

```bash
--singleEnd --readPaths 'path/to/data/' --designfile 'path/to/designfile/design.csv' --comparefile 'path/to/designfile/compare.txt'
```

Example:
>group_Endo_vs_group_ES

### `--aligners`, `--peakCalling_mode`, `--peakMerged_mode`, `--expression_analysis_mode`, `--methylation_analysis_mode`

```config
  // Setting main parameters of analysis mode
  aligners = "star"   // "star" OR "bwa" OR "tophat2" OR "hisat2" OR "none"
  peakCalling_mode = "independence" // "group" OR "independence"
  peakMerged_mode = "rank" // "rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc"
  expression_analysis_mode = "DESeq2" // "DESeq2" OR "edgeR" OR "none"
  methylation_analysis_mode = "QNB" // "MATK" OR "QNB" OR "Wilcox-test" OR "MeTDiff" OR "edgeR" OR "DESeq2"
```

If you prefer, you can specify the mode of meripseqpipe when you run the pipeline.

```bash
  --aligners star \
  --peakCalling_mode independence \
  --peakMerged_mode rank \
  --expression_analysis_mode DESeq2 \
  --methylation_analysis_mode QNB
```

## Reference genomes

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### `--genome` (using iGenomes)

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

- Human
  - `--genome GRCh37`
- Mouse
  - `--genome GRCm38`
- _Drosophila_
  - `--genome BDGP6`
- _S. cerevisiae_
  - `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

<!-- TODO nf-core: Update reference genome example according to what is needed -->

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

<!-- TODO nf-core: Describe reference path flags -->
### `--fasta`

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to Fasta reference]'
```

### `--gtf`

If you prefer, you can specify the full path to your reference genome annotation when you run the pipeline:

```bash
--gtf '[path to Annotation reference(.gtf)]'
```

### `--tophat2_index`, `hisat2_index`, `--bwa_index`, `--star_index`

If you prefer, you can specify the full path to your aligners index of reference genome when you run the pipeline. You only need to choose one for running pipeline which depend on `--aligner [star/bwa/hisat2/tophat2]`:

```bash
  --tophat2_index "/Path/to/Tophat2Index/*"
  --hisat2_index "/Path/to//Hisat2Index/*"
  --bwa_index "/Path/to/BWAIndex/*"
  --star_index "/Path/to/starindex"
```

### `--rRNA_fasta`

If you want to fiter the reads of rRNA, you need to offer the fasta of rRNA. For example:

```bash
--rRNA_fasta '[path to rRNA Fasta reference]'
```

### `--igenomesIgnore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

## Skipping steps

### `--skip_metpeak`, `--skip_macs2`, `--skip_matk`, `--skip_meyer`

The pipeline contains four tools for PeakCalling. Sometimes, it may not be desirable to run all of them.
The following options make this easy:

- `--skip_metpeak` -          Skip MeTPeak
- `--skip_macs2` -            Skip MACS2
- `--skip_matk` -             Skip MATK
- `--skip_meyer` -            Skip Meyer

### `--skip_fastp`, `--skip_fastqc`, `--skip_rseqc`, `--skip_createbedgraph`

The pipeline contains multiple tools for QC. Sometimes, it may not be desirable to run all of them.
The following options make this easy:

- `--skip_fastp` -             Skip MeTPeak
- `--skip_fastqc` -            Skip MACS2
- `--skip_rseqc` -             Skip MATK
- `--skip_createbedgraph` -    Skip Meyer

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use the `-awsbatch` profile and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region to run your job in. Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.

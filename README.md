# Transcriptome Mapping Pipeline
 https://github.com/palfalvi/rnaseq
 
 This pipeline is a basic workflow to QC and quantify raw RNA-seq reads with the following steps:
 
### 1, Raw fastq quality check 
Utilizing [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) software.

### 2, Index generation
This depends on the quantification mode selected. Currently [`salmon`](https://combine-lab.github.io/salmon/), [`kallisto`](https://pachterlab.github.io/kallisto/) and [`star`](https://github.com/alexdobin/STAR) are supported.

### 3, Mapping and quantification
This is connected to the index generation. 'salmon' and 'kallisto' natively quantifies reads, while the 'star' mode utilizes [`featureCounts`](http://subread.sourceforge.net/)


## Dependencies

As this is a nextflow project dependent on conda environments, the only 2 things you need before running are nextflow and conda:

### [Install nextflow](https://www.nextflow.io/docs/latest/getstarted.html)

In your `bin` directory, or in some places you can access from your `$PATH`, copy and execute the following:

```
wget -qO- https://get.nextflow.io | bash
```

Done. 

### [Install miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

Please download the miniconda installer into your `bin` directory form the [miniconda repository](https://docs.conda.io/en/latest/miniconda.html#linux-installers) and execute the follwing :

```
bash Miniconda3-latest-Linux-x86_64.sh
```

Done.

## How to run

First, you need to decide which mode you would like to run. The default and recommended way is `salmon`, but `kallisto` and `star` are also available. Just change specify with the `--mode` flag. e.g `--mode star`.


Then, you need to know if you have pair-end or single-end dataset. The pipeline runs by default on pair end data, but you can specify single end mode with the `--single` flag.

**If you run in single end mode with kallisto, please also specify `--fragment_length` and  `--fragment_sd`, which can be calculated from your BioAnalyzer file of the final libraries.**

Finally, you need to have a reference set. For `salmon` and `kallisto`, this should be a fasta file of transcripts (not genes!) specified with the `--transcriptome` flag. In the case of `star`, please specify the genome and a corresponding GTF annotation file with `--genome` and `--gtf` flags.

### Example commands

#### Salmon mapping with paire end reads.

```
nextflow run palfalvi/rnaseq --reads /path/to/reads/*R{1,2}.fastq --transcriptome transcripts.fasta
```

#### Kallisto with single end reads.

```
nextflow run palfalvi/rnaseq --mode kallisto --reads /path/to/reads/*.fastq --transcriptome transcripts.fasta --fragment_length 300 --fragment_sd 30 
```

#### STAR with pair end reads

```
nextflow run palfalvi/rnaseq --mode star --reads /path/to/reads/*R{1,2}.fastq --genome genome.fastq --gtf annotation.gtf
```





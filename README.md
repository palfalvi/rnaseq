# Transcriptome Mapping Pipeline
 https://github.com/palfalvi/rnaseq

 This pipeline is a basic workflow to QC and quantify raw RNA-seq reads with the following steps:

### 1, Raw fastq quality check and filtering
Utilizing [`fastp`](https://github.com/OpenGene/fastp) software, quality values are recorded and low quality reads removed. `fastp` also recognizes and trims adaptors.

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

* If it did not work on NIBB-BIAS5, please try on bias5-db.

### SRA accessions [currently not supported. Coming soon.]

If you wish to use SRA ids directly, you can provide with `--sra SRP043510` instead of `--reads`. This feature, however uses the [NCBI Esearch API](https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch), which requires an NCBI API Key in your environment. To get an API Key, follow [this link](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

After you get your API key, you need to set it into your environment. In your favorite editor (e.g. nano or vim) open `~/.bashrc` and in the end of the file, insert the following:
`export NCBI_API_KEY=0123456789abcdef` then save. You need to run `source ~/.bashrc` or log out and in again to make it work. From that point, you do not need to modify or rerun this part.

## How to run

First, you need to decide which mode you would like to run. The default and recommended way is `salmon`, but `kallisto` and `star` are also available. Just change specify with the `--mode` flag. e.g `--mode star`.


Then, you need to know if you have pair-end or single-end dataset. The pipeline runs by default on pair end data, but you can specify single end mode with the `--single` flag.

**If you run in single end mode with kallisto, please also specify `--fragment_length` and  `--fragment_sd`, which can be calculated from your BioAnalyzer file of the final libraries.**

Finally, you need to have a reference set. For `salmon` and `kallisto`, this should be a fasta file of transcripts (not genes!) specified with the `--transcriptome` flag. In the case of `star`, please specify the genome and a corresponding GTF annotation file with `--genome` and `--gtf` flags.

### Example commands

#### Salmon mapping with pair end reads.

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

#### Salmon with no QC and trimming and save the index files for later use. Also, save the results into `salmon_output` directory.

```
nextflow run palfalvi/rnaseq --reads /path/to/reads/*R{1,2}.fastq --transcriptome transcripts.fasta --save_index --skip_qc --out salmon_output
```

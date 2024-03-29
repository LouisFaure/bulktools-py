[![Install and Test](https://github.com/LouisFaure/bulktools-py/actions/workflows/test.yml/badge.svg)](https://github.com/LouisFaure/bulktools-py/actions/workflows/test.yml)

# bulktools-py

Analysis pipeline for bulk data. From fastq to anndata containing count matrix.

## Installation

```
conda create -n bulktools -c bioconda -c defaults star subread python=3.8 -y
conda activate bulktools
pip install git+https://github.com/LouisFaure/bulktools-py.git
```

## Usage

To work by default, the data structure should look like this:

```
├── fastq/ <-- the first element of the filename before '_' is the sample name
│   ├── sample1_S1_R1_L001.fastq.gz
│   └── sample2_S2_R1_L001.fastq.gz
├── genes.gtf <-- genome annotation
├── star_index/ <-- index folder created using STAR tool
│   ├── ...
```

### Command

```bash
bt -s star_index -g genes.gtf
```


## Full usage

```
usage: bt [-h] [--fq_path FQ_PATH] [--bam_path BAM_PATH] [--star_ref STAR_REF] [--gtf GTF]
          [--n_threads N_THREADS] [--adata_out ADATA_OUT]
          [cleanup]

Full bulk pipeline, from fastq to adata count matrix!
Performs the following: fastq -STAR-> bam -featureCounts-> anndata.h5ad

positional arguments:
  cleanup               remove temporary folders and files.

optional arguments:
  -h, --help            show this help message and exit
  --fq_path FQ_PATH, -f FQ_PATH
                        Path for input fastq files (relative, default: fastq).
  --bam_path BAM_PATH, -b BAM_PATH
                        Path for aligned BAMs (default: aligned).
  --star_ref STAR_REF, -s STAR_REF
                        STAR index path.
  --gtf GTF, -g GTF     GTF file path for featureCounts.
  --n_threads N_THREADS, -n N_THREADS
                        Total number of threads to use for both STAR and featureCounts.
  --adata_out ADATA_OUT, -o ADATA_OUT
                        Path for the adata output (relative, default: adata_bulk_star.h5ad).
```


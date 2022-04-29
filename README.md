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

```
usage: bt [-h] [--fq_path FQ_PATH] [--bam_path BAM_PATH] [--star_ref STAR_REF] [--gtf GTF] [--n_threads N_THREADS]
          [--mem MEM] [--adata_out ADATA_OUT]

Full bulk pipeline, from fastq to adata count matrix!
Performs the following: fastq -STAR-> bam -featureCounts-> anndata.h5ad

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
  --mem MEM, -m MEM     how much Gb to pass to --limitBAMsortRAM for STAR alignment (in Gb, default 10).
  --adata_out ADATA_OUT, -o ADATA_OUT
                        Path for the adata output (relative, default: adata_bulk_star.h5ad).
```

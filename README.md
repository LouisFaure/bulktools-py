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
usage: bt [-h] [--fq_path FQ_PATH] [--star_ref STAR_REF] [--gtf GTF] [--n_threads N_THREADS] [--adata_out ADATA_OUT]

Full bulk pipeline, from fastq to adata count matrix!
Performs the following: fastq -STAR-> bam -featureCounts-> anndata.h5ad

optional arguments:
  -h, --help            show this help message and exit
  --fq_path FQ_PATH     Path for input fastq files (relative, default: fastq).
  --star_ref STAR_REF   STAR index path.
  --gtf GTF             GTF file path for featureCounts.
  --n_threads N_THREADS, -n N_THREADS
                        number of threads per fastq file, for both STAR and featureCounts.
  --adata_out ADATA_OUT
                        Path for the adata output (relative, default: adata_bulk_star.h5ad).
```

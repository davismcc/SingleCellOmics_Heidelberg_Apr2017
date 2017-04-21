# Single Cell 'Omics: Analysis of single-cell methylation data

The final session of the course will cover pre-processing and basic analysis of
single-cell bisulfite sequencing data. We will assay two cell types, probably 16
cells in total.

## Goals
Two main goals:
1. Methylation profiles define cell type (i.e. cells will cluster apart by e.g. PCA)                                                             
2. Context specificity of methylation variance. E.g. in mouse ES cells, CGIs are homogenous (and low in methylation), repeat elements are homogenously high and active enhancer elements are heterogeneous. This is interesting because the enhancer elements are cell type specific and thus some variation in the methylation levels here implies plasticity in cell identity which could be important for lineage formation.

## First step

Clone or download this repository so that you have the necessary code, data and materials to hand.

If you're familiar with `git`:
```
git clone https://github.com/davismcc/SingleCellOmics_Heidelberg_Apr2017.git
```

If not, you can download a zip file of the repository by clicking the green "Clone or download" button above.

## Outline:

1. We will use `BISMARK` for alignments and methylation calling. For details, see this [protocol paper](http://www.nature.com/nprot/journal/v12/n3/full/nprot.2016.187.html).
2. QC (also see protocol paper)
    1.  Negative controls should not align
    1. bisulfite conversion efficiency (assessed using CHH methylation from bismark reports) should be >95%
    1. mapping efficiency (from bismark reports)  >10% (30-40% is normal here but may end up lower  in these practicals)
    1. number of CpG sites covered (I use 1M unique positions but this will depend on seq depth so maybe just exclude outliers)
3. Preprocessing and annotation
    *  Quantify methylation over regions of interest (promoters, gene bodies, enhancers, repeats, CpG islands).
        1. mean methylation rate (each covered position counts once – i.e. do not give extra weight to positions with >1 read)
        2. also record the coverage (number of CpG sites that were covered in the that cell at that locus) for the purpose of assigning weights to each cell in downstream analyses
4. Analysis
    1. Mean methylation by feature / cell type
    1. Variation by feature / cell type
    1. Dimension reduction    
    1. Clustering

We will manage the data processing and analysis "pipeline" using [snakemake](http://snakemake.readthedocs.io/en/stable/).

## Data

The aim will be for you to analyze the data you generate during the course in Heidelberg.

However, in case that data is unavailable for any reason and to have an alternative dataset that is processed and ready for analysis, we also have access to a small dataset from Stephen Clark and colleagues at the Babraham Institute, Cambridge. This dataset consists of 15 cells from mouse embryos.

1. Raw `fastq` files are available at this [link](https://www.dropbox.com/sh/1wy3gw7fpil73dd/AADIOGvbsYNdt45KnaHahmqqa?dl=0) (6GB; password required, which will be shared on the course Slack channel). Only if you want to work from raw `fastq` files (substantial computation needed) and have a high-bandwidth connection, download the files at the link and save to `data/fastq`.
2. Merged `Bismark` files are available at this [link](https://www.dropbox.com/sh/b3v55pdkkimo13s/AAA4gH-6uCxMqFSbFM72rwLna?dl=0) (76MB; password required). Download and copy these to `data/bismark/merged`.
3. Summarized, annotated methylation results that we will use for analysis are available in the results folder of this repository (we will generate these ourselves during the course).


## Software requirements:
* `R` >=3.3.0 with packages:
    * From CRAN: `tidyverse`, `data.table`, `docopt`
    * From Bioconductor: `scater`, `scran`, `GenomicRanges`, `SC3`, `pcaMethods`
* `RStudio`
* `Python` >=3.4 with packages: `snakemake`
* [`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), which requires [Cutadapt](https://github.com/marcelm/cutadapt/)
* [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [`Bismark`](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
* [`FastQC`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [`MultiQC`](http://multiqc.info/)
* [`MethylQA`](http://methylqa.sourceforge.net/index.php)

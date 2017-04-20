(I’m not sure which cell lines we will use at the moment but this shouldn’t affect the goals of the practical - I haven’t decided yet whether to use ES cells which we know are interesting methylation-wise or take a punt with the HEK and K562 lines they are using for RNA-seq / protein practicals. Also, I’m not sure how much time you will have? What I outline below would probably take a full day to teach, but I guess you could zip through some of it as a demonstration?)



1.       Goals etc.

a.       Experimental design: two cell types, ~16 cells each

b.       Two main goals:
   i.      Methylation profiles define cell type (i.e. cells will cluster apart by e.g. PCA)                                                             ii.      Context specificity of methylation variance. E.g. in mouse ES cells, CGIs are homogenous (and low in methylation), repeat elements are homogenously high and active enhancer elements are heterogeneous. This is interesting because the enhancer elements are cell type specific and thus some variation in the methylation levels here implies plasticity in cell identity which could be important for lineage formation.

2.       Processing:

a.       For alignments and methylation calling, see our protocol paper http://www.nature.com/nprot/journal/v12/n3/full/nprot.2016.187.html (this could be done in advance in order to concentrate on more interesting stuff)

b.       QC (also see protocol paper)
    i.      Negative controls should not align
    ii.      bisulfite conversion efficiency (assessed using CHH methylation from bismark reports) should be >95%
    iii.      mapping efficiency (from bismark reports)  >10% (30-40% is normal here but may end up lower  in these practicals)
    iv.      number of CpG sites covered (I use 1M unique positions but this will depend on seq depth so maybe just exclude outliers)

c.       Preprocessing – you could ask Ricard for some example code (I could also send you some but you’d probably prefer to read Ricard’s)
    i.      Quantify methylation over regions of interest (promoters, gene bodies, enhancers, repeats...).

        1.       mean methylation rate (each covered position counts once – i.e. do not give extra weight to positions with >1 read)

        2.       also record the coverage (number of CpG sites that were covered in the that cell at that locus) for the purpose of assigning weights to each cell in downstream analyses

d.       Analysis
    i.      Mean methylation by feature / cell type
    ii.      Variation by feature / cell type
    iii.      Clustering

Software requirements:
* R >=3.3.0 with packages:
    From CRAN: tidyverse, data.table
    From Bioconductor: scater, scran, GenomicRanges, SC3
* RStudio
* Python >=3.4 with packages: ipython, jupyter, scipy, numpy, snakemake
* FastQC
* MultiQC
* Trim Galore!
* Bowtie2
* Bismark
* MethylQA

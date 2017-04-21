"
Merge Bismark methylation output files

Usage:
merge.R --input_folder <DATA_IN> --output_folder <DATA_OUT> [-hv]

Options:
-i --input_folder <DATA_IN>    path to input file directory
-o --output_folder <DATA_OUT>  path to directory for output files
-h --help                      show this
-v --version                   print version and stop

This merges output bismark files where two files are produced for each sample.

## Input ##
(1) BISMARK
Bismark output files are expecting to contain columns for chromosome, position
methylation rate
#    chr   pos        rate
#    19    3152031    100
#    19    3152424    0

## Output ##
(1) a tmp folder with a bunch of tsv files with the preprocessed samples and
annotations.
#   For example:
#     tmp/[SAMPLE]_[FEATURE].tsv
(2) all.tsv: all annotations and samples in one dataframe
# output format=1
#   sample    id    anno    rate    weight
#   3289STDY6312493    super_enhancers_100    super_enhancers    42    31
#   3289STDY6312493    super_enhancers_1001    super_enhancers    0    1
#   3289STDY6312493    super_enhancers_1002    super_enhancers    0    2

This program requires the R packages 'data.table', 'docopt', 'purrr' and
'stringr'.

Ricard Argeluet and Davis McCarthy
April 2017

" -> doc

library(data.table)
library(purrr)

main <- function(indir, outdir) {
    filenames <- list.files(indir,pattern="(cov.gz)$")
    ## samples <- unique(substr(filenames,1,nchar(filenames)-10))
    samples <- unique(sapply(strsplit(filenames,split="_"),"[[",2))

    for (i in 1:length(samples)) {
        print(sprintf("%s (%d/%d)", samples[i], i, length(samples)))
        fname.in <- sprintf("%s/%s",indir,filenames[grep(samples[i],filenames)])

        if (length(fname.in) == 2) {
            dat1 <- fread(sprintf("zcat %s",fname.in[1]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            dat2 <- fread(sprintf("zcat %s",fname.in[2]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
            colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
            dat <- rbind(dat1,dat2) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)

                                        # In some cases, samples have four files (WHY?)
        } else if (length(fname.in) == 4) {
            dat1 <- fread(sprintf("zcat %s",fname.in[1]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            dat2 <- fread(sprintf("zcat %s",fname.in[2]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            dat3 <- fread(sprintf("zcat %s",fname.in[3]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            dat4 <- fread(sprintf("zcat %s",fname.in[4]), sep="\t", header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
            colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
            colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
            colnames(dat3) <- c("chr","pos","met_reads","nonmet_reads")
            colnames(dat4) <- c("chr","pos","met_reads","nonmet_reads")
            dat <- rbind(dat1,dat2,dat3,dat4) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)
        } else {
            stop("error")
        }
        fwrite(dat, file=sprintf("%s/%s.tsv",outdir,samples[i]), sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
    }
    system(sprintf("gzip -f %s/*.tsv",io$outdir))
}

opts <- docopt::docopt(doc, version = "version 0.0.1\n")

cat("working directory: ", getwd(), "\n")
cat("input folder: ", opts$input_folder, "\n")
cat("output folder: ", opts$output_folder, "\n")

## Run main function
main(opts$input_folder, opts$output_folder)

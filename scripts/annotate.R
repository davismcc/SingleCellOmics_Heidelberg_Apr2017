"
Overlap and annotate Bismark methylation output files with genomic features

Usage:
annotate.R --input_folder <DATA_IN> --anno_folder <ANNO_IN> --output_folder <DATA_OUT> [-hv]

Options:
-i --input_folder <DATA_IN>    path to input file directory
-a --anno_folder <ANNO_IN>     path to directory for output files
-o --output_folder <DATA_OUT>  path to directory for output files
-h --help                      show this
-v --version                   print version and stop

This script overlaps the output bismark files (individual CpG sites) with
genomic features such as promoters, gene bodies, etc. defined in BED files.

- Preprocessing of annotations: collect all CpG sites
- Preprocessing of samples: collect all CpG sites from mm10
- Annotate samples with the preprocessed annotations in BED files

## Input ##
(1) BISMARK
Bismark output files are expecting to contain columns for chromosome, position
methylation rate
#    chr   pos        rate
#    19    3152031    100
#    19    3152424    0

(1) ANNO
Genomic feature annotation files in BED6 format
#   1    3531624    3531843    *    CGI_1    CGI
#   1    3670619    3671074    *    CGI_2    CGI
#   1    3671654    3672156    *    CGI_3    CGI

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

options(warn = -1)
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(stringr))

main <- function(io, opts) {
    # Load samples to be kept
    opts$samples <- dir(io$input_folder)
    stopifnot(!anyDuplicated(opts$samples))
    ## Start processing ##
    cat("\nProcessing methylation samples with the following options:\n")
    cat(sprintf("- Input folder for annotation: %s\n", io$anno_folder))
    cat(sprintf("- Input folder for methylation files: %s\n", io$input_folder))
    cat(sprintf("- Output folder: %s\n", io$output_folder))
    cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse = " ")))
    cat(sprintf("- Valid chromosomes: %s\n",
                paste(opts$chr_list, collapse = " ")))
    cat(sprintf("- Samples: %s\n", paste(opts$chr_list, collapse = " ")))
    cat("\n")
    ## annotations
    anno_list <- list()
    for (i in opts$anno) {
        # Read annotation file
        anno.file <- sprintf("%s/%s.bed", io$anno_folder, i)
        dat_anno <- fread(anno.file, sep = "\t", header = FALSE,
                            select = c(1, 2, 3, 4, 5), verbose = FALSE)
        colnames(dat_anno) <- c("chr", "start", "end", "strand", "id")
        anno_list[[i]] <- dat_anno %>% .[chr %in% opts$chr_list, ] %>%
                            setkey(chr, start, end)
    }
    names(anno_list) <- opts$anno

    ## Preprocess and annotate the samples ##
    ## Create ouput temporary folder
    dir.create(file.path(io$output_folder, "tmp"), recursive = TRUE,
               showWarnings = FALSE)
    ## Run over samples
    pb <- txtProgressBar(min = 1, max = length(opts$samples),
                         style = 3)
    counter <- 1
    for (i in opts$samples) {
        setTxtProgressBar(pb, counter)
        counter <- counter + 1
        sample <- i
        samples_processed <- list.files(file.path(io$output_folder, "tmp"))
        ## Read and parse raw methylation data
        file_in <- file.path(io$input_folder, sample)
        dat_sample <- fread(sprintf("zcat < %s", file_in),
            sep = "\t", header = TRUE, verbose = FALSE, showProgress = FALSE)
        colnames(dat_sample) <- c("chr", "pos", "rate")
        # Remove weird chromosomes and add 'start' and 'end' columns for overlap
        dat_sample <- dat_sample[chr %in% opts$chr_list, ] %>%
                .[, c("start", "end") := list(pos, pos)] %>%
                .[, pos := NULL] %>% setkey(chr, start, end)
        # Overlap data with annotations
        for (anno in opts$anno) {
            fname.out <- sprintf("%s/tmp/%s_%s.tsv",
                                    io$output_folder, sample, anno)
            # Overlap
            ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch = 0) %>%
                .[, "i.end" := NULL]
            names(ov)[names(ov) == "i.start"] <- "pos"
            # Calculate methylation status for each region in the annotation by
            # summarising over all CG sites
            out <- ov[, c("sample", "anno") := list(sample, anno)] %>%
            .[, .(rate = round(mean(rate)), weight = .N),
                keyby = .(sample, id, anno)]
            # Store and save results
            fwrite(out, fname.out, quote = FALSE, sep = "\t",
                    col.names = TRUE, row.names = FALSE)
        }
    }
    close(pb)
    # Compress output files
    cat("Compressing...\n")
    system(sprintf("gzip -f %s/tmp/*.tsv", io$output_folder))
    # Concatenate everything and save it
    cat("Annotations finished, combining results...\n")
    files <- list.files(file.path(io$output_folder, "tmp"), pattern = ".tsv.gz")
    foo <- lapply(files,
        function(f) fread(sprintf("zcat < %s/tmp/%s", io$output_folder, f))) %>%
        rbindlist
    write.table(foo, sprintf("%s/all.tsv", io$output_folder), quote = FALSE,
                sep = "\t", col.names = TRUE, row.names = FALSE)
    system(sprintf("gzip -f %s/all.tsv", io$output_folder))
    message("Done!")
}

#####################
## Define options ##
####################
io <- docopt::docopt(doc, version = "version 0.0.1\n")

cat("working directory: ", getwd(), "\n")
cat("input folder: ", io$input_folder, "\n")
cat("anno folder: ", io$anno_folder, "\n")
cat("output folder: ", io$output_folder, "\n")

## Options ##
opts <- list()
# Define annotations to analyse
opts$annos <- c("active_enhancers", "CGI", "genebody", "IAP", "promoters")
# Define valid chromosomes
opts$chr_list <- c(1:19, "X", "Y")

## Run main function
main(io, opts)

## io <- list()
## io$anno_folder <- "annotation"
## io$input_folder <- "scbssseq_data"
## io$output_folder <- "results"

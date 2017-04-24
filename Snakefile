"""
Snakefile for single-cell methylation analysis

Davis McCarthy
EMBL-EBI
April 2017

Run on a cluster with a command like:
snakemake --jobs 1000 --latency-wait 30 --cluster 'bsub  -R "rusage[mem=32000]" -M 32000 -o ./snake_logs -e ./snake_logs'
"""

import glob
import os
#from subprocess import run
import pandas as pd
import re

TEST = True

if TEST:
    SAMPLES_LONG = glob.glob('data/fastq/test/*.fastq.gz')
    SAMPLES = [os.path.basename(w).replace('lane[123]+_', '') for w in SAMPLES_LONG]
    SAMPLES = [w.replace('.fastq.gz', '') for w in SAMPLES]
    SAMPLES_LONG = [os.path.basename(w).replace('.fastq.gz', '') for w in SAMPLES_LONG]
    SAMPLES_MERGE = [w.replace('_R1', '').replace('_R2', '') for w in SAMPLES]
else:
    SAMPLES_LONG = glob.glob('data/fastq/*.fastq.gz')
    SAMPLES = [w.replace('lane[123]+_', '') for w in SAMPLES_LONG]
    SAMPLES = [w.replace('_[ATCG]+_.*.fastq.gz', '') for w in SAMPLES]
    SAMPLES_LONG = [os.path.basename(w).replace('.fastq.gz', '') for w in SAMPLES_LONG]
    SAMPLES_MERGE = [w.replace('_R1', '').replace('_R2', '') for w in SAMPLES]

fastqc_html_reports = expand('reports/fastqc/{sample}_fastqc.html', sample = SAMPLES_LONG)

print(SAMPLES_MERGE)

rule all:
    input:
        'results/all.tsv.gz',
        'reports/multiqc/multiqc_report.html',
        #expand('data/bismark/methyl/{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz', sample = SAMPLES)


rule fastqc_reports:
    input:
        'data/fastq/test/{sample}.fastq.gz'
        ## if test = False, remove test/ from path above
    output:
        'reports/fastqc/{sample}_fastqc.html'
    params:
        output_dir="reports/fastqc/"
    shell:
        '/Users/davis/src/FastQC.app/Contents/MacOS/fastqc -o {params.output_dir} {input}'


rule trim_fastq:
    input:
        'data/fastq/test/{sample}.fastq.gz'
        ## if test = False, remove test/ from path above
    output:
        temp('{sample}_trimmed.fq.gz'),
        temp('{sample}.fastq.gz_trimming_report.txt')
    log:
        "logs/trim_fastq/{sample}.log"
    shell:
        'trim_galore --gzip --non_directional --rrbs '
        '{input} '


rule fastqc_reports_trimmed:
    input:
        '{sample}_trimmed.fq.gz'
    output:
        temp('reports/fastqc/{sample}_trimmed_fastqc.html')
    params:
        output_dir="reports/fastqc/"
    shell:
        '/Users/davis/src/FastQC.app/Contents/MacOS/fastqc -o {params.output_dir} {input}'


rule bismark_prepare_genome:
    input:
        'genome'
    output:
        'genome/Bisulfite_Genome'
    shell:
        'bismark_genome_preparation {input}'


rule bismark:
    input:
        '{sample}_trimmed.fq.gz',
        'genome/Bisulfite_Genome'
    output:
        temp('data/bismark/raw/{sample}_trimmed_bismark_bt2.bam') ## CHECK BISMARK OUTPUT
    shell:
        'bismark --non_directional --genome genome -o data/bismark/raw '
        '{input}'


rule bismark_dedup:
    input:
        'data/bismark/raw/{sample}_trimmed_bismark_bt2.bam'
    output:
        temp('data/bismark/raw/{sample}_trimmed_bismark_bt2.deduplicated.bam') ## CHECK BISMARK OUTPUT
    shell:
        'deduplicate_bismark --bam {input}'


rule bismark_methylation:
    input:
        'data/bismark/raw/{sample}_trimmed_bismark_bt2.deduplicated.bam'
    output:
        'data/bismark/methyl/{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz'
    shell:
        'bismark_methylation_extractor --gzip --bedGraph '
        '-o data/bismark/methyl {input}'


rule multiqc:
    input:
        fastqc_html_reports,
        expand('data/bismark/methyl/{sample}_trimmed_bismark_bt2.deduplicated.bismark.cov.gz', sample = SAMPLES),
        expand('{sample}_trimmed.fq.gz', sample = SAMPLES),
        expand('{sample}.fastq.gz_trimming_report.txt', sample = SAMPLES)
    output:
       'reports/multiqc/multiqc_report.html'
    shell:
        'multiqc --force --filename {output} '
        'reports/fastqc ./  '
        'data/bismark/methyl data/bismark/raw '


rule merge_methylation:
    input:
        expand('data/bismark/methyl/{sample}_R1_trimmed_bismark_bt2.deduplicated.bismark.cov.gz', sample = SAMPLES_MERGE),
        expand('data/bismark/methyl/{sample}_R2_trimmed_bismark_bt2.deduplicated.bismark.cov.gz', sample = SAMPLES_MERGE)
    output:
        expand('data/bismark/merged/{sample}.tsv.gz', sample = SAMPLES_MERGE)
    params:
        indir = 'data/bismark/methyl',
        outdir = 'data/bismark/merged'
    shell:
        'RScript scripts/merge.R -i {params.indir} -o {params.outdir}'


rule annotate_methylation:
    input:
        expand('data/bismark/merged/{sample}.tsv.gz', sample = SAMPLES_MERGE)
    output:
        'results/all.tsv.gz'
    params:
        indir = 'data/bismark/merged',
        annodir = 'annotation',
        outdir = 'results'
    shell:
        'RScript scripts/annotate.R -i {params.indir} -a {params.annodir} -o {params.outdir}'

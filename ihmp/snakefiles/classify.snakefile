'''
Author: Taylor Reiter & Phillip Brooks
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A Snakemake workflow to compute and classify sourmash signatures
Date: Sun Jan 28 2018
Run: snakemake --use-conda --use-singularity
Latest modification: Mon Feb 5
'''
##--------------------------------------------------------------------------------------##
## Variables declaration
## Declaring some variables
## (SAMPLES)
##--------------------------------------------------------------------------------------##

SAMPLES = ['HSMA33OT',
'HSM67VI9',
'MSMA26AV',]

rule all:
    input:
#        expand('output/classification/sourmash/{sample}.scaled10k.k51.gather.matches.csv',
#               sample=SAMPLES),
#        expand('outputs/classification/sourmash/{sample}.scaled10k.k51.lca.gather.matches.csv',
#               sample=SAMPLES),
#        expand('outputs/classification/sourmash/{sample}.scaled10k.k51.sig',
#               sample=SAMPLES)     
        expand('outputs/sigmaps/{sample}.scaled2k.sig'),
        expand('outputs/sigmaps/{sample}.sigmap')     

import os
                             
# Association between output files and source links
links = {
        'refseq-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/refseq-d2-k51.tar.gz',
        'genbank-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/genbank-d2-k51.tar.gz'}

# Make this association accessible via a function of wildcards
def chainfile2link(wildcards):
    return links[wildcards.chainfile]

rule download_databases:
    output:
        # We inform snakemake what this rule will generate
        os.path.join('inputs/databases/', '{chainfile}')
    message:
        '--- Downloading Data.'
    params:
        # using a function of wildcards in params
        link = chainfile2link,
        out_dir = 'inputs/database'
    shell:
        '''
        wget {params.link}
        tar xf refseq-d2-k51.tar.gz -C {params.out_dir}
        tar xf genbank-d2-k51.tar.gz -C {params.out_dir}
        '''

rule calculate_signatures:
    input:
        'inputs/data/{sample}.fastq.gz',
    output:
        sig = 'outputs/sigmaps/{sample}.scaled2k.sig',
        mapping = 'outputs/sigmaps/{sample}.mapping'
    message:
        '--- Compute sourmash signatures with quality trimmed data with sourmash'
    conda: 'env.yml'
    singularity:
        'docker://quay.io/biocontainers/sourmash:2.0.0a3--py36_0'
    log:
        'outputs/classification/sourmash/{sample}_compute.log'
    benchmark:
        'benchmarks/{sample}.compute.benchmark.txt'
    shell:
        '''
        sourmash compute {input} \
        --scaled 2000 \
        -k 21,31,51 \
        --track-abundance \
        --hash-to-reads {output.mapping}
        -o {output.sig}
        '''

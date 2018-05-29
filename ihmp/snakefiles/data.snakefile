'''
Author: Taylor Reiter & Phillip Brooks
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A simple Snakemake workflow to download the iHMP data
Date: Sun Feb 16 2018
Run: snakemake --use-conda --use-singularity
Latest modification:
'''
##--------------------------------------------------------------------------------------##
## Variables declaration
## Declaring some variables
## (links)
##--------------------------------------------------------------------------------------##

import os

rule all:
    input:
        expand(os.path.join("inputs/data/", {chainfile}), chainfile=LINKS.keys())
        
# Association between output files and source links
LINKS = { 
    'HSMA33OT.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/HSMA33OT.fastq.gz',
    'HSM67VI9.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/HSM67VI9.fastq.gz',
    'MSMA26AV.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/MSMA26AV.fastq.gz'}

SAMPLES = ['HSMA33OT',
'HSM67VI9',
'MSMA26AV']

# Make this association accessible via a function of wildcards
def chainfile2link(wildcards):
    return LINKS[wildcards.chainfile]
    

rule download_samples:
    output:
        # We inform snakemake what this rule will generate
        os.path.join("inputs/data/", "{chainfile}")
    params:
        # using a function of wildcards in params
        link = chainfile2link
    shell:
        """
        wget {params.link} -O {output}
        """
        

# Make association between db links and output file
DB_LINKS = {
        'refseq-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/refseq-d2-k51.tar.gz',
        'genbank-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/genbank-d2-k51.tar.gz'}

# Make this association accessible via a function of wildcards
def chainfile2link_db(wildcards):
    return DB_LINKS[wildcards.chainfile_db]

rule download_databases:
    output:
        # We inform snakemake what this rule will generate
        os.path.join('inputs/databases/', '{chainfile_db}')
    message:
        '--- Downloading Data.'
    params:
        # using a function of wildcards in params
        link = chainfile2link_db,
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
    #singularity:
    #    'docker://quay.io/biocontainers/sourmash:2.0.0a3--py36_0'
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

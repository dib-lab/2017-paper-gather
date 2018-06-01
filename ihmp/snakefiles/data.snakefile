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

# Association between output files and source links
links = { 
    'HSMA33OT.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/HSMA33OT.fastq.gz',
    'HSM67VI9.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/HSM67VI9.fastq.gz',
    'MSMA26AV.fastq.gz' : 'https://ibdmdb.org/tunnel/static/HMP2/WGS/1750/MSMA26AV.fastq.gz'}

SAMPLES = ['HSMA33OT',
'HSM67VI9',
'MSMA26AV']

# rule all:
#     input:
#         # expand generates the list of the final files we want
#         expand(os.path.join("inputs/data/", "{chainfile}"), chainfile=links.keys())
        
# Make association between output and source link accessible via a function of wildcards
def chainfile2link(wildcards):
    return links[wildcards.chainfile]
    

rule download_samples:
    output:
        os.path.join("inputs/data/", "{chainfile}")
    message:
        '--- Download pre-quality trimmed data from ibdmdb.org. As stated on the website, data is "filtered for quality and error checked for completeness."'
    params:
        # use a function of wildcards in params
        link = chainfile2link
    shell:'''
    wget {params.link} -O {output}
    '''

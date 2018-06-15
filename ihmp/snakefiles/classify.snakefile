'''
Author: Taylor Reiter & Phillip Brooks
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A Snakemake workflow to compute and classify sourmash signatures
Date: Sun Jan 28 2018
Run: snakemake --use-conda --use-singularity
Latest modification: Mon Feb 5
'''

SAMPLES = ['HSMA33OT',
'HSM67VI9',
'MSMA26AV']

import os
                             
# Association between output files and source links
links = {'genbank-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/genbank-d2-k51.tar.gz'}
        #'refseq-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/refseq-d2-k51.tar.gz',

# Make this association accessible via a function of wildcards
def chainfile2link(wildcards):
    return links[wildcards.chainfile]

rule download_databases:
    output:
        os.path.join('inputs/databases/', '{chainfile}')
    message:
        '--- download and unpack gather databases.'
    params:
        # using a function of wildcards in params
        link = chainfile2link,
        out_dir = 'inputs/databases'
    shell:
        '''
        wget {params.link}
        tar xf genbank-d2-k51.tar.gz -C {params.out_dir}
        #tar xf refseq-d2-k51.tar.gz -C {params.out_dir}
        '''

subworkflow data_snakefile:
    workdir: "."
    snakefile: "snakefiles/data.snakefile"

rule kmer_trim_data:
    input: data_snakefile('inputs/data/{sample}.fastq.gz')
    output: 'outputs/quality/{sample}.fq.gz'
    message: '--- kmer trim reads'
    conda: 'env.yml'
    log: 'outputs/quality/{sample}_trim.log'
    benchmark: 'benchmarks/{sample}.trim.benchmark.txt'
    shell:'''
    trim-low-abund.py -C 4 -Z 18 -V -M 4e9 --gzip -o {output} {input}
    '''

rule calculate_signatures:
    input: 'outputs/quality/{sample}.fq.gz'
    output:
        sig = 'outputs/sigs/{sample}.scaled2k.sig',
    message: '--- Compute sourmash signatures using quality trimmed data.'
    conda: 'env.yml'
    #singularity: 'docker://quay.io/biocontainers/sourmash:2.0.0a3--py36_0'
    log: 'outputs/sigs/{sample}_compute.log'
    benchmark: 'benchmarks/{sample}.compute.benchmark.txt'
    shell:'''
    sourmash compute -o {output.sig} --scaled 2000 -k 21,31,51 --track-abundance {input}
    '''

rule gather:
    input: 
        sig = 'outputs/sigs/{sample}.scaled2k.sig',
        db = 'inputs/databases/genbank-d2-k51.sbt.json'
    output:
        gather = 'outputs/gather/genbank/{sample}.gather',
        matches = 'outputs/gather/matches/{sample}.matches', 
        un = 'outputs/gather/unassigned/{sample}.un'
    message: '--- Classify signatures with gather.'
    conda: 'env.yml'
    log: 'outputs/gather/{sample}_gather.log'
    benchmark: 'benchmarks/{sample}.gather.benchmark.txt'
    shell:'''
    sourmash gather -o {output.gather} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 51 --ignore-abundance {input.sig} {input.db}
    '''
     

rule find_gather_genome_matches:
    input:
        gather = 'outputs/gather/genbank/{sample}.gather',
        db = 'inputs/databases/genbank-d2-k51.sbt.json'
    output: 'outputs/gather_genome_matches/{sample}.txt'
    message: '--- Grab names of genome fasta files that were output as matches by gather'
    log: 'outputs/gather_genome_matches/{sample}_genome_matches.log'
    benchmark: 'benchmarks/{sample}_grab_gather.benchmark.txt'
    run:
        import pandas as pd
        import os
        import glob
        from sourmash import signature 
        
        gather = pd.read_csv(str(input.gather))
        
        genome_files = list()
        for md5 in gather['md5']:
            for file in glob.glob(f'inputs/databases/.sbt.genbank-d2-k51/{md5}'):
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[0]
                genome_files.append('/data/databases/' + loaded_sig.d['filename']) # this requires the databases to be installed at /data/databases directory. We downloaded genbank with ncbi-genome-downloader. 
        
        with open(str(output), 'w') as file_handler:
            for i in genome_files:
                file_handler.write("{}\n".format(i))


rule subtract_genomes_from_reads:
    input: 
        genome_files = 'outputs/gather_genome_matches/{sample}.txt',
        reads = data_snakefile('inputs/data/{sample}.fastq.gz')
    output: 'outputs/subtracts/{sample}.fastq.gz.donut.fa'
    message: '--- Subtract genomes found by gather from reads'
    log: 'outputs/subtracts/{sample}.subtract.log'
    benchmark: 'benchmarks/{sample}_subtract_gather_genomes.benchmark.txt'
    shell:'''
    python scripts/make_donut.py -k 51 --query {input.reads} --subtract `cat {input.genome_files} | tr '\n' '\ '`
    mv {wildcards.sample}.fastq.gz.donut.fa {output}
    '''

rule assemble_subtracts:
    input: 'outputs/subtracts/{sample}.fastq.gz.donut.fa'
    output: 'outputs/megahit/{sample}.contigs.fa'
    message: '--- Assemble reads from unassigned hashes with MegaHit'
    log: 'outputs/megahit/{sample}.megahit.log'
    benchmark: 'benchmarks/{sample}_megahit.benchmark.txt'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 200 --out-dir {wildcards.sample} --out-prefix {wildcards.sample} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.sample}/{wildcards.sample}.contigs.fa {params.output_folder}/{wildcards.sample}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.sample}
    '''

rule prokka_subtract_assemblies:
    input: 'outputs/megahit/{sample}.contigs.fa'
    output: 'outputs/prokka/{sample}.faa'
    message: '--- Annotate megahit assemblies'
    log: 'outputs/prokka/{sample}.prokka.log'
    benchmark: 'benchmarks/{sample}_prokka.benchmark.txt'
    conda: 'env.yml'
    params: output_folder = 'outputs/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.sample} --metagenome --force --locustag {wildcards.sample}
    touch {output}
    '''

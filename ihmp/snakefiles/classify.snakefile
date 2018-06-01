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
links = {
        'refseq-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/refseq-d2-k51.tar.gz'} #,
        #'genbank-d2-k51.sbt.json' : 'https://s3-us-west-2.amazonaws.com/sourmash-databases/genbank-d2-k51.tar.gz'}

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
        tar xf refseq-d2-k51.tar.gz -C {params.out_dir}
        #tar xf genbank-d2-k51.tar.gz -C {params.out_dir}
        '''

subworkflow data_snakefile:
    workdir: "."
    snakefile: "snakefiles/data.snakefile"
    
rule calculate_signatures:
    input:
        files = data_snakefile('inputs/data/{sample}.fastq.gz')
    output:
        sig = 'outputs/sigmaps/{sample}.scaled2k.sig',
        sigmap = 'outputs/sigmaps/{sample}.sigmap'
    message:
        '--- Compute sourmash signatures using quality trimmed data.'
    conda: 'env.yml'
    #singularity:
    #    'docker://quay.io/biocontainers/sourmash:2.0.0a3--py36_0'
    log:
        'outputs/sigmaps/{sample}_compute.log'
    benchmark:
        'benchmarks/{sample}.compute.benchmark.txt'
    shell:
        '''
        sourmash compute -o {output.sig} --scaled 2000 -k 21,31,51 --track-abundance --hash-to-reads {output.sigmap} {input.files}
        '''

rule gather:
    input: 
        sig = 'outputs/sigmaps/{sample}.scaled2k.sig',
        db = 'inputs/databases/refseq-d2-k51.sbt.json'
    output:
        gather = 'outputs/gather/refseq/{sample}.gather',
        matches = 'outputs/gather/matches/{sample}.matches', 
        un = 'outputs/gather/unassigned/{sample}.un'
    message:
        '--- Classify signatures with gather.'
    conda: 'env.yml'
    log:
        'outputs/gather/{sample}_gather.log'
    benchmark:
        'benchmarks/{sample}.gather.benchmark.txt'
    shell:'''
    sourmash gather -o {output.gather} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 51 --ignore-abundance {input.sig} {input.db}
    '''

rule get_reads_of_un_hashes:
    input: 
        sigmap='outputs/sigmaps/{sample}.sigmap',
        un='outputs/gather/unassigned/{sample}.un'
    output: 'outputs/unassigned/{sample}.un.txt'
    message: '--- Map unassigned hashes back to reads' 
    log: 'outputs/unassigned/{sample}.log'
    benchmark: 'benchmarks/{sample}.sigmap.benchmark.txt'
    run:
        import json
        from sourmash import signature
        
        # import file with hash mappings
        with open(input.sigmap) as json_data:
            sigmap = json.load(json_data)
        
        # import the sourmash signature
        sigfp = open(input.un, 'rt')
        siglist = list(signature.load_signatures(sigfp))
        loaded_sig = siglist[0]
        
        lst = []
        for min in loaded_sig.minhash.get_mins():
            tmp = sigmap.get(str(min))
            lst.append(tmp)
        
        # filter out Nones
        lst = [x for x in lst if x is not None]
        # grab only the read name from the list 
        read_lst = [[t[0] for t in l] for l in lst]
        # flatten the list
        flat = [y for x in read_lst for y in x]
        
        # write to a file, one read per line
        with open(str(output), 'w') as file_handler:
            for i in flat:
                file_handler.write("{}\n".format(i))

rule grab_reads_of_un_hashes:
    input: 
        fq = data_snakefile('inputs/data/{sample}.fastq.gz'),
        readnames = 'outputs/unassigned/{sample}.un.txt'
    output: 'outputs/unassigned/{sample}.un.fq'
    message: '--- Grab reads from unassigned hashes'
    log: 'outputs/unassigned/{sample}_seqtk.log}'
    benchmark: 'benchmarks/{sample}.seqtk.benchmark.txt'
    conda: 'env.yml'
    shell: '''
    seqtk subseq {input.fq} {input.readnames} > {output}
    '''

rule assemble_unassigned:
    input: 'outputs/unassigned/{sample}.un.fq'
    output: 'outputs/megahit/{sample}.contigs.fa'
    message: '--- Assemble reads from unassigned hashes with MegaHit'
    log: 'outputs/megahit/{sample}.megahit.log'
    benchmark: 'benchmarks/{sample}_megahit.benchmark.txt'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.sample} --out-prefix {wildcards.sample} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.sample}/{wildcards.sample}.contigs.fa {params.output_folder}/{wildcards.sample}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.sample}
    '''

rule prokka_unassigned_assemblies:
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
     

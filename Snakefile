
'''
Download and build databases
Add prokka, mingap


# Requirements
Bash: Samtools, bbtools, bowtie2, fastqc, megahit, dimond, krakenhll, prodigal, picard, awk
Python >= 3.5: snakemake, psutil, multiqc, 

'''

import os


## CONFIG ##
configfile: "Config"

## DEFINE WORKDIR ##
WORKDIR=workflow.basedir+'/'

##### DEFS #####


##### SETUP #####
# Gather files to be processed.
file_ids = glob_wildcards(WORKDIR+"data/sequencing/{id}_R{read}.fastq{gz}")
# Grab the GLDS number.
DS_NUM=config['glds_num']
    
##### RULES #####

rule setup:
    shell:
        '''
        cd {WORKDIR}
        mkdir -p logs data/trim report/plots report/fastqc data/metadata/isa_files
        '''
        
rule help:
    run:
        print('''
        ++==============================================++
        || WGS pipeline for metagenomics GLDS datasets  ||
        ||             Paired End Edition               ||
        ++==============================================++
        
        
        Please see the Config file for set up of the pipeline.
        ''')

rule fastqc:
    input: 
        'data/{type}/{clean}/{id}_{read}.fastq.gz'
    output:
        'report/fastqc/{type}/{clean}/{id}_{read}_fastqc.html'
    shell:
        '''
        fastqc {input} -o report/fastqc/{wildcards.type}/{wildcards.clean}/
        '''
        
rule multiqc:
    input:
        decan=expand('report/fastqc/decontaminate/clean/{IDS}_R{READ}-clean_fastqc.html', IDS=file_ids.id, READ=file_ids.read),
    output:
        'report/multiqc_report.html'
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''       
               
ruleorder: bbduk_pe > bbduk_se
               
rule bbduk_pe:
    input: 
        read1='data/sequencing/{id}_R1.fastq.gz',
        read2='data/sequencing/{id}_R2.fastq.gz'
    output:
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz',
        read1_fail='data/trim/fail/{id}_R1-failed.fastq.gz',
        read2_fail='data/trim/fail/{id}_R2-failed.fastq.gz'
    log:
        stats='logs/bbduk/{id}/stats.txt',
        refstats='logs/bbduk/{id}/refstats.txt',
        stderr='logs/bbduk/{id}/stderr.txt'
    benchmark: 'logs/bbduk/{id}/benchmark.tsv'
    threads: config['duk_threads']
    shell:
        '''
        bbduk.sh in={input.read1} in2={input.read2} out={output.read1} out2={output.read2} outm={output.read1_fail} outm2={output.read2_fail}\
        stats={log.stats} refstats={log.refstats} ref={config[duk_ref]} ktrim={config[duk_ktrim]} k={config[duk_k]} mink={config[duk_mink]}\
        qtrim={config[duk_qtrim]} trimq={config[duk_trimq]} cardinalityout=t mbq={config[duk_minavgquality]} maxns={config[duk_maxns]}\
        zl={config[duk_ziplevel]} tossjunk=t overwrite=t threads={threads} 2> {log.stderr}
        '''
        
rule bbduk_se:
    input: 
        read1='data/sequncing/{id}_R1.fastq.gz'
    output:
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read1_fail='data/trim/fail/{id}_R1-failed.fastq.gz'
    log:
        stats='logs/bbduk/{id}/stats.txt',
        refstats='logs/bbduk/{id}/refstats.txt',
        stderr='logs/bbduk/{id}/stderr.txt'
    benchmark: 'logs/bbduk/{id}/benchmark.tsv'
    threads: config['duk_threads']
    shell:
        '''
        bbduk.sh in={input.read1} out={output.read1} outm={output.read1_fail}\
        stats={log.stats} refstats={log.refstats} ref={config[duk_ref]} ktrim={config[duk_ktrim]} k={config[duk_k]} mink={config[duk_mink]} \ 
        qtrim={config[duk_qtrim]} trimq={config[duk_trimq]} cardinalityout=t mbq={config[duk_minavgquality]} maxns={config[duk_maxns]} \
        zl={config[duk_ziplevel]} tossjunk=t overwrite=t 2> {log.stderr}        
        '''
        
ruleorder: contam_bowtie_pe > contam_bowtie_se
  
rule contam_bowtie_pe:
    input: 
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz'
    output:
        sam='data/decontaminate/sam/{id}.sam',
        read1_clean='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read2_clean='data/decontaminate/clean/{id}_R2-clean.fastq.gz',
        read1_contam='data/decontaminate/dirty/{id}_R1-dirty.fastq.gz',
        read2_contam='data/decontaminate/dirty/{id}_R2-dirty.fastq.gz'
    log:
        stderr='logs/contam_bowtie/{id}/stderr.txt',
        stdout='logs/contam_bowtie/{id}/stdout.txt'
    benchmark: 'logs/contam_bowtie/{id}/benchmark.tsv'
    threads: config['contam_threads']
    shell:
        '''
        bowtie2 -x {config[contam_indexes]} -1 {input.read1} -2 {input.read2} -S {output.sam} --al-conc-gz data/decontaminate/dirty/{wildcards.id}_R%-dirty.fastq.gz --un-conc-gz data/decontaminate/clean/{wildcards.id}_R%-clean.fastq.gz --{config[contam_mode]} --no-discordant --threads {threads} 2> {log.stderr} > {log.stdout} 
        '''

rule contam_bowtie_se:
    input: 
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz'
    output:
        sam='data/decontaminate/sam/{id}.sam',
        read1_clean='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read1_contam='data/decontaminate/dirty/{id}_R1-dirty.fastq.gz'
    log:
        stderr='logs/contam_bowtie/{id}/stderr.txt',
        stdout='logs/contam_bowtie/{id}/stdout.txt'
    benchmark: 'logs/contam_bowtie/{id}/benchmark.tsv'
    threads: config['contam_threads']
    shell:
        '''
        bowtie2 -x {config[contam_indexes]} -1 {input.read1} -S {output.sam} --al-conc-gz data/decontaminate/dirty/{wildcards.id}_R%-dirty.fastq.gz --un-conc-gz data/decontaminate/clean/{wildcards.id}_R%-clean.fastq.gz --{config[contam_mode]} --no-discordant --threads {threads}  2> {log.stderr} > {log.stdout}
        '''

ruleorder: megahit_pe > megahit_se

rule megahit_pe:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read2='data/decontaminate/clean/{id}_R2-clean.fastq.gz'
    output: 
        contig='data/assembly/{id}/{id}.contigs.fa'
    log:
        stderr='logs/assembly/{id}/stderr.txt',
        stdout='logs/assembly/{id}/stdout.txt'
    benchmark: 'logs/assembly/{id}/benchmark.tsv'
    threads: config['megahit_threads']
    shell:
        '''
        megahit -1 {input.read1} -2 {input.read2} -o data/assembly/{wildcards.id} --out-prefix {wildcards.id} --continue -t {threads} 2> {log.stderr} > {log.stdout}
        '''
      
rule megahit_se:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz'
    output: 
        contig='data/assembly/{id}/{id}.contigs.fa'
    log:
        stderr='logs/assembly/{id}/stderr.txt',
        stdout='logs/assembly/{id}/stdout.txt'
    benchmark: 'logs/assembly/{id}/benchmark.tsv'
    threads: config['megahit_threads']
    shell:
        '''
        megahit -1 {input.read1} -o data/assembly/{wildcards.id} --out-prefix {wildcards.id} --continue -t {threads} 2> {log.stderr} > {log.stdout}
        '''

rule prodigal_metagenome:
    input:
        contigs='data/assembly/{id}/{id}.contigs.fa'
    output:
        protiens='data/meta_gene_calls/{id}/{id}_protiens.fasta',
        genes='data/meta_gene_calls/{id}/{id}_genes.fasta',
        gff='data/meta_gene_calls/{id}/{id}_metagenome.gff',
        start='data/meta_gene_calls/{id}/{id}_prodigal.start'
    shell:
        '''
        prodigal -i {input.contigs} -a {output.protiens} -d {output.genes} -f gff -o {output.gff} -p meta                    
        '''

rule remap_bowtie_build:
    input:
        contig='data/assembly/{id}/{id}.contigs.fa'
    output: 
        outdir=directory('data/remap/{id}/index/')
        
    log:
        stderr='logs/bowtie_build/{id}/stderr.txt',
        stdout='logs/bowtie_build/{id}/stdout.txt'
    benchmark: 'logs/bowtie_build/{id}/benchmark.tsv'
    shell:
        '''
        bowtie2-build {input.contig} {output.outdir}{wildcards.id} 2> {log.stderr} > {log.stdout}
        '''

ruleorder: remap_bowtie_pe > remap_bowtie_se        
    
rule remap_bowtie_pe:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read2='data/decontaminate/clean/{id}_R2-clean.fastq.gz',
        index='data/remap/{id}/index/'
    output:
        sam=temporary('data/remap/{id}/{id}.sam')
    log:
        stderr='logs/bowtie_remap/{id}/stderr.txt',
        stdout='logs/bowtie_remap/{id}/stdout.txt'
    benchmark: 'logs/bowtie_remap/{id}/benchmark.tsv'
    threads: config['remap_threads']
    shell:
        '''
        bowtie2 -x {input.index}{wildcards.id} -U {input.read1},{input.read2} --{config[remap_mode]} --threads {threads} -S {output.sam} 2> {log.stderr} > {log.stdout}
        '''

rule remap_bowtie_se:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        index='data/remap/{id}/index/'
    output:
        sam=temporary('data/remap/{id}/{id}.sam')
    log:
        stderr='logs/bowtie_remap/{id}/stderr.txt',
        stdout='logs/bowtie_remap/{id}/stdout.txt'
    benchmark: 'logs/bowtie_remap/{id}/benchmark.tsv'
    threads: config['remap_threads']
    shell:
        '''
        bowtie2 -x {input.index}{wildcards.id} -U {input.read1} --{config[remap_mode]} --threads {threads} -S {output.sam} 2> {log.stderr} > {log.stdout}
        '''


rule remap_samtools:
    input: 
        ref='data/assembly/{id}/{id}.contigs.fa',
        sam='data/remap/{id}/{id}.sam'
    output:
        bam='data/remap/{id}/{id}.bam'
    shell:
        '''
        samtools faidx {input.ref} -o data/remap/{wildcards.id}.fai
        samtools view -bt data/remap/{wildcards.id}.fai {input.sam} > {output.bam}.1
        samtools sort {output.bam}.1 > {output.bam}
        samtools index {output.bam}
        rm {output.bam}.1
        '''

rule picard_mark_dup:
    input:
        bam='data/remap/{id}/{id}.bam'
    output:
        bam='data/remap/{id}/{id}_md.bam'
    threads: config['picard_threads']
    log:
        metrics='logs/picard/{id}/{id}.metrics',
        stderr='logs/picard/{id}/stderr.txt',
        stdout='logs/picard/{id}/stdout.txt'
    benchmark: 'logs/picard/{id}/benchmark.tsv'
    shell:
        '''
         java -Xms1g -Xmx6g -XX:ParallelGCThreads={threads} -XX:MaxPermSize=1g -XX:+CMSClassUnloadingEnabled \
             -jar {config[picard_jar]} MarkDuplicates INPUT={input.bam} OUTPUT={output.bam}.tmp1 \
             METRICS_FILE={log.metrics} \
             AS=TRUE \
             VALIDATION_STRINGENCY=LENIENT \
             MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
             REMOVE_DUPLICATES=TRUE \
             2> {log.stderr} > {log.stdout}
         samtools sort {output.bam}.tmp1 > {output.bam}
         samtools index {output.bam}
         rm {output.bam}.tmp1
         '''

rule coverage_and_abundance_calc:
    input:
        bam='data/remap/{id}/{id}_md.bam'
    output:
        cov='data/remap/{id}/{id}_coverage.txt',
        abund='data/remap/{id}/{id}_abundance.txt'
    shell:
        '''
        pileup.sh in={input.bam}  out={output.cov}
        awk '{{print $1"\t"$5}}' {output.cov} | grep -v '^#' > {output.abund}
        '''

rule maxbin:
    input:
        contig='data/assembly/{id}/{id}.contigs.fa',
        abund='data/remap/{id}/{id}_abundance.txt'
    output:
        base=dynamic('data/binning/{id}/{id}.{bin}.fasta')
#    log:
#        stderr='logs/maxbin/{id}/stderr.txt',
#        stdout='logs/maxbin/{id}/stdout.txt'
#    benchmark: 'logs/maxbin/{id}/benchmark.tsv'
    threads: config['maxbin_threads']
    shell:
        '''
        {config[maxbin_dir]}/run_MaxBin.pl -contig {input.contig} -out data/binning/{wildcards.id}/{wildcards.id} -abund {input.abund} -thread {threads}
        '''

rule prodigal_draft_genome:
    input:
        contigs='data/binning/{id}/{id}.{bin}.fasta'
    output:
        contigs='data/draft_genome/{id}/{bin}/{id}_{bin}_contigs.fasta',
        protiens='data/draft_genome/{id}/{bin}/{id}_{bin}_protiens.fasta',
        genes='data/draft_genome/{id}/{bin}/{id}_{bin}_genes.fasta',
        gff='data/draft_genome/{id}/{bin}/{id}_{bin}_genome.gff'
    shell:
        '''
        ln -s $(pwd)/{input.contigs} $(pwd)/{output.contigs}
        prodigal -i {input.contigs} -a {output.protiens} -d {output.genes} -f gff -o {output.gff} -p single                    
        '''

        
rule do_annotation:
    input: dynamic('data/draft_genome/{id}/{bin}/{id}_{bin}_genome.gff')
    output: touch('data/draft_genome/{id}/annotation.done')

rule kraken2_build:
    output:
        directory(config['kraken_ref_dir'])
    threads: config['kraken_threads']
    run:
        shell('kraken2-build --db {config[kraken_ref_dir]} --download-taxonomy --threads {threads}')
        for lib in config['kraken_bld_libraries']:
            shell('kraken2-build --db {config[kraken_ref_dir]} --download-library {lib} --threads {threads}')
        shell('kraken2-build --db {config[kraken_ref_dir]} --build --kmer-len {config[kraken_bld_kmer_len]} --minimizer-len {config[kraken_bld_minimizer_len]} --minimizer-spaces {config[kraken_bld_minimizer_spaces]} --threads {threads}')
        

ruleorder: kraken2_pe > kraken2_se     

rule kraken2_pe:
    input:
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read2='data/decontaminate/clean/{id}_R2-clean.fastq.gz',
        index=directory(config['kraken_ref_dir'])
    output:
        unclass1='data/metagenome_taxa_assignment/{id}/{id}_unclassified_1.fq',
        unclass2='data/metagenome_taxa_assignment/{id}/{id}_unclassified_2.fq',
        class1='data/metagenome_taxa_assignment/{id}/{id}_classified_1.fq',
        class2='data/metagenome_taxa_assignment/{id}/{id}_classified_2.fq',
        taxa='data/metagenome_taxa_assignment/{id}/{id}_kraken_taxa.txt',
        report='data/metagenome_taxa_assignment/{id}/{id}_kraken_report.txt'   
    log:
        stderr='logs/kraken/{id}/stderr.txt',
        stdout='logs/kraken/{id}/stdout.txt'
    benchmark: 'logs/kraken/{id}/benchmark.tsv'
    threads: config['kraken_threads']
    shell:
        '''
        kraken2 --db {input.index} --threads {threads} \
        --unclassified-out data/metagenome_taxa_assignment/{wildcards.id}/{wildcards.id}_unclassified\
        --classified-out data/metagenome_taxa_assignment/{wildcards.id}/{wildcards.id}_classified \
        --output {output.taxa} --report {output.report} --confidence {config[kraken_confidence]}\
        --paired --use-names --gzip-compressed {input.read1} {input.read2}
        '''
        
rule kraken2_se:
    input:
        read='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        index=directory(config['kraken_ref_dir'])
    output:
        unclass='data/metagenome_taxa_assignment/{id}/{id}_unclassified.fasta.gz',
        classed='data/metagenome_taxa_assignment/{id}/{id}_classified_R2.fasta.gz',
        taxa='data/metagenome_taxa_assignment/{id}/{id}_kraken_taxa.txt',
        report='data/metagenome_taxa_assignment/{id}/{id}_kraken_report.txt' 
    log:
        stderr='logs/kraken/{id}/stderr.txt',
        stdout='logs/kraken/{id}/stdout.txt'
    benchmark: 'logs/kraken/{id}/benchmark.tsv'
    threads: config['kraken_threads']
    shell:
        '''
        kraken2 --db {input.index} --threads {threads} \
        --unclassified-out {output.unclass}
        --classified-out {output.classed} \
        --output {output.taxa} --report {output.report} --confidence {config[kraken_percision]}
        --paired --use-names --gzip-compressed {input.read}
        '''

rule grab_uniprot:
    output:
        ref=ancient(config['uniprot_ref_dir']+'uniref100.fasta.gz'),
        go=ancient(config['uniprot_ref_dir']+'goa_uniprot_all.gaf.gz'),
        meta=ancient(config['uniprot_ref_dir']+'goa_uniprot_all.gpi.gz')
    shell:
        '''
        cd {config['uniprot_ref_dir']}
        wget -N ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        wget -N ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.release_note
        wget -N ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
        wget -N ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz
        '''

rule diamond_build: 
    input:
        uniref=config['uniprot_ref_dir']+'uniref100.fasta.gz'
    output:
        directory(config['diamond_ref_dir']+'uniref100'),
        config['diamond_ref_dir']+'/uniref100/uniref100.dmnd'         
    shell:
        '''
        cd config['diamond_ref_dir']
        mkdir -p uniref100
        cd uniref100
        diamond makedb --in {input.uniref} --db uniref100 -v --log 
        '''

rule diamond:
    input:
    output:
    shell:
        '''
        diamond blastx --db config['diamond_ref_dir']+'/uniref100/uniref100.dmnd' --out {output.table} --outfmt 6 --query {input.query} -max-target-seqs 1 --strand both --compress 1 --sensitive
        '''

rule taxa_assignment:
    input:
        expand('data/metagenome_taxa_assignment/{id}/{id}_classified.fasta', id=file_ids.id)
    
rule assembly:
    input: expand('data/draft_genome/{id}/annotation.done', id=file_ids.id)
    output: touch("assembly.done")

rule single:
    input:
        'data/draft_genome/{id}/annotation.done',
        expand('report/fastqc/decontaminate/clean/{{id}}_R{READ}-clean_fastqc.html', READ=file_ids.read),
        'data/metagenome_taxa_assignment/{id}/{id}_kraken_report.txt',
        
    output:
        touch('single_{id}.done')

rule all:
    input:
        qc='report/multiqc_report.html',
        assembly='assembly.done',
        
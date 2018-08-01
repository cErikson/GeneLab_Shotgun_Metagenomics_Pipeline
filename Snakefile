
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
        trim=expand('report/fastqc/trim/clean/{IDS}_R{READ}-trimmed_fastqc.html', IDS=file_ids.id, READ=file_ids.read),
    output:
        'report/multiqc_report.html'
    shell:
        '''
        multiqc -d -f -ip -o report/ report/
        '''       
               
ruleorder: bbduk_pe > bbduk_se
               
rule bbduk_pe:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz',
        read2='data/decontaminate/clean/{id}_R2-clean.fastq.gz',
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
        tbo={config[duk_tbo]} forcetrimleft={config[duk_trimleft]} forcetrimright={config[duk_trimright]}\
        qin={config[duk_qin]} zl={config[duk_ziplevel]} tossjunk=t overwrite=t threads={threads} 2> {log.stderr}
        '''
        
rule bbduk_se:
    input: 
        read1='data/decontaminate/clean/{id}_R1-clean.fastq.gz'
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
        tbo={config[duk_tbo]}\
        qin={config[duk_qin]} zl={config[duk_ziplevel]} tossjunk=t overwrite=t 2> {log.stderr}        
        '''

rule contam_bowtie_build:
    input:
        masked=config['contam_fa']
    output:
        index=directory(config['contam_index'])
    threads:config['bowtie_build_threads']
    shell:
        '''
        bowtie2-build {input.masked} {output.index}/{config[contam_name]} --threads {threads}
        '''
 
ruleorder: contam_bowtie_pe > contam_bowtie_se
  
rule contam_bowtie_pe:
    input: 
        read1='data/sequencing/{id}_R1.fastq.gz',
        read2='data/sequencing/{id}_R2.fastq.gz',
        index=directory(config['contam_index'])
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
        bowtie2 -x {input.index}/{config[contam_name]} -1 {input.read1} -2 {input.read2} -S {output.sam} --al-conc-gz data/decontaminate/dirty/{wildcards.id}_R%-dirty.fastq.gz --un-conc-gz data/decontaminate/clean/{wildcards.id}_R%-clean.fastq.gz --{config[contam_mode]} --no-discordant --threads {threads} 2> {log.stderr} > {log.stdout} 
        '''

rule contam_bowtie_se:
    input: 
        read1='data/sequencing/{id}_R1.fastq.gz',
        index=directory(config['contam_index'])
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
        bowtie2 -x {input.index}/{config[contam_name]} -1 {input.read1} -S {output.sam} --al-conc-gz data/decontaminate/dirty/{wildcards.id}_R%-dirty.fastq.gz --un-conc-gz data/decontaminate/clean/{wildcards.id}_R%-clean.fastq.gz --{config[contam_mode]} --no-discordant --threads {threads}  2> {log.stderr} > {log.stdout}
        '''

ruleorder: megahit_pe > megahit_se

rule megahit_pe:
    input: 
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz'
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
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz'
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
    benchmark: 'logs/prodigal/{id}/benchmark.tsv'
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
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz',
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
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
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
        samtools faidx {input.ref} 
        samtools view -bt data/remap/{wildcards.id}/{wildcards.id}.fai {input.sam} > {output.bam}.1
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
    threads: config['maxbin_threads']
    shell:
        '''
        mkdir -p logs/maxbin/{wildcards.id}/ logs/maxbin/{wildcards.id}/
        {config[maxbin_dir]}/run_MaxBin.pl -contig {input.contig} -out data/binning/{wildcards.id}/{wildcards.id} -abund {input.abund} -thread {threads} 2> logs/maxbin/{wildcards.id}/stderr.txt > logs/maxbin/{wildcards.id}/stdout.txt
        '''
        
rule kraken2_id:
    input:
        read='data/binning/{id}/{id}.{bin}.fasta',
        index=config['kraken_ref_dir']+'/hash.k2d'
    output:
        taxa='data/binning/{id}/taxa_id/{id}.{bin}_kraken_taxa.txt',
        report='data/binning/{id}/taxa_id/{id}.{bin}_kraken_report.txt' 
    threads: config['kraken_threads']
    shell:
        '''
        kraken2 --db $(dirname {input.index}) --threads {threads} \
        --output {output.taxa} --report {output.report} --confidence {config[kraken_id_confidence]} \
        --use-names {input.read}
        '''

rule genome_ids:
    input:
        report=dynamic('data/binning/{id}/taxa_id/{id}.{bin}_kraken_report.txt')
    output:
        calls='data/binning/{id}/taxa_id/{id}_genome-calls.tsv'
    params:
        cutoff=float(config['id_cutoff'])
    run:  
        import sys
        out_file=open(output.calls, 'w+')
        for file in input.report:
            COUNT=1; RANK=3; TAXA=4; NAME=5; 
            assign=(0, '', '', '')
            classed=0
            king=''
            with open(file, 'r') as fh:
                cell=fh.readline().strip().split('\t')
                unclass= int(cell[COUNT])
                cell=fh.readline().strip().split('\t')
                while cell != ['']:
                    if cell[RANK]=='R':
                        classed = int(cell[COUNT])
                    if cell[RANK] == 'D':
                        king=cell[NAME]
                    if int(cell[COUNT]) > classed//2:
                        assign=(cell[COUNT], cell[RANK], cell[TAXA], cell[NAME])
                    else:
                        break
                    cell=fh.readline().strip().split('\t')
                if classed==0 or classed/unclass <= params.cutoff:
                    assign=(0, 'U', 0, 'unclassified')
                    classed=1
                out_file.write('{f}\t{t}\t{r}__{n}\t{o}\t{c}\t{k}\n'.format(f=file.rsplit('/',1)[1][:-18]+'.fasta', t=assign[2], r=assign[1], n=assign[3].strip(), o=int(assign[0])/(classed+unclass), c=int(assign[0])/classed, k=king))
        out_file.close()

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
    input: 
        dynamic('data/draft_genome/{id}/{bin}/{id}_{bin}_genome.gff'),
        'data/binning/{id}/taxa_id/{id}_genome-calls.tsv'
    output: touch('data/draft_genome/{id}/annotation.done')

rule kraken2_build:
    output:
        config['kraken_ref_dir']+'/hash.k2d'
    threads: config['kraken_threads']
    run:
        shell('kraken2-build --db {config[kraken_ref_dir]} --download-taxonomy --threads {threads}')
        for lib in config['kraken_bld_libraries']:
            shell('kraken2-build --db {config[kraken_ref_dir]} --download-library {lib} --threads {threads}')
        shell('kraken2-build --db {config[kraken_ref_dir]} --build --kmer-len {config[kraken_bld_kmer_len]} --minimizer-len {config[kraken_bld_minimizer_len]} --minimizer-spaces {config[kraken_bld_minimizer_spaces]} --threads {threads}')
        

ruleorder: kraken2_pe > kraken2_se     

rule kraken2_pe:
    input:
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz',
        index=config['kraken_ref_dir']+'/hash.k2d'
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
        kraken2 --db $(dirname {input.index}) --threads {threads} \
        --unclassified-out data/metagenome_taxa_assignment/{wildcards.id}/{wildcards.id}_unclassified\
        --classified-out data/metagenome_taxa_assignment/{wildcards.id}/{wildcards.id}_classified \
        --output {output.taxa} --report {output.report} --confidence {config[kraken_confidence]}\
        --paired --use-names --gzip-compressed {input.read1} {input.read2}
        '''
        
rule kraken2_se:
    input:
        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
        index=config['kraken_ref_dir']+'/hash.k2d'
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
        kraken2 --db $(dirname {input.index}) --threads {threads} \
        --unclassified-out {output.unclass}
        --classified-out {output.classed} \
        --output {output.taxa} --report {output.report} --confidence {config[kraken_percision]}
        --use-names --gzip-compressed {input.read}
        '''

rule braken_build:
    input:
        db=config['kraken_ref_dir']+'/hash.k2d'
    output:
        dist=config['kraken_ref_dir']+'/database_kmer_distr_'+str(config['braken_readlen'])+'mers.txt'
    threads: config['kraken_threads']
    shell:
        '''
        kraken2 --db $(dirname {input.db}) --threads={threads}  --out {config[kraken_ref_dir]}/database_align.kraken $( find -L {config[kraken_ref_dir]}/library -name "*.fna" -o -name "*.fa" -o -name "*.fasta" | tr '\n' ' ' ) 
        perl {config[braken_dir]}/src/count-kmer-abundances.pl --db=$(dirname {input.db}) --threads={threads} --read-length={config[braken_readlen]} {config[kraken_ref_dir]}/database_align.kraken > {config[kraken_ref_dir]}/database{config[braken_readlen]}mers.kraken_cnts 
        python {config[braken_dir]}/src/generate_kmer_distribution.py -i {config[kraken_ref_dir]}/database{config[braken_readlen]}mers.kraken_cnts -o {config[kraken_ref_dir]}/database_kmer_distr_{config[braken_readlen]}mers.txt
        '''

rule braken:
    input:
        report='data/metagenome_taxa_assignment/{id}/{id}_kraken_report.txt',
        dist=config['kraken_ref_dir']+'/database_kmer_distr_'+str(config['braken_readlen'])+'mers.txt'
    output:
        abund='data/metagenome_taxa_assignment/{id}/{id}_braken_abund.txt'
    shell:
        '''
        python {config[braken_dir]}/src/est_abundance.py -i {input.report} -k {input.dist} -l {config[braken_level]} -t {config[braken_cutoff]} -o {output.abund}
        '''
        
rule combine_braken:
    input:
        abund=expand('data/metagenome_taxa_assignment/{id}/{id}_braken_abund.txt', id=file_ids.id)
    output:
        combine='data/{DS_NUM}_metagenomics_braken-abundances.txt'
    shell:
        '''
        python {config[braken_dir]}/analysis_scripts/combine_bracken_outputs.py --files {input.abund} --names $(basename {input.abund} -a -s _braken_abund.txt | tr '\n' ' ') -o {output.combine}
        '''

rule convert2fa:
    input:
        'data/trim/clean/{id}_{read}-trimmed.fastq.gz'
    output:
        temporary('data/metagenome_function_assignment/tmp/{id}/{id}_{read}-trimmed.fasta')
    shell:
        '''
        if [[ $(zcat {input} | head -n 5 | sed -n '5p') == @* ]]
        then zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}
        else echo 'The FASTQ does not meet the assumption that a record is 4 lines' 
        fi
        '''

rule fraggenescan:
    input:
        'data/metagenome_function_assignment/tmp/{id}/{id}_{read}-trimmed.fasta'
    output:
        'data/metagenome_function_assignment/frag/{id}_{read}.faa'
    log:
        stdout='logs/frag/{id}/{id}_{read}.stdout',
        stderr='logs/frag/{id}/{id}_{read}.stderr'
    benchmark: 'logs/frag/{id}/{id}_{read}_benchmark.tsv'
    threads: config['frag_threads']
    shell:
        '''
        {config[fraggenescan_dir]}/run_FragGeneScan.pl -genome=$(pwd)/{input} -out=$(dirname {output})/{wildcards.id}_{wildcards.read} -complete=0 -train={config[frag_train]} -thread=8 > {log.stdout} 2> {log.stderr}
        sed -i '/\*/d' {output}
        '''
        
rule chuckify:
    input:
        'data/metagenome_function_assignment/frag/{id}_{read}.faa'
    output:
        temporary(dynamic('data/metagenome_function_assignment/tmp/{id}/{id}_{read}.{chunk}.faa'))
    params:
        
    shell:
        '''
        rm -f data/metagenome_function_assignment/tmp/{wildcards.id}/{wildcards.id}_{wildcards.read}.*.faa
        awk 'BEGIN {{n_seq=0;}} /^>/ {{if(n_seq%100000==0){{file=sprintf("data/metagenome_function_assignment/tmp/{wildcards.id}/{wildcards.id}_{wildcards.read}.%d.faa",n_seq);}} print >> file; n_seq++; next;}} {{ print >> file; }}' < {input}
        '''
        
rule interprotscan:
    input:
        'data/metagenome_function_assignment/tmp/{id}/{id}_{read}.{chunk}.faa'
    output:
        temporary('data/metagenome_function_assignment/tmp/{id}/{id}_{read}.{chunk}.tsv')
    log:
        stdout='logs/inter/{id}/{id}_{read}.{chunk}.stdout',
        stderr='logs/inter/{id}/{id}_{read}.{chunk}.stderr'
    benchmark: 'logs/inter/{id}/{id}_{read}.{chunk}.benchmark.tsv'
    threads: config['inter_threads']
    shell:
        '''
        interproscan.sh -appl {config[inter_appl]} -i {input} -f tsv -o {output} --pa --goterms -dp --cpu {threads} -dp -dra
        '''

rule cat_interproscan:
    input:
        dynamic('data/metagenome_function_assignment/tmp/{id}/{id}_{read}.{chunk}.tsv')
    output:
        'data/metagenome_function_assignment/inter/{id}_{read}.tsv'
    shell:
        '''
        cat {input} > {output}
        '''

#rule grab_uniprot:
#    output:
#        ref=ancient(config['uniprot_ref_dir']+config['uniprot_ftp'].rsplit('/',1)[1]),
#        go=ancient(config['uniprot_ref_dir']+'goa_uniprot_all.gaf'),
#        meta=ancient(config['uniprot_ref_dir']+'goa_uniprot_all.gpi.gz')
#    threads:8
#    shell:
#        '''
#        cd {config[uniprot_ref_dir]}
#        wget -N {config[uniprot_ftp]}
#        wget -N ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
#        wget -N ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpi.gz
#        '''

#rule diamond_build: 
#    input:
#        uniref=config['uniprot_ref_dir']+config['uniprot_ftp'].rsplit('/',1)[1]
#    output:
#        config['diamond_ref_db']+'/'+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'/'+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'.dmnd'         
#    params:
#        name=config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]
#    shell:
#        '''
#        cd {config[diamond_ref_db]}
#        mkdir -p {params.name}
#        cd {params.name}
#       diamond makedb --in {input.uniref} --db {params.name} -v --log 
#        '''

#rule diamond_pe:
#    input:
#        read1='data/trim/clean/{id}_R1-trimmed.fastq.gz',
#        read2='data/trim/clean/{id}_R2-trimmed.fastq.gz',
#        index=config['diamond_ref_db']+'/'+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'/'+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'.dmnd'
#    output:
#        sam='data/metagenome_function_assignment/diamond/{id}_diamond.sam.gz'
#    threads: config ['dmnd_threads']
#    shell:
#        '''
#        diamond blastx --db {input.index} --outfmt 101 --query {input.read1} --max-target-seqs 1 --strand both --sensitive --threads {threads} | gzip -c > {output.sam}
#        diamond blastx --db {input.index} --outfmt 101 --query {input.read2} --max-target-seqs 1 --strand both --sensitive --threads {threads} | gzip -c >> {output.sam}
#        '''

#rule get_prot_len:
#    input:
#        prot= config['uniprot_ref_dir']+config['uniprot_ftp'].rsplit('/',1)[1]
#    output:
#        lens= config['uniprot_ref_dir']+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'.len'
#    shell:
#        '''
#        cat {input.prot} | awk 'BEGIN{SEQ=0}; /^>/{head=$0}; /^[ATCGN]/{print(head "\t" length($0))};' | sort -k 2g > {output.len}
#        '''
        
#rule get_library_len:
#    input: 
#        expand('data/decontaminate/clean/{id}_R1-clean.fastq.gz', id=file_ids.id)
#    output:
#        len=temporary('data/metagenome_function_assignment/diamond/lib-len')
#    shell:
#        '''
#        wc -l {input.read1} > {output.len}
#        '''
        
#rule count_go:
#    input:
#        dmnd='data/metagenome_function_assignment/diamond/{id}_diamond.sam.gz',
#        gaf=config['uniprot_ref_dir']+'goa_uniprot_all.gaf'
#    output:
#        gohits=temporary('data/metagenome_function_assignment/diamond/{id}_gohits'),
#        goanno='data/metagenome_function_assignment/diamond/{id}_go-annotations.tsv'        
#    shell:
#        '''
#        join -t $'\t' -e 'NA' -1 2 -2 2 <(zcat {input.gaf} | sed '/^!/ d' | awk '{{if ($1 == "UniProtKB") {{print $0}}}}') <(zcat {input.dmnd} | awk '{{if ($3 != "*") {{print $0}}}}' | cut -f 3 | cut -d \| -f 2 | sort | uniq -c | sort -k 2 | awk '{{$1=$1;printf("%s\t%s\n", $1, $2)}}') > {output.goanno}
#        awk '{{printf ("%s\t%s\n", $NF, $4)}}' | sort -k 1g | awk '{{i[$2]+=$1}} END{{for(x in i){{printf ("%s\t%s\n" , i[x], x)}}}}' | sort -k 2 < {output.goanno} > {output.gohits}
#        '''        

#rule agrregate_go:
#    input:
#        prot_len= config['uniprot_ref_dir']+config['uniprot_ftp'].rsplit('/',1)[1].split('.')[0]+'.len',
#        lib_len='data/metagenome_function_assignment/diamond/{id}_lib-len',
#        go_hits=expand('data/metagenome_function_assignment/diamond/{id}_gohits', id=file_ids.id)
#    output:
#        raw='data/metagenome_function_assignment/agrregated-go-hits_raw.tsv',
#        scale='data/metagenome_function_assignment/agrregated-go-hits_scale.tsv'
#    script:
#        'snake/scripts/agrregate_go.R'
    
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
        'data/metagenome_taxa_assignment/{id}/{id}_braken_abund.txt',
        expand('data/metagenome_function_assignment/inter/{{id}}/{{id}}_R{read}.tsv', read=file_ids.read),
    output:
        touch('single_{id}.done')

rule all:
    input:
        expand('data/metagenome_function_assignment/inter/{id}/{id}_R{read}.tsv', id=file_ids.id, read=file_ids.read),
        qc='report/multiqc_report.html',
        assembly='assembly.done',
        abundance='data/{pre}_metagenomics_braken-abundances.txt'.format(pre=DS_NUM)

    

# WGS metagenomics pipeline 
```
                                                          _____
   ______________________________________________________/  O  \___/
 _/ WGS Metagenomics Snakemake Pipeline for GeneLab        ____/   \
< _/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/_\_/__/ 
```
# Requirements
Bash: Samtools, bbtools, bowtie2, fastqc, megahit, diamond, kraken2, braken, prodigal, picard, awk, maxbin
Python >= 3.5: snakemake, psutil, multiqc

## Assumptions 
- Reads have been demuliplexed. Any barcode used for pcr/optical duplicates has been used.



## Setup
- Run `snakemake setup'
- Place Reads in `data/sequencing/`
- Edit the `Config` File.
- run snakemake with one of the targets. `snakemake {report/multiqc_report.html, abundence.done, assembly.done, function.done, all, 
single_<samp_id>.done} -j <cores> -kpr`

## Output
```

+-- assembly
|   +-- 5632_S4_L005
|       +--{sample}.contigs.fa		# assembled metagenome
|       +-- ... 					# files associated with megahit assembly
+-- binning			# binning of contigs in to genomes 
|   +-- {sample}
|       +-- {sample}.001.fasta		# Genome bin
|       +-- {sample}.002.fasta
|       +-- {sample}.{bin}.fasta
|       +-- {sample}.marker			# Marker genestats
|       +-- {sample}.summary		# summery statistics
|       +-- taxa_id
|           +-- {sample}_kraken_report.txt	# genome kraken report
|           +-- {sample}_kraken_taxa.txt	# per contig taxa assignment
|           +-- {sample}_genome-calls.tsv	# summery of genome taxa assignment based on least common ansestor who wins a majority vote.
+-- decontaminate		# Contaminate removal
|   +-- clean
|   |   +-- {sample}_{Read}-clean.fastq.gz		# reads that did not align to referance.
|   +-- dirty
|   |   +-- {sample}_{Read}.fastq.gz			# reads that aligned to the referance.
|   +-- sam
|       +-- {sample}.sam
+-- draft_genome			# Annoation of binned contigs
|   +-- {sample}
|       +-- {bin}
|       |   +-- {sample}_{bin}_contigs.fasta 	# symlink to bin 
|       |   +-- {sample}_{bin}_genes.fasta		# genes fron the annotation
|       |   +-- {sample}_{bin}_genome.gff		# GFF of the annoation
|       |   +-- {sample}_{bin}_protiens.fasta	# protiens from the annotation
|       +-- ...
|           +-- ...
+-- metadata
|   +-- ...
+-- metagenome_function_assignment
|   +-- frag			# Output of fraggenescan
|   |   +-- {sample}
|   |       +-- {sample}_{Read}.faa		# amino acid sequence 
|   |       +-- {sample}_{read}.ffn		# nuclotide sequence
|   |       +-- {sample}_{read}.gff		# gff of ORF
|   |       +-- {sample}_{read].out		# report
|   +-- inter
|   |   +-- {id}
|   |       +-- {id}_{read}....
|   +-- diamond
|       +-- {id}_{read}_diamond.daa		# Diamond alignment archives
+-- metagenome_taxa_assignment
|   +-- {study id}_metagenomics_braken-abundances.txt
|   +-- {sample}
|       +-- {sample}_braken_abund.txt
|       +-- {sample}_kraken_report.txt
+-- remap
|   +-- {sample}
|       +-- {sample}_abundance.txt
|       +-- {sample}.bam
|       +-- {sample}_coverage.txt
|       +-- {sample}_md.bam
|       +-- index
|           +-- {Bowtie2 index}
+-- sequencing
|   +-- {id}_{read}.fastq.gz
+-- tmp
|   +-- ....
+-- trim
    +-- clean
    |   +-- {sample}_{read}-trimmed.fastq.gz
    +-- fail
        +-- {sample}_{read}-failed.fastq.gz
```

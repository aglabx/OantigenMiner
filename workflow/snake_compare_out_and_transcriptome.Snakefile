import os
from pathlib import Path

data = Path('data_2')
results = data / Path('results')
the_results = results / 'final'

GENOME = "E.coli_genome"

assembly = data / f"{GENOME}.fna"
transcriptome = data / f'{GENOME}.cdna.all.fa'

transcriptome_comparison = data / 'transcriptome_comparison'
os.makedirs(transcriptome_comparison, exist_ok=True)

maxthreads = 10

rule all:
    input:
        transcriptome_comparison / 'intersection.bed'

rule align:
    input:
        genome = assembly,
        transcript = transcriptome
    output:
        sam = transcriptome_comparison / 'alignment.sam'
    conda:
        'envs/minimap2.yml'
    threads:
        maxthreads
    shell:
        """
        minimap2 -a {input.genome} {input.transcript} > {output.sam} 
        """

rule sam_to_bam:
    input:
        sam = transcriptome_comparison / 'alignment.sam'
    output:
        bam = transcriptome_comparison / 'alignment.bam'
    conda:
        'envs/samtools.yml'
    threads:
        maxthreads
    shell:
        """
        samtools view -@ {threads} -S -b {input.sam} > {output.bam}
        """

rule sort_bam:
    input:
        sam = transcriptome_comparison / 'alignment.bam'
    output:
        bam = transcriptome_comparison / 'sorted_alignment.bam'
    conda:
        'envs/samtools.yml'
    threads:
        maxthreads
    shell:
        """
        samtools sort -@ {threads} {input.sam} -o {output.bam}
        """

rule index_bam:
    input:
        bam = transcriptome_comparison / 'sorted_alignment.bam'
    output:
        transcriptome_comparison / 'sorted_alignment.bam.bai'
    conda:
        'envs/samtools.yml'
    threads:
        maxthreads
    shell:
        """
        samtools index {input.bam} -@ {threads}
        """

rule bam_to_bedgraph:
    input:
        bam = transcriptome_comparison / 'sorted_alignment.bam',
        bai = transcriptome_comparison / 'sorted_alignment.bam.bai'
    output:
        bedgraph = transcriptome_comparison / 'sorted_alignment.bedgraph'
    conda:
        'envs/deeptools.yml'
    threads:
        maxthreads
    shell:
        """
        bamCoverage --bam {input.bam} -o {output.bedgraph} -of bedgraph
        """

rule bedgraph_to_bed:
    input:
        bedgraph = transcriptome_comparison / 'sorted_alignment.bedgraph'
    output:
        bed = transcriptome_comparison / 'sorted_alignment.bed'
    shell:
        """
        cat {input.bedgraph} | awk '{{print $1 "\t" $2 "\t" $3}}' > {output.bed}
        """


rule get_target_tsv:
    input:
        gff = the_results / 'operons_annotation.gff3'
    output:
        transcriptome_comparison / 'o_ant_operons.tsv'
    params:
        script_path = 'scripts/extract_operons_from_operonmapper.py',
        out_dir = transcriptome_comparison / 'o_ant_operons'
    shell:
        """
        pip install genomenotebook
        python3 {params.script_path} -i {input.gff} -o {params.out_dir} 
        """

rule tsv_to_bed:
    input:
        tsv = transcriptome_comparison / 'o_ant_operons.tsv'
    output:
        bed = transcriptome_comparison / 'o_ant_operons.bed'
    shell:
        """
        cat {input.tsv} | cut -f 5,7,8 | grep -v chr > {output.bed}
        """

rule intersect:
    input:
        ref = transcriptome_comparison / 'sorted_alignment.bed',
        ops = transcriptome_comparison / 'o_ant_operons.bed'
    output:
        bed = transcriptome_comparison / 'intersection.bed'
    conda:
        'envs/bedtools.yml'
    shell:
        """
        bedtools intersect -a {input.ref} -b {input.ops} > {output.bed}
        """

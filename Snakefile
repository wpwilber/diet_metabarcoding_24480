PRIMER1_F = "GATATCCGTTGCCGAGAGTC"
PRIMER1_R = "CCGAAGGCGTCAAGGAACAC"

PRIMER2_F = "GGGCAATCCTGAGCCAA"
PRIMER2_R = "CCATTGAGTCTCTGCACCTATC"

# Auto-discover sample names from fastq directory
SAMPLES = glob_wildcards("fastq/{sample}_R1_001.fastq.gz").sample

rule all:
    input:
        # QC for raw samples
        expand("qc_raw/{sample}_R1_001_fastqc.html", sample=SAMPLES),
        expand("qc_raw/{sample}_R2_001_fastqc.html", sample=SAMPLES),
        # demux QC for amp1, amp2, unnasigned
        expand("qc_demux/{sample}_amp1_R1_fastqc.html", sample=SAMPLES),
        expand("qc_demux/{sample}_amp1_R2_fastqc.html", sample=SAMPLES),
        expand("qc_demux/{sample}_amp2_R1_fastqc.html", sample=SAMPLES),
        expand("qc_demux/{sample}_amp2_R2_fastqc.html", sample=SAMPLES),
        expand("qc_demux/{sample}_unassigned_R1_fastqc.html", sample=SAMPLES),
        expand("qc_demux/{sample}_unassigned_R2_fastqc.html", sample=SAMPLES),
        # demuxed samples
        expand("demux/{sample}_amp1_R1.fastq.gz", sample=SAMPLES),
        expand("demux/{sample}_amp1_R2.fastq.gz", sample=SAMPLES),
        expand("demux/{sample}_amp2_R1.fastq.gz", sample=SAMPLES),
        expand("demux/{sample}_amp2_R2.fastq.gz", sample=SAMPLES)

rule fastqc_raw:
    conda: "envs/environment.yaml"
    input:
        r1 = "fastq/{sample}_R1_001.fastq.gz",
        r2 = "fastq/{sample}_R2_001.fastq.gz"
    output:
        "qc_raw/{sample}_R1_001_fastqc.html",
        "qc_raw/{sample}_R1_001_fastqc.zip",
        "qc_raw/{sample}_R2_001_fastqc.html",
        "qc_raw/{sample}_R2_001_fastqc.zip"

    threads: 2
    shell:
        r"""
        mkdir -p qc_raw
        fastqc --threads {threads} --outdir qc_raw {input.r1} {input.r2}
        """

rule demux_by_primer:
    conda: "envs/environment.yaml"
    input:
        r1 = "fastq/{sample}_R1_001.fastq.gz",
        r2 = "fastq/{sample}_R2_001.fastq.gz"
    output:
        amp1_r1 = "demux/{sample}_amp1_R1.fastq.gz",
        amp1_r2 = "demux/{sample}_amp1_R2.fastq.gz",
        amp2_r1 = "demux/{sample}_amp2_R1.fastq.gz",
        amp2_r2 = "demux/{sample}_amp2_R2.fastq.gz",
        un_r1   = "demux/{sample}_unassigned_R1.fastq.gz",
        un_r2   = "demux/{sample}_unassigned_R2.fastq.gz"
    threads: 8
    shell:
        r"""
        cutadapt \
            -j {threads} \
            -e 0.10 \
            --action=none \
            \
            -g amp1=^{PRIMER1_F} \
            -g amp2=^{PRIMER2_F} \
            \
            --pair-filter=any \
            \
            --untrimmed-output        {output.un_r1} \
            --untrimmed-paired-output {output.un_r2} \
            \
            -o demux/{wildcards.sample}_{{name}}_R1.fastq.gz \
            -p demux/{wildcards.sample}_{{name}}_R2.fastq.gz \
            {input.r1} {input.r2}
        """

rule fastqc_demux:
    conda: "envs/environment.yaml"
    input:
        "demux/{sample}_{amp}_R1.fastq.gz",
        "demux/{sample}_{amp}_R2.fastq.gz"
    output:
        "qc_demux/{sample}_{amp}_R1_fastqc.html",
        "qc_demux/{sample}_{amp}_R1_fastqc.zip",
        "qc_demux/{sample}_{amp}_R2_fastqc.html",
        "qc_demux/{sample}_{amp}_R2_fastqc.zip"
    threads: 2
    shell:
        r"""
        mkdir -p qc_demux
        fastqc --threads {threads} --outdir qc_demux {input}
        """

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
        expand("demux/{sample}_amp2_R2.fastq.gz", sample=SAMPLES),
        # read counts
        "read_counts/read_counts.tsv",

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

rule count_reads:
    input:
        raw = expand("fastq/{sample}_R{read}_001.fastq.gz",
                     sample=SAMPLES, read=[1,2]),
        demux = expand("demux/{sample}_{amp}_R{read}.fastq.gz",
                       sample=SAMPLES,
                       amp=["amp1", "amp2", "unassigned"],
                       read=[1,2])
    output:
        "read_counts/read_counts.tsv"
    threads: 2
    shell:
        r"""
        mkdir -p read_counts

        echo -e "file\treads" > {output}

        # Count reads in raw FASTQs
        for f in fastq/*.fastq.gz; do
            n=$(gzip -dc "$f" | wc -l)
            echo -e "$(basename "$f")\t$((n/4))" >> {output}
        done

        # Count reads in demultiplexed FASTQs
        for f in demux/*.fastq.gz; do
            n=$(gzip -dc "$f" | wc -l)
            echo -e "$(basename "$f")\t$((n/4))" >> {output}
        done
        """

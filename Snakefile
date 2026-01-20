# Define primers

ITS1_F = "GATATCCGTTGCCGAGAGTC"
ITS1_R = "CCGAAGGCGTCAAGGAACAC"

trnL_F = "GGGCAATCCTGAGCCAA"
trnL_R = "CCATTGAGTCTCTGCACCTATC"

# Auto-discover sample names from fastq directory
SAMPLES = glob_wildcards("fastq/{sample}_R1_001.fastq.gz").sample

# Define snakemake inputs.
rule all:
    input:
        # demuxed samples
        expand("trim_clean_qc/demux/{sample}_{amp}_R1.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/demux/{sample}_{amp}_R2.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # trimmed reads from cutadapt
        expand("trim_clean_qc/trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/trimmed/{sample}_{amp}_R2.primertrim.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/trimmed_reports/{sample}_{amp}_cutadapt.txt",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # cleaned reads from fastp
        expand("trim_clean_qc/cleaned/{sample}_{amp}_R1.cleaned.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/cleaned/{sample}_{amp}_R2.cleaned.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/cleaned_reports/{sample}_{amp}_fastp.html",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/cleaned_reports/{sample}_{amp}_fastp.json",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # length filtered files
        expand("trim_clean_qc/length_filtered/{sample}_{amp}_R1.lenfilt.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/length_filtered/{sample}_{amp}_R2.lenfilt.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

       # final QC
        expand("trim_clean_qc/qc_final/{sample}_{amp}_R1.lenfilt_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/qc_final/{sample}_{amp}_R2.lenfilt_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # read retention table
        "trim_clean_qc/read_summary/read_retention.tsv",

        # DADA2 output
        "dada2/trnL_done", 
        "dada2/ITS1_done"

# This rule separates the fastq files by target amplicon (ITS1 or trnL). In this pipeline I refer to this step as "demuxing" even though it does not fit a precise definition of demultiplexing. This rule leaves a small percentage of unassigned reads that often very nearly match the expected primer sequence. It might be worth relaxing e or dropping ^ to see if we can include a few more reads, but it's pretty marginal.

rule demux_by_primer:
    conda: "envs/environment.yaml"
    input:
        r1 = "fastq/{sample}_R1_001.fastq.gz",
        r2 = "fastq/{sample}_R2_001.fastq.gz"
    output:
        ITS1_r1 = "trim_clean_qc/demux/{sample}_ITS1_R1.fastq.gz",
        ITS1_r2 = "trim_clean_qc/demux/{sample}_ITS1_R2.fastq.gz",
        trnL_r1 = "trim_clean_qc/demux/{sample}_trnL_R1.fastq.gz",
        trnL_r2 = "trim_clean_qc/demux/{sample}_trnL_R2.fastq.gz",
        un_r1   = "trim_clean_qc/demux/{sample}_unassigned_R1.fastq.gz",
        un_r2   = "trim_clean_qc/demux/{sample}_unassigned_R2.fastq.gz"
    threads: 8
    shell:
        r"""
        cutadapt \
            -j {threads} \
            -e 0.10 \
            --action=none \
            \
            -g ITS1=^{ITS1_F} \
            -g trnL=^{trnL_F} \
            \
            --pair-filter=any \
            \
            --untrimmed-output        {output.un_r1} \
            --untrimmed-paired-output {output.un_r2} \
            \
            -o trim_clean_qc/demux/{wildcards.sample}_{{name}}_R1.fastq.gz \
            -p trim_clean_qc/demux/{wildcards.sample}_{{name}}_R2.fastq.gz \
            {input.r1} {input.r2}
        """

# This rule trims primers from demuxed fastq files and applies an N base filter. Zero tolerance N filtering can remove a lot of pairs from the small number of samples I have run so far. One possible compromise here could be to truncate low quality tails and then run N filtering. I need to look at the distribution of Ns in demuxed samples to determine if this could be useful. I'll look into this later.

rule trim_primers_cutadapt:
    conda: "envs/environment.yaml"
    input:
        r1 = "trim_clean_qc/demux/{sample}_{amp}_R1.fastq.gz",
        r2 = "trim_clean_qc/demux/{sample}_{amp}_R2.fastq.gz"
    output:
        r1_trim = "trim_clean_qc/trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2_trim = "trim_clean_qc/trimmed/{sample}_{amp}_R2.primertrim.fastq.gz",
        report = "trim_clean_qc/trimmed_reports/{sample}_{amp}_cutadapt.txt"
    threads: 4
    params:
        F = lambda wc: {"ITS1": ITS1_F, "trnL": trnL_F}[wc.amp],
        R = lambda wc: {"ITS1": ITS1_R, "trnL": trnL_R}[wc.amp]
    shell:
        r"""
        mkdir -p trim_clean_qc/trimmed trim_clean_qc/trimmed_reports

        cutadapt \
            -j {threads} \
            -g ^{params.F} \
            -G ^{params.R} \
            --max-n 0 \
            -o {output.r1_trim} \
            -p {output.r2_trim} \
            {input.r1} {input.r2} \
            > {output.report} 2>&1
        """

# This ruke trims nextera adapter tails with fastp. I'm trying this because I'm having a hard time with explicitly defining adapter tails in fastp, but I prefer its explicit primer trimming. This is just an experiment for now to see if I can get better tail trimming by combining trimming tools. 

# I'm getting good results trimming tails with fastp this way so I am going to leave the pipeline in tact this way for now. One more step I'm interested in implementing is a maxium length filter, which fastp doesn't explicity support. Looking at the fastp report, my insert size distribution looks really good but I'm getting some leftover debris that is long and low quality towards the tail. I believe this is just junk and a max length filter makes sense here, maybe 80 for ITS1 and 60 for trnL.

rule trim_adapters_fastp:
    conda: "envs/environment.yaml"
    input:
        r1 = "trim_clean_qc/trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2 = "trim_clean_qc/trimmed/{sample}_{amp}_R2.primertrim.fastq.gz"
    output:
        r1_clean = "trim_clean_qc/cleaned/{sample}_{amp}_R1.cleaned.fastq.gz",
        r2_clean = "trim_clean_qc/cleaned/{sample}_{amp}_R2.cleaned.fastq.gz",
        html     = "trim_clean_qc/cleaned_reports/{sample}_{amp}_fastp.html",
        json     = "trim_clean_qc/cleaned_reports/{sample}_{amp}_fastp.json"
    params:
        min_len = lambda wc: {"ITS1": 50, "trnL": 35}[wc.amp]
    threads: 4
    shell:
        r"""
        mkdir -p trim_clean_qc/cleaned trim_clean_qc/cleaned_reports

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1_clean} \
            -O {output.r2_clean} \
            --detect_adapter_for_pe \
            -q 30 \
            -l {params.min_len} \
            -w {threads} \
            --html {output.html} \
            --json {output.json}
        """

# The following rule applies a max length filter to the fastp output. This discards junk sequences that were not trimmed to the expect size range.

rule filter_length_python:
    input:
        r1 = "trim_clean_qc/cleaned/{sample}_{amp}_R1.cleaned.fastq.gz",
        r2 = "trim_clean_qc/cleaned/{sample}_{amp}_R2.cleaned.fastq.gz"
    output:
        r1_filt = "trim_clean_qc/length_filtered/{sample}_{amp}_R1.lenfilt.fastq.gz",
        r2_filt = "trim_clean_qc/length_filtered/{sample}_{amp}_R2.lenfilt.fastq.gz"
    params:
        max_len = lambda wc: {"ITS1": 90, "trnL": 70}[wc.amp]
    threads: 2
    shell:
        r"""
        mkdir -p trim_clean_qc/length_filtered

        python3 - << 'EOF'
import gzip

max_len = {params.max_len}

def filter_fastq_pair(r1_in, r2_in, r1_out, r2_out):
    with gzip.open(r1_in, "rt") as f1, gzip.open(r2_in, "rt") as f2, \
         gzip.open(r1_out, "wt") as o1, gzip.open(r2_out, "wt") as o2:

        while True:
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break

            s1 = f1.readline().rstrip()
            s2 = f2.readline().rstrip()

            p1 = f1.readline()
            p2 = f2.readline()

            q1 = f1.readline().rstrip()
            q2 = f2.readline().rstrip()

            # Keep the pair only if BOTH reads pass
            if len(s1) <= max_len and len(s2) <= max_len:
                o1.write(h1); o1.write(s1 + "\n"); o1.write(p1); o1.write(q1 + "\n")
                o2.write(h2); o2.write(s2 + "\n"); o2.write(p2); o2.write(q2 + "\n")

filter_fastq_pair("{input.r1}", "{input.r2}", "{output.r1_filt}", "{output.r2_filt}")
EOF
        """

# Run fastqc on the final length filtered files.

rule fastqc_final:
    conda: "envs/environment.yaml"
    input:
        r1 = "trim_clean_qc/length_filtered/{sample}_{amp}_R1.lenfilt.fastq.gz",
        r2 = "trim_clean_qc/length_filtered/{sample}_{amp}_R2.lenfilt.fastq.gz"
    output:
        "trim_clean_qc/qc_final/{sample}_{amp}_R1.lenfilt_fastqc.html",
        "trim_clean_qc/qc_final/{sample}_{amp}_R1.lenfilt_fastqc.zip",
        "trim_clean_qc/qc_final/{sample}_{amp}_R2.lenfilt_fastqc.html",
        "trim_clean_qc/qc_final/{sample}_{amp}_R2.lenfilt_fastqc.zip"
    threads: 2
    shell:
        r"""
        mkdir -p trim_clean_qc/qc_final
        fastqc --threads {threads} --outdir trim_clean_qc/qc_final {input.r1} {input.r2}
        """

# Rule for generating a read retention table across pipeline steps.

rule summarize_read_retention:
    input:
        raw = expand("fastq/{sample}_R1_001.fastq.gz", sample=SAMPLES),
        unassigned = expand("trim_clean_qc/demux/{sample}_unassigned_R1.fastq.gz", sample=SAMPLES),
        its_demux = expand("trim_clean_qc/demux/{sample}_ITS1_R1.fastq.gz", sample=SAMPLES),
        its_trim = expand("trim_clean_qc/trimmed/{sample}_ITS1_R1.primertrim.fastq.gz", sample=SAMPLES),
        its_clean = expand("trim_clean_qc/cleaned/{sample}_ITS1_R1.cleaned.fastq.gz", sample=SAMPLES),
        its_filt = expand("trim_clean_qc/length_filtered/{sample}_ITS1_R1.lenfilt.fastq.gz", sample=SAMPLES),
        trnl_demux = expand("trim_clean_qc/demux/{sample}_trnL_R1.fastq.gz", sample=SAMPLES),
        trnl_trim = expand("trim_clean_qc/trimmed/{sample}_trnL_R1.primertrim.fastq.gz", sample=SAMPLES),
        trnl_clean = expand("trim_clean_qc/cleaned/{sample}_trnL_R1.cleaned.fastq.gz", sample=SAMPLES),
        trnl_filt = expand("trim_clean_qc/length_filtered/{sample}_trnL_R1.lenfilt.fastq.gz", sample=SAMPLES)
    output:
        "trim_clean_qc/read_summary/read_retention.tsv"
    run:
        import gzip, os
        os.makedirs("trim_clean_qc/read_summary", exist_ok=True)

        def count_reads(path):
            with gzip.open(path, "rt") as f:
                return sum(1 for _ in f) // 4

        with open(output[0], "w") as out:
            out.write(
                "sample\traw\tunassigned\tITS_demuxed\tITS_trimmed\tITS_cleaned\tITS_filtered\t"
                "trnL_demuxed\ttrnL_trimmed\ttrnL_cleaned\ttrnL_filtered\n"
            )

            for sample in SAMPLES:
                raw = count_reads(f"fastq/{sample}_R1_001.fastq.gz")
                unassigned = count_reads(f"trim_clean_qc/demux/{sample}_unassigned_R1.fastq.gz")
                its_demux = count_reads(f"trim_clean_qc/demux/{sample}_ITS1_R1.fastq.gz")
                its_trim = count_reads(f"trim_clean_qc/trimmed/{sample}_ITS1_R1.primertrim.fastq.gz")
                its_clean = count_reads(f"trim_clean_qc/cleaned/{sample}_ITS1_R1.cleaned.fastq.gz")
                its_filt = count_reads(f"trim_clean_qc/length_filtered/{sample}_ITS1_R1.lenfilt.fastq.gz")
                trnl_demux = count_reads(f"trim_clean_qc/demux/{sample}_trnL_R1.fastq.gz")
                trnl_trim = count_reads(f"trim_clean_qc/trimmed/{sample}_trnL_R1.primertrim.fastq.gz")
                trnl_clean = count_reads(f"trim_clean_qc/cleaned/{sample}_trnL_R1.cleaned.fastq.gz")
                trnl_filt = count_reads(f"trim_clean_qc/length_filtered/{sample}_trnL_R1.lenfilt.fastq.gz")

                out.write(
                    f"{sample}\t{raw}\t{unassigned}\t{its_demux}\t{its_trim}\t{its_clean}\t{its_filt}\t"
                    f"{trnl_demux}\t{trnl_trim}\t{trnl_clean}\t{trnl_filt}\n"
                )

# The following rule runs DADA2 sample inference on trnL samples.

rule dada2_trnl:
    input:
        expand("trim_clean_qc/length_filtered/{sample}_trnL_R1.lenfilt.fastq.gz", sample=SAMPLES),
        expand("trim_clean_qc/length_filtered/{sample}_trnL_R2.lenfilt.fastq.gz", sample=SAMPLES)
    output:
        "dada2/trnL_done"
    threads: 8
    shell:
        r"""
        mkdir -p dada2
        Rscript scripts/run_dada2_trnl.R
        touch dada2/trnL_done
        """

# The following rule runs DADA2 sample inference on ITS1 samples.

rule dada2_its1:
    input:
        expand("trim_clean_qc/length_filtered/{sample}_ITS1_R1.lenfilt.fastq.gz", sample=SAMPLES),
        expand("trim_clean_qc/length_filtered/{sample}_ITS1_R2.lenfilt.fastq.gz", sample=SAMPLES)
    output:
        "dada2/ITS1_done"
    threads: 8
    shell:
        r"""
        mkdir -p dada2
        Rscript scripts/run_dada2_ITS1.R
        touch dada2/ITS1_done
        """


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
        # demux QC
        expand("trim_clean_qc/qc_demux/{sample}_{amp}_R1_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL", "unassigned"]),
        expand("trim_clean_qc/qc_demux/{sample}_{amp}_R2_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL", "unassigned"]),

        # demuxed samples
        expand("trim_clean_qc/demux/{sample}_{amp}_R1.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trim_clean_qc/demux/{sample}_{amp}_R2.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # read counts
        "trim_clean_qc/read_counts/read_counts.tsv",

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

        # primer QC summary
        "trim_clean_qc/primer_qc/primer_qc_summary.tsv",

        # read retention table
        "trim_clean_qc/read_summary/read_retention.tsv"


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

# This rule runs fastqc on the demultiplexed samples.

rule fastqc_demux:
    conda: "envs/environment.yaml"
    input:
        "trim_clean_qc/demux/{sample}_{amp}_R1.fastq.gz",
        "trim_clean_qc/demux/{sample}_{amp}_R2.fastq.gz"
    output:
        "trim_clean_qc/qc_demux/{sample}_{amp}_R1_fastqc.html",
        "trim_clean_qc/qc_demux/{sample}_{amp}_R1_fastqc.zip",
        "trim_clean_qc/qc_demux/{sample}_{amp}_R2_fastqc.html",
        "trim_clean_qc/qc_demux/{sample}_{amp}_R2_fastqc.zip"
    threads: 2
    shell:
        r"""
        mkdir -p trim_clean_qc/qc_demux
        fastqc --threads {threads} --outdir trim_clean_qc/qc_demux {input}
        """

# This rule counts the number of reads in each fastq after demuxing. See the output file read_counts.tsv.

rule count_reads:
    input:
        raw = expand("fastq/{sample}_R{read}_001.fastq.gz",
                     sample=SAMPLES, read=[1,2]),
        demux = expand("trim_clean_qc/demux/{sample}_{amp}_R{read}.fastq.gz",
                       sample=SAMPLES,
                       amp=["ITS1", "trnL", "unassigned"],
                       read=[1,2])
    output:
        "trim_clean_qc/read_counts/read_counts.tsv"
    threads: 2
    shell:
        r"""
        mkdir -p trim_clean_qc/read_counts

        echo -e "file\treads" > {output}

        for f in fastq/*.fastq.gz; do
            n=$(gzip -dc "$f" | wc -l)
            echo -e "$(basename "$f")\t$((n/4))" >> {output}
        done

        for f in trim_clean_qc/demux/*.fastq.gz; do
            n=$(gzip -dc "$f" | wc -l)
            echo -e "$(basename "$f")\t$((n/4))" >> {output}
        done
        """

# This rule trims primers from demuxed fastq files and applies an N base filter. Zero tolerance N filtering can remove a lot of pairs from the small number of samples I have run so far. One possible compromise here could be to truncate low quality tails and then run N filtering. I need to look at the distribution of Ns in demuxed samples to determine if this could be useful. I'll look into this later.

rule trim_primers:
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
            -l 50 \
            -w {threads} \
            --html {output.html} \
            --json {output.json}
        """

# The following rule applies a max length filter to the fastp output. This discards junk sequences that were not trimmed to the expect size range.

rule filter_length:
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

def filter_fastq(infile, outfile):
    with gzip.open(infile, "rt") as fin, gzip.open(outfile, "wt") as fout:
        while True:
            h = fin.readline()
            if not h:
                break
            seq = fin.readline().rstrip()
            plus = fin.readline()
            qual = fin.readline().rstrip()

            if len(seq) <= max_len:
                fout.write(h)
                fout.write(seq + "\n")
                fout.write(plus)
                fout.write(qual + "\n")

filter_fastq("{input.r1}", "{output.r1_filt}")
filter_fastq("{input.r2}", "{output.r2_filt}")
EOF
        """

# Run a final fastqc on the final length filtered files.

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

# This rule generates a QC report from cutadapts reporting generated in the rule above.

rule summarize_primer_qc:
    input:
        reports = expand(
            "trim_clean_qc/trimmed_reports/{sample}_{amp}_cutadapt.txt",
            sample=SAMPLES,
            amp=["ITS1", "trnL"]
        )
    output:
        "trim_clean_qc/primer_qc/primer_qc_summary.tsv"
    run:
        import re
        import os
        from pathlib import Path

        os.makedirs("trim_clean_qc/primer_qc", exist_ok=True)

        # Header for the summary table
        with open(output[0], "w") as out:
            out.write(
                "sample\tamplicon\ttotal_pairs\tR1_primer_pct\tR2_primer_pct\t"
                "N_filtered_pct\tretained_pct\tperfect_matches\tmismatch1\tmismatch2\ttruncated\n"
            )

            # Parse each cutadapt report
            for report in input.reports:
                fname = Path(report).name                      # e.g. 23412_S1_ITS1_cutadapt.txt
                base = fname.replace("_cutadapt.txt", "")      # 23412_S1_ITS1
                sample, amp = base.rsplit("_", 1)              # sample = 23412_S1, amp = ITS1

                with open(report) as f:
                    text = f.read()

                # Extract key metrics using regex
                def extract(pattern, cast=float):
                    m = re.search(pattern, text)
                    return cast(m.group(1)) if m else None

                total_pairs = extract(r"Total read pairs processed:\s+([\d,]+)", lambda x: int(x.replace(",", "")))
                r1_pct = extract(r"Read 1 with adapter:\s+[\d,]+\s+\(([\d\.]+)%\)")
                r2_pct = extract(r"Read 2 with adapter:\s+[\d,]+\s+\(([\d\.]+)%\)")
                n_filtered = extract(r"Pairs with too many N:\s+[\d,]+\s+\(([\d\.]+)%\)")
                retained = extract(r"Pairs written \(passing filters\):\s+[\d,]+\s+\(([\d\.]+)%\)")

                # Primer mismatch table (first read only)
                # Get primer length and max.err for first read (Adapter 1)
                m = re.search(
                    r"=== First read: Adapter 1 ===.*?Length:\s+(\d+);.*?No\. of allowed errors:\s+(\d+)",
                    text,
                    flags=re.S,
                )
                if m:
                    primer_len = int(m.group(1))
                    max_err = int(m.group(2))
                else:
                    primer_len = None
                    max_err = None

                perfect = mismatch1 = mismatch2 = None

                if primer_len is not None and max_err is not None:
                    # Find the row in the "Overview of removed sequences" table
                    # that corresponds to the full primer length
                    row_match = re.search(
                        rf"\n{primer_len}\s+([0-9,]+)\s+[0-9\.]+\s+{max_err}\s+([0-9,\s]+)",
                        text,
                    )
                    if row_match:
                        # counts_str contains the error counts columns (space-separated)
                        counts_str = row_match.group(2).strip()
                        counts = [
                            int(x.replace(",", ""))
                            for x in counts_str.split()
                        ]
                        # First column = perfect matches (0 mismatches)
                        perfect = counts[0] if len(counts) > 0 else 0
                        # Second column = 1 mismatch (if present)
                        mismatch1 = counts[1] if len(counts) > 1 else 0
                        # Third column = 2 mismatches (if present)
                        mismatch2 = counts[2] if len(counts) > 2 else 0
                    else:
                        perfect = mismatch1 = mismatch2 = 0

                # Truncated primer matches (length 18 or 19 or 21)
                truncated = 0
                for length in ["18", "19", "21"]:
                    m = re.search(rf"\n{length}\s+([\d,]+)", text)
                    if m:
                        truncated += int(m.group(1).replace(",", ""))

                # Write row
                out.write(
                    f"{sample}\t{amp}\t{total_pairs}\t{r1_pct}\t{r2_pct}\t"
                    f"{n_filtered}\t{retained}\t{perfect}\t{mismatch1}\t{mismatch2}\t{truncated}\n"
                )

# Rule for generation a read retention table across pipeline steps.
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


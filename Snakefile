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
        expand("qc_demux/{sample}_{amp}_R1_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL", "unassigned"]),
        expand("qc_demux/{sample}_{amp}_R2_fastqc.html",
               sample=SAMPLES, amp=["ITS1", "trnL", "unassigned"]),

        # demuxed samples
        expand("demux/{sample}_{amp}_R1.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("demux/{sample}_{amp}_R2.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # read counts
        "read_counts/read_counts.tsv",

        # trimmed reads from cutadapt
        expand("trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trimmed/{sample}_{amp}_R2.primertrim.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("trimmed_reports/{sample}_{amp}_cutadapt.txt",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # cleaned reads from fastp
        expand("cleaned/{sample}_{amp}_R1.cleaned.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("cleaned/{sample}_{amp}_R2.cleaned.fastq.gz",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("cleaned_reports/{sample}_{amp}_fastp.html",
               sample=SAMPLES, amp=["ITS1", "trnL"]),
        expand("cleaned_reports/{sample}_{amp}_fastp.json",
               sample=SAMPLES, amp=["ITS1", "trnL"]),

        # primer QC summary
        "primer_qc/primer_qc_summary.tsv"


# This rule separates the fastq files by target amplicon (ITS1 or trnL). In this pipeline I refer to this step as "demuxing" even though it does not fit a precise definition of demultiplexing. This rule leaves a small percentage of unassigned reads that often very nearly match the expected primer sequence. It might be worth relaxing e or dropping ^ to see if we can include a few more reads, but it's pretty marginal.

rule demux_by_primer:
    conda: "envs/environment.yaml"
    input:
        r1 = "fastq/{sample}_R1_001.fastq.gz",
        r2 = "fastq/{sample}_R2_001.fastq.gz"
    output:
        ITS1_r1 = "demux/{sample}_ITS1_R1.fastq.gz",
        ITS1_r2 = "demux/{sample}_ITS1_R2.fastq.gz",
        trnL_r1 = "demux/{sample}_trnL_R1.fastq.gz",
        trnL_r2 = "demux/{sample}_trnL_R2.fastq.gz",
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
            -g ITS1=^{ITS1_F} \
            -g trnL=^{trnL_F} \
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

# This rule runs fastqc on the demultiplexed samples.

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

# This rule counts the number of reads in each fastq after demuxing. See the output file read_counts.tsv.

rule count_reads:
    input:
        raw = expand("fastq/{sample}_R{read}_001.fastq.gz",
                     sample=SAMPLES, read=[1,2]),
        demux = expand("demux/{sample}_{amp}_R{read}.fastq.gz",
                       sample=SAMPLES,
                       amp=["ITS1", "trnL", "unassigned"],
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

# This rule trims primers from demuxed fastq files and applies an N base filter. Zero tolerance N filtering can remove a lot of pairs from the small number of samples I have run so far. One possible compromise here could be to truncate low quality tails and then run N filtering. I need to look at the distribution of Ns in demuxed samples to determine if this could be useful. I'll look into this later.

rule trim_primers:
    conda: "envs/environment.yaml"
    input:
        r1 = "demux/{sample}_{amp}_R1.fastq.gz",
        r2 = "demux/{sample}_{amp}_R2.fastq.gz"
    output:
        r1_trim = "trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2_trim = "trimmed/{sample}_{amp}_R2.primertrim.fastq.gz",
        report = "trimmed_reports/{sample}_{amp}_cutadapt.txt"
    threads: 4
    params:
        F = lambda wc: {"ITS1": ITS1_F, "trnL": trnL_F}[wc.amp],
        R = lambda wc: {"ITS1": ITS1_R, "trnL": trnL_R}[wc.amp]
    shell:
        r"""
        mkdir -p trimmed trimmed_reports

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
        r1 = "trimmed/{sample}_{amp}_R1.primertrim.fastq.gz",
        r2 = "trimmed/{sample}_{amp}_R2.primertrim.fastq.gz"
    output:
        r1_clean = "cleaned/{sample}_{amp}_R1.cleaned.fastq.gz",
        r2_clean = "cleaned/{sample}_{amp}_R2.cleaned.fastq.gz",
        html     = "cleaned_reports/{sample}_{amp}_fastp.html",
        json     = "cleaned_reports/{sample}_{amp}_fastp.json"
    threads: 4
    shell:
        r"""
        mkdir -p cleaned cleaned_reports

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

# This rule generates a QC report from cutadapts reporting generated in the rule above.

rule summarize_primer_qc:
    input:
        reports = expand(
            "trimmed_reports/{sample}_{amp}_cutadapt.txt",
            sample=SAMPLES,
            amp=["ITS1", "trnL"]
        )
    output:
        "primer_qc/primer_qc_summary.tsv"
    run:
        import re
        import os
        from pathlib import Path

        os.makedirs("primer_qc", exist_ok=True)

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



rule fastqc_untrimmed:
    input:
        lambda wc: get_fastq(wc)[wc.read]
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_fastqc.zip"
    log:
        f"{LOGDIR}/fastqc_untrimmed/{{sample}}-{{unit}}-{{read}}.log"
    threads: get_resource("fastqc","threads")
    resources:
        mem_mb = get_resource("fastqc","mem_mb"),
        runtime = get_resource("fastqc","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}-{{read}}.fastqc.txt"
    wrapper:
        "v3.5.0/bio/fastqc"

rule fastqc_trimmed:
    input:
        lambda wc: get_trimmed_reads_qc(wc)[wc.read]
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_trimmed_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}-{{unit}}-{{read}}_trimmed_fastqc.zip"
    log:
        f"{LOGDIR}/fastqc_trimmed/{{sample}}-{{unit}}-{{read}}.log"
    threads: get_resource("fastqc","threads")
    resources:
        mem_mb = get_resource("fastqc","mem_mb"),
        runtime = get_resource("fastqc","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}-{{read}}.trimmed.fastqc.txt"
    wrapper:
        "v3.5.0/bio/fastqc"

rule samtools_stats:
    input:
        f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/qc/samtools-stats/{{sample}}-{{unit}}.txt"
    log:
        f"{LOGDIR}/samtools-stats/{{sample}}-{{unit}}.log"
    threads: get_resource("samtools_stats","threads")
    resources:
        mem_mb = get_resource("samtools_stats","mem_mb"),
        runtime = get_resource("samtools_stats","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.samtools_stats.txt"
    wrapper:
        "v3.5.0/bio/samtools/stats"

rule genome_dict:
    input:
        genome=config["ref"]["genome"],
        fai=f"{config['ref']['genome']}.fai"
    output:
        dict=re.sub(r"\.fa","",os.path.splitext(config["ref"]["genome"])[0]) + ".dict"
    log:
        f"{LOGDIR}/picard/genome_dict/CreateSequenceDictionary.log"
    params:
        java_opts = "-XX:ParallelGCThreads={}".format(get_resource("genome_dict","threads"))
    resources:
        mem_mb = get_resource("genome_dict","mem_mb"),
        runtime = get_resource("genome_dict","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/genome_dict.txt"
    wrapper:
        "v3.5.0/bio/picard/createsequencedictionary"

if "restrict_regions" in config["processing"]:
    rule bed_to_interval:
        input:
            bed=config["processing"]["restrict_regions"],
            dict=re.sub(r"\.fa.*","",os.path.splitext(config["ref"]["genome"])[0]) + ".dict"
        output:
            temp(f"{OUTDIR}/regions.intervals")
        log:
            f"{LOGDIR}/picard/bed_to_interval/bedtointervals.log"
        params:
            extra = "",
            java_opts = "-XX:ParallelGCThreads={}".format(get_resource("bed_to_interval","threads"))
        resources:
            mem_mb =  get_resource("bed_to_interval","mem_mb"),
            runtime = get_resource("bed_to_interval","runtime")
        benchmark:
            f"{LOGDIR}/benchmarks/bed_to_interval.txt"
        wrapper:
            "v3.5.0/bio/picard/bedtointervallist"

    rule picard_collect_hs_metrics:
        input:
            bam=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam",
            reference=config["ref"]["genome"],
            bait_intervals=f"{OUTDIR}/regions.intervals",
            target_intervals=f"{OUTDIR}/regions.intervals"
        output:
            f"{OUTDIR}/qc/picard/{{sample}}-{{unit}}.txt"
        log:
            f"{LOGDIR}/picard_collect_hs_metrics/{{sample}}-{{unit}}.log"
        params:
            extra = "",
            java_opts = "-XX:ParallelGCThreads={}".format(get_resource("picard_collect_hs_metrics","threads"))
        resources:
            mem_mb = get_resource("picard_collect_hs_metrics","mem_mb"),
            runtime = get_resource("picard_collect_hs_metrics","runtime")
        benchmark:
            f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.picard_collect_hs_metrics.txt"
        wrapper:
            "v3.5.0/bio/picard/collecthsmetrics"

rule multiqc:
    input:
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_fastqc.zip", r=["r1"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) == "nan")],
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_fastqc.zip", r=["r1","r2"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) != "nan")],
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_trimmed_fastqc.zip", r=["r1"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) == "nan")],
         [expand(f"{OUTDIR}/qc/fastqc/{row.sample}-{row.unit}-{{r}}_trimmed_fastqc.zip", r=["r1","r2"]) for row in units.itertuples() if (str(getattr(row, 'fq2')) != "nan")],
         expand(f"{OUTDIR}/qc/samtools-stats/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/dedup/{{u.sample}}-{{u.unit}}.metrics.txt", u=units.itertuples()),
         expand(f"{OUTDIR}/qc/picard/{{u.sample}}-{{u.unit}}.txt", u=units.itertuples()) if config["processing"].get("restrict_regions") else [],
         #expand(f"{OUTDIR}/snpeff/{{u.group}}.csv", u=samples.itertuples())
    output:
        report(f"{OUTDIR}/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    params:
        extra = "--config res/config/multiqc_config.yaml",
        use_input_files_only=True
    log:
        f"{LOGDIR}/multiqc.log"
    resources:
        threads = get_resource("multiqc","threads"),
        mem_mb = get_resource("multiqc","mem_mb"),
        runtime = get_resource("multiqc","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/multiqc.txt"
    wrapper:
        "v3.5.0/bio/multiqc"

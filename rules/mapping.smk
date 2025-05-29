rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    log:
        f"{LOGDIR}/trimmomatic/{{sample}}-{{unit}}.log"
    threads: get_resource("trim_reads","threads")
    resources:
        mem_mb = get_resource("trim_reads","mem_mb"),
        runtime = get_resource("trim_reads","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.trim_reads_se.txt"
    wrapper:
        "v3.5.0/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.fastq.gz"),
        r2=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.fastq.gz"),
        r1_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.unpaired.fastq.gz"),
        trimlog=temp(f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.trimlog.txt")
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        f"{LOGDIR}/trimmomatic/{{sample}}-{{unit}}.log"
    threads: get_resource("trim_reads","threads")
    resources:
        mem_mb = get_resource("trim_reads","mem_mb"),
        runtime = get_resource("trim_reads","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.trim_reads_pe.txt"
    wrapper:
        "v3.5.0/bio/trimmomatic/pe"

idx_cmd = "bwa index {input} > {log.out} 2> {log.err}"
rule bwa_idx_genome:
    shadow:"shallow"
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".0123",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".amb",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".ann",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".bwt.2bit.64",
        f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome'])+".pac"
    threads: get_resource("bwa_idx_genome","threads")
    resources:
        mem_mb = get_resource("bwa_idx_genome","mem_mb"),
        runtime = get_resource("bwa_idx_genome","runtime")
    log:
        f"{LOGDIR}/bwa_idx_genome/bwa_idx_genome.log"
    params:
        prefix=config['ref']['genome_idx']+os.path.basename(config['ref']['genome'])
    benchmark:
        f"{LOGDIR}/benchmarks/bwa_idx_genome.txt"
    wrapper:
        "v3.5.0/bio/bwa-mem2/index"

rule map_reads:
    input:
        r1=f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.1.fastq.gz",
        r2=f"{OUTDIR}/trimmed/{{sample}}-{{unit}}.2.fastq.gz",
        idx=multiext(f"{config['ref']['genome_idx']}"+os.path.basename(config['ref']['genome']), ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        fasta=f"{config['ref']['genome']}"    
    output:
        temp(f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam")
    log:
        f"{LOGDIR}/bwa_mem/{{sample}}-{{unit}}.log"
    params:
        extra=get_read_group,
    shadow: "shallow"
    threads: get_resource("map_reads","threads")
    resources:
        mem_mb = get_resource("map_reads","mem_mb"),
        runtime = get_resource("map_reads","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.map_reads.txt"
    conda:
        "../envs/bwa-picard.yml"
    shell:
        (
            "bwa-mem2 mem {params.extra} -t {threads} {input.fasta} {input.r1} {input.r2} 2> {log} | samtools sort -@ {threads} -O bam -o {output}"
        )

rule mark_duplicates:
    input:
        bams=f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    output:
        bam=temp(f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"),
        metrics=f"{OUTDIR}/qc/dedup/{{sample}}-{{unit}}.metrics.txt"
    log:
        f"{LOGDIR}/picard/dedup/{{sample}}-{{unit}}.log"
    params:
        extra=config["params"]["picard"]["MarkDuplicates"],
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("mark_duplicates","threads"))
    resources:
        mem_mb = get_resource("mark_duplicates","mem_mb"),
        runtime = get_resource("mark_duplicates","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.mark_duplicates.txt"
    wrapper:
        "v3.5.0/bio/picard/markduplicates"

checkpoint genome_faidx:
    input:
        config["ref"]["genome"]
    output:
        f"{config['ref']['genome']}.fai"
    log:
        f"{LOGDIR}/genome_faidx/genome_faidx.log"
    threads: get_resource("genome_faidx","threads")
    resources:
        mem_mb = get_resource("genome_faidx","mem_mb"),
        runtime = get_resource("genome_faidx","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/genome_faidx.txt"
    wrapper:
        "v3.5.0/bio/samtools/faidx"

rule obtain_recal_table:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        ref_idx=f"{config['ref']['genome']}.fai",
        dict=re.sub(r"\.fa.*","",os.path.splitext(config["ref"]["genome"])[0]) + ".dict",
        known=config["ref"]["known_variants"],
        known_idx=f"{config['ref']['known_variants']}.tbi"
    output:
        recal_table=f"{OUTDIR}/recal/{{sample}}-{{unit}}.grp"
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("obtain_recal_table","threads"))
    log:
        f"{LOGDIR}/gatk/bqsr/{{sample}}-{{unit}}.log"
    resources:
        mem_mb = get_resource("obtain_recal_table","mem_mb"),
        runtime = get_resource("obtain_recal_table","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.obtain_recal_table.txt"
    wrapper:
        "v3.5.0/bio/gatk/baserecalibrator"

rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        ref_idx=f"{config['ref']['genome']}.fai",
        recal_table=f"{OUTDIR}/recal/{{sample}}-{{unit}}.grp"
    output:
        bam=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bam",
        bai=f"{OUTDIR}/recal/{{sample}}-{{unit}}.bai"
    params:
        extra="",
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("recalibrate_base_qualities","threads"))
    log:
        f"{LOGDIR}/gatk/gatk_applybqsr/{{sample}}-{{unit}}.log"
    resources:
        mem_mb = get_resource("recalibrate_base_qualities","mem_mb"),
        runtime = get_resource("recalibrate_base_qualities","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.recalibrate_base_qualities.txt"
    wrapper:
        "v3.5.0/bio/gatk/applybqsr"

rule samtools_index:
    input:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bam"
    output:
        f"{OUTDIR}/dedup/{{sample}}-{{unit}}.bai"
    threads: get_resource("samtools_index","threads")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    resources:
        mem_mb = get_resource("samtools_index","mem_mb"),
        runtime = get_resource("samtools_index","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.samtools_index.txt"
    wrapper:
        "v3.5.0/bio/samtools/index"

rule samtools_index_sorted:
    input:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bam"
    output:
        f"{OUTDIR}/mapped/{{sample}}-{{unit}}.sorted.bai"
    threads: get_resource("samtools_index","threads")
    log:
        f"{LOGDIR}/samtools/index/{{sample}}-{{unit}}.log"
    resources:
        mem_mb = get_resource("samtools_index","mem_mb"),
        runtime = get_resource("samtools_index","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}-{{unit}}.samtools_index_sorted.txt"
    wrapper:
        "v3.5.0/bio/samtools/index"

rule index_known_variants:
    input:
        f"{config['ref']['known_variants']}"
    output:
        f"{config['ref']['known_variants']}.tbi"
    params:
        extra=""
    log:
        f"{LOGDIR}/gatk/index_known_variants.log"
    resources:
        threads = get_resource("index_known_variants","threads"),
        mem_mb = get_resource("index_known_variants","mem_mb"),
        runtime = get_resource("index_known_variants","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/index_known_variants.txt"
    wrapper:
        "v3.5.0/bio/bcftools/index"

rule snpeff_download:
    output:
        directory(f"{config['annotation']['snpeff']['cache_directory']}/{config['annotation']['snpeff']['database']}")
    log:
        f"{LOGDIR}/snpeff/download/snpeff_download.log"
    params:
        reference=f"{config['annotation']['snpeff']['database']}"
    resources:
        threads = get_resource("snpeff_download","threads"),
        mem_mb = get_resource("snpeff_download","mem_mb"),
        runtime = get_resource("snpeff_download","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/snpeff_download.txt"
    wrapper:
        "v3.5.0/bio/snpeff/download"

rule snpeff:
    input:
        calls = f"{OUTDIR}/filtered/{{group}}.vcf.gz",
        db = f"{config['annotation']['snpeff']['cache_directory']}/{config['annotation']['snpeff']['database']}"
    output:
        calls=report(f"{OUTDIR}/annotated/{{group}}.snpeff.vcf.gz",
        caption="../report/vcf.rst", category="Calls"),
        csvstats=f"{OUTDIR}/snpeff/{{group}}.csv",
        stats="snpeff/{{group}}.html"
    log:
        f"{LOGDIR}/snpeff/{{group}}.snpeff.log"
    params:
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("snpeff","threads")),
        extra=""
    resources:
        mem_mb = get_resource("snpeff","mem_mb"),
        runtime = get_resource("snpeff","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.snpeff.txt"
    wrapper:
        "v3.5.0/bio/snpeff/annotate"

rule get_vep_cache:
    output:
        directory(f"{config['annotation']['vep']['cache_directory']}/cache")
    params:
        species=f"{config['annotation']['vep']['species']}",
        build=f"{config['annotation']['vep']['assembly']}",
        release=f"{config['annotation']['vep']['cache_version']}"
    log:
        f"{LOGDIR}/vep/get_vep_cache.log"
    resources:
        threads = get_resource("get_vep_cache","threads"),
        mem_mb = get_resource("get_vep_cache","mem_mb"),
        runtime = get_resource("get_vep_cache","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/get_vep_cache.tsv"
    wrapper:
        "v3.5.0/bio/vep/cache"

rule download_vep_plugins:
    output:
        directory(f"{config['annotation']['vep']['cache_directory']}/plugins")
    params:
        release=f"{config['annotation']['vep']['cache_version']}"
    log:
        f"{LOGDIR}/vep/download_vep_plugins.log"
    resources:
        threads = get_resource("download_vep_plugins","threads"),
        mem_mb = get_resource("download_vep_plugins","mem_mb"),
        runtime = get_resource("download_vep_plugins","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/download_vep_plugin.tsv"
    wrapper:
        "v3.5.0/bio/vep/plugins"

rule vep_gatk:
    input:
        calls=f"{OUTDIR}/filtered/{{group}}.vcf.gz",
        cache=f"{config['annotation']['vep']['cache_directory']}/cache",
        plugins=f"{config['annotation']['vep']['cache_directory']}/plugins"
    output:
        calls=f"{OUTDIR}/annotated/{{group}}.vep.vcf.gz",
        stats=f"{OUTDIR}/annotated/{{group}}.vep.vcf.gz_summary.html"
    params:
        plugins=f"{config['annotation']['vep']['plugins']}",
        extra=f"{config['annotation']['vep']['extra']}"
    log:
        f"{LOGDIR}/vep/{{group}}.gatk_vep.log"
    threads: get_resource("vep_gatk","threads")
    resources:
        mem_mb = get_resource("vep_gatk","mem_mb"),
        runtime = get_resource("vep_gatk","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.vep_gatk.txt"
    wrapper:
        "v3.5.0/bio/vep/annotate"

rule vep_mutect:
    input:
        calls=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered_custom_minAF.vcf.gz",
        cache=f"{config['annotation']['vep']['cache_directory']}/cache",
        plugins=f"{config['annotation']['vep']['cache_directory']}/plugins"
    output:
        tsv=f"{OUTDIR}/annotated/{{sample}}_mutect.vep.tsv",
        stats=f"{OUTDIR}/annotated/{{sample}}_mutect.vep.tsv_summary.html"
    params:
        cache_version = config["annotation"]["vep"]["cache_version"],
        species = config["annotation"]["vep"]["species"],
        assembly = config["annotation"]["vep"]["assembly"],
        plugins=f"{config['annotation']['vep']['plugins']}",
        extra=f"{config['annotation']['vep']['extra']}"
    log:
        f"{LOGDIR}/vep/mutect_{{sample}}_vep.log"
    threads: get_resource("vep_mutect","threads")
    resources:
        mem_mb = get_resource("vep_mutect","mem_mb"),
        runtime = get_resource("vep_mutect","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.vep_mutect.txt"
    conda:
        "../envs/vep.yaml"
    shell:
        """
        vep \
            --input_file {input.calls} \
            --output_file {output.tsv} \
            --format vcf \
            --tab \
            --force_overwrite \
            --offline \
            --cache \
            --cache_version {params.cache_version} \
            --dir_cache {input.cache} \
            --species {params.species} \
            --assembly {params.assembly} \
            --dir_plugins {input.plugins} \
            {params.extra} \
            --stats_file {output.stats} \
            > {log} 2>&1
        """

rule filter_tumoralgene_variants:
    input:
        vcf="annotated/{sample}_vep.tsv",
        genes="res/Intogen.txt"
    output:
        tsv="annotated/{sample}_tumoralgene_vep.tsv"
    log:
        "log/filter_tumoralgene/{sample}.log"
    shell:
        """
        python scripts/filter_tumor_genes.py {input.vcf} {input.genes} {output.tsv} > {log} 2>&1
        """


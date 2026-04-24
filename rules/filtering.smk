def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/genotyped/{{group}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz")
    params:
        extra=get_vartype_arg,
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("select_calls","threads"))
    log:
        f"{LOGDIR}/gatk/selectvariants/{{group}}.{{vartype}}.log"
    resources:
        mem_mb = get_resource("select_calls","mem_mb"),
        runtime = get_resource("select_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.{{vartype}}.select_calls.txt"
    wrapper:
        "v3.5.0/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz"
    output:
        vcf=temp(f"{OUTDIR}/filtered/{{group}}.{{vartype}}.hardfiltered.vcf.gz")
    params:
        filters=get_filter,
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("hard_filter_calls","threads"))
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{group}}.{{vartype}}.log"
    resources:
        mem_mb = get_resource("hard_filter_calls","mem_mb"),
        runtime = get_resource("hard_filter_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.{{vartype}}.hard_filter_calls.txt"
    wrapper:
        "v3.5.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.vcf.gz",
        ref=config["ref"]["genome"],
        hapmap=config["params"]["gatk"]["VariantRecalibrator"]["hapmap"],
        omni=config["params"]["gatk"]["VariantRecalibrator"]["omni"],
        g1k=config["params"]["gatk"]["VariantRecalibrator"]["g1k"],
        dbsnp=config["params"]["gatk"]["VariantRecalibrator"]["dbsnp"],
        aux=config["params"]["gatk"]["VariantRecalibrator"]["aux"]
    output:
        vcf=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.recalibrated.vcf.gz",
        tranches=f"{OUTDIR}/filtered/{{group}}.{{vartype}}.tranches"
    params:
        mode=lambda wc: "SNP"
                if wc.vartype == "snvs"
                else "INDEL",
        resources=config["params"]["gatk"]["VariantRecalibrator"]["parameters"],
        annotation=config["params"]["gatk"]["VariantRecalibrator"]["annotation"],
        extra=config["params"]["gatk"]["VariantRecalibrator"]["extra"],
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("recalibrate_calls","threads"))
    log:
        f"{LOGDIR}/gatk/variantrecalibrator/{{group}}.{{vartype}}.log"
    resources:
        mem_mb = get_resource("recalibrate_calls","mem_mb"),
        runtime = get_resource("recalibrate_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.{{vartype}}.recalibrate_calls.txt"
    wrapper:
        "v3.5.0/bio/gatk/variantrecalibrator"

rule merge_calls:
    input:
        vcfs=lambda wc: expand(f"{OUTDIR}/filtered/{wc.group}.{{vartype}}.{{filtertype}}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf=f"{OUTDIR}/filtered/{{group}}.vcf.gz"
    params:
        extra = "",
        java_opts = "-XX:ParallelGCThreads={}".format(get_resource("merge_calls","threads"))
    log:
        f"{LOGDIR}/picard/{{group}}.merge-filtered.log"
    resources:
        mem_mb = get_resource("merge_calls","mem_mb"),
        runtime = get_resource("merge_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{group}}.merge_calls.txt"
    wrapper:
        "v3.5.0/bio/picard/mergevcfs"

rule learn_read_orientation_model:
    input:
        f1r2=f"{OUTDIR}/mutect/{{sample}}.f1r2.tar.gz"
    output:
        rom=f"{OUTDIR}/mutect/{{sample}}_read_orientation_model.tar.gz"
    log:
        f"{LOGDIR}/gatk/read_orientation_model.{{sample}}.log"
    params:
        extra="",
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("learn_read_orientation_model", "threads"))
    resources:
        mem_mb = get_resource("learn_read_orientation_model","mem_mb"),
        runtime = get_resource("learn_read_orientation_model","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.learn_read_orientation_model.txt"
    wrapper:
        "v3.5.0/bio/gatk/learnreadorientationmodel"

rule filter_mutect_calls:
    input:
        vcf=f"{OUTDIR}/mutect/{{sample}}.vcf.gz",
        ref=config["ref"]["genome"],
        f1r2=f"{OUTDIR}/mutect/{{sample}}_read_orientation_model.tar.gz" if config["filtering"]["mutect"]["lrom"] else []
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel.vcf.gz"
    params:
        extra=config["params"]["gatk"]["mutect"],
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("filter_mutect_calls", "threads"))
    log:
        f"{LOGDIR}/gatk/mutect_filter.{{sample}}.log"
    resources:
        mem_mb = get_resource("filter_mutect_calls","mem_mb"),
        runtime = get_resource("filter_mutect_calls","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.filter_mutect_calls.txt"
    wrapper:
        "v3.5.0/bio/gatk/filtermutectcalls"

rule filter_mutect_2:
    input:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel.vcf.gz",
        ref=config["ref"]["genome"]
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered.vcf.gz"
    params:
        filters={"DPfilter": config["filtering"]["mutect"]["depth"]},
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("filter_mutect_2","threads"))
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{sample}}_mutect.log"
    resources:
        mem_mb = get_resource("filter_mutect_2","mem_mb"),
        runtime = get_resource("filter_mutect_2","runtime")
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.filter_mutect_2.txt"
    wrapper:
        "v3.5.0/bio/gatk/variantfiltration"

rule filter_mutect_custom:
    input:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered.vcf.gz",
        ref=config["ref"]["genome"],
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered_custom.vcf.gz",
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{sample}}_mutect_custom.log",
    params:
        extra="--exclude-filtered",
        java_opts="-XX:ParallelGCThreads={}".format(get_resource("filter_mutect_custom","threads"))
    resources:
        mem_mb=get_resource("filter_mutect_custom","mem_mb"),
        runtime=get_resource("filter_mutect_custom","runtime"),
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.filter_mutect_custom.txt",
    wrapper:
        "v3.5.0/bio/gatk/selectvariants"

def get_min_af_filter(wildcards):
    threshold = config["filtering"]["bcftools_filtering"]["min_AF_threshold"]
    control = samples.loc[wildcards.sample, "control"]

    if control == wildcards.sample:
        af_index = "0:0"   # tumor-only
    else:
        af_index = "1:0"   # tumor-normal

    return f'-i "FORMAT/AF[{af_index}] >= {threshold}"'

rule filter_minAF:
    input:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered_custom.vcf.gz",
    output:
        vcf=f"{OUTDIR}/mutect_filter/{{sample}}_passlabel_filtered_custom_minAF.vcf.gz",
    log:
        f"{LOGDIR}/gatk/variantfiltration/{{sample}}_mutect_custom_minAF.log",
    params:
        extra = get_min_af_filter,
    resources:
        mem_mb=get_resource("filter_AF","mem_mb"),
        runtime=get_resource("filter_AF","runtime"),
    benchmark:
        f"{LOGDIR}/benchmarks/{{sample}}.filter_mutect_custom_minAF.txt",
    wrapper:
        "v3.5.0/bio/bcftools/view"

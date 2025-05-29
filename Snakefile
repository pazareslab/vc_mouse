import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
min_version("8.6.0")

singularity: "docker://continuumio/miniconda3:4.6.14"

report: "report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]

samples = pd.read_csv(config["samples"],sep="\t").set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"],sep="\t", dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

include: "rules/common.smk"

##### Warnings #####
class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)

if config['ref']['known_variants'].endswith(".vcf") == True:
    warning(f"ERROR: Known variants file is not compressed. It has to be provided compressed in bgzip format.")
    sys.exit(1)

warning("Check if your genome contains alternative contigs and configure the execution accordingly in the ref section of the config file.")

warning("After each rule is executed, it generates a benchmark report. You can find them in the benchmarks directory within the logs directory. Please, adjust the parameters in the config file accordingly for subsequent executions.")

##### Target rules #####

rule all:
    input:
        #["{OUTDIR}/annotated/{group}.snpeff.vcf.gz".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        #["{OUTDIR}/annotated/{group}.vep.vcf.gz".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        ["{OUTDIR}/annotated/{sample}_mutect.vep.tsv".format(OUTDIR=OUTDIR,sample=getattr(row, 'sample')) for row in samples.itertuples() if (getattr(row, 'control') != "-")],
        f"{OUTDIR}/qc/multiqc.html",
        #[f"{OUTDIR}/plots/{{group}}.depths.svg".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()],
        #[f"{OUTDIR}/plots/{{group}}.allele-freqs.svg".format(OUTDIR=OUTDIR,group=getattr(row, 'group')) for row in samples.itertuples()]


##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"

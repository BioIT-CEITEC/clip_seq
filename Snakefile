import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
#GLOBAL_REF_PATH = "/mnt/references"
#GLOBAL_TMPD_PATH = "tmp"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as BR_*

##### Config processing #####

config = BR.load_organism()

sample_tab = BR.load_sample()

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),

##### Target rules #####
rule all:
    input: "CLIPseq_analysis_report.html"

##### Modules #####

include: "rules/CLIP-seq.smk"
# include: "rules/prepare_reference.smk"


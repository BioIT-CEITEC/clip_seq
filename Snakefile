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

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["reference"] = "GRCh38-p10"
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")


##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),

##### Target rules #####
rule all:
    input: "CLIPseq_analysis_report.html"

##### Modules #####

include: "rules/CLIP-seq.smk"
# include: "rules/prepare_reference.smk"


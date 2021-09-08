import math
import subprocess
import json
import re
import os.path
import pandas as pd
from snakemake.utils import R
from snakemake.utils import report
from os.path import split
from helper import define_variable

include: "prepare_reference.smk"



###########################################
# DEFINITION OF VARIABLES
#

cfg = pd.DataFrame(config)

REF_DIR = define_variable(cfg, "REF_DIR")
PROJECT_NAME = define_variable(cfg, "PROJECT_NAME")
PROJECT_DIR = define_variable(cfg, "PROJECT_DIR")
INPUTS_DIR = define_variable(cfg, "INPUTS_DIR")
ADIR = define_variable(cfg, "ANALYSIS_DIR")
SAMPLE = define_variable(cfg, "FULL_NAME")
SAMPLE_REP = define_variable(cfg, "FULL_NAME_REP")
SMK_HOME = define_variable(cfg, "SMK_HOME")




##########################################
# FINAL RULE
#

rule all:
    input:  report = ADIR + "/" + PROJECT_NAME + "_finished"
    
def specify_inputs_for_final_report(wildcards):
    sample_list = []
    for x in set(cfg.condition.tolist()):
        if len(cfg.loc[cfg.condition == x, "tag"].tolist()) > 1:
            sample_list += [x+"_"+y for y in cfg.loc[cfg.condition == x, "tag"].tolist()+["merged"]]
        else:
            sample_list += [x]
    # print(sample_list)
    inputs = []
    for caller in config["callers"].split(";"):
        inputs.append(["RCAS/from_"+caller+"/"+x+"/"+x+".uniq_reads.no_dups.RCAS_report.html" for x in sample_list])
        inputs.append(["annotated_beds/from_"+caller+"/"+x+"/"+x+".uniq_reads.no_dups.annotated.bed" for x in sample_list])
        if config['keep_dups']:
            inputs.append(["RCAS/from_"+caller+"/"+x+"/"+x+".uniq_reads.keep_dups.RCAS_report.html" for x in sample_list])
            inputs.append(["annotated_beds/from_"+caller+"/"+x+"/"+x+".uniq_reads.keep_dups.annotated.bed" for x in sample_list])
        if config['multimapped']:
            inputs.append(["RCAS/from_"+caller+"/"+x+"/"+x+".all_reads.no_dups.RCAS_report.html" for x in sample_list])
            inputs.append(["annotated_beds/from_"+caller+"/"+x+"/"+x+".all_reads.no_dups.annotated.bed" for x in sample_list])
            if config['keep_dups']:
                inputs.append(["RCAS/from_"+caller+"/"+x+"/"+x+".all_reads.keep_dups.RCAS_report.html" for x in sample_list])
                inputs.append(["annotated_beds/from_"+caller+"/"+x+"/"+x+".all_reads.keep_dups.annotated.bed" for x in sample_list])
    print(inputs)
    return inputs
    
rule final_report:
    input:  unpack(specify_inputs_for_final_report),
    output: report = ADIR + "/" + PROJECT_NAME + "_finished",
    shell: "touch {output.report}"
    #shell: ""
    
    
#########################################
# Peak calling
#

def annotate_peaks_input(wildcards):
    inputs = dict()
    if wildcards.caller == "pureClip":
        inputs['bed'] = "pureClip/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".binding_regions.bed"
    elif wildcards.caller == "CLAM":
        inputs['bed'] = "CLAM/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+"/narrow_peak.permutation.processed.bed"
    elif wildcards.caller == "macs2":
        inputs['bed'] = "macs2/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".peaks.narrowPeak"
    inputs['gtf'] = expand("{ref-dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0]
    return inputs

rule annotate_peaks:
    input:  unpack(annotate_peaks_input),
    output: bed = "annotated_beds/from_{caller}/{sample}/{sample}.{multi}.{dups}.annotated.bed",
    log:    run = "annotated_beds/from_{caller}/{sample}/{sample}.{multi}.{dups}.annotate_peaks.log",
    resources: mem=10 if config["organism"] == "homsap" else 5
    params: rscript = workflow.basedir+"/../scripts/annotate_bed_file.R",
            feat_type=config["feat_type"],
            annotate_by = config["annotate_by"],
    conda:  "../wrappers/annotate_peaks/env.yaml"
    script: "../wrappers/annotate_peaks/script.py"
    

def post_process_by_RCAS_input(wildcards):
    inputs = dict()
    if wildcards.caller == "pureClip":
        inputs['bed'] = "pureClip/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".binding_regions.bed"
    elif wildcards.caller == "CLAM":
        inputs['bed'] = "CLAM/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+"/narrow_peak.permutation.processed.bed"
    elif wildcards.caller == "macs2":
        inputs['bed'] = "macs2/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".peaks.narrowPeak"
    inputs['gtf'] = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0]
    inputs['msigdb'] = expand("{ref_dir}/other/MSigDB_for_RCAS/c2.all.v7.1.entrez.gmt", ref_dir=reference_directory,ref=config["reference"])[0]
    return inputs

rule post_process_by_RCAS:
    input:  unpack(post_process_by_RCAS_input),
    output: html    = "RCAS/from_{caller}/{sample}/{sample}.{multi}.{dups}.RCAS_report.html",
            tmp_bed = "RCAS/from_{caller}/{sample}/{sample}.{multi}.{dups}.input.bed",
    log:    run     = "RCAS/from_{caller}/{sample}/{sample}.{multi}.{dups}.post_process_by_RCAS.log",
    params: organism = config["organism"],
            dir = "RCAS/from_{caller}/{sample}/",
            html= "RCAS/from_{caller}/{sample}/{sample}.{multi}.{dups}.input.bed.RCAS.report.html",
            rscript= workflow.basedir+"/../wraps/CLIP-seq_analysis/post_process_by_RCAS/RCAS_script.R",
    conda:  "../wrappers/post_process_by_RCAS/env.yaml"
    script: "../wrappers/post_process_by_RCAS/script.py"
    
    
rule CLAM_postprocess:
    input:  bed = "CLAM/{name}.{multi}.{dups}/narrow_peak.permutation.bed",
    output: bed = "CLAM/{name}.{multi}.{dups}/narrow_peak.permutation.processed.bed",
    log:    run = "CLAM/{name}.{multi}.{dups}.CLAM_postprocess.log",
    shell: "awk -F '\\t' 'BEGIN {{OFS = FS}}{{split($4,n,\",\"); split(n[1],v,\":\"); split($5,s,\",\"); $4=v[1]\":\"v[2]\":\"s[1]; $5=v[3]; split($6,t,\",\"); $6=t[1]; print $0}}' {input.bed} > {output.bed} 2>> {log.run}" 
        # """
        # CMD="awk -F '\\t' 'BEGIN {{OFS = FS}}{{split($4,n,\",\"); split(n[1],v,\":\"); split($5,s,\",\"); $4=v[1]\":\"v[2]\":\"s[1]; $5=v[3]; split($6,t,\",\"); $6=t[1]; print $0}}' {input.bed} > {output.bed} 2>> {log.run}"
        # echo "## COMMAND: $CMD" > {log.run}
        # eval "$CMD"
        # """
    
    
rule call_CLAM:
    input:  uniq_bam = "CLAM/{name}.{multi}.{dups}/unique.sorted.bam",
            mult_bam ="CLAM/{name}.{multi}.{dups}/realigned.sorted.bam",
            gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: bed = "CLAM/{name}.{multi}.{dups}/narrow_peak.permutation.bed",
    log:    run = "CLAM/{name}.{multi}.{dups}.call_CLAM.log",
    threads: 20
    params: strand =  config["strandness"], # strandness
            qval_cutoff = config["qval_cutoff"], # adjusted p-values cutoff [float: 0-1]
            merge_size = config["merge_size"], # Select best peak within this size [integer]
            extend_peak = config["extend_peak"], # Extend peak to this size if less than it [integer]
            dir = "CLAM/{name}.{multi}.{dups}/",
    conda:  "../wrappers/call_CLAM/env.yaml",
    script: "../wrappers/call_CLAM/script.py"
    
rule call_CLAM_preprocess:
    input:  bam = "bams/{name}.{multi}.{dups}.bam",
            bai = "bams/{name}.{multi}.{dups}.bam.bai",
    output: uniq_bam = "CLAM/{name}.{multi}.{dups}/unique.sorted.bam",
            mult_bam = "CLAM/{name}.{multi}.{dups}/realigned.sorted.bam",
    log:    run = "CLAM/{name}.{multi}.{dups}.call_CLAM_preprocess.log",
    threads: 1
    resources: mem=50
    params: strand =  config["strandness"], # strandness
            max_multi_hits = config["max_multi_hits"], # maximum hits allowed for multi-mapped reads [integer]
            read_tagger = config["read_tagger"], # read tagger method, 'median' for read center, 'start' for read start site ['median', 'start']
            dir = "CLAM/{name}.{multi}.{dups}/",
    conda:  "../wrappers/call_CLAM/env.yaml", # it would be the exact same env as in rule call_CLAM but as it was tough to create one we will use the same
    script: "../wrappers/call_CLAM_preprocess/script.py"
    

rule call_macs2:
    input:  bam = "bams/{name}.{multi}.{dups}.bam",
            bai = "bams/{name}.{multi}.{dups}.bam.bai",
            chrs= expand("{ref_dir}/seq/{ref}.chrom.sizes", ref_dir=reference_directory,ref=config["reference"])[0],
    output: trt_bdg = "macs2/{name}/{name}.{multi}.{dups}.bdg",
            trt_bwg = "macs2/{name}/{name}.{multi}.{dups}.bigWig",
            # ctl_bdg = ADIR+"/macs2/{name}/{name}.{multi}.{dups}.control.bdg",
            # ctl_bwg = ADIR+"/macs2/{name}/{name}.{multi}.{dups}.control.bigWig",
            xls_tab = "macs2/{name}/{name}.{multi}.{dups}.peaks.xls",
            sum_tab = "macs2/{name}/{name}.{multi}.{dups}.summits.bed",
            nar_tab = "macs2/{name}/{name}.{multi}.{dups}.peaks.narrowPeak",
    log:    run = "macs2/{name}.{multi}.{dups}.call_macs2.log",
    threads: 1
    params: trt_bdg = "macs2/{name}/{name}.{multi}.{dups}_treat_pileup.bdg",
            trt_bwg = "macs2/{name}/{name}.{multi}.{dups}_treat_pileup.bigWig",
            # ctl_bdg = ADIR+"/macs2/{name}/{name}.{multi}.{dups}_control_lambda.bdg",
            # ctl_bwg = ADIR+"/macs2/{name}/{name}.{multi}.{dups}_control_lambda.bigWig",
            xls_tab = "macs2/{name}/{name}.{multi}.{dups}_peaks.xls",
            sum_tab = "macs2/{name}/{name}.{multi}.{dups}_summits.bed",
            nar_tab = "macs2/{name}/{name}.{multi}.{dups}_peaks.narrowPeak",
            effective_GS = config["effective_GS"], #effective genome size stored in DB or references
            # frag_len = "unk", #sequencing fragment length
            qval_cutof = config["qval_cutof"], #q-value cutof
            dir = "macs2/{name}/",
            name= "{name}.{multi}.{dups}",
            temp= "/mnt/ssd/ssd_1/tmp/",
    conda:  "../wrappers/call_macs2/env.yaml"
    script: "../wrappers/call_macs2/script.py"


rule call_pureClip:
    input:  bam = "bams/{name}.{multi}.{dups}.bam",
            bai = "bams/{name}.{multi}.{dups}.bam.bai",
            gen = expand("{ref_dir}/seq/{ref}.fasta.gz", ref_dir=reference_directory,ref=config["reference"])[0],
    output: bed = "pureClip/{name}/{name}.{multi}.{dups}.crosslink_sites.bed",
            bed2= "pureClip/{name}/{name}.{multi}.{dups}.binding_regions.bed",
    log:    run = "pureClip/{name}/{name}.{multi}.{dups}.call_pureClip.log",
    threads: 10
    conda:  "../wrappers/call_pureClip/env.yaml"
    script: "../wrappers/call_pureClip/script.py"


# ##########################################
# # QC AFTER MAPPING
# #
# 
# rule inspect_bam_coverage:
#     input:  bam = ADIR+"/bams/{full_name}.bam",
#     output: bw = ADIR+"/bam_QC/{full_name}.bam_cov.bigWig",
#             bw2= ADIR+"/bam_QC/{full_name}.bam_cov.no_dups.bigWig",
#     log:    run = ADIR+"/bam_QC/{full_name}.inspect_bam_coverage.log",
#     threads: 5
#     params: effective_GS = 12000000, #TODO: should be project pramater for effective genome size stored in DB or references
#             frag_len = 450, #TODO: should be project parameter for sequencing fragment length
#     conda:  "../wraps/CLIP-seq_analysis/inspect_bam_coverage/env.yaml"
#     script:  "../wraps/CLIP-seq_analysis/inspect_bam_coverage/script.py"
# 
# rule multi_bam_summary:
#     input:  bam = [ADIR+"/bams/"+x+".bam" for x in cfg.full_name.tolist()],
#     output: plot = ADIR+"/bam_QC/correlation_heatmap.pdf",
#             table =ADIR+"/bam_QC/correlation_heatmap.tsv",
#             plot2 =ADIR+"/bam_QC/correlation_heatmap.no_dups.pdf",
#             table2=ADIR+"/bam_QC/correlation_heatmap.no_dups.tsv",
#             finger = ADIR+"/bam_QC/fingerprint.pdf",
#             finger2 =ADIR+"/bam_QC/fingerprint.no_dups.pdf",
#     log:    run = ADIR+"/bam_QC/multi_bam_summary.log",
#     threads: 20
#     params: matrix = ADIR+"/bam_QC/multi_bam_summary.mtx",
#             matrix2= ADIR+"/bam_QC/multi_bam_summary.no_dups.mtx",
#             corr_method = "spearman", #TODO: should be project pramater ['spearman', 'pearson']
#     conda:  "../wraps/CLIP-seq_analysis/multi_bam_summary/env.yaml"
#     script:  "../wraps/CLIP-seq_analysis/multi_bam_summary/script.py"
#     # shell:
#     #   """
#     #   multiBamSummary bins -bs 10000 --smartLabels -p {threads} -o {params.matrix} -b {input.bam} > {log.run} 2>&1
#     #   plotCorrelation -in {params.matrix} --whatToPlot heatmap --corMethod {params.corr_method} -o {output.plot} --outFileCorMatrix {output.table} --plotNumbers >> {log.run} 2>&1
#     #   """


##########################################
# PREPARE INPUT
#

def process_bams_input(wildcards):
    if not "_merged" in wildcards.full_name:
        return {
            'bam': "mapped/"+wildcards.full_name+".bam",
        }
    else:
        return {
            'bam': ["mapped/"+x+".bam" for x in cfg.loc[cfg.condition == re.sub('_merged$', '', wildcards.full_name), "full_name"].tolist()],
        }
        
rule process_bams:
    input:  unpack(process_bams_input),
    output: bam = "bams/{full_name}.{multi}.{dups}.bam",
            bai = "bams/{full_name}.{multi}.{dups}.bam.bai",
    log:    "bams/{full_name}.{multi}.{dups}.process_bams.log"
    threads: 5,
    params: quality_cutof = config["quality_cutof"],
            tmp_bam = "bams/{full_name}.{multi}.{dups}.tmp.bam",
    conda:  "../wrappers/process_bams/env.yaml"
    script: "../wrappers/process_bams/script.py"


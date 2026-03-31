import math
import subprocess
import json
import re
import os.path
import pandas as pd
from snakemake.utils import R
from snakemake.utils import report
from os.path import split


##########################################
# FINAL RULE
#
    
def specify_inputs_for_final_report(wildcards):
    inputs = list()
    callers = config["callers"].split(";")
    pureclip_w_ctl = ['G3BP1_oxid_rep1_VS_G3BP1_ctrl_rep1', 'G3BP1_oxid_rep2_VS_G3BP1_ctrl_rep2']
    pureclip_multi_sample = ['G3BP1_oxid', 'G3BP1_ctrl']
    inputs.append(expand("results/RCAS/{name}/{name}.uniq_reads.no_dups.from_{caller}/{name}.uniq_reads.no_dups.from_{caller}.RCAS_report.html", caller=callers, name=sample_tab.sample_name))
    inputs.append(expand("results/annotated_beds/{name}/{name}.uniq_reads.no_dups.from_{caller}.annotated.bed", caller=callers, name=sample_tab.sample_name))
    inputs.append(expand("results/pureClip_w_ctl/{name}/{name}.uniq_reads.no_dups.binding_regions.bed", name=pureclip_w_ctl))
    inputs.append(expand("results/pureClip/{name}/{name}.uniq_reads.no_dups.binding_regions.bed", name=pureclip_multi_sample))
    if config['keep_dups']:
        inputs.append(expand("results/RCAS/{name}/{name}.uniq_reads.keep_dups.from_{caller}/{name}.uniq_reads.keep_dups.from_{caller}.RCAS_report.html", caller=callers, name=sample_tab.sample_name))
        inputs.append(expand("results/annotated_beds/{name}/{name}.uniq_reads.keep_dups.from_{caller}.annotated.bed", caller=callers, name=sample_tab.sample_name))
        inputs.append(expand("results/pureClip_w_ctl/{name}/{name}.uniq_reads.keep_dups.binding_regions.bed", name=pureclip_w_ctl))
        inputs.append(expand("results/pureClip/{name}/{name}.uniq_reads.keep_dups.binding_regions.bed", name=pureclip_multi_sample))
    if config['multimapped']:
        inputs.append(expand("results/RCAS/{name}/{name}.all_reads.no_dups.from_{caller}/{name}.all_reads.no_dups.from_{caller}.RCAS_report.html", caller=callers, name=sample_tab.sample_name))
        inputs.append(expand("results/annotated_beds/{name}/{name}.all_reads.no_dups.from_{caller}.annotated.bed", caller=callers, name=sample_tab.sample_name))
        inputs.append(expand("results/pureClip_w_ctl/{name}/{name}.all_reads.no_dups.binding_regions.bed", name=pureclip_w_ctl))
        inputs.append(expand("results/pureClip/{name}/{name}.all_reads.no_dups.binding_regions.bed", name=pureclip_multi_sample))
        if config['keep_dups']:
            inputs.append(expand("results/RCAS/{name}/{name}.all_reads.keep_dups.from_{caller}/{name}.all_reads.keep_dups.from_{caller}.RCAS_report.html", caller=callers, name=sample_tab.sample_name))
            inputs.append(expand("results/annotated_beds/{name}/{name}.all_reads.keep_dups.from_{caller}.annotated.bed", caller=callers, name=sample_tab.sample_name))
            inputs.append(expand("results/pureClip_w_ctl/{name}/{name}.all_reads.keep_dups.binding_regions.bed", name=pureclip_w_ctl))
            inputs.append(expand("results/pureClip/{name}/{name}.all_reads.keep_dups.binding_regions.bed", name=pureclip_multi_sample))
    #print(inputs)
    return inputs
    
rule final_report:
    input:  unpack(specify_inputs_for_final_report),
    output: report = "CLIPseq_analysis_report.html",
    shell: "touch {output.report}"
    #shell: ""
    
    
#########################################
# Peak calling
#

def annotate_peaks_input(wildcards):
    inputs = dict()
    if wildcards.caller == "pureClip":
        inputs['bed'] = "results/pureClip/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".binding_regions.bed"
    elif wildcards.caller == "CLAM":
        inputs['bed'] = "results/CLAM/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+"/narrow_peak.permutation.processed.bed"
    elif wildcards.caller == "macs2":
        inputs['bed'] = "results/macs2/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".peaks.narrowPeak"
    inputs['gtf'] = expand("{ref_dir}/annot/{ref}.gtf", ref_dir=reference_directory, ref=config["reference"])[0]
    return inputs

rule annotate_peaks:
    input:  unpack(annotate_peaks_input),
    output: bed = "results/annotated_beds/{sample}/{sample}.{multi}.{dups}.from_{caller}.annotated.bed",
    log:    run = "logs/{sample}/{sample}.{multi}.{dups}.annotate_peaks.from_{caller}.log",
    resources: mem=10 if config["organism"] == "homo_sapiens" else 5
    params: rscript = workflow.basedir+"/wrappers/annotate_peaks/annotate_peaks.R",
            feat_type=config["feat_type"],
            annotate_by = config["annotate_by"],
    conda:  "../wrappers/annotate_peaks/env.yaml"
    script: "../wrappers/annotate_peaks/script.py"
    

def post_process_by_RCAS_input(wildcards):
    inputs = dict()
    if wildcards.caller == "pureClip":
        inputs['bed'] = "results/pureClip/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".binding_regions.bed"
    elif wildcards.caller == "CLAM":
        inputs['bed'] = "results/CLAM/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+"/narrow_peak.permutation.processed.bed"
    elif wildcards.caller == "macs2":
        inputs['bed'] = "results/macs2/"+wildcards.sample+"/"+wildcards.sample+"."+wildcards.multi+"."+wildcards.dups+".peaks.narrowPeak"
    inputs['gtf'] = expand("{ref_dir}/annot/{ref}.gtf", ref_dir=reference_directory, ref=config["reference"])[0]
    #inputs['msigdb'] = expand("{ref_dir}/other/MSigDB_for_RCAS/c2.all.v7.1.entrez.gmt", ref_dir=reference_directory, ref=config["reference"])[0]
    return inputs

rule post_process_by_RCAS:
    input:  unpack(post_process_by_RCAS_input),
    output: html    = "results/RCAS/{sample}/{sample}.{multi}.{dups}.from_{caller}/{sample}.{multi}.{dups}.from_{caller}.RCAS_report.html",
            bed = "results/RCAS/{sample}/{sample}.{multi}.{dups}.from_{caller}/{sample}.{multi}.{dups}.from_{caller}.input.bed",
    log:    run     = "logs/{sample}/{sample}.{multi}.{dups}.post_process_by_RCAS.from_{caller}.log",
    params: organism = config["organism"],
            dir = "results/RCAS/{sample}/{sample}.{multi}.{dups}.from_{caller}/",
            html= "results/RCAS/{sample}/{sample}.{multi}.{dups}.from_{caller}/{sample}.{multi}.{dups}.from_{caller}.input.bed.RCAS.report.html",
            rscript= workflow.basedir+"/wrappers/post_process_by_RCAS/RCAS_script.R",
    conda:  "../wrappers/post_process_by_RCAS/env.yaml"
    script: "../wrappers/post_process_by_RCAS/script.py"
    
    
rule call_CLAM_postprocess:
    input:  bed = "results/CLAM/{name}/{name}.{multi}.{dups}/narrow_peak.permutation.bed",
    output: bed = "results/CLAM/{name}/{name}.{multi}.{dups}/narrow_peak.permutation.processed.bed",
    log:    run = "logs/{name}/{name}.{multi}.{dups}.CLAM_postprocess.log",
    shell: "awk -F '\\t' 'BEGIN {{OFS = FS}}{{split($4,n,\",\"); split(n[1],v,\":\"); split($5,s,\",\"); $4=v[1]\":\"v[2]\":\"s[1]; $5=v[3]; split($6,t,\",\"); $6=t[1]; print $0}}' {input.bed} > {output.bed} 2>> {log.run}" 
    
    
rule call_CLAM:
    input:  uniq_bam = "results/CLAM/{name}/{name}.{multi}.{dups}/unique.sorted.bam",
            mult_bam ="results/CLAM/{name}/{name}.{multi}.{dups}/realigned.sorted.bam",
            gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
    output: bed = "results/CLAM/{name}/{name}.{multi}.{dups}/narrow_peak.permutation.bed",
    log:    run = "logs/{name}/{name}.{multi}.{dups}.call_CLAM.log",
    threads: 20
    params: strand =  config["strandness"], # strandness
            qval_cutoff = config["qval_cutof"], # adjusted p-values cutoff [float: 0-1]
            merge_size = config["merge_size"], # Select best peak within this size [integer]
            extend_peak = config["extend_peak"], # Extend peak to this size if less than it [integer]
            dir = "results/CLAM/{name}/{name}.{multi}.{dups}/",
    conda:  "../wrappers/call_CLAM/env.yaml",
    script: "../wrappers/call_CLAM/script.py"
    
    
rule call_CLAM_preprocess:
    input:  bam = "mapped/{name}.{multi}.{dups}.bam",
            bai = "mapped/{name}.{multi}.{dups}.bam.bai",
            tool_ok = "results/CLAM/CLAM_installed",
    output: uniq_bam = "results/CLAM/{name}/{name}.{multi}.{dups}/unique.sorted.bam",
            mult_bam = "results/CLAM/{name}/{name}.{multi}.{dups}/realigned.sorted.bam",
    log:    run = "logs/{name}/{name}.{multi}.{dups}.call_CLAM_preprocess.log",
    threads: 1
    resources: mem=50
    params: strand =  config["strandness"], # strandness
            max_multi_hits = config["max_multi_hits"], # maximum hits allowed for multi-mapped reads [integer]
            read_tagger = config["read_tagger"], # read tagger method, 'median' for read center, 'start' for read start site ['median', 'start']
            dir = "results/CLAM/{name}/{name}.{multi}.{dups}/",
    conda:  "../wrappers/call_CLAM/env.yaml", # it would be the exact same env as in rule call_CLAM but as it was tough to create one we will use the same
    script: "../wrappers/call_CLAM_preprocess/script.py"
    
rule install_CLAM:
    output: "results/CLAM/CLAM_installed",
    log:    "logs/install_CLAM.log",
    conda:  "../wrappers/call_CLAM/env.yaml",
    script: "../wrappers/install_CLAM/script.py"
    

rule call_macs2:
    input:  bam = "mapped/{name}.{multi}.{dups}.bam",
            bai = "mapped/{name}.{multi}.{dups}.bam.bai",
            chrs= expand("{ref_dir}/seq/{ref}.chrom.sizes", ref_dir=reference_directory, ref=config["reference"])[0],
    output: trt_bdg = "results/macs2/{name}/{name}.{multi}.{dups}.bdg",
            trt_bwg = "results/macs2/{name}/{name}.{multi}.{dups}.bigWig",
            # ctl_bdg = ADIR+"/results/macs2/{name}/{name}.{multi}.{dups}.control.bdg",
            # ctl_bwg = ADIR+"/results/macs2/{name}/{name}.{multi}.{dups}.control.bigWig",
            xls_tab = "results/macs2/{name}/{name}.{multi}.{dups}.peaks.xls",
            sum_tab = "results/macs2/{name}/{name}.{multi}.{dups}.summits.bed",
            nar_tab = "results/macs2/{name}/{name}.{multi}.{dups}.peaks.narrowPeak",
    log:    run = "logs/{name}/{name}.{multi}.{dups}.call_macs2.log",
    threads: 1
    params: trt_bdg = "results/macs2/{name}/{name}.{multi}.{dups}_treat_pileup.bdg",
            trt_bwg = "results/macs2/{name}/{name}.{multi}.{dups}_treat_pileup.bigWig",
            # ctl_bdg = ADIR+"/results/macs2/{name}/{name}.{multi}.{dups}_control_lambda.bdg",
            # ctl_bwg = ADIR+"/results/macs2/{name}/{name}.{multi}.{dups}_control_lambda.bigWig",
            xls_tab = "results/macs2/{name}/{name}.{multi}.{dups}_peaks.xls",
            sum_tab = "results/macs2/{name}/{name}.{multi}.{dups}_summits.bed",
            nar_tab = "results/macs2/{name}/{name}.{multi}.{dups}_peaks.narrowPeak",
            effective_GS = config["effective_GS"], #effective genome size stored in DB or references
            # frag_len = "unk", #sequencing fragment length
            qval_cutof = config["qval_cutof"], #q-value cutof
            dir = "results/macs2/{name}/",
            name= "{name}.{multi}.{dups}",
            temp= GLOBAL_TMPD_PATH,
    conda:  "../wrappers/call_macs2/env.yaml"
    script: "../wrappers/call_macs2/script.py"


# rule call_pureClip:
#     input:  bam = "mapped/{name}.{multi}.{dups}.bam",
#             bai = "mapped/{name}.{multi}.{dups}.bam.bai",
#             gen = expand("{ref_dir}/seq/{ref}.fasta.gz", ref_dir=reference_directory, ref=config["reference"])[0],
#     output: bed = "results/pureClip/{name}/{name}.{multi}.{dups}.crosslink_sites.bed",
#             bed2= "results/pureClip/{name}/{name}.{multi}.{dups}.binding_regions.bed",
#     log:    run = "logs/{name}/{name}.{multi}.{dups}.call_pureClip.log",
#     threads: 10
#     conda:  "../wrappers/call_pureClip/env.yaml"
#     script: "../wrappers/call_pureClip/script.py"
    
def call_pureClip_inputs(wcs):
  inputs = {}
  sample=set([i for i in sample_tab.sample_name if i.startswith(wcs.name)])
  inputs['bam'] = expand("mapped/{name}.{{multi}}.{{dups}}.bam", name=sample)
  inputs['bai'] = expand("mapped/{name}.{{multi}}.{{dups}}.bam.bai", name=sample)
  inputs['gen'] = expand("{ref_dir}/seq/{ref}.fasta.gz", ref_dir=reference_directory, ref=config["reference"])[0]
  print(inputs)
  return inputs

rule call_pureClip:
    input:  unpack(call_pureClip_inputs),
    output: bed = "results/pureClip/{name}/{name}.{multi}.{dups}.crosslink_sites.bed",
            bed2= "results/pureClip/{name}/{name}.{multi}.{dups}.binding_regions.bed",
    log:    run = "logs/{name}/{name}.{multi}.{dups}.call_pureClip.log",
    threads: 10
    conda:  "../wrappers/call_pureClip/env.yaml"
    script: "../wrappers/call_pureClip/script.py"
    
rule call_pureClip_with_ctrl:
    input:  bam = "mapped/{name}.{multi}.{dups}.bam",
            bai = "mapped/{name}.{multi}.{dups}.bam.bai",
            cbam = "mapped/{cname}.{multi}.{dups}.bam",
            cbai = "mapped/{cname}.{multi}.{dups}.bam.bai",
            gen = expand("{ref_dir}/seq/{ref}.fasta.gz", ref_dir=reference_directory, ref=config["reference"])[0],
    output: bed = "results/pureClip_w_ctl/{name}_VS_{cname}/{name}_VS_{cname}.{multi}.{dups}.crosslink_sites.bed",
            bed2= "results/pureClip_w_ctl/{name}_VS_{cname}/{name}_VS_{cname}.{multi}.{dups}.binding_regions.bed",
    log:    run = "logs/{name}_VS_{cname}/{name}_VS_{cname}.{multi}.{dups}.call_pureClip_with_ctrl.log",
    threads: 10
    conda:  "../wrappers/call_pureClip_with_ctrl/env.yaml"
    script: "../wrappers/call_pureClip_with_ctrl/script.py"


# ##########################################
# # QC AFTER MAPPING
# #
# 
# rule inspect_bam_coverage:
#     input:  bam = ADIR+"/mapped/{full_name}.bam",
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
#     input:  bam = [ADIR+"/mapped/"+x+".bam" for x in cfg.full_name.tolist()],
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

# def process_bams_input(wildcards):
#     if not "_merged" in wildcards.full_name:
#         return {
#             'bam': "mapped/"+wildcards.full_name+".bam",
#         }
#     else:
#         return {
#             'bam': ["mapped/"+x+".bam" for x in cfg.loc[cfg.condition == re.sub('_merged$', '', wildcards.full_name), "full_name"].tolist()],
#         }
        
rule process_bams:
    input:  bam = "mapped/{sample}.bam",
    output: bam = "mapped/{sample}.{multi}.{dups}.bam",
            bai = "mapped/{sample}.{multi}.{dups}.bam.bai",
    log:    run = "logs/{sample}/{sample}.{multi}.{dups}.process_bams.log"
    threads: 5,
    params: quality_cutof = config["quality_cutof"],
            tmp_bam = "mapped/{sample}.{multi}.{dups}.tmp.bam",
    conda:  "../wrappers/process_bams/env.yaml"
    script: "../wrappers/process_bams/script.py"


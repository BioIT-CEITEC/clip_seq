#########################################
# wrapper for rule: process_bams
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

SAMTOOLS = "samtools"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: process_bams \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

extra_params = ""
if snakemake.wildcards.multi == "uniq_reads":
    extra_params += " -q "+str(snakemake.params.quality_cutof)
if snakemake.wildcards.dups == "no_dups":
    extra_params += " -F 1024"

if isinstance(snakemake.input.bam, list):
    command = SAMTOOLS+" merge -@ "+str(snakemake.threads)+" "+snakemake.params.tmp_bam+" "+" ".join(snakemake.input.bam)+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = SAMTOOLS+" view -@ "+str(snakemake.threads)+" "+extra_params+" -b "+snakemake.params.tmp_bam+" > "+snakemake.output.bam+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    command = "rm "+snakemake.params.tmp_bam+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

else:
    # samtools view -b -@ {threads} -q {params.quality_cutof} {input.bam} > {output.bam} 2> {log.run} && samtools index {output.bam} >> {log.run} 2>&1
    command = SAMTOOLS+" view -@ "+str(snakemake.threads)+" "+extra_params+" -b "+snakemake.input.bam+" > "+snakemake.output.bam+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

command = SAMTOOLS+" index "+snakemake.output.bam+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

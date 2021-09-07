#########################################
# wrapper for rule: call_pureClip
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

SAMTOOLS = "samtools"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_pureClip \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen(SAMTOOLS+" --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

version = str(subprocess.Popen("pureclip --version 2>&1 | grep \"pureclip\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

# ${PURECLIP} -i ${sample} -bai ${sample}.bai -g ${SCRATCH}/${GENOME} -nt ${THREADS} -iv "1;2;3;" -oa -o $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.crosslink_sites.bed -or $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.binding_regions.bed &> $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.binding_regions.log
command = "pureclip -i "+snakemake.input.bam+" -bai "+snakemake.input.bai+" -g "+snakemake.input.gen[0]+" -nt "+str(snakemake.threads)+" -iv '1;2;3;' -o "+snakemake.output.bed+" -or "+snakemake.output.bed2+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

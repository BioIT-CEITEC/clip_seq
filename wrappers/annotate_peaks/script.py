#########################################
# wrapper for rule: annotate_peaks
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

SAMTOOLS = "samtools"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: annotate_peaks \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("R --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
    command = "Rscript "+snakemake.params.rscript+" "+snakemake.input.bed+" "+snakemake.input.gtf+" "+snakemake.output.bed+" "+snakemake.params.feat_type+" "+snakemake.params.annotate_by+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
else:
    with open(snakemake.log.run, 'at') as f:
        f.write("## NOTE: "+snakemake.input.bed+" is empty\n")
  
    command = "touch "+snakemake.output.bed+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)


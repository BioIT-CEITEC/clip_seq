#########################################
# wrapper for rule: post_process_by_RCAS
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

SAMTOOLS = "samtools"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: post_process_by_RCAS \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("R --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

if snakemake.params.organism == "homsap":
    try:
        command = "Rscript "+snakemake.params.rscript+" "+snakemake.input.bed+" "+snakemake.input.gtf+" "+snakemake.params.dir+" "+snakemake.output.tmp_bed+" "+snakemake.input.msigdb+" >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        
        with open(snakemake.log.run, 'at') as f:
            f.write("## ERROR: "+message+"\n")
        
        command = "touch "+snakemake.params.html+" >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
else:
    # other species must be added in RCAS_script.R
    raise ValueError("Only Human is currently supported, other species must be added!")

command = "mv "+snakemake.params.html+" "+snakemake.output.html+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

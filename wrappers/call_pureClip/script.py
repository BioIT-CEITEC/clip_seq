#########################################
# wrapper for rule: call_pureClip
#########################################
import subprocess
import os
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log.run)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: call_pureClip \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

try:
    # ${PURECLIP} -i ${sample} -bai ${sample}.bai -g ${SCRATCH}/${GENOME} -nt ${THREADS} -iv "1;2;3;" -oa -o $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.crosslink_sites.bed -or $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.binding_regions.bed &> $SCRATCH/pureclip/`basename ${sample%%.*}`/`basename ${sample%%.*}`.PureCLIP.binding_regions.log
    command = "$(which time) pureclip -i "+snakemake.input.bam+" -bai "+snakemake.input.bai+" -g "+snakemake.input.gen+" -nt "+str(snakemake.threads)+" -iv '1;2;3;' -o "+snakemake.output.bed+" -or "+snakemake.output.bed2+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
except Exception as ex:
    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
    message = template.format(type(ex).__name__, ex.args)
    print(message)
    
    with open(log_filename, 'at') as f:
        f.write("## ERROR: "+message+"\n")
    
    command = "touch "+snakemake.output.bed+" "+snakemake.output.bed2+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

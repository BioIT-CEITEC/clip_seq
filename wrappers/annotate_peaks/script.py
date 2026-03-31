#########################################
# wrapper for rule: annotate_peaks
#########################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log.run)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: annotate_peaks \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
    command = "$(which time) Rscript "+snakemake.params.rscript+\
              " "+snakemake.input.bed+" "+snakemake.input.gtf+\
              " "+snakemake.output.bed+" "+snakemake.params.feat_type+\
              " "+snakemake.params.annotate_by+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
else:
    with open(log_filename, 'at') as f:
        f.write("## NOTE: "+snakemake.input.bed+" is empty\n")
  
    command = "touch "+snakemake.output.bed+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)


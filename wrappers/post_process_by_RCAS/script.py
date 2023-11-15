#########################################
# wrapper for rule: post_process_by_RCAS
#########################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log.run)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: post_process_by_RCAS \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if snakemake.params.organism == "homo_sapiens":
  if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
    try:
        command = "Rscript "+snakemake.params.rscript+\
                  " "+snakemake.input.bed+\
                  " "+snakemake.input.gtf+\
                  " "+snakemake.params.dir+\
                  " "+snakemake.output.tmp_bed+\
                  " >> "+log_filename+" 2>&1"
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

        command = "touch "+snakemake.params.html+" >> "+log_filename+" 2>&1"
        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
  else:
    with open(log_filename, 'at') as f:
        f.write("## NOTE: "+snakemake.input.bed+" is empty\n")

    command = "touch "+snakemake.params.html+" "+snakemake.output.tmp_bed+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
else:
    # other species must be added in RCAS_script.R
    raise ValueError("Only Human is currently supported, other species must be added!")

command = "mv "+snakemake.params.html+" "+snakemake.output.html+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

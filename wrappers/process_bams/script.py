#########################################
# wrapper for rule: process_bams
#########################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: process_bams \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

extra_params = ""
if snakemake.wildcards.multi == "uniq_reads":
    extra_params += " -q "+str(snakemake.params.quality_cutof)
if snakemake.wildcards.dups == "no_dups":
    extra_params += " -F 1024"

if isinstance(snakemake.input.bam, list):
    command = "samtools merge -@ "+str(snakemake.threads)+" "+snakemake.params.tmp_bam+" "+" ".join(snakemake.input.bam)+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "samtools view -@ "+str(snakemake.threads)+" "+extra_params+" -b "+snakemake.params.tmp_bam+" > "+snakemake.output.bam+" 2>> "+log_filename
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    command = "rm "+snakemake.params.tmp_bam+" >> "+log_filename+" 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

else:
    # samtools view -b -@ {threads} -q {params.quality_cutof} {input.bam} > {output.bam} 2> {log.run} && samtools index {output.bam} >> {log.run} 2>&1
    command = "samtools view -@ "+str(snakemake.threads)+" "+extra_params+" -b "+snakemake.input.bam+" > "+snakemake.output.bam+" 2>> "+log_filename
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

command = "samtools index "+snakemake.output.bam+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

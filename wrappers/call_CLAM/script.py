#########################################
# wrapper for rule: call_CLAM
#########################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

############################
# Be aware that CLAM is installed dirrectly on our machine bcf-server and so far is not on conda
############################

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: call_CLAM \n##\n")
f.close()

shell("which CLAM")

version = str(subprocess.Popen("CLAM --version 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

if snakemake.params.strand == "none":
  strand = "--unstranded"
else:
  strand = ""
# $CLAM permutation_callpeak -i unique.sorted.bam realigned.sorted.bam -o $SCRATCH/clam/${i%.bam*} --gtf $SCRATCH/$GTF -p $THREADS --random-state 777 --qval-cutoff 0.005 --unstranded
command = "CLAM permutation_callpeak -i "+snakemake.input.uniq_bam+" "+snakemake.input.mult_bam+" -o "+snakemake.params.dir+" --gtf "+snakemake.input.gtf+" -p "+str(snakemake.threads)+" --qval-cutoff "+str(snakemake.params.qval_cutoff)+" --merge-size "+str(snakemake.params.merge_size)+" --extend "+str(snakemake.params.extend_peak)+" --random-state 777 "+strand+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

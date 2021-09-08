#########################################
# wrapper for rule: call_CLAM_preprocess
#########################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

############################
# Be aware that CLAM is installed dirrectly on our machine bcf-server and so far is not on conda
############################

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: call_CLAM_preprocess \n##\n")
f.close()

shell("which CLAM")

version = str(subprocess.Popen("CLAM --version 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

# $CLAM realigner -i $i -o $SCRATCH/clam/${i%.bam*} --read-tagger-method median --max-multihits 100 --strandness none --winsize 50
command = "CLAM realigner -i "+snakemake.input.bam+" -o "+snakemake.params.dir+" --read-tagger-method "+snakemake.params.read_tagger+" --max-multihits "+str(snakemake.params.max_multi_hits)+" --strandness "+snakemake.params.strand+" --winsize 50 >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

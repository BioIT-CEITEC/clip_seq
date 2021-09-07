#########################################
# wrapper for rule: call_CLAM
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

############################
# Be aware that CLAM is installed dirrectly on our machine bcf-server and so far is not on conda
############################

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_CLAM \n##\n")
f.close()

shell.executable("/bin/bash")

shell("which CLAM")

version = str(subprocess.Popen("CLAM --version 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

# $CLAM realigner -i $i -o $SCRATCH/clam/${i%.bam*} --read-tagger-method median --max-multihits 100 --strandness none --winsize 50
# command = "CLAM realigner -i "+snakemake.input.bam+" -o "+snakemake.params.dir+" --read-tagger-method "+snakemake.params.read_tagger+" --max-multihits "+str(snakemake.params.max_multi_hits)+" --strandness "+snakemake.params.strand+" --winsize 50 >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

if snakemake.params.strand == "none":
  strand = "--unstranded"
else:
  strand = ""
# $CLAM permutation_callpeak -i unique.sorted.bam realigned.sorted.bam -o $SCRATCH/clam/${i%.bam*} --gtf $SCRATCH/$GTF -p $THREADS --random-state 777 --qval-cutoff 0.005 --unstranded
command = "CLAM permutation_callpeak -i "+snakemake.input.uniq_bam+" "+snakemake.input.mult_bam+" -o "+snakemake.params.dir+" --gtf "+snakemake.input.gtf[0]+" -p "+str(snakemake.threads)+" --qval-cutoff "+str(snakemake.params.qval_cutoff)+" --merge-size "+str(snakemake.params.merge_size)+" --extend "+str(snakemake.params.extend_peak)+" --random-state 777 "+strand+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

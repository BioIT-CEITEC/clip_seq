#########################################
# wrapper for rule: call_macs2
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_macs2 \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("macs2 --version 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

version = str(subprocess.Popen("bedGraphToBigWig 2>&1 | grep \"bedGraphToBigWig v\" | cut -f1 -d- ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

if snakemake.wildcards.dups == "no_dups":
   keep_dups = "1"
else:
   keep_dups = "all"
   
nomodel = "--nolambda --nomodel"
# if str(snakemake.params.frag_len) == "unk":
#   nomodel = ""
# else:
#   nomodel = "--nomodel --extsize "+str(snakemake.params.frag_len)

command = "mkdir -p "+snakemake.params.dir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# macs2 callpeak -t $i -c ${CTL} -g ${GSIZE} -s ${RL} --outdir ${i%.bam} --name ${i%.bam} --nomodel --extsize ${FL} --bdg --tempdir ./tmp -q 0.05 2>&1 | tee ${i%.bam}.macs2.log
command = "macs2 callpeak -t "+snakemake.input.bam+" --keep-dup "+keep_dups+" -g "+str(snakemake.params.effective_GS)+" --outdir "+snakemake.params.dir+" --name "+snakemake.params.name+" "+nomodel+" --bdg --tempdir "+snakemake.params.temp+" -q "+str(snakemake.params.qval_cutof)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.trt_bdg+" "+snakemake.output.trt_bdg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.xls_tab+" "+snakemake.output.xls_tab+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.sum_tab+" "+snakemake.output.sum_tab+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.nar_tab+" "+snakemake.output.nar_tab+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "bedGraphToBigWig "+snakemake.output.trt_bdg+" "+snakemake.input.chrs[0]+" "+snakemake.output.trt_bwg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

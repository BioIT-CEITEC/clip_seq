#########################################
# wrapper for rule: call_macs2
#########################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log.run)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: call_macs2 \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()


keep_dups = "all"
# if snakemake.wildcards.dups == "no_dups":
#    keep_dups = "1"

nomodel = "--nolambda --nomodel"
# if str(snakemake.params.frag_len) == "unk":
#   nomodel = ""
# else:
#   nomodel = "--nomodel --extsize "+str(snakemake.params.frag_len)

command = "mkdir -p "+snakemake.params.dir+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# macs2 callpeak -t $i -c ${CTL} -g ${GSIZE} -s ${RL} --outdir ${i%.bam} --name ${i%.bam} --nomodel --extsize ${FL} --bdg --tempdir ./tmp -q 0.05 2>&1 | tee ${i%.bam}.macs2.log
command = "$(which time) macs2 callpeak -t "+snakemake.input.bam+\
          " --keep-dup "+keep_dups+" -g "+str(snakemake.params.effective_GS)+\
          " --outdir "+snakemake.params.dir+" --name "+snakemake.params.name+\
          " "+nomodel+" --bdg --tempdir "+snakemake.params.temp+\
          " -q "+str(snakemake.params.qval_cutof)+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.trt_bdg+" "+snakemake.output.trt_bdg+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.xls_tab+" "+snakemake.output.xls_tab+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.sum_tab+" "+snakemake.output.sum_tab+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.nar_tab+" "+snakemake.output.nar_tab+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) bedGraphToBigWig "+snakemake.output.trt_bdg+" "+snakemake.input.chrs+" "+snakemake.output.trt_bwg+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

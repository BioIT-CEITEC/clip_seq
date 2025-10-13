library(data.table)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(scales)

# setwd("/mnt/ssd/ssd_1/snakemake/stage557_Anil20220203.ChIP_peak_calling/ChIP-seq_general_analysis/")
# setwd("/mnt/ssd/ssd_1/sequia/5232__chipseq_analysis__peak_calling_test__220810")
# homer_file="macs2/PCF11_ChIP/PCF11_ChIP.no_dups.peaks.all.annotated_by_HOMER.tsv"
# homer_file="macs2/PHRF1_ChIP/PHRF1_ChIP.no_dups.peaks.all.annotated_by_HOMER.tsv"
# homer_file="results/MSPC_peaks/p53wt_360/p53wt_360.no_dups.peaks.all.annotated.new.tsv"
# fdr_cutof = 0.05
# best = 10
#args = c("results/MSPC_peaks/Mut_H33_Fl/Mut_H33_Fl.no_dups.peaks.all.annotated.tsv", "0.05", "10")

args = commandArgs(trailingOnly=TRUE)
print(args)
homer_file = args[1]
fdr_cutof = as.numeric(args[2])
best = as.numeric(args[3])

prefix = sub("\\.tsv$","",homer_file)
homer = fread(homer_file, sep = '\t')
homer[,ann_type:=gsub(" \\(.*","",Annotation)]
if(homer[is.na(Distance_to_TSS),.N]>0) {
  homer[is.na(Distance_to_TSS),ann_type:="unannotated"]
  order_ann_types = c("promoter-TSS", "exon", "intron", "TTS", "Intergenic", "unannotated")
} else {
  order_ann_types = c("promoter-TSS", "exon", "intron", "TTS", "Intergenic")
}
homer[,FDR:=10^-qvalue]
max_chr_plot = 25
homer[,new_chr:=chr]
if(length(unique(homer$chr))>max_chr_plot) {
  order_chr_types = c(homer[,.N,by=chr][order(-N)][1:(max_chr_plot-1),chr],"other")
  homer[!chr %in% order_chr_types, new_chr:="other"]
} else {
  order_chr_types = homer[,.N,by=chr][order(-N),chr]
}

if(Inf %in% range(homer$score)) {
  y_scale = breaks_extended(n = 9)(homer[score != Inf, range(score)])
  y_scale_labels = sapply(y_scale, function(x){paste0(x,"\n(",signif(10^(-0.1*x),2),")")})
  y_scale_step = y_scale[length(y_scale)]-y_scale[length(y_scale)-1]
  y_scale = c(y_scale, Inf)
  y_scale_labels = c(y_scale_labels, "Inf\n(0)")
} else {
  y_scale = breaks_extended(n = 10)(range(homer$score))
  y_scale_labels = sapply(y_scale, function(x){paste0(x,"\n(",signif(10^(-0.1*x),2),")")})
}

# annstat_basic=fread(file="macs2/PCF11_ChIP/PCF11_ChIP.no_dups.peaks.all.ann_stats.tsv",nrows = 5)
# annstat_ext=fread(cmd="tail -n +7 macs2/PCF11_ChIP/PCF11_ChIP.no_dups.peaks.all.ann_stats.tsv")

# plot HOMER peak score and FDR against the distance from nearest TSS (line not points)
nTSS_plot = ggplot(homer[!is.na(Distance_to_TSS)], aes(x=Distance_to_TSS, y=score, color=FDR<fdr_cutof))+
  geom_segment(aes(x=Distance_to_TSS-((end-start)*0.5), y=score, xend=Distance_to_TSS+((end-start)*0.5), yend=score), alpha=.33)+
  scale_color_manual(values = c("red", "black"))+
  theme(legend.position = "top")+
  # geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed")+
  scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
  geom_text_repel(aes(Distance_to_TSS, score, label = Gene_Name), data = homer[!is.na(Distance_to_TSS)][order(FDR)][1:min(best,.N)], max.overlaps = best)+
  coord_cartesian(xlim=c(-1e4,1e4))+
  labs(title="Plot of peaks distance to nearest TSS",
       color=paste0("FDR<",fdr_cutof), 
       x="distance to TSS: upstream (-), downstream (+)", 
       y="peak score by MACS2 (peak FDR)")

pdf(paste0(prefix,".dist_to_TSS.pdf"), width = 12, height = 9)
print(nTSS_plot)
dev.off()

# plot HOMER peak score and FDR against the distance from nearest TSS (line not points) and zoomed to smaller distance
nTSS_plot_2 = ggplot(homer[!is.na(Distance_to_TSS)], aes(x=Distance_to_TSS, y=score, color=FDR<fdr_cutof))+
  geom_segment(aes(x=Distance_to_TSS-((end-start)*0.5), y=score, xend=Distance_to_TSS+((end-start)*0.5), yend=score), alpha=.33)+
  scale_color_manual(values = c("red", "black"))+
  theme(legend.position = "top")+
  # geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed")+
  scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
  geom_text_repel(aes(Distance_to_TSS, score, label = Gene_Name), data = homer[!is.na(Distance_to_TSS)][order(FDR)][1:min(best,.N)], max.overlaps = best)+
  coord_cartesian(xlim=c(-1e3,1e3))+
  labs(title="Plot of peaks distance to nearest TSS",
       color=paste0("FDR<",fdr_cutof), 
       x="distance to TSS: upstream (-), downstream (+)", 
       y="peak score by MACS2 (peak FDR)")
pdf(paste0(prefix,".dist_to_TSS_2.pdf"), width = 12, height = 9)
print(nTSS_plot_2)
dev.off()

# # plot FDR against the distance from nearest TSS (line not points)
# ggplot(homer[!is.na(Distance_to_TSS)], aes(x=Distance_to_TSS, y=FDR, color=FDR<fdr_cutof))+
#   geom_segment(aes(x=Distance_to_TSS-((end-start)*0.5), y=FDR, xend=Distance_to_TSS+((end-start)*0.5), yend=FDR), alpha=.33)+
#   scale_color_manual(values = c("red", "black"))+
#   theme(legend.position = "top")+
#   # geom_hline(yintercept = fdr_cutof, linetype="dashed")+
#   scale_y_continuous(trans = trans_reverser('log10'))+
#   geom_text_repel(aes(Distance_to_TSS, FDR, label = Gene_Name), data = homer[!is.na(Distance_to_TSS)][order(FDR)][1:best], max.overlaps = best)+
#   coord_cartesian(xlim=c(-1e4,1e4))+
#   labs(title="Plot of peaks distance to nearest TSS", 
#        color=paste0("FDR<",fdr_cutof),
#        x="Distance to TSS: upstream (-), downstream (+)", 
#        y="peaks FDR (in log10 scale)")

# plot HOMER score and FDR of all peaks against the annotation type (TSS, exon, intron, etc)
annot_types_plot = ggplot(homer, aes(x=ann_type, y=score))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  # geom_jitter(aes(color=FDR<fdr_cutof), width = .4, alpha = .25)+
  # scale_color_manual(values = c("red", "black"))+
  geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed", colour = "red")+
  scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
  theme(legend.position = "top")+
  scale_x_discrete(
    limits=order_ann_types,
    breaks=order_ann_types,
    # labels=c("TSS", "Exon", "Intron", "TTS", "Intergenic"))+
    # labels=c(paste0("TSS\n#",homer[ann_type=="promoter-TSS",.N],"\n(",homer[, round(sum(ann_type=="promoter-TSS")/.N*100,2)],"%)"), 
    #          paste0("Exon\n#",homer[ann_type=="exon",.N],"\n(",homer[, round(sum(ann_type=="exon")/.N*100,2)],"%)"), 
    #          paste0("Intron\n#",homer[ann_type=="intron",.N],"\n(",homer[, round(sum(ann_type=="intron")/.N*100,2)],"%)"), 
    #          paste0("TTS\n#",homer[ann_type=="TTS",.N],"\n(",homer[, round(sum(ann_type=="TTS")/.N*100,2)],"%)"), 
    #          paste0("Intergenic\n#",homer[ann_type=="Intergenic",.N],"\n(",homer[, round(sum(ann_type=="Intergenic")/.N*100,2)],"%)")))+
    labels=homer[, paste0(.BY,"\n#",.N,"\n(",round(.N/homer[,.N]*100,2),"%)"), keyby=.(ann_type)][order_ann_types, V1])+
  labs(title="Plot of (all) peaks distribution according to annotation types",
       subtitle=paste0("filter on FDR < ",fdr_cutof), 
       x ="genomic annotation type", 
       y ="peak score by MACS2 (peak FDR)")
pdf(paste0(prefix,".annot_types.pdf"), width = 12, height = 9)
print(annot_types_plot)
dev.off()

# plot HOMER score and FDR of filtered peaks against the annotation type (TSS, exon, intron, etc)
if(homer[FDR<fdr_cutof,.N] > 0) {
  annot_types_plot_2 = ggplot(homer[FDR<fdr_cutof], aes(x=ann_type, y=score))+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
    # geom_jitter(aes(color=FDR<fdr_cutof), width = .4, alpha = .25)+
    # scale_color_manual(values = c("red", "black"))+
    geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed", colour = "red")+
    scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
    theme(legend.position = "top")+
    scale_x_discrete(
      limits=order_ann_types,
      breaks=order_ann_types,
      labels=homer[FDR<fdr_cutof, paste0(.BY,"\n#",.N,"\n(",round(.N/homer[FDR<fdr_cutof,.N]*100,2),"%)"), keyby=.(ann_type)][order_ann_types, V1])+
    labs(title="Plot of (filtered) peaks distribution according to annotation types", 
         subtitle=paste0("filter on FDR < ",fdr_cutof),
         x ="genomic annotation type", 
         y ="peak score by MACS2 (peak FDR)")
} else {
  annot_types_plot_2 = ggplot() +
    theme_void() +
    geom_text(aes(0,0,label=paste0('Not enough data after using filter FDR<',fdr_cutof,'!'))) +
    xlab(NULL)
}
pdf(paste0(prefix,".annot_types_filt.pdf"), width = 12, height = 9)
print(annot_types_plot_2)
dev.off()

# plot HOMER score and FDR of all peaks against the chromosomes/contigs
chr_plot = ggplot(homer, aes(x=new_chr, y=score, color=FDR<fdr_cutof))+
  geom_point(data=homer[,.(.N,score,FDR),by=.(new_chr)][N<2,.(new_chr,score,FDR)])+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed", colour = "red")+
  scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
  scale_color_manual(values = c("red", "black"))+
  theme(legend.position = "top")+
  scale_x_discrete(
    limits=order_chr_types,
    breaks=order_chr_types,
    labels=homer[, paste0(.BY,"\n#",.N,"\n(",round(.N/homer[,.N]*100,2),"%)"), keyby=.(new_chr)][order_chr_types, V1])+
  labs(title="Plot of (all) peaks distribution according to chromosomes/contigs",
       subtitle=paste0("filter on FDR < ",fdr_cutof), 
       x ="chromosome/contig", 
       y ="peak score by MACS2 (peak FDR)")
pdf(paste0(prefix,".chr_types.pdf"), width = 12, height = 9)
print(chr_plot)
dev.off()

order_chr_types = order_chr_types[which(order_chr_types %in% homer[FDR<fdr_cutof,unique(new_chr)])]
# plot HOMER score and FDR of filtered peaks against the chromosomes/contigs
if(homer[FDR<fdr_cutof,.N] > 0) {
  chr_plot_2 = ggplot(homer[FDR<fdr_cutof], aes(x=new_chr, y=score))+
    geom_point(data=homer[FDR<fdr_cutof,.(.N,score),by=.(new_chr)][N<2,.(new_chr,score)])+
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
    # geom_jitter(aes(color=FDR<fdr_cutof), width = .4, alpha = .25)+
    # scale_color_manual(values = c("red", "black"))+
    geom_hline(yintercept = -10*log10(fdr_cutof), linetype="dashed", colour = "red")+
    scale_y_continuous(breaks=y_scale, labels=y_scale_labels)+
    theme(legend.position = "top")+
    scale_x_discrete(
      limits=order_chr_types,
      breaks=order_chr_types,
      labels=homer[FDR<fdr_cutof, paste0(.BY,"\n#",.N,"\n(",round(.N/homer[FDR<fdr_cutof,.N]*100,2),"%)"), keyby=.(new_chr)][order_chr_types, V1])+
    labs(title="Plot of (filtered) peaks distribution according to chromosomes/contigs",
         subtitle=paste0("filter on FDR < ",fdr_cutof), 
         x ="chromosome/contig", 
         y ="peak score by MACS2 (peak FDR)")
} else {
  chr_plot_2 = ggplot() +
    theme_void() +
    geom_text(aes(0,0,label=paste0('Not enough data after using filter FDR<',fdr_cutof,'!'))) +
    xlab(NULL)
}
pdf(paste0(prefix,".chr_types_filt.pdf"), width = 12, height = 9)
print(chr_plot_2)
dev.off()

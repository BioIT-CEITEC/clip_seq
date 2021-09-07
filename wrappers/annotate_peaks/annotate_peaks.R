library(data.table)
library(rtracklayer)
library(stringr)

args <- commandArgs(trailingOnly = T)
bed_file <- args[1]
gtf_file <- args[2]
out_bed  <- args[3]
feat_type <- args[4]
identifier<- str_split(args[5], ",")[[1]]  # must be comma separated list of valid GTF TAGS (e.g., gene_name, gene_id)

comp_tab <- fread(bed_file, sep = "\t", select=1:3)
names(comp_tab) <- c("chr","start","end")

ref <- as.data.table(rtracklayer::import(gtf_file))[type == feat_type, c("seqnames","start","end",identifier), with=F]

overlapped <- comp_tab[ref,on=c("chr==seqnames"),allow.cartesian=T][start+1<=i.end & end>=i.start]
res <- overlapped[, paste0(unique(get(identifier[1])),collapse = ","), by=.(chr,start,end)]
setnames(res, "V1", identifier[1])
for(idf in identifier[-1]) {
  res <- merge(res,overlapped[, paste0(unique(get(idf)),collapse = ","), by=.(chr,start,end)], by=c("chr","start","end"))
  setnames(res, "V1", idf)
}
res <- merge(res, comp_tab, by = c("chr","start","end"), all.y = T)
res[is.na(res)] <- ""

fwrite(res, out_bed, sep="\t", row.names = F, col.names = F, quote = F)

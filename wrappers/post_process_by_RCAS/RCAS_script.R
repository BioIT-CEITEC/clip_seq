# source('http://bioconductor.org/biocLite.R')
# biocLite('RCAS')
library(RCAS)
library(data.table)

# For other species (supported are mouse, fly and worm) a genome file must be downloaded as a library (see line below) and genomeVersion parameter must be changed in runReport
if (!require('BSgenome.Hsapiens.UCSC.hg38')) BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')

args = commandArgs(trailingOnly=TRUE)
# 
# bed = "results/CLAM/ZNF-FCLx/ZNF-FCLx.all_reads.keep_dups/narrow_peak.permutation.processed.bed"
# gtf = "/mnt/nfs/shared/CFBioinformatics/references_backup/homo_sapiens/GRCh38-p10/annot/GRCh38-p10.gtf"
# dir = "results/RCAS/from_CLAM/ZNF-FCLx/"
# tmp_bed = "results/RCAS/from_CLAM/ZNF-FCLx/ZNF-FCLx.all_reads.keep_dups.input.bed"
# msigdb = "/mnt/nfs/shared/CFBioinformatics/references_backup/homo_sapiens/GRCh38-p10/other/MSigDB_for_RCAS/c2.all.v7.1.entrez.gmt"

bed <- args[1]
gtf <- args[2]
dir <- args[3]
tmp_bed <- args[4]

bed_tab <- fread(bed, sep = "\t", select=1:6)
names(bed_tab) <- c("chr","start","end","name","score","strand")
class(bed_tab[["chr"]]) <- "character"
fwrite(bed_tab, tmp_bed, sep = "\t", quote = F, row.names = F, col.names = F)

setwd(dir)

runReport( queryFilePath = basename(tmp_bed),
           gffFilePath = gtf,
           printProcessedTables = T,
           genomeVersion = "hg38",
           motifAnalysis = T)


# queryRegions <- importBed(filePath = tmp_bed, sampleN = 10000)
# # setwd("/mnt/nfs/shared/999993-Bioda/projects/a96_dragana_clipseq/results/pub1/peak_calling/pureclip/oa/SRR1688578/")
# # bed <- fread("SRR1688578.PureCLIP.crosslink_sites.bed",sep = "\t")
# bed$V7 <- NULL
# fwrite(bed,"SRR1688578.PureCLIP.crosslink_sites-V7.bed",sep = "\t",quote = F,row.names = F,col.names = F)
# queryRegions <- importBed(filePath = "SRR1688578.PureCLIP.crosslink_sites-V7.bed", sampleN = 10000)
# gff <- importGtf(filePath = "/mnt/nfs/home/408320/000000-My_Documents/VM-home/annotation/ensembl91/Homo_sapiens.GRCh38.91.gtf.gz")
# 
# #Querying the annotation file
# overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))
# 
# #To find out the distribution of the query regions across gene types:
# biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
# df <- overlaps[,length(unique(overlappingQuery)), by = biotype_col]
# colnames(df) <- c("feature", "count")
# df$percent <- round(df$count / length(queryRegions) * 100, 1)
# df <- df[order(count, decreasing = TRUE)]
# 
# ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
#   geom_bar(stat = 'identity', aes(fill = feature)) + 
#   geom_label(aes(y = percent + 0.5), label = df$count) + 
#   labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
#   theme_bw(base_size = 14) + 
#   theme(axis.text.x = element_text(angle = 90))
# 
# #GTF files contain some annotation features (e.g. exons, transcripts) that are usually explicitly defined, however, some transcript features 
# #such as introns, exon-intron boundaries, promoter regions are only implicitly defined. Such implicit features can be extracted from a GTF file 
# #using makeTxDb family of functions from the GenomicFeatures library.
# 
# #First we create a list of GRanges objects, where each list element contains all the available coordinates of transcript features such 
# #as transcripts, exons, introns, 5’/3’ UTRs, exon-intron boundaries, and promoter regions.
# 
# txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
# 
# 
# #To have a global overview of the distribution of query regions across gene features, we can use the summarizeQueryRegions function. 
# #If a given query region does not overlap with any of the given coordinates of the transcript features, it is categorized under NoFeatures.
# 
# summary <- summarizeQueryRegions(queryRegions = queryRegions, 
#                                  txdbFeatures = txdbFeatures)
# 
# df <- data.frame(summary)
# df$percent <- round((df$count / length(queryRegions)), 3) * 100
# df$feature <- rownames(df)
# ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
#   geom_bar(stat = 'identity', aes(fill = feature)) + 
#   geom_label(aes(y = percent + 3), label = df$count) + 
#   labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
#   theme_bw(base_size = 14) + 
#   theme(axis.text.x = element_text(angle = 90))
# 
# 
# #Obtaining a table of overlap counts between query regions and genes
# #To find out which genes overlap with how many queries and categorise overlaps by transcript features;
# #we use getTargetedGenesTable function, which returns a data.frame object.
# 
# dt <- getTargetedGenesTable(queryRegions = queryRegions, 
#                             txdbFeatures = txdbFeatures)
# dt <- dt[order(transcripts, decreasing = TRUE)]
# 
# knitr::kable(dt[1:10,])
# 
# #Profiling the coverage of query regions across transcript features
# #Coverage profile of query regions at feature boundaries
# #It may be useful to look at the distribution of query regions at the boundaries of transcript features. For instance, it may be important 
# #to see the relative signal at transcript ends (transcription start sites versus transcription end sites). Or, it may be important to see 
# #how the signal is distributed at exon boundaries, which may give an idea about the regulation of the transcript. Here we demonstrate 
# #how to get such signal distributions at transcription start/end sites. The same approach can be done for any other collection of transcript 
# #features (exons, introns, promoters, UTRs etc.)
# #this analysis shows pretty much nothing
# cvgF <- getFeatureBoundaryCoverage(queryRegions = queryRegions, 
#                                    featureCoords = txdbFeatures$transcripts, 
#                                    flankSize = 1000, 
#                                    boundaryType = 'fiveprime', 
#                                    sampleN = 1000000)
# cvgT <- getFeatureBoundaryCoverage(queryRegions = queryRegions, 
#                                    featureCoords = txdbFeatures$transcripts, 
#                                    flankSize = 1000, 
#                                    boundaryType = 'threeprime', 
#                                    sampleN = 1000000)
# 
# cvgF$boundary <- 'fiveprime'
# cvgT$boundary <- 'threeprime'
# 
# df <- rbind(cvgF, cvgT)
# 
# ggplot2::ggplot(df, aes(x = bases, y = meanCoverage)) + 
#   geom_ribbon(fill = 'lightgreen', 
#               aes(ymin = meanCoverage - standardError * 1.96, 
#                   ymax = meanCoverage + standardError * 1.96)) + 
#   geom_line(color = 'black') + 
#   facet_grid( ~ boundary) + theme_bw(base_size = 14) 
# 
# 
# 
# #Coverage profile of query regions for all transcript features
# #Coverage profiles can be obtained for a single type of transcript feature or a list of transcript features. Here we demonstrate how to get 
# #coverage profile of query regions across all available transcript features. It might be a good idea to use sampleN parameter to randomly 
# #downsample the target regions to speed up the calculations.
# 
# cvgList <- calculateCoverageProfileList(queryRegions = queryRegions, 
#                                         targetRegionsList = txdbFeatures,
#                                         sampleN = 100000)
# 
# ggplot2::ggplot(cvgList, aes(x = bins, y = meanCoverage)) + 
#   geom_ribbon(fill = 'lightgreen', 
#               aes(ymin = meanCoverage - standardError * 1.96, 
#                   ymax = meanCoverage + standardError * 1.96)) + 
#   geom_line(color = 'black') + theme_bw(base_size = 14) +
#   facet_wrap( ~ feature, ncol = 3)
# 
# 
# 
# 
# #GO term analysis
# #Biological processes enriched among targeted genes
# #RCAS can perform GO term enrichment analysis to find out enriched functions in genes that overlap the query regions. Below is demonstrated 
# #how to get biological processes terms (‘BP’) enriched in the genes that overlap query regions and the top 10 GO terms with most fold change 
# #increase relative to the background are provided.
# 
# #get all genes from the GTF data
# backgroundGenes <- unique(gff$gene_id)
# #get genes that overlap query regions
# targetedGenes <- unique(overlaps$gene_id)
# 
# #run TopGO
# goBP <- runTopGO(ontology = 'BP', 
#                  species = 'human', 
#                  backgroundGenes = backgroundGenes, 
#                  targetedGenes = targetedGenes)
# 
# goBP <- goBP[order(goBP$foldEnrichment, decreasing = TRUE),]
# rownames(goBP) <- goBP$GO.ID
# goBP <- subset(goBP, select = -c(Annotated,classicFisher, bh, GO.ID))
# 
# knitr::kable(goBP[1:10,])
# 
# 
# 
# 
# #MSIGDB gene sets enriched among targeted genes
# #RCAS can use gene sets from Molecular Signatures Database and calculate gene set enrichment analysis (GSEA) to find out which gene sets 
# #are enriched among the genes targeted by the query regions.
# 
# geneSets <- parseMsigdb('/mnt/nfs/home/408320/000000-My_Documents/VM-home/annotation/human_msigdb_c2.all.v5.1.entrez.gmt')
# 
# resultsGSEA <- runGSEA(geneSetList = geneSets,
#                        backgroundGenes = backgroundGenes, 
#                        targetedGenes = targetedGenes)
# 
# knitr::kable(x = resultsGSEA[1:10,])

##########################################################################
##########################################################################
# Project: Emilio's mir51 project
# Script purpose: analyze the time series quant-seq data 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct 22 15:56:58 2019
##########################################################################
##########################################################################
library("openxlsx")
require('DESeq2')
RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"


### data verision and analysis version
version.Data = 'Quantseq_R8521'
version.analysis = paste0("_", version.Data, "_20191022")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE

### Directories to save results
#design.file = "../exp_design/NGS_Samples_Philipp_mRNAseq_all.xlsx"
#dataDir = "../data/"

resDir = paste0("../results/", version.Data, "/")
tabDir =  paste0(resDir, "tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


########################################################
########################################################
# Section : Load the processed data from Philipp's analysis
# because Emilio's R8521 were sequenced together with Philipp's data
########################################################
########################################################
Counts.to.Use = "UMI"
#QC.for.cpm = FALSE
#EDA.with.normalized.table = FALSE
#Add.miRNA.Targets = TRUE
### specify parameters for DESEeq2 and pairwise comparisons
require(DESeq2)
source(RNA_QCfunctions)
source(RNAfunctions)

load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI_All_incl_Paula_Emilio_Quantseq_R8043_R8521_20190926.Rdata'))

outDir = paste0(resDir, "time_sereis_N2.vs.mutant")
if(!dir.exists(outDir)) dir.create(outDir)


## select Emilio's data
## manually modify design matrix
## extract umi count table
kk = which(design$SampleID >= 100003)
design = design[kk, ]
design$condition[which(design$strain=="MLC1800")] = 'mir51.mutant'
design = design[, -which(colnames(design) == 'strain')]
design$stage = sapply(design$stage, function(x) gsub(' ', '.', x))

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
}else{
  if(Counts.to.Use == "UMI"){
    all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene


##########################################
# start the QC and DESeq2 analysis
##########################################
samples.sels = c(1:nrow(design))
lowlyExpressed.readCount.threshold = 10

design$condTime = paste0(design$stage, "_", design$condition)

pdfname = paste0(outDir, "/Data_QC_wt_mutant", version.analysis, "_", Counts.to.Use, ".pdf")
pdf(pdfname, width = 12, height = 10)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, c(1, 3, 2)])

#dev.off()

##  start DE analysis
dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ stage + condition )
dds$group <- factor(paste0(dds$stage, "_", dds$condition))
design(dds) <- ~ group
#dds$condition = relevel(dds$condition, "wt")

dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
dds <- estimateSizeFactors(dds)

cpm = fpm(dds, robust = TRUE)

colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")

dds = estimateDispersions(dds, fitType = "parametric")

plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

#res = cpm
res.ii = results(dds, contrast=c("group", '190.cell_mir51.mutant', '190.cell_wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_mutant.vs.wt_190.cell")
res = data.frame(res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("group", 'Bean_mir51.mutant', 'Bean_wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_mutant.vs.wt_Bean")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("group", '2.3.fold_mir51.mutant', '2.3.fold_wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_mutant.vs.wt_2.3.fold")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res = data.frame(res, res.ii[, c(2, 5, 6)])


xx = data.frame(cpm, res)

write.csv(xx,
          file = paste0(outDir, "/GeneAll_wt_mutant_", Counts.to.Use, "_", version.analysis, ".csv"), 
          row.names = TRUE)


dev.off()


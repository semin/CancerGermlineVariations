rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(2)

require(sqldf)
require(gdata)
require(Hmisc)
require(gtools)
require(annotate)
require(snpStats)
require(data.table)
require(MatrixEQTL)
require(org.Hs.eg.db)
require(GenomicRanges)
require(VariantAnnotation)


##
## Initialize global variables
##
if (!is.na(pmatch("darwin", version$os))) {
    rootDir = "/Volumes/orchestra"
} else {
    rootDir = "/home/sl279"
}

chrs = c(1:22, "X")
chrs = factor(chrs, level=chrs)
cchrs = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs = factor(cchrs, level=cchrs)

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
gdacDir = file.path(germDir, "GDAC")
tcgaDir = file.path(germDir, "TCGA")
eqtlDir = file.path(germDir, "EQTL")
tableDir = file.path(germDir, "Table")
figDir = file.path(germDir, "Figure")
rdataDir = file.path(germDir, "RData")

gdacStdDataDate = "2015_02_04"
gdacStdDataDateNoUb = gsub("_", "", gdacStdDataDate)
gdacStdDataDir = file.path(gdacDir, sprintf("stddata__%s", gdacStdDataDate))

gdacAnalysesDate = "2014_10_17"
gdacAnalysesDateNoUb = gsub("_", "", gdacAnalysesDate)
gdacAnalysesDir = file.path(gdacDir, sprintf("analyses__%s", gdacAnalysesDate))

baseFontSize = 12
args = commandArgs(TRUE)
canType = args[1]

tcgaCanDir = file.path(tcgaDir, canType)
dir.create(tcgaCanDir, showWarnings = F, recursive = T)

imageFile = file.path(rdataDir, sprintf("GCC-Germline-Postprocess_Broad_GDAC_Gene_Expression_Data-%s.RData", canType))
#load(imageFile)
#save.image(imageFile)


##
## Load RefSeq gene info table (ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz)
##
geneInfoFile = file.path(germDir, "Entrez/gene_info_human.txt")
geneInfoDt = fread(geneInfoFile, header = F)
geneInfoDt = geneInfoDt[, c(2:5), with = F]
setnames(geneInfoDt, c("GeneID", "Symbol", "LocusTag", "Synonyms"))
geneInfoExtDt = geneInfoDt[, list(Synonym = unlist(strsplit(Synonyms, "|", fixed = T), use.names = FALSE)), by = "GeneID,Symbol,LocusTag"]
geneInfoExtDf = as.data.frame(geneInfoExtDt)


##
## Read TCGA Broad GDAC data
##
ncanType = ifelse(canType == "CRC", "COADREAD", canType)
gdacStdDataCanDir = file.path(gdacStdDataDir, ncanType, gdacStdDataDateNoUb)
gdacAnalysesCanDir = file.path(gdacAnalysesDir, ncanType, gdacAnalysesDateNoUb)
canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  


##
## Read gene expression data
##

## GA data
if (ncanType == "STAD") {
    geneExpGaFile = file.path(gdacStdDataCanDir,
                              sprintf("gdac.broadinstitute.org_%s.Merge_rnaseq__illuminaga_rnaseq__bcgsc_ca__Level_3__gene_expression__data.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                              sprintf("%s.rnaseq__illuminaga_rnaseq__bcgsc_ca__Level_3__gene_expression__data.data.txt", ncanType))
} else {
    geneExpGaFile = file.path(gdacStdDataCanDir,
                              sprintf("gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                              sprintf("%s.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", ncanType))
}

if (file.exists(geneExpGaFile)) {
    geneExpGaLines = readLines(geneExpGaFile)
    geneExpGaLines = geneExpGaLines[-2]
    geneExpGaDf = read.delim(textConnection(geneExpGaLines), header = TRUE, as.is = TRUE)
    colnames(geneExpGaDf) = c(colnames(geneExpGaDf)[1], as.vector(sapply(colnames(geneExpGaDf[,2:ncol(geneExpGaDf)]), function(x) { paste(strsplit(x, "\\.")[[1]][1:4], collapse="_") })))
    gaSampleIds = grep("TCGA", colnames(geneExpGaDf), value = TRUE)
} else {
    geneExpGaDf = data.frame()
    gaSampleIds = c()
}

## HiSeq data
if (ncanType == "STAD") {
    geneExpHiseqFile = file.path(gdacStdDataCanDir,
                                 sprintf("gdac.broadinstitute.org_%s.Merge_rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                                 sprintf("%s.rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.data.txt", ncanType))
} else {
    geneExpHiseqFile = file.path(gdacStdDataCanDir,
                                 sprintf("gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                                 sprintf("%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", ncanType))
}

if (file.exists(geneExpHiseqFile)) {
    geneExpHiseqLines = readLines(geneExpHiseqFile)
    geneExpHiseqLines = geneExpHiseqLines[-2]
    geneExpHiseqDf = read.delim(textConnection(geneExpHiseqLines), header = TRUE, as.is = TRUE)
    if (ncanType == "STAD") { geneExpHiseqDf = geneExpHiseqDf[, c(1, seq(2, ncol(geneExpHiseqDf), by = 3))] }
    colnames(geneExpHiseqDf) = c(colnames(geneExpHiseqDf)[1], as.vector(sapply(colnames(geneExpHiseqDf[,2:ncol(geneExpHiseqDf)]), function(x) { paste(strsplit(x, "\\.")[[1]][1:4], collapse="_") })))
    hiseqSampleIds = grep("TCGA", colnames(geneExpHiseqDf), value = TRUE)
} else {
    stop("Cannot find ", geneExpHiseqFile)
}

## Merge GA and HiSeq data
sharedGeneExpSampleIds = intersect(gaSampleIds, hiseqSampleIds)
if (nrow(geneExpGaDf) > 0) {
    geneExpMergedDf = cbind(geneExpGaDf[,-c(1, which(colnames(geneExpGaDf) %in% sharedGeneExpSampleIds))], geneExpHiseqDf[,-1])
    rownames(geneExpMergedDf) = geneExpGaDf[,1]
} else {
    geneExpMergedDf = geneExpHiseqDf[,-1]
    rownames(geneExpMergedDf) = geneExpHiseqDf[,1]
}

## primary tumors only in expression data!
if (ncanType == "SKCM") {
    geneExpPrimDf = geneExpMergedDf[, which(as.integer(gsub("TCGA_\\S{2}_\\S{4}_(\\S{2})\\S{1}", "\\1", colnames(geneExpMergedDf))) == 6)]
} else {
    geneExpPrimDf = geneExpMergedDf[, which(as.integer(gsub("TCGA_\\S{2}_\\S{4}_(\\S{2})\\S{1}", "\\1", colnames(geneExpMergedDf))) == 1)]
}
colnames(geneExpPrimDf) = gsub(".*(TCGA_\\S{2}_\\S{4}).*", "\\1", colnames(geneExpPrimDf), perl = T)

## Update gene names using Entrez gene ID
geneExpPrimModDf = data.frame()
for (i in 1:nrow(geneExpPrimDf)) {
    rname = rownames(geneExpPrimDf)[i]
    rSymbol = strsplit(rname, "|", fixed = T)[[1]][1]
    rEntrez = strsplit(rname, "|", fixed = T)[[1]][2]
    if (ncanType == "STAD") {
        nSymbol = rSymbol
    } else {
        rEntrez2Symbol = getSYMBOL(rEntrez, data='org.Hs.eg')
        nSymbol = ifelse(is.na(rEntrez2Symbol), rSymbol, rEntrez2Symbol)
    }
    if (nSymbol != "?") { 
        cat(sprintf("%s|%s\n", rname, nSymbol))
        #geneExpPrimModDf = rbind(geneExpPrimModDf, data.frame(Gene = rname, geneExpPrimDf[i,], row.names = nSymbol, stringsAsFactors = F))
        geneExpPrimModDf = rbind(geneExpPrimModDf, data.frame(Gene = nSymbol, geneExpPrimDf[i,], row.names = nSymbol, stringsAsFactors = F))
    }
}

if (ncanType == "STAD") {
    geneNameFreqDf = data.frame(symbol = geneExpPrimModDf$Gene)
    geneNameFreqGrpDf = sqldf('SELECT symbol, COUNT(*) AS cnt FROM geneNameFreqDf GROUP BY symbol')
    geneNamesRddDf = subset(geneNameFreqGrpDf, cnt > 1)
    geneNamesRdd = as.character(unique(geneNamesRddDf$symbol))
    geneExpPrimModDf = subset(geneExpPrimModDf, Gene %nin% geneNamesRdd)
}

geneExpPrimFile = file.path(tcgaCanDir, sprintf("Gene_Expression-Tumor-%s.new.txt", canType))
write.table(geneExpPrimModDf, geneExpPrimFile, sep = "\t", na = "NA", row.names = T, col.names = T, quote = F)
save.image(imageFile)


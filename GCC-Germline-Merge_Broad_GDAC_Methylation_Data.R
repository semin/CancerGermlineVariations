rm(list=ls())

##
## Load libraries
##
require(sqldf)
require(gdata)
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

chrs = c(1:22)
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

args = commandArgs(TRUE)
canType = args[1]

tcgaCanDir = file.path(tcgaDir, canType)
dir.create(tcgaCanDir, showWarnings = F, recursive = T)
imageFile = file.path(rdataDir, sprintf("GCC-Germline-Postprocess_Broad_GDAC_Methylation_Data-%s.RData", canType))


##
## Read TCGA Broad GDAC data
##
ncanType = ifelse(canType == "CRC", "COADREAD", canType)
gdacStdDataCanDir = file.path(gdacStdDataDir, ncanType, gdacStdDataDateNoUb)
gdacAnalysesCanDir = file.path(gdacAnalysesDir, ncanType, gdacAnalysesDateNoUb)
canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  


##
## Read methylation data for promoter regions
##
methyl27File = file.path(gdacStdDataCanDir,
                         sprintf("gdac.broadinstitute.org_%s.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                         sprintf("%s.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.promoter.txt", ncanType))
if (file.exists(methyl27File)) {
    methyl27Df = read.delim(methyl27File, header = T, as.is = T)
    if (ncanType == "SKCM") {
        methyl27PrimDf = cbind(methyl27Df$Gene,
                               methyl27Df[, 2:ncol(methyl27Df)][, which(as.integer(gsub(".*TCGA\\.\\S{2}\\.\\S{4}\\.(\\S{2})\\S{1}.*", "\\1", colnames(methyl27Df)[2:ncol(methyl27Df)], perl = T)) == 6)])
    } else {
        methyl27PrimDf = cbind(methyl27Df$Gene,
                               methyl27Df[, 2:ncol(methyl27Df)][, which(as.integer(gsub(".*TCGA\\.\\S{2}\\.\\S{4}\\.(\\S{2})\\S{1}.*", "\\1", colnames(methyl27Df)[2:ncol(methyl27Df)], perl = T)) == 1)])
    }
    colnames(methyl27PrimDf)[2:ncol(methyl27PrimDf)] = gsub(".*(TCGA\\.\\S{2}\\.\\S{4}).*", "\\1", colnames(methyl27PrimDf)[2:ncol(methyl27PrimDf)], perl = T)
    colnames(methyl27PrimDf)[2:ncol(methyl27PrimDf)] = gsub(".", "_", colnames(methyl27PrimDf)[2:ncol(methyl27PrimDf)], fixed = T)
    colnames(methyl27PrimDf)[1] = "Gene"
    methyl27PatientIds = colnames(methyl27PrimDf)[2:ncol(methyl27PrimDf)]
}

methyl450File = file.path(gdacStdDataCanDir,
                          sprintf("gdac.broadinstitute.org_%s.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                          sprintf("%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.promoter.txt", ncanType))
methyl450Df = read.delim(methyl450File, header = T, as.is = T)
colnames(methyl450Df) = gsub(".", "_", colnames(methyl450Df), fixed = T)
if (ncanType == "SKCM") {
    methyl450PrimDf = cbind(methyl450Df$Gene,
                            methyl450Df[, 2:ncol(methyl450Df)][, which(as.integer(gsub(".*TCGA_\\S{2}_\\S{4}_(\\S{2})\\S{1}.*", "\\1", colnames(methyl450Df)[2:ncol(methyl450Df)], perl = T)) == 6)])
} else {
    methyl450PrimDf = cbind(methyl450Df$Gene,
                            methyl450Df[, 2:ncol(methyl450Df)][, which(as.integer(gsub(".*TCGA_\\S{2}_\\S{4}_(\\S{2})\\S{1}.*", "\\1", colnames(methyl450Df)[2:ncol(methyl450Df)], perl = T)) == 1)])
}
colnames(methyl450PrimDf)[2:ncol(methyl450PrimDf)] = gsub(".*(TCGA_\\S{2}_\\S{4}).*", "\\1", colnames(methyl450PrimDf)[2:ncol(methyl450PrimDf)], perl = T)
colnames(methyl450PrimDf)[1] = "Gene"
methyl450PatientIds = colnames(methyl450PrimDf)[2:ncol(methyl450PrimDf)]

if (file.exists(methyl27File)) {
    sharedMethylPatientIds = intersect(methyl27PatientIds, methyl450PatientIds)
    if (length(sharedMethylPatientIds) > 0) {
        methyl27PrimDf2 = methyl27PrimDf[, -c(which(colnames(methyl27PrimDf) %in% sharedMethylPatientIds))]
    } else {
        methyl27PrimDf2 = methyl27PrimDf
    }
    methylMergedDf = sqldf('SELECT m4.*, m2.* FROM methyl450PrimDf AS m4 LEFT JOIN methyl27PrimDf2 AS m2 ON m4.Gene = m2.Gene')
    methylMergedDf = methylMergedDf[, -which(colnames(methylMergedDf) == "Gene")[2]]
} else {
    methylMergedDf = methyl450PrimDf
}

rownames(methylMergedDf) = methylMergedDf[,1]

methylMergedFile = file.path(tcgaCanDir, sprintf("Methylation-Tumor-%s.new.txt", canType))
write.table(methylMergedDf, methylMergedFile, sep = "\t", na = "NA", row.names = F, col.names = T, quote = F)
save.image(imageFile)


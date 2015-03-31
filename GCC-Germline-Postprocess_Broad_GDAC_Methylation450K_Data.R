rm(list=ls())

require('doMC')
registerDoMC(1)

require(tools)
require(data.table)

if (!is.na(pmatch("darwin", version$os))) {
        rootDir = "/Volumes/orchestra"
} else {
        rootDir = "/home/sl279"
}

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
gdacDir = file.path(germDir, "GDAC")
tcgaDir = file.path(germDir, "TCGA")

gdacStdDataDate = "2015_02_04"
gdacStdDataDateNoUb = gsub("_", "", gdacStdDataDate)
gdacStdDataDir = file.path(gdacDir, sprintf("stddata__%s", gdacStdDataDate))

args = commandArgs(TRUE)
canType = args[1]
ncanType = ifelse(canType == "CRC", "COADREAD", canType)
gdacStdDataCanDir = file.path(gdacStdDataDir, ncanType, gdacStdDataDateNoUb)

methylTcga450File = file.path(gdacStdDataCanDir,
                              sprintf("gdac.broadinstitute.org_%s.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                              sprintf("%s.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", ncanType))

if (!file.exists(methylTcga450File)) {
    cat(sprintf("Cannot find %s\n", methylTcga450File))
    quit(save = "no")
}

methylTcga450ModForGeneAvgFile = paste(file_path_sans_ext(methylTcga450File), ".promoter.txt", sep = "")

if (file.exists(methylTcga450ModForGeneAvgFile)) {
  cat(sprintf("%s exists!\n", methylTcga450ModForGeneAvgFile))
}

print(sprintf("Loading %s ...", methylTcga450File))
methylTcga450Dt = fread(methylTcga450File, header = T, stringsAsFactors=F)
methylTcga450ModDt = methylTcga450Dt[, c(1, seq(2,ncol(methylTcga450Dt),by=4)), with = F]
methylTcga450ModDt = methylTcga450ModDt[-1,]
setnames(methylTcga450ModDt, colnames(methylTcga450ModDt)[1], "Probe")

# Load Illumina Human Methylation 450K probe set information
methyl450PromoterExpandedFile = file.path(germDir, "Illumina/HumanMethylation450_15017482_v1-2.expanded.promoter.csv")
print(sprintf("Loading %s ...", methyl450PromoterExpandedFile))
methyl450PromoterExpandedDf = read.csv(methyl450PromoterExpandedFile, header = T, as.is = T)
genes = sort(unique(methyl450PromoterExpandedDf$UCSC_RefGene_Name))

methylTcga450ModForGeneAvgDf <- foreach (i=1:length(genes), .combine=rbind) %dopar% {
  gene = genes[i]
  print(gene)
  probeIdsForGene = unique(subset(methyl450PromoterExpandedDf, UCSC_RefGene_Name == gene)$IlmnID)
  methylTcga450ModForGeneDt = subset(methylTcga450ModDt, Probe %in% probeIdsForGene)
  methylTcga450ModForGeneDf = as.data.frame(methylTcga450ModForGeneDt)
  rownames(methylTcga450ModForGeneDf) = methylTcga450ModForGeneDf[,1]
  methylTcga450ModForGeneDf = methylTcga450ModForGeneDf[,-1]
  methylTcga450ModForGeneAvg = colMeans(data.matrix(methylTcga450ModForGeneDf), na.rm = T)
  #methylTcga450ModForGeneAvgDf = rbind(methylTcga450ModForGeneAvgDf, t(c(gene, methylTcga450ModForGeneAvg)))
  t(c(gene, methylTcga450ModForGeneAvg))
}

colnames(methylTcga450ModForGeneAvgDf)[1] = "Gene"
methylTcga450ModForGeneAvgDf[,1] = as.character(methylTcga450ModForGeneAvgDf[,1])
for (i in 2:ncol(methylTcga450ModForGeneAvgDf)) { 
  methylTcga450ModForGeneAvgDf[,i] = as.numeric(as.character(methylTcga450ModForGeneAvgDf[,i]))
}

methylTcga450ModForGeneAvgFile = paste(file_path_sans_ext(methylTcga450File), ".promoter.txt", sep = "")
print(sprintf("Saving %s ...", methylTcga450ModForGeneAvgFile))
write.table(methylTcga450ModForGeneAvgDf, methylTcga450ModForGeneAvgFile, quote = F, col.names = T, row.names = F, sep = "\t")


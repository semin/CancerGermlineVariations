rm(list=ls())

require(doMC)
registerDoMC(2)
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

methylTcga27File = file.path(gdacStdDataCanDir,
                             sprintf("gdac.broadinstitute.org_%s.Merge_methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                             sprintf("%s.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", ncanType))

if (!file.exists(methylTcga27File)) {
  cat(sprintf("Cannot find %s\n", methylTcga27File))
  quit(save = "no")  
}

print(sprintf("Loading %s ...", methylTcga27File))
methylTcga27Dt = fread(methylTcga27File, header = T, stringsAsFactors=F)
methylTcga27ModDt = methylTcga27Dt[, c(1, seq(2,ncol(methylTcga27Dt),by=4)), with = F]
methylTcga27ModDt = methylTcga27ModDt[-1,]
setnames(methylTcga27ModDt, colnames(methylTcga27ModDt)[1], "Probe")

# Load Illumina Human Methylation 27K probe set information
methyl27PromoterExpandedFile = file.path(germDir, "Illumina/illumina_humanmethylation27_content.promoter.txt")
print(sprintf("Loading %s ...", methyl27PromoterExpandedFile))
methyl27PromoterExpandedDf = read.delim(methyl27PromoterExpandedFile, header = T, as.is = T)
genes = sort(unique(methyl27PromoterExpandedDf$Symbol), na.last = NA)

methylTcga27ModForGeneAvgDf <- foreach (i=1:length(genes), .combine=rbind) %dopar% {
  gene = genes[i]
  print(gene)
  probeIdsForGene = unique(subset(methyl27PromoterExpandedDf, Symbol == gene)$Name)
  methylTcga27ModForGeneDt = subset(methylTcga27ModDt, Probe %in% probeIdsForGene)
  methylTcga27ModForGeneDf = as.data.frame(methylTcga27ModForGeneDt)
  rownames(methylTcga27ModForGeneDf) = methylTcga27ModForGeneDf[,1]
  methylTcga27ModForGeneDf = methylTcga27ModForGeneDf[,-1]
  methylTcga27ModForGeneAvg = colMeans(data.matrix(methylTcga27ModForGeneDf), na.rm = T)
  #methylTcga27ModForGeneAvgDf = rbind(methylTcga27ModForGeneAvgDf, t(c(gene, methylTcga27ModForGeneAvg)))
  t(c(gene, methylTcga27ModForGeneAvg))
}

colnames(methylTcga27ModForGeneAvgDf)[1] = "Gene"
methylTcga27ModForGeneAvgDf[methylTcga27ModForGeneAvgDf == "NaN"] = NA
#colnames(methylTcga27ModForGeneAvgDf)[2:ncol(methylTcga27ModDt)] = colnames(methylTcga27ModDt)[2:ncol(methylTcga27ModDt)]

methylTcga27ModForGeneAvgFile = paste(file_path_sans_ext(methylTcga27File), ".promoter.txt", sep = "")
print(sprintf("Saving %s ...", methylTcga27ModForGeneAvgFile))
write.table(methylTcga27ModForGeneAvgDf, methylTcga27ModForGeneAvgFile, quote = F, col.names = T, row.names = F, sep = "\t")


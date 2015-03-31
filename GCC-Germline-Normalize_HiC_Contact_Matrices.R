args = commandArgs(TRUE)
binSize = as.numeric(args[1])
rawObsFile = args[2]
krNormFile = args[3]
krExpFile = args[4]
outFile = args[5]

rawObsDf = read.delim(rawObsFile, header = F)
colnames(rawObsDf) = c("i", "j", "c")

krNormDf = read.delim(krNormFile, header = F)
colnames(krNormDf) = c("n") 

rawObsDf$nc = rawObsDf$c / (krNormDf[rawObsDf$i / binSize + 1,] * krNormDf[rawObsDf$j / binSize + 1,])

krExpDf = read.delim(krExpFile, header = F)
colnames(krExpDf) = c("e") 

if (basename(krExpFile) != "NA") {
    rawObsDf$oe = rawObsDf$nc / krExpDf[(abs(rawObsDf$j - rawObsDf$i) / binSize + 1),]
    rawObsDf$loe = log2(rawObsDf$oe)
}

save(rawObsDf, file = outFile)

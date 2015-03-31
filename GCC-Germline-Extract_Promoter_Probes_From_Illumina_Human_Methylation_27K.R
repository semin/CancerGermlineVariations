rm(list=ls())

require(data.table)
require(annotate)
require(org.Hs.eg.db)
require(ChIPpeakAnno)

# Load TSS table
data(TSS.human.NCBI36)
tssHumanNcbi36 = TSS.human.NCBI36
tssHumanNcbi36$feature = rownames(TSS.human.NCBI36)
tssHumanNcbi36 = addGeneIDs(tssHumanNcbi36, "org.Hs.eg.db", c("symbol"))
tssHumanNcbi36Df = as.data.frame(tssHumanNcbi36)

# Load RefSeq ID to gene name mapping table
geneInfoFile = file.path("/groups/kucherlapati/GCC/Germline/Entrez/gene_info_human.txt")
geneInfoDt = fread(geneInfoFile, header = F)
geneInfoDt = geneInfoDt[, c(2:5), with = F]
setnames(geneInfoDt, c("GeneID", "Symbol", "LocusTag", "Synonyms"))
geneInfoExtDt = geneInfoDt[, list(Synonym = unlist(strsplit(Synonyms, "|", fixed = T), use.names = FALSE)), by = "GeneID,Symbol,LocusTag"]
geneInfoExtDf = as.data.frame(geneInfoExtDt)

# Extract and expand promoter probes from Illumina Human Methylation 27K probe set
methyl27ProbeFile = file.path("/groups/kucherlapati/GCC/Germline/Illumina/illumina_humanmethylation27_content.txt")
methyl27ProbeDf = read.delim(methyl27ProbeFile, header = T, as.is = T)

# Update TSS_Coordinate
for (i in 1:nrow(methyl27ProbeDf)) {
    geneId = strsplit(methyl27ProbeDf[i, "Gene_ID"], ":", fixed = T)[[1]][2]
    print(geneId)
    geneId2Symbol = geneInfoExtDf[geneInfoExtDf$GeneID == as.numeric(geneId),][1,]$Symbol
    if (is.na(geneId2Symbol)) {
        entrez2Symbol = getSYMBOL(geneId, data='org.Hs.eg')
        if (!is.na(entrez2Symbol)) {
            methyl27ProbeDf[i, "Symbol"] = entrez2Symbol
        } else {
            alias2Symbol = geneInfoExtDf[geneInfoExtDf$Synonym == methyl27ProbeDf[i, "Symbol"],][1,]$Symbol
            if (!is.na(alias2Symbol)) {
                methyl27ProbeDf[i, "Symbol"] = alias2Symbol
            }
        }
    } else {
        methyl27ProbeDf[i, "Symbol"] = geneId2Symbol
    }
    if (is.na(methyl27ProbeDf[i, "TSS_Coordinate"])) {
        if (methyl27ProbeDf[i, "Symbol"] == "ZIM2") {
            # ZIM2 (hg18), Position: chr19:61,977,740-62,043,887 Size: 66,148 Total Exon Count: 11 Strand: -
            methyl27ProbeDf[i, "TSS_Coordinate"] = 62043887
            next
        }
        if (methyl27ProbeDf[i, "Symbol"] == "GAGE7B") { 
            methyl27ProbeDf[i, "Symbol"] = "GAGE7"
        }
        methyl27ProbeTmpDf = methyl27ProbeDf[(methyl27ProbeDf$Symbol == methyl27ProbeDf[i, "Symbol"]) &
                                             (!is.na(methyl27ProbeDf$TSS_Coordinate)),]
        if (nrow(methyl27ProbeTmpDf) > 0) {
            methyl27ProbeDf[i, "TSS_Coordinate"] = methyl27ProbeTmpDf[1, "TSS_Coordinate"]
        } else {
            tssHumanNcbi36TmpDf = subset(tssHumanNcbi36Df, grepl(methyl27ProbeDf[i, "Symbol"], symbol))
            if (nrow(tssHumanNcbi36TmpDf) > 0) {
                methyl27ProbeDf[i, "TSS_Coordinate"] = ifelse(tssHumanNcbi36TmpDf[1, "strand"] == 1,
                                                              tssHumanNcbi36TmpDf[1, "start"],
                                                              tssHumanNcbi36TmpDf[1, "end"])
            }
        }
    }
}

methyl27ProbeDf$TSS_Distance = as.numeric(methyl27ProbeDf$MapInfo) - as.numeric(methyl27ProbeDf$TSS_Coordinate)
methyl27ProbeDf$TSS_Distance = sapply(1:nrow(methyl27ProbeDf),
                                      function(x) { ifelse(methyl27ProbeDf[x,]$Gene_Strand == "-",
                                                           -methyl27ProbeDf[x,]$TSS_Distance,
                                                           methyl27ProbeDf[x,]$TSS_Distance) })

subset(methyl27ProbeDf, Symbol == "A1CF")
subset(methyl27ProbeDf, Symbol == "ABR")

#methyl27PromoterProbeDf = subset(methyl27ProbeDf, TSS_Distance >= -1500 & TSS_Distance <= 500)
methyl27PromoterProbeFile = file.path("/groups/kucherlapati/GCC/Germline/Illumina/illumina_humanmethylation27_content.promoter.txt3")
write.table(methyl27ProbeDf, methyl27PromoterProbeFile, quote = F, col.names = T, row.names = F, sep = "\t")

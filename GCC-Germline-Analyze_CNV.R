options(stringsAsFactors = FALSE)

require(ggplot2)
require(GenomicRanges)

chrs    = c(1:22)
chrs    = factor(chrs, level=chrs)
cchrs   = sapply(chrs, function(x) paste("chr", x, sep=""))
cchrs   = factor(cchrs, level=cchrs)

# contruct patient ID to cancer type mapping table
gdacDir     = file.path("/groups/kucherlapati/GCC/Germline/GDAC/CNV_Analysis")
gisticDir   = file.path("/groups/kucherlapati/GCC/Germline/CNV")
cancerDirs  = Sys.glob(file.path(gisticDir, "*"))
cancerDirs  = cancerDirs[grep("PANCAN", cancerDirs, invert = TRUE)]

my.t.test.p.value <- function(...) { 
      obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

my.wilcox.test.p.value <- function(...) { 
      obj<-try(wilcox.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else return(obj$p.value) 
} 

# Collect RNA expression data
date                = 20141017
sampleListDf        = data.frame()
rnaLog2FoldChangeLt = list()

for (cancerDir in cancerDirs) {
    # collect sample list for each cancer type
    cancerType = basename(cancerDir)

    if (cancerType == "CRC") cancerType = "COADREAD"

    print(sprintf("[%s] collecting mRNA expression data", cancerType))
    segFile = file.path(cancerDir, "segmentationfile.txt")
    segDf = read.delim(segFile, header = TRUE, as.is = TRUE)
    sampleIDs = as.vector(sapply(unique(segDf[,1]), function(x) gsub("-", "_", strsplit(x, "---")[[1]][1])))
    sampleListDf = rbind(sampleListDf, data.frame(sampleId = sampleIDs, cancerType = rep(cancerType, length(sampleIDs))))

    # collect RNA fold changes for each cancer type
    gdacRnaSeqRsemL3File = file.path(gdacDir,
                                     sprintf("gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.%d00.0.0", cancerType, date),
                                     sprintf("%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", cancerType))
    gdacRnaSeqRpkmL3File = file.path(gdacDir,
                                     sprintf("gdac.broadinstitute.org_%s.Merge_rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.Level_3.%d00.0.0", cancerType, date),
                                     sprintf("%s.rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.data.txt", cancerType))
    if (file.exists(gdacRnaSeqRsemL3File)) {
        gdacRnaSeqRsemL3Df = read.delim(gdacRnaSeqRsemL3File, header = FALSE, as.is = TRUE, skip = 2)
        headerLine = readLines(gdacRnaSeqRsemL3File, 1)
        headers = gsub("-", "_", as.vector(unlist(strsplit(headerLine, '\t'))))
        colnames(gdacRnaSeqRsemL3Df) = headers
        tumorSamples = as.vector(sapply(colnames(gdacRnaSeqRsemL3Df), function(x) { ifelse(as.integer(substr(strsplit(x, "_")[[1]][4], 1, 2)) < 10, x, NA) }))
        tumorSamples = tumorSamples[!is.na(tumorSamples)]
        gdacRnaSeqRsemL3TumorDf = gdacRnaSeqRsemL3Df[, tumorSamples]
        rownames(gdacRnaSeqRsemL3TumorDf) = gdacRnaSeqRsemL3Df[, 1]
        gdacRnaSeqFoldChangeDf = log2(gdacRnaSeqRsemL3TumorDf + 1) - log2(rowMeans(gdacRnaSeqRsemL3TumorDf + 1))
        rnaLog2FoldChangeLt[[cancerType]] = gdacRnaSeqFoldChangeDf
    } else if (file.exists(gdacRnaSeqRpkmL3File)) {
        gdacRnaSeqRpkmL3Df = read.delim(gdacRnaSeqRpkmL3File, header = FALSE, as.is = TRUE, skip = 2)
        headerLine = readLines(gdacRnaSeqRpkmL3File, 1)
        headers = gsub("-", "_", as.vector(unlist(strsplit(headerLine, '\t'))))
        colnames(gdacRnaSeqRpkmL3Df) = headers
        gdacRnaSeqRpkmL3Df = gdacRnaSeqRpkmL3Df[, c(1, seq(4, ncol(gdacRnaSeqRpkmL3Df), 3))]
        tumorSamples = as.vector(sapply(colnames(gdacRnaSeqRpkmL3Df), function(x) { ifelse(as.integer(substr(strsplit(x, "_")[[1]][4], 1, 2)) < 10, x, NA) }))
        tumorSamples = tumorSamples[!is.na(tumorSamples)]
        gdacRnaSeqRpkmL3TumorDf = gdacRnaSeqRpkmL3Df[, tumorSamples]
        rownames(gdacRnaSeqRpkmL3TumorDf) = gdacRnaSeqRpkmL3Df[,1]
        gdacRnaSeqFoldChangeDf = log2(gdacRnaSeqRpkmL3TumorDf + 1) - log2(rowMeans(gdacRnaSeqRpkmL3TumorDf + 1))
        rnaLog2FoldChangeLt[[cancerType]] = gdacRnaSeqFoldChangeDf
    } else {
        warning(sprintf("[%s] Cannot find RNA-seq data", cancerType))
        next
    }
}


# Collect clinical information (patient.bcrpatientbarcode, patient.gender, patient.race, patient.ethnicity
clinDf = data.frame()
clinFields = c("Hybridization REF", "gender", "race", "ethnicity")
for (cancerDir in cancerDirs) {
    # collect sample list for each cancer type
    cancerType = basename(cancerDir)
    if (cancerType == "CRC") cancerType = "COADREAD"
    print(sprintf("[%s] collecting clinical information", cancerType))
    clinFile    = file.path(gdacDir, sprintf("gdac.broadinstitute.org_%s.Clinical_Pick_Tier1.Level_4.%d00.0.0/%s.clin.merged.picked.txt", cancerType, date, cancerType))
    clinTmpDf   = read.delim(clinFile, header = FALSE, as.is = TRUE)
    clinTmpDf   = clinTmpDf[clinTmpDf[,1] %in% clinFields,]
    clinDf      = rbind(clinDf, data.frame(patientId    = gsub("-", "_", toupper(as.vector(unlist(clinTmpDf[clinTmpDf[, 1] == "Hybridization REF", 2:ncol(clinTmpDf)])))),
                                           gender       = as.vector(unlist(clinTmpDf[clinTmpDf[, 1] == "gender", 2:ncol(clinTmpDf)])),
                                           race         = as.vector(unlist(clinTmpDf[clinTmpDf[, 1] == "race", 2:ncol(clinTmpDf)])),
                                           ethnicity    = as.vector(unlist(clinTmpDf[clinTmpDf[, 1] == "ethnicity", 2:ncol(clinTmpDf)]))))
}

nonHispanicWhitePatients = subset(clinDf, race == "white" & ethnicity == "not hispanic or latino")

# Collect DGV data
dgvFile             = file.path("/groups/park/semin/BiO/Research/GCC/Germline/DGV/GRCh37_hg19_variants_2013-07-23.txt")
dgvDf               = read.delim(dgvFile, header = TRUE, as.is = TRUE)
colnames(dgvDf)[2]  = "space"
dgvDf$space         = factor(dgvDf$space, level=chrs)
dgvGr               = as(as(dgvDf, "RangedData"), "GRanges")

# Pop code
dgvPopCodeFile = file.path("/groups/park/semin/BiO/Research/GCC/Germline/DGV/Population_Code.tsv")
dgvPopCodeDf = read.delim(dgvPopCodeFile, header = FALSE, as.is = TRUE, comment.char = "#")
colnames(dgvPopCodeDf) = c("popCode", "supPopCode")
# PMID: 23128226
dgv1kgMatchFile = file.path("/groups/park/semin/BiO/Research/GCC/Germline/DGV/1000_Genomes_Sample_Population_Matching_Table.tsv")
dgv1kgMatchDf   = read.delim(dgv1kgMatchFile, header = FALSE, as.is = TRUE)
colnames(dgv1kgMatchDf) = c("sampleId", "popCode")
dgv1kgMatchDf$supPopCode = sapply(dgv1kgMatchDf$popCode, function(x) dgvPopCodeDf[dgvPopCodeDf$popCode == x,]$supPopCode)
dgv1kgAfrDf     = subset(dgv1kgMatchDf, supPopCode == "AFR")
dgv1kgAmrDf     = subset(dgv1kgMatchDf, supPopCode == "AMR")
dgv1kgAsnDf     = subset(dgv1kgMatchDf, supPopCode == "ASN")
dgv1kgEurDf     = subset(dgv1kgMatchDf, supPopCode == "EUR")
dgv1kgSanDf     = subset(dgv1kgMatchDf, supPopCode == "SAN")

# PMID: 20811451
dgvHapMap3MatchFile = file.path("/groups/park/semin/BiO/Research/GCC/Germline/DGV/HapMap3_Sample_Population_Matching_Table.tsv")
dgvHapMap3MatchDf   = read.delim(dgvHapMap3MatchFile, header = FALSE, as.is = TRUE)
colnames(dgvHapMap3MatchDf) = c("sampleId", "popCode")
dgvHapMap3MatchDf$supPopCode = sapply(dgvHapMap3MatchDf$popCode, function(x) dgvPopCodeDf[dgvPopCodeDf$popCode == x,]$supPopCode)
dgvHapMap3MatchDf$sampleId = gsub("GM", "NA", dgvHapMap3MatchDf$sampleId, fixed = TRUE)
dgvHapMap3AfrDf     = subset(dgvHapMap3MatchDf, supPopCode == "AFR")
dgvHapMap3AmrDf     = subset(dgvHapMap3MatchDf, supPopCode == "AMR")
dgvHapMap3AsnDf     = subset(dgvHapMap3MatchDf, supPopCode == "ASN")
dgvHapMap3EurDf     = subset(dgvHapMap3MatchDf, supPopCode == "EUR")
dgvHapMap3SanDf     = subset(dgvHapMap3MatchDf, supPopCode == "SAN")

# Check frequency of certain peaks in the pop
gisticPanCanDir                 = file.path(gisticDir, "PANCAN/Results")
gisticAllLesionsFile            = file.path(gisticPanCanDir, "all_lesions.conf_99.txt")
gisticAllLesionsDf              = read.delim(gisticAllLesionsFile, header = TRUE, as.is = TRUE)
colnames(gisticAllLesionsDf)    = gsub("\\.", "_", colnames(gisticAllLesionsDf))
colnames(gisticAllLesionsDf)    = gsub("___L3_B1000_NS_NBICseq", "", colnames(gisticAllLesionsDf))

# Read amplification/deletion peaks information and check distribution
ampGenesFile    = file.path(gisticPanCanDir, "amp_genes.conf_99.txt")
ampGenesDf      = read.delim(ampGenesFile, header = TRUE, as.is = TRUE)
delGenesFile    = file.path(gisticPanCanDir, "del_genes.conf_99.txt")
delGenesDf      = read.delim(delGenesFile, header = TRUE, as.is = TRUE)

# Check overlapping between germline CNVs and known DGV entries
peakStatsDf     = data.frame()
cancerTypes     = unique(sampleListDf$cancerType)

for (peakType in c("Amp", "Del")) {
    if (peakType == "Amp") {
        peaksDf = ampGenesDf
    } else {
        peaksDf = delGenesDf
    }
    # Amp:
    # i = 67: SNX16
    # i = 139: GOLPH3
    # Del:
    # i = 5: GSTM1
    # i = 339: APOBEC
    for (i in 2:(ncol(peaksDf)-1)) {
        cytoband        = gsub("^X((\\d+|X|Y)(p|q)\\d+)\\.(\\S+)", "\\1.\\2", colnames(peaksDf)[i], perl = TRUE)
        qvalue          = peaksDf[1, i]
        peakBoundary    = peaksDf[3, i]
        peakChr         = gsub("chr", "", strsplit(peakBoundary, ":")[[1]][1])
        peakPosStart    = as.integer(strsplit(strsplit(peakBoundary, ":")[[1]][2], "-")[[1]][1])
        peakPosStop     = as.integer(strsplit(strsplit(peakBoundary, ":")[[1]][2], "-")[[1]][2])
        peakGr          = GRanges(seqnames = peakChr, ranges = IRanges(peakPosStart, peakPosStop), strand = "*")
        genesInPeak     = peaksDf[4:nrow(peaksDf), i]
        genesInPeak     = genesInPeak[genesInPeak != ""]
        dgvPeakOl       = findOverlaps(peakGr, dgvGr, maxgap = 0L, minoverlap = 1L, select = "all")
        dgvPeakOlWidth  = width(pintersect(peakGr[queryHits(dgvPeakOl)], dgvGr[subjectHits(dgvPeakOl)]))

        print(sprintf(">>> Processing a %s peak, %s (%d/%d) ... ", peakType, peakBoundary, i - 1, ncol(peaksDf) - 1))

        cnvPeak                     = gisticAllLesionsDf[grep(peakBoundary, gisticAllLesionsDf$Wide_Peak_Limits, fixed = TRUE),][1, 10:(ncol(gisticAllLesionsDf) - 1)]
        affectedSamples             = colnames(cnvPeak[, cnvPeak > 0])
        unaffectedSamples           = colnames(cnvPeak[, cnvPeak < 1])
        affectedPatients            = as.vector(unlist(sapply(affectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
        unaffectedPatients          = as.vector(unlist(sapply(unaffectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
        affectedWhitePatients       = affectedPatients[affectedPatients %in% nonHispanicWhitePatients$patientId]
        unaffectedWhitePatients     = unaffectedPatients[unaffectedPatients %in% nonHispanicWhitePatients$patientId]
        frcAffectedPatients         = length(affectedPatients) / (length(affectedPatients) + length(unaffectedPatients))
        frcAffectedWhitePatients    = length(affectedWhitePatients) / (length(affectedWhitePatients) + length(unaffectedWhitePatients))

        if (length(dgvPeakOl) > 0) {
            for (j in 1:length(dgvPeakOl)) {
                hit = dgvGr[subjectHits(dgvPeakOl)[j]]
                if ((peakType == "Amp" & !is.na(hit$observedgains) & !is.na(hit$samplesize)) |
                    (peakType == "Del" & !is.na(hit$observedlosses) & !is.na(hit$samplesize))) {
                    if (((peakType == "Amp" & hit$observedgains > 0) | (peakType == "Del" & hit$observedlosses > 0)) & hit$samplesize > 100 & width(hit) > 1000) {

                        sampleIds = strsplit(hit$samples, ",", fixed = TRUE)[[1]]

                        if (peakType == "Amp") {
                            contTable       = as.table(cbind(c(length(affectedPatients), length(unaffectedPatients)),
                                                             c(hit$observedgains, hit$samplesize - hit$observedgains)))
                            contTableWhite  = as.table(cbind(c(length(affectedWhitePatients), length(unaffectedWhitePatients)),
                                                             c(hit$observedgains, hit$samplesize - hit$observedgains)))
                        } else {
                            contTable       = as.table(cbind(c(length(affectedPatients), length(unaffectedPatients)),
                                                             c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                            contTableWhite  = as.table(cbind(c(length(affectedWhitePatients), length(unaffectedWhitePatients)),
                                                             c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                        }

                        chisqTest       = chisq.test(contTable)
                        chisqTestWhite  = chisq.test(contTableWhite)
                        fisherTest      = fisher.test(contTable)
                        fisherTestWhite = fisher.test(contTableWhite)

                        if (hit$pubmedid == 23128226 | hit$pubmedid == 20811451) {
                            if (hit$pubmedid == 23128226) {
                                dgvEurDf    = dgv1kgEurDf
                                dgvAfrDf    = dgv1kgAfrDf
                                dgvAmrDf    = dgv1kgAmrDf
                                dgvAsnDf    = dgv1kgAsnDf
                                dgvSanDf    = dgv1kgSanDf
                            } else {
                                dgvEurDf    = dgvHapMap3EurDf
                                dgvAfrDf    = dgvHapMap3AfrDf
                                dgvAmrDf    = dgvHapMap3AmrDf
                                dgvAsnDf    = dgvHapMap3AsnDf
                                dgvSanDf    = dgvHapMap3SanDf
                            }

                            eurSamples  = sampleIds[sampleIds %in% dgvEurDf$sampleId]
                            afrSamples  = sampleIds[sampleIds %in% dgvAfrDf$sampleId]
                            amrSamples  = sampleIds[sampleIds %in% dgvAmrDf$sampleId]
                            asnSamples  = sampleIds[sampleIds %in% dgvAsnDf$sampleId]
                            sanSamples  = sampleIds[sampleIds %in% dgvSanDf$sampleId]

                            if (peakType == "Amp") {
                                contTableWhiteEur = as.table(cbind(c(length(affectedWhitePatients), length(unaffectedWhitePatients)),
                                                                    c(length(eurSamples), nrow(dgvEurDf) - length(eurSamples))))
                            } else {
                                contTableWhiteEur = as.table(cbind(c(length(affectedWhitePatients), length(unaffectedWhitePatients)),
                                                                    c(length(eurSamples), nrow(dgvEurDf) - length(eurSamples))))
                            }

                            chisqTestWhiteEur   = chisq.test(contTableWhiteEur)
                            fisherTestWhiteEur  = fisher.test(contTableWhiteEur)

                            print(sprintf("!!! Found an overlapping %s (chr%s:%d-%d, query: %6.2f%%, subject: %6.2f%%). Total: %6.2f%% (%d/%d), EUR: %6.2f%% (%d/%d) in %s", 
                                          peakType, seqnames(hit), start(hit), end(hit),
                                          100 * dgvPeakOlWidth[j] / width(peakGr),
                                          100 * dgvPeakOlWidth[j] / width(hit),
                                          100 * ifelse(peakType == "Amp", hit$observedgains, hit$observedlosses) / hit$samplesize,
                                          ifelse(peakType == "Amp", hit$observedgains, hit$observedlosses), hit$samplesize,
                                          100 * length(eurSamples) / nrow(dgvEurDf), length(eurSamples), nrow(dgvEurDf),
                                          hit$reference))

                            peakStatsDf         = rbind(peakStatsDf, data.frame(eventType                   = peakType,
                                                                                peakBoundary                = peakBoundary,
                                                                                peakOverlapPct              = round(100 * dgvPeakOlWidth[j] / width(peakGr), 2),
                                                                                cytoband                    = cytoband,
                                                                                qvalue                      = qvalue,
                                                                                genesInPeak                 = paste(genesInPeak, collapse = ","),
                                                                                cancerType                  = "PANCAN12",
                                                                                numAffectedPatients         = length(affectedPatients),
                                                                                numUnaffectedPatients       = length(unaffectedPatients),
                                                                                pctAffectedPatients         = round(100 * frcAffectedPatients, 2),
                                                                                numAffectedWhitePatients    = length(affectedWhitePatients),
                                                                                numUnaffectedWhitePatients  = length(unaffectedWhitePatients),
                                                                                pctAffectedWhitePatients    = round(100 * frcAffectedWhitePatients, 2),
                                                                                dgvAccession                = hit$variantaccession,
                                                                                dgvChr                      = seqnames(hit),
                                                                                dgvPosStart                 = start(hit),
                                                                                dgvPosStop                  = end(hit),
                                                                                dgvOverlapPct               = round(100 * dgvPeakOlWidth[j] / width(hit), 2),
                                                                                dgvVarType                  = hit$varianttype,
                                                                                dgvVarSubType               = hit$variantsubtype,
                                                                                dgvReference                = hit$reference,
                                                                                dgvPubmedId                 = hit$pubmedid,
                                                                                dgvSampleSize               = hit$samplesize,
                                                                                dgvObsGains                 = hit$observedgains,
                                                                                dgvObsLosses                = hit$observedlosses,
                                                                                dgvPctGains                 = round(100 * hit$observedgains / hit$samplesize, 2),
                                                                                dgvPctLosses                = round(100 * hit$observedlosses / hit$samplesize, 2),
                                                                                dgvEurSampleSize            = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, nrow(dgvEurDf)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, nrow(dgvEurDf))),
                                                                                dgvEurObsGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, length(eurSamples)),
                                                                                                                     NA),
                                                                                dgvEurObsLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, length(eurSamples)),
                                                                                                                     NA),
                                                                                dgvEurPctGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(100 * length(eurSamples) / nrow(dgvEurDf), 2)),
                                                                                                                     NA),
                                                                                dgvEurPctLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(100 * length(eurSamples) / nrow(dgvEurDf), 2)),
                                                                                                                     NA),
                                                                                dgvAfrSampleSize            = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, nrow(dgvAfrDf)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, nrow(dgvAfrDf))),
                                                                                dgvAfrObsGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, length(afrSamples)),
                                                                                                                     NA),
                                                                                dgvAfrObsLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, length(afrSamples)),
                                                                                                                     NA),
                                                                                dgvAfrPctGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(100 * length(afrSamples) / nrow(dgvAfrDf), 2)),
                                                                                                                     NA),
                                                                                dgvAfrPctLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(100 * length(afrSamples) / nrow(dgvAfrDf), 2)),
                                                                                                                     NA),
                                                                                dgvAmrSampleSize            = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, nrow(dgvAmrDf)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, nrow(dgvAmrDf))),
                                                                                dgvAmrObsGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, length(amrSamples)),
                                                                                                                     NA),
                                                                                dgvAmrObsLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, length(amrSamples)),
                                                                                                                     NA),
                                                                                dgvAmrPctGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(100 * length(amrSamples) / nrow(dgvAmrDf), 2)),
                                                                                                                     NA),
                                                                                dgvAmrPctLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(100 * length(amrSamples) / nrow(dgvAmrDf), 2)),
                                                                                                                     NA),
                                                                                dgvAsnSampleSize            = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, nrow(dgvAsnDf)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, nrow(dgvAsnDf))),
                                                                                dgvAsnObsGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, length(asnSamples)),
                                                                                                                     NA),
                                                                                dgvAsnObsLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, length(asnSamples)),
                                                                                                                     NA),
                                                                                dgvAsnPctGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(100 * length(asnSamples) / nrow(dgvAsnDf), 2)),
                                                                                                                     NA),
                                                                                dgvAsnPctLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(100 * length(asnSamples) / nrow(dgvAsnDf), 2)),
                                                                                                                     NA),
                                                                                dgvSanSampleSize            = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, nrow(dgvSanDf)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, nrow(dgvSanDf))),
                                                                                dgvSanObsGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, length(sanSamples)),
                                                                                                                     NA),
                                                                                dgvSanObsLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, length(sanSamples)),
                                                                                                                     NA),
                                                                                dgvSanPctGains              = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(100 * length(sanSamples) / nrow(dgvSanDf), 2)),
                                                                                                                     NA),
                                                                                dgvSanPctLosses             = ifelse(peakType == "Del",
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(100 * length(sanSamples) / nrow(dgvSanDf), 2)),
                                                                                                                     NA),
                                                                                dgvGenes                    = hit$genes,
                                                                                chisqTestPvalue             = round(chisqTest$p.value, 2),
                                                                                chisqTestWhitePvalue        = round(chisqTestWhite$p.value, 2),
                                                                                chisqTestWhiteEurPvalue     = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(chisqTestWhiteEur$p.value, 2)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(chisqTestWhiteEur$p.value, 2))),
                                                                                fisherTestPvalue            = round(fisherTest$p.value, 2),
                                                                                fisherTestWhitePvalue       = round(fisherTestWhite$p.value, 2),
                                                                                fisherTestWhiteEurPvalue    = ifelse(peakType == "Amp",
                                                                                                                     ifelse(hit$observedlosses > 0, NA, round(fisherTestWhiteEur$p.value, 2)),
                                                                                                                     ifelse(hit$observedgains > 0, NA, round(fisherTestWhiteEur$p.value, 2)))))


                            for (cancerType in cancerTypes) {
                                cancerTypeSamples                   = as.vector(sampleListDf[sampleListDf$cancerType == cancerType,]$sampleId)
                                cancerTypeAffectedSamples           = cancerTypeSamples[cancerTypeSamples %in% affectedSamples]
                                cancerTypeUnaffectedSamples         = cancerTypeSamples[cancerTypeSamples %in% unaffectedSamples]
                                cancerTypeAffectedPatients          = as.vector(unlist(sapply(cancerTypeAffectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
                                cancerTypeUnaffectedPatients        = as.vector(unlist(sapply(cancerTypeUnaffectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
                                cancerTypeAffectedWhitePatients     = cancerTypeAffectedPatients[cancerTypeAffectedPatients %in% nonHispanicWhitePatients$patientId]
                                cancerTypeUnaffectedWhitePatients   = cancerTypeUnaffectedPatients[cancerTypeUnaffectedPatients %in% nonHispanicWhitePatients$patientId]
                                cancerTypeFrcAffectedPatients       = length(cancerTypeAffectedPatients) / (length(cancerTypeAffectedPatients) + length(cancerTypeUnaffectedPatients))
                                cancerTypeFrcAffectedWhitePatients  = length(cancerTypeAffectedWhitePatients) / (length(cancerTypeAffectedWhitePatients) + length(cancerTypeUnaffectedWhitePatients))

                                if (peakType == "Amp") {
                                    cancerTypeContTable           = as.table(cbind(c(length(cancerTypeAffectedPatients), length(cancerTypeUnaffectedPatients)),
                                                                         c(hit$observedgains, hit$samplesize - hit$observedgains)))
                                    cancerTypeContTableWhite      = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(hit$observedgains, hit$samplesize - hit$observedgains)))
                                    cancerTypeContTableWhiteEur   = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(length(eurSamples), nrow(dgvEurDf) - length(eurSamples))))
                                } else {
                                    cancerTypeContTable           = as.table(cbind(c(length(cancerTypeAffectedPatients), length(cancerTypeUnaffectedPatients)),
                                                                         c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                                    cancerTypeContTableWhite      = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                                    cancerTypeContTableWhiteEur   = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(length(eurSamples), nrow(dgvEurDf) - length(eurSamples))))
                                }

                                cancerTypeChisqTest                 = chisq.test(cancerTypeContTable)
                                cancerTypeChisqTestWhite            = chisq.test(cancerTypeContTableWhite)
                                cancerTypeChisqTestWhiteEur         = chisq.test(cancerTypeContTableWhiteEur)
                                cancerTypeFisherTest                = fisher.test(cancerTypeContTable)
                                cancerTypeFisherTestWhite           = fisher.test(cancerTypeContTableWhite)
                                cancerTypeFisherTestWhiteEur        = fisher.test(cancerTypeContTableWhiteEur)

                                peakStatsDf                         = rbind(peakStatsDf, data.frame(eventType                   = peakType,
                                                                                                    peakBoundary                = peakBoundary,
                                                                                                    peakOverlapPct              = round(100 * dgvPeakOlWidth[j] / width(peakGr), 2),
                                                                                                    cytoband                    = cytoband,
                                                                                                    qvalue                      = qvalue,
                                                                                                    genesInPeak                 = paste(genesInPeak, collapse = ","),
                                                                                                    cancerType                  = cancerType,
                                                                                                    numAffectedPatients         = length(cancerTypeAffectedPatients),
                                                                                                    numUnaffectedPatients       = length(cancerTypeUnaffectedPatients),
                                                                                                    pctAffectedPatients         = round(100 * cancerTypeFrcAffectedPatients, 2),
                                                                                                    numAffectedWhitePatients    = length(cancerTypeAffectedWhitePatients),
                                                                                                    numUnaffectedWhitePatients  = length(cancerTypeUnaffectedWhitePatients),
                                                                                                    pctAffectedWhitePatients    = round(100 * cancerTypeFrcAffectedWhitePatients, 2),
                                                                                                    dgvAccession                = hit$variantaccession,
                                                                                                    dgvChr                      = seqnames(hit),
                                                                                                    dgvPosStart                 = start(hit),
                                                                                                    dgvPosStop                  = end(hit),
                                                                                                    dgvOverlapPct               = round(100 * dgvPeakOlWidth[j] / width(hit), 2),
                                                                                                    dgvVarType                  = hit$varianttype,
                                                                                                    dgvVarSubType               = hit$variantsubtype,
                                                                                                    dgvReference                = hit$reference,
                                                                                                    dgvPubmedId                 = hit$pubmedid,
                                                                                                    dgvSampleSize               = hit$samplesize,
                                                                                                    dgvObsGains                 = hit$observedgains,
                                                                                                    dgvObsLosses                = hit$observedlosses,
                                                                                                    dgvPctGains                 = round(100 * hit$observedgains / hit$samplesize, 2),
                                                                                                    dgvPctLosses                = round(100 * hit$observedlosses / hit$samplesize, 2),
                                                                                                    dgvEurSampleSize            = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, nrow(dgvEurDf)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, nrow(dgvEurDf))),
                                                                                                    dgvEurObsGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, length(eurSamples)),
                                                                                                                                         NA),
                                                                                                    dgvEurObsLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, length(eurSamples)),
                                                                                                                                         NA),
                                                                                                    dgvEurPctGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(100 * length(eurSamples) / nrow(dgvEurDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvEurPctLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(100 * length(eurSamples) / nrow(dgvEurDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAfrSampleSize            = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, nrow(dgvAfrDf)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, nrow(dgvAfrDf))),
                                                                                                    dgvAfrObsGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, length(afrSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAfrObsLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, length(afrSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAfrPctGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(100 * length(afrSamples) / nrow(dgvAfrDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAfrPctLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(100 * length(afrSamples) / nrow(dgvAfrDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAmrSampleSize            = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, nrow(dgvAmrDf)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, nrow(dgvAmrDf))),
                                                                                                    dgvAmrObsGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, length(amrSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAmrObsLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, length(amrSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAmrPctGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(100 * length(amrSamples) / nrow(dgvAmrDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAmrPctLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(100 * length(amrSamples) / nrow(dgvAmrDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAsnSampleSize            = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, nrow(dgvAsnDf)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, nrow(dgvAsnDf))),
                                                                                                    dgvAsnObsGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, length(asnSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAsnObsLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, length(asnSamples)),
                                                                                                                                         NA),
                                                                                                    dgvAsnPctGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(100 * length(asnSamples) / nrow(dgvAsnDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvAsnPctLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(100 * length(asnSamples) / nrow(dgvAsnDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvSanSampleSize            = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, nrow(dgvSanDf)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, nrow(dgvSanDf))),
                                                                                                    dgvSanObsGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, length(sanSamples)),
                                                                                                                                         NA),
                                                                                                    dgvSanObsLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, length(sanSamples)),
                                                                                                                                         NA),
                                                                                                    dgvSanPctGains              = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(100 * length(sanSamples) / nrow(dgvSanDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvSanPctLosses             = ifelse(peakType == "Del",
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(100 * length(sanSamples) / nrow(dgvSanDf), 2)),
                                                                                                                                         NA),
                                                                                                    dgvGenes                    = hit$genes,
                                                                                                    chisqTestPvalue             = round(cancerTypeChisqTest$p.value, 2),
                                                                                                    chisqTestWhitePvalue        = round(cancerTypeChisqTestWhite$p.value, 2),
                                                                                                    chisqTestWhiteEurPvalue     = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(cancerTypeChisqTestWhiteEur$p.value, 2)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(cancerTypeChisqTestWhiteEur$p.value, 2))),
                                                                                                    fisherTestPvalue            = round(fisherTest$p.value, 2),
                                                                                                    fisherTestWhitePvalue       = round(fisherTestWhite$p.value, 2),
                                                                                                    fisherTestWhiteEurPvalue    = ifelse(peakType == "Amp",
                                                                                                                                         ifelse(hit$observedlosses > 0, NA, round(cancerTypeFisherTestWhiteEur$p.value, 2)),
                                                                                                                                         ifelse(hit$observedgains > 0, NA, round(cancerTypeFisherTestWhiteEur$p.value, 2)))))
                            }
                        } else {
                            print(sprintf("!!! Found an overlapping %s (chr%s:%d-%d, query: %6.2f%%, subject: %6.2f%%). Total: %6.2f%% (%d/%d) in %s", 
                                          peakType, seqnames(hit), start(hit), end(hit),
                                          100 * dgvPeakOlWidth[j] / width(peakGr),
                                          100 * dgvPeakOlWidth[j] / width(hit),
                                          100 *ifelse(peakType == "Amp", hit$observedgains, hit$observedlosses) / hit$samplesize,
                                          ifelse(peakType == "Amp", hit$observedgains, hit$observedlosses), hit$samplesize,
                                          hit$reference))

                            peakStatsDf = rbind(peakStatsDf, data.frame(eventType                   = peakType,
                                                                        peakBoundary                = peakBoundary,
                                                                        peakOverlapPct              = round(100 * dgvPeakOlWidth[j] / width(peakGr), 2),
                                                                        cytoband                    = cytoband,
                                                                        qvalue                      = qvalue,
                                                                        genesInPeak                 = paste(genesInPeak, collapse = ","),
                                                                        cancerType                  = "PANCAN12",
                                                                        numAffectedPatients         = length(affectedPatients),
                                                                        numUnaffectedPatients       = length(unaffectedPatients),
                                                                        pctAffectedPatients         = round(100 * frcAffectedPatients, 2),
                                                                        numAffectedWhitePatients    = length(affectedWhitePatients),
                                                                        numUnaffectedWhitePatients  = length(unaffectedWhitePatients),
                                                                        pctAffectedWhitePatients    = round(100 * frcAffectedWhitePatients, 2),
                                                                        dgvAccession                = hit$variantaccession,
                                                                        dgvChr                      = seqnames(hit),
                                                                        dgvPosStart                 = start(hit),
                                                                        dgvPosStop                  = end(hit),
                                                                        dgvOverlapPct               = round(100 * dgvPeakOlWidth[j] / width(hit), 2),
                                                                        dgvVarType                  = hit$varianttype,
                                                                        dgvVarSubType               = hit$variantsubtype,
                                                                        dgvReference                = hit$reference,
                                                                        dgvPubmedId                 = hit$pubmedid,
                                                                        dgvSampleSize               = hit$samplesize,
                                                                        dgvObsGains                 = hit$observedgains,
                                                                        dgvObsLosses                = hit$observedlosses,
                                                                        dgvPctGains                 = round(100 * hit$observedgains / hit$samplesize, 2),
                                                                        dgvPctLosses                = round(100 * hit$observedlosses / hit$samplesize, 2),
                                                                        dgvEurSampleSize            = NA,
                                                                        dgvEurObsGains              = NA,
                                                                        dgvEurObsLosses             = NA,
                                                                        dgvEurPctGains              = NA,
                                                                        dgvEurPctLosses             = NA,
                                                                        dgvAfrSampleSize            = NA,
                                                                        dgvAfrObsGains              = NA,
                                                                        dgvAfrObsLosses             = NA,
                                                                        dgvAfrPctGains              = NA,
                                                                        dgvAfrPctLosses             = NA,
                                                                        dgvAmrSampleSize            = NA,
                                                                        dgvAmrObsGains              = NA,
                                                                        dgvAmrObsLosses             = NA,
                                                                        dgvAmrPctGains              = NA,
                                                                        dgvAmrPctLosses             = NA,
                                                                        dgvAsnSampleSize            = NA,
                                                                        dgvAsnObsGains              = NA,
                                                                        dgvAsnObsLosses             = NA,
                                                                        dgvAsnPctGains              = NA,
                                                                        dgvAsnPctLosses             = NA,
                                                                        dgvSanSampleSize            = NA,
                                                                        dgvSanObsGains              = NA,
                                                                        dgvSanObsLosses             = NA,
                                                                        dgvSanPctGains              = NA,
                                                                        dgvSanPctLosses             = NA,
                                                                        dgvGenes                    = hit$genes,
                                                                        chisqTestPvalue             = round(chisqTest$p.value, 2),
                                                                        chisqTestWhitePvalue        = round(chisqTestWhite$p.value, 2),
                                                                        chisqTestWhiteEurPvalue     = NA,
                                                                        fisherTestPvalue            = round(fisherTest$p.value, 2),
                                                                        fisherTestWhitePvalue       = round(fisherTestWhite$p.value, 2),
                                                                        fisherTestWhiteEurPvalue    = NA))

                            for (cancerType in cancerTypes) {
                                cancerTypeSamples                   = as.vector(sampleListDf[sampleListDf$cancerType == cancerType,]$sampleId)
                                cancerTypeAffectedSamples           = cancerTypeSamples[cancerTypeSamples %in% affectedSamples]
                                cancerTypeUnaffectedSamples         = cancerTypeSamples[cancerTypeSamples %in% unaffectedSamples]
                                cancerTypeAffectedPatients          = as.vector(unlist(sapply(cancerTypeAffectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
                                cancerTypeUnaffectedPatients        = as.vector(unlist(sapply(cancerTypeUnaffectedSamples, function(x) paste(strsplit(x, "_")[[1]][1:3], collapse = "_"))))
                                cancerTypeAffectedWhitePatients     = cancerTypeAffectedPatients[cancerTypeAffectedPatients %in% nonHispanicWhitePatients$patientId]
                                cancerTypeUnaffectedWhitePatients   = cancerTypeUnaffectedPatients[cancerTypeUnaffectedPatients %in% nonHispanicWhitePatients$patientId]
                                cancerTypeFrcAffectedPatients       = length(cancerTypeAffectedPatients) / (length(cancerTypeAffectedPatients) + length(cancerTypeUnaffectedPatients))
                                cancerTypeFrcAffectedWhitePatients  = length(cancerTypeAffectedWhitePatients) / (length(cancerTypeAffectedWhitePatients) + length(cancerTypeUnaffectedWhitePatients))

                                if (peakType == "Amp") {
                                    cancerTypeContTable         = as.table(cbind(c(length(cancerTypeAffectedPatients), length(cancerTypeUnaffectedPatients)),
                                                                         c(hit$observedgains, hit$samplesize - hit$observedgains)))
                                    cancerTypeContTableWhite    = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(hit$observedgains, hit$samplesize - hit$observedgains)))
                                } else {
                                    cancerTypeContTable         = as.table(cbind(c(length(cancerTypeAffectedPatients), length(cancerTypeUnaffectedPatients)),
                                                                         c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                                    cancerTypeContTableWhite    = as.table(cbind(c(length(cancerTypeAffectedWhitePatients), length(cancerTypeUnaffectedWhitePatients)),
                                                                         c(hit$observedlosses, hit$samplesize - hit$observedlosses)))
                                }

                                cancerTypeChisqTest         = chisq.test(cancerTypeContTable)
                                cancerTypeChisqTestWhite    = chisq.test(cancerTypeContTableWhite)
                                cancerTypeFisherTest        = fisher.test(cancerTypeContTable)
                                cancerTypeFisherTestWhite   = fisher.test(cancerTypeContTableWhite)

                                peakStatsDf = rbind(peakStatsDf, data.frame(eventType                   = peakType,
                                                                            peakBoundary                = peakBoundary,
                                                                            peakOverlapPct              = round(100 * dgvPeakOlWidth[j] / width(peakGr), 2),
                                                                            cytoband                    = cytoband,
                                                                            qvalue                      = qvalue,
                                                                            genesInPeak                 = paste(genesInPeak, collapse = ","),
                                                                            cancerType                  = cancerType,
                                                                            numAffectedPatients         = length(cancerTypeAffectedPatients),
                                                                            numUnaffectedPatients       = length(cancerTypeUnaffectedPatients),
                                                                            pctAffectedPatients         = round(100 * cancerTypeFrcAffectedPatients, 2),
                                                                            numAffectedWhitePatients    = length(cancerTypeAffectedWhitePatients),
                                                                            numUnaffectedWhitePatients  = length(cancerTypeUnaffectedWhitePatients),
                                                                            pctAffectedWhitePatients    = round(100 * cancerTypeFrcAffectedWhitePatients, 2),
                                                                            dgvAccession                = hit$variantaccession,
                                                                            dgvChr                      = seqnames(hit),
                                                                            dgvPosStart                 = start(hit),
                                                                            dgvPosStop                  = end(hit),
                                                                            dgvOverlapPct               = round(100 * dgvPeakOlWidth[j] / width(hit), 2),
                                                                            dgvVarType                  = hit$varianttype,
                                                                            dgvVarSubType               = hit$variantsubtype,
                                                                            dgvReference                = hit$reference,
                                                                            dgvPubmedId                 = hit$pubmedid,
                                                                            dgvSampleSize               = hit$samplesize,
                                                                            dgvObsGains                 = hit$observedgains,
                                                                            dgvObsLosses                = hit$observedlosses,
                                                                            dgvPctGains                 = round(100 * hit$observedgains / hit$samplesize, 2),
                                                                            dgvPctLosses                = round(100 * hit$observedlosses / hit$samplesize, 2),
                                                                            dgvEurSampleSize            = NA,
                                                                            dgvEurObsGains              = NA,
                                                                            dgvEurObsLosses             = NA,
                                                                            dgvEurPctGains              = NA,
                                                                            dgvEurPctLosses             = NA,
                                                                            dgvAfrSampleSize            = NA,
                                                                            dgvAfrObsGains              = NA,
                                                                            dgvAfrObsLosses             = NA,
                                                                            dgvAfrPctGains              = NA,
                                                                            dgvAfrPctLosses             = NA,
                                                                            dgvAmrSampleSize            = NA,
                                                                            dgvAmrObsGains              = NA,
                                                                            dgvAmrObsLosses             = NA,
                                                                            dgvAmrPctGains              = NA,
                                                                            dgvAmrPctLosses             = NA,
                                                                            dgvAsnSampleSize            = NA,
                                                                            dgvAsnObsGains              = NA,
                                                                            dgvAsnObsLosses             = NA,
                                                                            dgvAsnPctGains              = NA,
                                                                            dgvAsnPctLosses             = NA,
                                                                            dgvSanSampleSize            = NA,
                                                                            dgvSanObsGains              = NA,
                                                                            dgvSanObsLosses             = NA,
                                                                            dgvSanPctGains              = NA,
                                                                            dgvSanPctLosses             = NA,
                                                                            dgvGenes                    = hit$genes,
                                                                            chisqTestPvalue             = round(chisqTest$p.value, 2),
                                                                            chisqTestWhitePvalue        = round(chisqTestWhite$p.value, 2),
                                                                            chisqTestWhiteEurPvalue     = NA,
                                                                            fisherTestPvalue            = round(fisherTest$p.value, 2),
                                                                            fisherTestWhitePvalue       = round(fisherTestWhite$p.value, 2),
                                                                            fisherTestWhiteEurPvalue    = NA))
                            }
                        }
                    }
                }
            }
        }

        affectedPatientsPattern = sprintf("(%s)", paste(affectedPatients, collapse="|"))
        unaffectedPatientsPattern = sprintf("(%s)", paste(unaffectedPatients, collapse="|"))

        if (length(affectedPatients) > 2 & length(unaffectedPatients) > 2) {
            for (gene in genesInPeak) {
                gene = gsub("(\\[|\\])", "", gene, perl = TRUE)
                geneRnaSeqFoldChange = c()
                for (cancerType in cancerTypes) {
                    rnaLog2FoldChangeDf = rnaLog2FoldChangeLt[[cancerType]]
                    if (!is.null(rnaLog2FoldChangeDf)) {
                        rnaLog2FoldChangeGeneDf = rnaLog2FoldChangeDf[grep(sprintf("^%s\\|\\S+", gene, perl = TRUE), rownames(rnaLog2FoldChangeDf)),]
                        if (nrow(rnaLog2FoldChangeGeneDf) > 0) {
                            geneRnaSeqFoldChange = c(geneRnaSeqFoldChange, rnaLog2FoldChangeGeneDf)
                        }
                    }
                }
                geneRnaSeqFoldChangeDf = as.data.frame(geneRnaSeqFoldChange)
                if (nrow(geneRnaSeqFoldChangeDf) > 0) {
                    rnaSeqFoldChangeAffectedPatients = as.vector(unlist(geneRnaSeqFoldChangeDf[, grep(affectedPatientsPattern, colnames(geneRnaSeqFoldChangeDf), perl = TRUE)]))
                    rnaSeqFoldChangeUnaffectedPatients = as.vector(unlist(geneRnaSeqFoldChangeDf[, grep(unaffectedPatientsPattern, colnames(geneRnaSeqFoldChangeDf), perl = TRUE)]))
                    if (length(rnaSeqFoldChangeAffectedPatients) > 2 & length(rnaSeqFoldChangeUnaffectedPatients) > 2) {
                        tTestResultP = my.t.test.p.value(rnaSeqFoldChangeAffectedPatients, rnaSeqFoldChangeUnaffectedPatients)
                        rTestResultP = my.wilcox.test.p.value(rnaSeqFoldChangeAffectedPatients, rnaSeqFoldChangeUnaffectedPatients)
                        rsemValuesDf = data.frame(group = rep(peakType, length(rnaSeqFoldChangeAffectedPatients)), rsem = rnaSeqFoldChangeAffectedPatients)
                        rsemValuesDf = rbind(rsemValuesDf, data.frame(group = rep(sprintf("No_%s", peakType), length(rnaSeqFoldChangeUnaffectedPatients)), rsem = rnaSeqFoldChangeUnaffectedPatients))
                        p = ggplot(rsemValuesDf, aes(x = factor(group), y = rsem)) + geom_boxplot() + scale_x_discrete("") + scale_y_continuous("log2 fold change") +
                        ggtitle(sprintf("%s\n(t-test: %f, wilcox.test: %f)\n", gene, tTestResultP, rTestResultP))
                        f = file.path("/groups/park/semin/BiO/Research/GCC/Germline/Figures", sprintf("mRNAExpressionFoldChange_%s_%s_%s.pdf", peakBoundary, gene, peakType))
                        ggsave(p, file = f)
                    }
                }
            }
        }
    }
}


peakStatsFile = file.path("/groups/park/semin/BiO/Research/GCC/Germline/Tables", "PANCAN-Peak_Stats.tsv")
write.table(peakStatsDf, file=peakStatsFile, sep="\t", row.names=FALSE, quote=FALSE, col.names = c("eventType", "peakBoundary", "peakOverlapPct", "cytoband", "qvalue", "genesInPeak", "cancerType", "numAffectedPatients", "numUnaffectedPatients", "pctAffectedPatients", "numAffectedWhitePatients", "numUnaffectedWhitePatients", "pctAffectedWhitePatients", "dgvAccession", "dgvChr", "dgvPosStart", "dgvPosStop", "dgvOverlapPct", "dgvVarType", "dgvVarSubType", "dgvReference", "dgvPubmedId", "dgvSampleSize", "dgvObsGains", "dgvObsLosses", "dgvPctGains", "dgvPctLosses", "dgvEurSampleSize", "dgvEurObsGains", "dgvEurObsLosses", "dgvEurPctGains", "dgvEurPctLosses", "dgvAfrSampleSize", "dgvAfrObsGains", "dgvAfrObsLosses", "dgvAfrPctGains", "dgvAfrPctLosses", "dgvAmrSampleSize", "dgvAmrObsGains", "dgvAmrObsLosses", "dgvAmrPctGains", "dgvAmrPctLosses", "dgvAsnSampleSize", "dgvAsnObsGains", "dgvAsnObsLosses", "dgvAsnPctGains", "dgvAsnPctLosses", "dgvSanSampleSize", "dgvSanObsGains", "dgvSanObsLosses", "dgvSanPctGains", "dgvSanPctLosses", "dgvGenes", "chisqTestPvalue", "chisqTestWhitePvalue", "chisqTestWhiteEurPvalue", "fisherTestPvalue", "fisherTestWhitePvalue", "fisherTestWhiteEurPvalue"))

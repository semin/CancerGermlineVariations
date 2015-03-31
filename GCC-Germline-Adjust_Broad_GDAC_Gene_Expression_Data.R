rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(2)

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

imageFile = file.path(rdataDir, sprintf("GCC-Germline-Adjust_Broad_GDAC_Gene_Expression_Data-%s.RData", canType))
#load(imageFile)


##
## Read bam file list file to get TCGA patient IDs
##
bamListFile = file.path(germDir, sprintf("VCF/CANCERS/%s/BamFileListForGATK-%s-Normal.list", canType, canType))
bamListDf = read.delim(bamListFile, header = F, as.is = T)
colnames(bamListDf)[1] = "File"
bamListDf$Patient_ID = sapply(bamListDf$File, function(x) { paste(strsplit(basename(x), "-", fixed = T)[[1]][1:3], collapse = "_") })
snpPatientIds = unique(bamListDf$Patient_ID)


##
## Read predicted ancestry information
##
ansFile = file.path(tableDir, sprintf("Continental_Ancestry_Prediction_For_TCGA_Individuals_With_Shared_SNPs-%s.txt", canType))
ansDf = read.delim(ansFile, header = T, as.is = T)
ansEurDf = subset(ansDf, predicted_super_pop_id == "EUR")
eurSnpPatientIds = unique(ansEurDf$tcga_patient_id)


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
## Read gene-level somatic SCNA file
## : The copy number values in the table are in units of (copy number -2), 
##   so that no amplification or deletion is 0, genes with amplifications have positive values, and genes with deletions are negative values.
##
if (ncanType == "SKCM") {
    scnaFile = file.path(gdacAnalysesCanDir,
                         sprintf("gdac.broadinstitute.org_%s-TM.CopyNumber_Gistic2.Level_4.%s00.0.0/all_data_by_genes.txt", ncanType, gdacAnalysesDateNoUb))
} else {
    scnaFile = file.path(gdacAnalysesCanDir,
                         sprintf("gdac.broadinstitute.org_%s-TP.CopyNumber_Gistic2.Level_4.%s00.0.0/all_data_by_genes.txt", ncanType, gdacAnalysesDateNoUb))
}
scnaDf = read.delim(scnaFile, header = T, as.is = T)
rownames(scnaDf) = scnaDf$Gene.Symbol
scnaDf = scnaDf[,4:ncol(scnaDf)]
scnaDf = scnaDf + 2
colnames(scnaDf) = sapply(colnames(scnaDf), function(x) { paste(strsplit(x, ".", fixed = T)[[1]][1:3], collapse = "_") })
scnaPatientIds = colnames(scnaDf)


##
## Read methylation data for promoter regions
##
methylMergedFile = file.path(tcgaCanDir, sprintf("Methylation-Tumor-%s.new.txt", canType))
methylMergedDf = read.delim(methylMergedFile, header = T, as.is = T)
rownames(methylMergedDf) = methylMergedDf[, "Gene"]


##
## Read gene expression data
##
geneExpMergedFile = file.path(tcgaCanDir,  sprintf("Gene_Expression-Tumor-%s.new.txt", canType))
geneExpMergedDf = read.delim(geneExpMergedFile, header = T, row.names = 1, as.is = T)


##
## Patient IDs falling in the intersetion between data sets
##
sharedSnpExpEurPatients = mixedsort(intersect(eurSnpPatientIds, colnames(geneExpMergedDf)))
sharedExpScnaMethylPatients = mixedsort(intersect(intersect(colnames(geneExpMergedDf), colnames(scnaDf)), colnames(methylMergedDf)))
sharedSnpExpScnaMethylEurPatients = mixedsort(intersect(eurSnpPatientIds, sharedExpScnaMethylPatients))


##
## Detect genes with low/missing expression level in EUR patients
##
geneExpMergedEurDf = geneExpMergedDf[, sharedSnpExpEurPatients]
genesToRemove = c()
for (i in 1:nrow(geneExpMergedEurDf)) {
    gene = rownames(geneExpMergedEurDf)[i]
    avgExp = mean(as.numeric(geneExpMergedEurDf[i,]))
    perNa = 100 * sum(is.na(geneExpMergedEurDf[i,])) / ncol(geneExpMergedEurDf)
    if (avgExp == 0) {
        warning(sprintf("%s has %f average expression level (RSEM)!", gene, avgExp))
        genesToRemove = c(genesToRemove, gene)
    }
    if (perNa > 90) {
        warning(sprintf("%s has more than %f% missing expression data (RSEM)!", gene, perNa))
        genesToRemove = c(genesToRemove, gene)
    }
}


##
## Create gene location file for MatrixEQTL using GAF hg19 from NCI (https://tcga-data.nci.nih.gov/docs/GAF)
##
gafFile = file.path(germDir, "GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf")
geneLocFile = file.path(dirname(gafFile), gsub("gaf", "gene_locations.txt", basename(gafFile), fixed = T))

if (file.exists(geneLocFile)) {
    geneLocSrtDf = read.delim(geneLocFile, header = T, as.is = T)
} else {
    gafDf = read.delim(gafFile, header = T, as.is = T)
    gafDf$GeneLocus1 = sapply(gafDf$GeneLocus, function(x) { strsplit(x, ";", fixed = T)[[1]][1] })
    gafDf$chromosome = sapply(gafDf$GeneLocus1, function(x) { strsplit(x, ":", fixed = T)[[1]][1] })
    gafDf$chr_start = as.numeric(sapply(gafDf$GeneLocus1, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][1] }))
    gafDf$chr_stop = as.numeric(sapply(gafDf$GeneLocus1, function(x) { strsplit(strsplit(x, ":", fixed = T)[[1]][2], "-", fixed = T)[[1]][2] }))
    gafGeneDf = subset(gafDf, FeatureType == "gene")[, c("FeatureID", "Gene", "chromosome", "chr_start", "chr_stop")]

    geneLocDf = data.frame()
    geneLocDf <- foreach (i = 1:nrow(geneExpMergedDf), .combine=rbind) %dopar% {
        ngene = rownames(geneExpMergedDf)[i]
        ggene = geneExpMergedDf[i, "Gene"]
        cat(sprintf("%d:%s:%s\n", i, ngene, ggene))
        gafGeneTmpDf = subset(gafGeneDf, Gene == ggene)[1,]
        if (nrow(gafGeneTmpDf) == 0) {
            warning(sprintf("%d:%s:%s has no start and end positions!", i, ngene, ggene))
            data.frame()
        } else {
            data.frame(Symbol = ngene, gafGeneTmpDf[, c("Gene", "chromosome", "chr_start", "chr_stop")])
        }
    }

    geneLocSrtDf = geneLocDf[, c(1, 3:5)]
    colnames(geneLocSrtDf) = c("geneid", "chr", "s1", "s2")
    geneLocSrtDf = geneLocSrtDf[order(geneLocSrtDf$chr, geneLocSrtDf$s1, geneLocSrtDf$s2),]
    geneLocSrtDf$geneid = as.character(geneLocSrtDf$geneid)
    write.table(geneLocSrtDf, geneLocFile, sep = "\t", na = "NA", row.names = F, col.names = T, quote = F)
}


##
## Store unadjusted gene expression data for EUR patients
##
geneExpMergedDf2 = data.frame(id = rownames(geneExpMergedDf), geneExpMergedDf[, sharedSnpExpEurPatients])
geneExpMergedDf3 = geneExpMergedDf2[intersect(geneLocSrtDf$geneid, rownames(geneExpMergedDf)),]
geneExpMergedEurFile = file.path(tcgaCanDir, sprintf("Gene_Expression-Tumor-EUR-%s.new.txt", canType))
write.table(geneExpMergedDf3, geneExpMergedEurFile, sep = "\t", row.names = F, col.names = T, quote = F)

geneLog2ExpMergedDf = cbind(geneExpMergedDf[,1], log2(geneExpMergedDf[, 2:ncol(geneExpMergedDf)] + 1))
geneLog2ExpMergedDf2 = data.frame(id = rownames(geneLog2ExpMergedDf), geneLog2ExpMergedDf[, sharedSnpExpEurPatients])
geneLog2ExpMergedDf3 = geneLog2ExpMergedDf2[intersect(geneLocSrtDf$geneid, rownames(geneExpMergedDf)),]
geneLog2ExpMergedFile = file.path(tcgaCanDir, sprintf("Log2_Gene_Expression-Tumor-EUR-%s.new.txt", canType))
write.table(geneLog2ExpMergedDf3, geneLog2ExpMergedFile, sep = "\t", row.names = F, col.names = T, quote = F)

expGenes = rownames(geneExpMergedDf)

##
## Compute residual expressions using SCNA and methylation data (inspired by Li et al. (2013, Cell)
## Ti = Sci + Mi + εi, εi = Gi + ui
##
geneExpRsdDf <- foreach(i=1:length(expGenes), .combine=rbind) %dopar% {
    expGene = expGenes[i]
    cat(sprintf("%s (%d/%d)\n", expGene, i, length(expGenes)))
    geneExpTmpDf = geneExpMergedDf[expGene, sharedExpScnaMethylPatients]
    scnaTmpDf = scnaDf[expGene, sharedExpScnaMethylPatients]
    if (nrow(scnaTmpDf) == 0) {
        nexpGene = geneInfoExtDf[geneInfoExtDf$Synonym == expGene, "Symbol"]
        scnaTmpDf = scnaDf[rownames(scnaDf) == nexpGene, sharedExpScnaMethylPatients]
    }
    methylTmpDf = methylMergedDf[expGene, sharedExpScnaMethylPatients]
    if ((nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
        (nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        resDf = data.frame(id = expGene, geneExpTmpDf)
        return(resDf)
    } else if ((nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
               !(nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneExpTmpDf),
                                     beta_value = unlist(methylTmpDf))
        expScnaMethylFit = lm(rsem ~ beta_value, data = expScnaMethylDf)
    } else if (!(nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
               (nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneExpTmpDf),
                                     cn = unlist(scnaTmpDf))
        expScnaMethylFit = lm(rsem ~ cn, data = expScnaMethylDf)
    } else {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneExpTmpDf),
                                     cn = unlist(scnaTmpDf),
                                     beta_value = unlist(methylTmpDf))
        expScnaMethylFit = lm(rsem ~ cn + beta_value, data = expScnaMethylDf)
    }
    residualsDf = as.data.frame(t(residuals(expScnaMethylFit)[sharedExpScnaMethylPatients]))
    colnames(residualsDf) = sharedExpScnaMethylPatients
    data.frame(id = expGene, residualsDf)
}

geneExpRsdEurDf = cbind(geneExpRsdDf$id, geneExpRsdDf[, 2:ncol(geneExpRsdDf)][, sharedSnpExpScnaMethylEurPatients])
colnames(geneExpRsdEurDf)[1] = "id"
rownames(geneExpRsdEurDf) = geneExpRsdEurDf$id
geneExpRsdEurDf = geneExpRsdEurDf[setdiff(intersect(geneLocSrtDf$geneid, rownames(geneExpMergedDf)), genesToRemove),]
geneExpRsdFile = file.path(tcgaCanDir, sprintf("Residual_Gene_Expression-Tumor-EUR-%s.new.txt", canType))
write.table(geneExpRsdEurDf, geneExpRsdFile, sep = "\t", row.names = F, col.names = T, quote = F)


##
## Compute residual log2(expressions + 1) using SCNA and methylation data (inspired by Li et al. (2013, Cell)
## Ti = Sci + Mi + εi, εi = Gi + ui
##
geneLog2ExpRsdDf <- foreach(i=1:length(expGenes), .combine=rbind) %dopar% {
    expGene = expGenes[i]
    cat(sprintf("%s (%d/%d)\n", expGene, i, length(expGenes)))
    geneLog2ExpTmpDf = geneLog2ExpMergedDf[rownames(geneLog2ExpMergedDf) == expGene, sharedExpScnaMethylPatients]
    scnaTmpDf = scnaDf[rownames(scnaDf) == expGene, sharedExpScnaMethylPatients]
    if (nrow(scnaTmpDf) == 0) {
        nexpGene = geneInfoExtDf[geneInfoExtDf$Synonym == expGene, "Symbol"]
        scnaTmpDf = scnaDf[rownames(scnaDf) == nexpGene, sharedExpScnaMethylPatients]
    }
    methylTmpDf = methylMergedDf[rownames(methylMergedDf) == expGene, sharedExpScnaMethylPatients]
    if (nrow(methylTmpDf) == 0) {
        nexpGene = geneInfoExtDf[geneInfoExtDf$Synonym == expGene, "Symbol"]
        methylTmpDf = methylMergedDf[rownames(methylMergedDf) == nexpGene, sharedExpScnaMethylPatients]
    }
    if ((nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
        (nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        resDf = data.frame(id = expGene, geneLog2ExpTmpDf)
        return(resDf)
    } else if ((nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
               !(nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneLog2ExpTmpDf),
                                     beta_value = unlist(methylTmpDf))
        expScnaMethylFit = lm(rsem ~ beta_value, data = expScnaMethylDf)
    } else if (!(nrow(scnaTmpDf) == 0 | all(is.na(scnaTmpDf))) &
               (nrow(methylTmpDf) == 0 | all(is.na(methylTmpDf)))) {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneLog2ExpTmpDf),
                                     cn = unlist(scnaTmpDf))
        expScnaMethylFit = lm(rsem ~ cn, data = expScnaMethylDf)
    } else {
        expScnaMethylDf = data.frame(gene = expGene,
                                     rsem = unlist(geneLog2ExpTmpDf),
                                     cn = unlist(scnaTmpDf),
                                     beta_value = unlist(methylTmpDf))
        expScnaMethylFit = lm(rsem ~ cn + beta_value, data = expScnaMethylDf)
    }
    residualsDf = as.data.frame(t(residuals(expScnaMethylFit)[sharedExpScnaMethylPatients]))
    colnames(residualsDf) = sharedExpScnaMethylPatients
    data.frame(id = expGene, residualsDf)
}

geneLog2ExpRsdEurDf = cbind(geneLog2ExpRsdDf$id, geneLog2ExpRsdDf[, 2:ncol(geneLog2ExpRsdDf)][, sharedSnpExpScnaMethylEurPatients])
colnames(geneLog2ExpRsdEurDf)[1] = "id"
rownames(geneLog2ExpRsdEurDf) = geneLog2ExpRsdEurDf$id
geneLog2ExpRsdEurDf = geneLog2ExpRsdEurDf[setdiff(intersect(geneLocSrtDf$geneid, rownames(geneExpMergedDf)), genesToRemove),]
geneLog2ExpRsdFile = file.path(tcgaCanDir, sprintf("Residual_Log2_Gene_Expression-Tumor-EUR-%s.new.txt", canType))
write.table(geneLog2ExpRsdEurDf, geneLog2ExpRsdFile, sep = "\t", row.names = F, col.names = T, quote = F)

save.image(imageFile)


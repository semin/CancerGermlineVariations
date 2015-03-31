rm(list=ls())

require(reshape2)
require(sqldf)
require(scales)
require(ggplot2)
require(RSQLite)
require(vegan)

my.t.test.p.value <- function(...) { 
  obj<-try(t.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else if (is.nan(obj$p.value)) return(NA) else return(obj$p.value) 
} 

my.wilcox.test.p.value <- function(...) { 
  obj<-try(wilcox.test(...), silent=TRUE) 
  if (is(obj, "try-error")) return(NA) else if (is.nan(obj$p.value)) return(NA) else return(obj$p.value) 
} 


if (!is.na(pmatch("darwin", version$os))) {
    rootDir = "/Volumes/orchestra"
} else if (!is.na(pmatch("mingw32", version$os))) {
    rootDir = "Z:"
} else {
    rootDir = "/home/sl279"
}

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
figDir = file.path(germDir, "Figure/FoldChange")
tableDir = file.path(germDir, "Table")

canTypes = c("CRC")
canType = canTypes[1]

sqldf()

for (canType in canTypes) {
    if (canType == "CRC") {
        ncanType = "COADREAD"
    } else {
        ncanType = canType
    }
    canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  
    sqlite3Db = file.path(canVcfDir, sprintf("%s.db", canType))
    sqldf(sprintf("ATTACH '%s' AS %s;", sqlite3Db, canType))

    # Read MutSig SMGs
    sigGenesFile = file.path(germDir, sprintf("GDAC/gdac.broadinstitute.org_%s-TP.MutSigNozzleReport1.5.Level_4.2014011500.0.0/%s-TP.sig_genes.txt", ncanType, ncanType))
    sigGenesDf = read.delim(sigGenesFile, header = T, as.is = T)
    sigGenes = head(sigGenesDf$gene, 10)

    # Read MAF file
    mafFile = file.path(germDir, sprintf("GDAC/gdac.broadinstitute.org_%s-TP.MutSigNozzleReport1.5.Level_4.2014011500.0.0/%s-TP.final_analysis_set.maf", ncanType, ncanType))
    mafDf = read.delim(mafFile, header = T, as.is = T)
    nrow(mafDf)
    unique(mafDf$Variant_Classification)
    mafNonsilSigDf = subset(mafDf, gene %in% sigGenes & Variant_Classification != "Silent")
    nrow(mafNonsilSigDf)

    mafPatientIds = unique(sapply(mafNonsilSigDf$Tumor_Sample_Barcode, function(x) { paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "_") }))

    ## Over-represented known SNPs
    #orSnps1Df = sqldf(sprintf("SELECT *
                            #FROM %s.snp_ext_annotations 
                            #WHERE (chrom != 'X' AND chrom != 'Y') AND
                                #(id != '.' AND rsid1 != '' AND all_af1 > tg_all_af1 AND tg_all_af1 < 0.05 AND tg_eur_af1 < 0.05) AND
                                #(fisher_adj_pval_all_ac1 < 0.01 AND fisher_adj_pval_all_ac1 != '') AND
                                #(chisq_adj_pval_all_ac1 < 0.01 AND chisq_adj_pval_all_ac1 != '')", canType), dbname = sqlite3Db)

    # Novel SNPs with AF > 0.05
    novel05Snps1Df = sqldf(sprintf("SELECT *
                                   FROM %s.snp_ext_annotations 
                                   WHERE (chrom != 'X' AND chrom != 'Y') AND
                                   (id == '.' AND rsid1 == '' AND tg_all_af1 == '' AND all_af1 >= 0.05)", canType), dbname = sqlite3Db)

    novel05SnpsSampleIds = grep("TCGA", colnames(novel05Snps1Df), value = T)
    novel05SnpsPatientIds = as.vector(unlist(sapply(novel05SnpsSampleIds, function(x) { paste(strsplit(x, "_", fixed = T)[[1]][1:3], collapse = "_") })))
    
    # Patient IDs shared by MAF and germline calls
    sharedPatientIds = intersect(mafPatientIds, novel05SnpsPatientIds)

    # Construct binary matrices for MAF and germline calls based on shared patient IDs
    mafNonsilSigDf$patientId = sapply(mafNonsilSigDf$Tumor_Sample_Barcode, function(x) { paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "_") })
    mafNonsilSigGrpDf = sqldf('SELECT Hugo_Symbol AS gene, patientId, COUNT(*) AS mutCnt FROM mafNonsilSigDf GROUP BY patientId')
    mafNonsilSigGrpCastDf = dcast(mafNonsilSigGrpDf, "gene ~ patientId", value.var = "mutCnt")
    rownames(mafNonsilSigGrpCastDf) = mafNonsilSigGrpCastDf$gene
    mafNonsilSigGrpCastDf = mafNonsilSigGrpCastDf[,-1]
    mafNonsilSigGrpCastDf[is.na(mafNonsilSigGrpCastDf)] = 0
    mafNonsilSigGrpCastDf[mafNonsilSigGrpCastDf > 0] = 1 
    head(mafNonsilSigGrpCastDf, 2)
    ncol(mafNonsilSigGrpCastDf)
    mafNonsilSigGrpCastSharedDf = mafNonsilSigGrpCastDf[,sharedPatientIds]
    head(mafNonsilSigGrpCastSharedDf)
    ncol(mafNonsilSigGrpCastSharedDf)
    nrow(mafNonsilSigGrpCastSharedDf)

    novel05Snps1GtsDf = novel05Snps1Df[, grep("TCGA", colnames(novel05Snps1Df))]
    ncol(novel05Snps1GtsDf)
    rownames(novel05Snps1GtsDf) = gsub(" ", "", apply(novel05Snps1Df[, c("chrom", "pos", "ref", "alt")], 1, paste, collapse = "_"))
    colnames(novel05Snps1GtsDf) = as.vector(unlist(sapply(colnames(novel05Snps1GtsDf), function(x) { paste(strsplit(x, "_", fixed = T)[[1]][1:3], collapse = "_") })))
    novel05Snps1GtsSharedDf = novel05Snps1GtsDf[, sharedPatientIds]
    ncol(novel05Snps1GtsSharedDf)
    novel05Snps1GtsSharedDf[novel05Snps1GtsSharedDf == "'0/0'"] = 0
    novel05Snps1GtsSharedDf[novel05Snps1GtsSharedDf == "'./.'"] = 0
    novel05Snps1GtsSharedDf[novel05Snps1GtsSharedDf != 0] = 1
    novel05Snps1GtsSharedMat = data.matrix(novel05Snps1GtsSharedDf)
    rownames(novel05Snps1GtsSharedMat) = rownames(novel05Snps1GtsSharedDf)
    novel05Snps1GtsSharedDf = as.data.frame(novel05Snps1GtsSharedMat)
    nrow(novel05Snps1GtsSharedDf)

    rSums = rowSums(novel05Snps1GtsSharedDf)
    novel05Snps1GtsSharedSubDf = novel05Snps1GtsSharedDf[rSums > 8,]


    mergedVarsDf = rbind(mafNonsilSigGrpCastSharedDf, novel05Snps1GtsSharedSubDf)
    jcDist = vegdist(mergedVarsDf, method="jaccard")

    for (j in 1:nrow(mafNonsilSigGrpCastSharedDf)) {
        gene = rownames(mafNonsilSigGrpCastSharedDf)[j]
        print(sprintf("Finding germline mutations mutually exclusive/co-occurring with %s ...", gene))
        for (k in 1:nrow(novel05Snps1GtsSharedSubDf)) {
            mat = rbind(mafNonsilSigGrpCastSharedDf[j,], novel05Snps1GtsSharedSubDf[k,])
            jc = vegdist(mat, method="jaccard")
        }
    }

    miBetweenSomaticAndGermlineFile = file.path(canVcfDir, "miBetweenSomaticAndGermline.txt")
    write.table(miBetweenSomaticAndGermlineDf, miBetweenSomaticAndGermlineFile, sep = "\t", quote = FALSE, row.names = FALSE)
}

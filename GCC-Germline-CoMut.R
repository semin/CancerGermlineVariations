rm(list=ls())
options(width=250)

# require(vegan)
require(reshape2)
require(sqldf)
require(scales)
require(ggplot2)

if (!is.na(pmatch("darwin", version$os))) {
#   rootDir = "/Volumes/transfer.orchestra.med.harvard.edu"
   rootDir = "/Volumes/orchestra"
} else if (!is.na(pmatch("mingw32", version$os))) {
  rootDir = "Z:"
} else {
  rootDir = "~"
}

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
canTypes = c("HNSC")
canType = canTypes[1]

for (canType in canTypes) {
  # Collpase individual mutations (SNVs and INDELs) into genes
  if (canType == "HNSC") {
    mafFile = file.path(germDir, "MAF/HNSC_pair_set_279_freeze_09252012.final_analysis_set.maf")
    mafDf = read.delim(mafFile, header = TRUE, as.is = TRUE)
    mafDf$patientId = sapply(mafDf$Tumor_Sample_Barcode, function(x) { paste(strsplit(x, "-", fixed = TRUE)[[1]][2:4], collapse = "_") })
  } else {
    mafFile = file.path(germDir, sprintf("MAF/%s_cleaned.maf", tolower(canType)))
    mafDf = read.delim(mafFile, header = TRUE, as.is = TRUE)
    mafDf$patientId = sapply(mafDf$Tumor_Sample_Barcode, function(x) { paste(strsplit(x, "-", fixed = TRUE)[[1]][1:3], collapse = "_") })
  }
  #   unique(mafDf$Variant_Classification)
  mafNosilDf = subset(mafDf, Variant_Classification != "RNA", Variant_Classification != "Silent")
  mafNosilGrpDf = sqldf('SELECT Hugo_Symbol AS gene, patientId, COUNT(*) AS cnt FROM mafNosilDf GROUP BY patientId, Hugo_Symbol')
  mafNosilCastDf = dcast(mafNosilGrpDf, "gene ~ patientId", value.var = "cnt")
  rownames(mafNosilCastDf) = mafNosilCastDf$gene
  mafNosilCastDf = mafNosilCastDf[,-1]
  
  # Collect gene-level somatic CNA information
  geneScnaFile = file.path(germDir, "GDAC", sprintf("gdac.broadinstitute.org_%s-TP.CopyNumber_Gistic2.Level_4.2013092300.0.0/all_data_by_genes.txt", toupper(canType)))
  geneScnaDf = read.delim(geneScnaFile, header = TRUE, as.is = TRUE)[, -c(2,3)]
  colnames(geneScnaDf) = c("gene", as.vector(unlist(sapply(colnames(geneScnaDf)[2:ncol(geneScnaDf)],
                                             function(x) { paste(strsplit(x, ".", fixed = TRUE)[[1]][1:3], collapse = "_") }))))
  geneScnaMeltDf = melt(geneScnaDf, id=c("gene"))
  colnames(geneScnaMeltDf)[2] = "patientId"
  colnames(geneScnaMeltDf)[3] = "cn"
  geneScnaCastDf = dcast(geneScnaMeltDf, "gene ~ patientId", value.var = "cn")
  rownames(geneScnaCastDf) = geneScnaCastDf$gene
  geneScnaCastDf = geneScnaCastDf[,-1]
  
  # Load MutSig results (for driver genes)
  mutSigGenesFile = file.path(germDir,
                              sprintf("GDAC/gdac.broadinstitute.org_%s-TP.MutSigNozzleReport1.5.Level_4.2013092300.0.0/%s-TP.sig_genes.txt", toupper(canType), toupper(canType)))
  mutSigGenesDf = read.delim(mutSigGenesFile, header = TRUE, as.is = TRUE)
  mutSigTop10Genes = head(mutSigGenesDf$gene, 10)

  # Collect germline variations
  germSnpSplicingVarFile = file.path(germDir, 
                                  sprintf("VCF/CANCERS/%s/Harvard_GCC_WGS-%s-Normal.vqsr.snp.common.maf5.annovar.splicing_variant_function", canType, canType))
  germSnpSplicingVarDf = read.delim(germSnpSplicingVarFile, header = TRUE, as.is = TRUE)
  germIndelSplicingVarFile = file.path(germDir, 
                                     sprintf("VCF/CANCERS/%s/Harvard_GCC_WGS-%s-Normal.vqsr.indel.common.maf5.annovar.splicing_variant_function", canType, canType))
  germIndelSplicingVarDf = read.delim(germIndelSplicingVarFile, header = TRUE, as.is = TRUE)
  germSplicingVarDf = rbind(germSnpSplicingVarDf, germIndelSplicingVarDf)
  colnames(germSplicingVarDf)[1] = "TYPE"
  colnames(germSplicingVarDf)[2] = "FEATURE"
  colnames(germSplicingVarDf)[3] = "CHROM"
  germSplicingVarDf$GENE = sapply(germSplicingVarDf$FEATURE, function(x) { gsub("\\s", "", strsplit(x, "(", fixed = TRUE)[[1]][1]) })
  
  germSnpExonicVarFile = file.path(germDir, 
                                sprintf("VCF/CANCERS/%s/Harvard_GCC_WGS-%s-Normal.vqsr.snp.common.maf5.annovar.exonic_variant_function", canType, canType))
  germSnpExonicVarDf = read.delim(germSnpExonicVarFile, header = FALSE, as.is = TRUE)
  germIndelExonicVarFile = file.path(germDir, 
                                   sprintf("VCF/CANCERS/%s/Harvard_GCC_WGS-%s-Normal.vqsr.indel.common.maf5.annovar.exonic_variant_function", canType, canType))
  germIndelExonicVarDf = read.delim(germIndelExonicVarFile, header = FALSE, as.is = TRUE)
  germExonicVarDf = rbind(germSnpExonicVarDf, germIndelExonicVarDf)
  germExonicVarDf = germExonicVarDf[,-1]
  germExonicVarDf$GENE = sapply(germExonicVarDf[,2], function(x) { gsub("\\s", "", strsplit(x, ":", fixed = TRUE)[[1]][1]) })
  
  colnames(germExonicVarDf) = colnames(germSplicingVarDf)
  unique(germExonicVarDf$TYPE)
#   germExonicNosilVarDf = subset(germExonicVarDf, TYPE != "unknown")

  germGenicVarDf = rbind(germExonicVarDf, germSplicingVarDf)
  germGenicVarDf = germGenicVarDf[,c(ncol(germGenicVarDf), 1:(ncol(germGenicVarDf)-1))]

#   germGenicVarTp53Df = subset(germGenicVarDf, GENE == "TP53")
#   germGenicVarTp53Df.2 = germGenicVarTp53Df[, c(9:ncol(germGenicVarTp53Df))]
#   germGenicVarTp53Df.2[germGenicVarTp53Df.2 == "."] = 0
#   germGenicVarTp53Mat = data.matrix(germGenicVarTp53Df.2)

  # Load HPV information
  if (canType == "HNSC") {
    hpvFile = file.path(germDir, "Pathogen", "HNSC_GDAC_FREEZE_Mar272013_load_file_HPVpos35_sampleIDs.txt")
    hpvDf = read.delim(hpvFile, header = TRUE, as.is = TRUE)
    hpvDf$patientId = sapply(hpvDf$sample_id, function(x) { paste(strsplit(x, "-", fixed = TRUE)[[1]][1:3], collapse = "_") })
  }


  

}

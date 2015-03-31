rm(list=ls())

require(sqldf)
require(RSQLite)

if (!is.na(pmatch("darwin", version$os))) {
    rootDir = "/Volumes/orchestra"
} else if (!is.na(pmatch("mingw32", version$os))) {
    rootDir = "Z:"
} else {
    rootDir = "/home/sl279"
}

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
#canTypes = c("LUAD")
canTypes = c("CRC")
canType = canTypes[1]

sqldf()

for (canType in canTypes) {
    canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  
    sqlite3Db = file.path(canVcfDir, sprintf("%s.db", canType))
    sqldf(sprintf("ATTACH '%s' AS %s;", sqlite3Db, canType))

    # Calculate FDRs and create a separate table
    colNames = sqldf("pragma table_info(snp_annotations)", dbname = sqlite3Db)$name
    pvalueCols = grep("adj", grep("pval", colNames, value=T), invert=T, value=T)

    # Fetch known rare SNPs with MAF < 0.05 (AF1 < 0.05 or AF > 0.95)
    snps05Df = sqldf(sprintf("SELECT *
                             FROM %s.snp_annotations 
                             WHERE  (chrom != 'X' AND chrom != 'Y') AND
                             (id != '.' OR rsid1 != '' OR dbsnp137 IS NOT NULL OR dbsnp138 IS NOT NULL) AND 
                             ((tg_all_af1 < 0.05 AND tg_asn_af1 < 0.05 AND tg_amr_af1 < 0.05 AND tg_afr_af1 < 0.05 AND tg_eur_af1 < 0.05) OR
                              (tg_all_af1 > 0.95 AND tg_asn_af1 > 0.95 AND tg_amr_af1 > 0.95 AND tg_afr_af1 > 0.95 AND tg_eur_af1 > 0.95))", canType), dbname = sqlite3Db)


    for (pvalueCol in pvalueCols) {
        adjPvalueCol = gsub("pval", "adj_pval", pvalueCol)
        snps05Df[, adjPvalueCol] = p.adjust(snps05Df[, pvalueCol], "fdr")
    }

    sqldf(sprintf("CREATE TABLE %s.rare_snps AS SELECT * FROM snps05Df ", canType), dbname = sqlite3Db)

    #pvaluesDf[is.na(pvaluesDf)] = ""

    #fdrFile = file.path(canVcfDir, "fdrs.tsv")
    #write.table(pvaluesDf, fdrFile, sep = "\t", quote = FALSE, row.names = FALSE)
    #read.csv.sql(fdrFile, sql = sprintf("CREATE TABLE %s.snp_fdrs AS SELECT * FROM file", canType), dbname = sqlite3Db, sep = "\t")

    #for (pvalueCol in pvalueCols) {
        #print(pvalueCol)
        #pvaluesDf = sqldf(sprintf("SELECT chrom, pos, ref, alt, %s FROM snp_annotations", pvalueCols), dbname = sqlite3Db)
        #pvaluesDf[pvaluesDf == ""] = NA
        #adjPvalueCol = gsub("pval", "adj_pval", pvalueCol)
        #pvaluesDf[, adjPvalueCol] = p.adjust(pvaluesDf[, pvalueCol], "fdr")
        #fdrFile = file.path(canVcfDir, "fdrs.tsv")
        #write.table(pvaluesDf, fdrFile, sep = "\t", quote = FALSE, row.names = FALSE)
    #}

    # Create a new table containing adjusted pvalues
    #sqldf(sprintf("CREATE TABLE %s.snp_ext_annotations AS
                  #SELECT a.*,
                  #p.chisq_adj_pval_all_ac1,
                  #p.chisq_adj_pval_all_ac2,
                  #p.chisq_adj_pval_all_ac3,
                  #p.fisher_adj_pval_all_ac1,
                  #p.fisher_adj_pval_all_ac2,
                  #p.fisher_adj_pval_all_ac3,
                  #p.chisq_adj_pval_eur_ac1,
                  #p.chisq_adj_pval_eur_ac2,
                  #p.chisq_adj_pval_eur_ac3,
                  #p.fisher_adj_pval_eur_ac1,
                  #p.fisher_adj_pval_eur_ac2,
                  #p.fisher_adj_pval_eur_ac3
                  #FROM %s.snp_annotations AS a
                  #LEFT JOIN %s.snp_fdrs AS p
                  #ON a.chrom = p.chrom AND a.pos = p.pos AND a.ref = p.ref AND a.alt = p.alt", canType, canType, canType),
                  #dbname = sqlite3Db)
}

rm(list=ls())

##
## Load libraries
##
require(snpStats)
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
args = commandArgs(TRUE)
canType = args[1]


##
## Read VCFs and create genotypes file for MatrixEQTL
##
.integerGt = function(x) {
    if (is.character(x)) {
        if (nchar(x) == 3) {
            alleles = unlist(strsplit(x, "[\\||/]", perl=T))
            return(ifelse(alleles[1] == "." | alleles[2] == ".",  NA, as.character(sum(as.integer(alleles)))))
        } else if (nchar(x) == 1) {
            return(ifelse(x == ".", NA, x))
        }
    } else {
        return(NA)
    }
}
.integerGtVec <- Vectorize(.integerGt, "x")


for (i in 1:length(chrs)) {
    chr = chrs[i]
    vcfFiles = Sys.glob(file.path(germDir, sprintf("VCF/CANCERS/%s/*/chrs/%s/Harvard_GCC_WGS-%s-Normal.dbsnp.*.vqsr.pass.an100.chr%s.*eur_af.biallelic.nodups.vcf.bgz", canType, chr, canType, chr)))
    for (vcfFile in vcfFiles) {
        vcfDir = dirname(vcfFile)
        cat(sprintf("Integer %s ...\n", vcfFile))
        varType = ifelse(grepl("indel", basename(vcfFile)), "indel", "snp")
        vcf = readVcf(vcfFile, "hg19")
        if (varType == "indel") {
            vcfGtMat = geno(vcf)$GT
            for (j in 1:ncol(vcfGtMat)) { 
                #cat(sprintf("Integer %s (%d/%d) genotypes ...\n", colnames(vcfGtMat)[j], j, ncol(vcfGtMat)))
                vcfGtMat[,j] = .integerGtVec(vcfGtMat[,j]) 
            }
            snpNumDf = as.data.frame(vcfGtMat)
        } else {
            snpMat = genotypeToSnpMatrix(vcf)
            snpNumDf = as.data.frame(t(as(snpMat$genotypes, "numeric")))
        }
        colnames(snpNumDf) = sapply(colnames(snpNumDf), function(x) { paste(strsplit(x, "-", fixed = T)[[1]][1:3], collapse = "_") })
        eqtlGtEurFile = file.path(vcfDir, gsub("vcf.bgz", "genotypes.txt", basename(vcfFile)))
        write.table(snpNumDf, eqtlGtEurFile, sep = "\t", row.names = T, col.names = T, quote = F)
        snpRd = rowData(vcf)
        eqtlGtLocFile = file.path(vcfDir, gsub("vcf.bgz", "genotype_locations.txt", basename(vcfFile)))
        write.table(data.frame(snp = names(snpRd), chr = paste("chr", seqnames(snpRd), sep = ""), pos = start(snpRd)),
                    eqtlGtLocFile, sep = "\t", row.names = F, col.names = T, quote = F)
    }
}


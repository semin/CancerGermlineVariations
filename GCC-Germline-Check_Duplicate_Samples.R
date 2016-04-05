rm(list=ls())

##
## Load libraries
##
require(doMC)
registerDoMC(4)
require(gtools)
require(VariantAnnotation)

##
## Initialize global variables
##
if (!is.na(pmatch("darwin", version$os))) {
    rootDir = "/Volumes/orchestra"
} else {
    rootDir = "/home/sl279"
}

.jaccardIdx = function(g1, g2) {
    if (length(g1) != length(g2)) {
        return(0)
    } else {
        cntEq = sum(g1 == g2)
        return(cntEq / length(g1))
    }
}

germDir = file.path("/groups/kucherlapati/GCC/Germline")
chr22VcfFile = file.path(germDir, "VCF/PANCAN/snp/chrs/22/Harvard_GCC-WGS-PANCAN.snp.vqsr.pass.an100.chr22.maf05.vcf.gz")
chr22Vcf = readVcf(chr22VcfFile, "hg19")
sampleIds = samples(header(chr22Vcf))
genotypes = geno(chr22Vcf)$GT[1:10000,]
nrow(genotypes)
ncol(genotypes)

genotypeSimsDf <- foreach(i=1:ncol(genotypes), .combine=rbind) %dopar% {
    print(i)
    genotypeSimsTmpDf = data.frame()
    for (j in i:ncol(genotypes)) {
        genotypeSimsTmpDf = rbind(genotypeSimsTmpDf,
                                  data.frame(index1 = i, 
                                             index2 = j,
                                             sample1 = sampleIds[i],
                                             sample2 = sampleIds[j],
                                             jaccard_similarity = .jaccardIdx(genotypes[,i], genotypes[,j])))
    }
    genotypeSimsTmpDf
}

length(sort(unique(genotypeSimsDf$index1)))

hist(genotypeSimsDf$jaccard_similarity)
genotypeSimsIdnDf = subset(genotypeSimsDf, jaccard_similarity == 1)
genotypeSimsIdnDf2 = subset(genotypeSimsIdnDf, index1 != index2)
genotypeSimsIdnDf2

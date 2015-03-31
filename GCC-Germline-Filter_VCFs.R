rm(list=ls())

##
## Load libraries
##
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

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
args = commandArgs(TRUE)
canType = args[1]


##
## Read VCFs and create genotypes file for prioritization
##
vcfLst = list()
vcfFiles = mixedsort(Sys.glob(file.path(germDir, sprintf("VCF/CANCERS/%s/*/chrs/*/*eur_af.vcf.gz", canType, canType))))

for (vcfFile in vcfFiles) {
    cat(sprintf("Reading %s ...\n", vcfFile))
    varType = ifelse(grepl("indel", basename(vcfFile)), "indel", "snp")
    chr = gsub(".*(chr\\S+?).*", "\\1", basename(vcfFile))
    vcf = readVcf(vcfFile, "hg19")

    # Multi-allelic, segmental duplication, simple repeats, and miscrosatelites overlapping variants were filtered out
    vcfFiltered = vcf[elementLengths(alt(vcf)) == 1 &
                      elementLengths(info(vcf)$genomicSuperDups) == 0 &
                      elementLengths(info(vcf)$simpleRepeat) == 0 &
                      elementLengths(info(vcf)$microsat) == 0]
    vcfFileterdFile = gsub("eur_af", "eur_af.biallelic.nodups", vcfFile, fixed = T)
    vcfFileterdFile = gsub(".gz", "", vcfFileterdFile, fixed = T)
    #cat(sprintf("Writing %s ...\n", vcfFileterdFile))
    #writeVcf(vcfFiltered, vcfFileterdFile, index = T)

    # Rare
    rareVcf = vcfFiltered[unlist(!is.na(info(vcfFiltered)$TG_ALL_AF) &
                                 info(vcfFiltered)$TG_ALL_AF < 0.01 &
                                 info(vcfFiltered)$TG_EAS_AF < 0.01 &
                                 info(vcfFiltered)$TG_EUR_AF < 0.01 &
                                 info(vcfFiltered)$TG_AFR_AF < 0.01 &
                                 info(vcfFiltered)$TG_AMR_AF < 0.01 &
                                 info(vcfFiltered)$TG_SAS_AF < 0.01 &
                                 (is.na(info(vcfFiltered)$TA_AF) | info(vcfFiltered)$TA_AF < 0.01) &
                                 (is.na(info(vcfFiltered)$EA_AF) | info(vcfFiltered)$EA_AF < 0.01) &
                                 (is.na(info(vcfFiltered)$AA_AF) | info(vcfFiltered)$AA_AF < 0.01) &
                                 (info(vcfFiltered)$EUR_AF >= 0.05 & info(vcfFiltered)$EUR_AF <= 0.95)
                                 )]

    rareVcfFile = gsub("nodups", "nodups.rare", vcfFileterdFile, fixed = T)
    cat(sprintf("Writing %s ...\n", rareVcfFile))
    writeVcf(rareVcf, rareVcfFile, index = T)

    # Novel
    novelVcf = vcfFiltered[unlist(is.na(info(vcfFiltered)$TG_ALL_AF) &
                                  is.na(info(vcfFiltered)$TA_AF) &
                                  (info(vcfFiltered)$EUR_AF >= 0.05 & info(vcfFiltered)$EUR_AF <= 0.95) &
                                  info(vcfFiltered)$DBSNP == "." &
                                  elementLengths(info(vcfFiltered)$genomicSuperDups) == 0 &
                                  elementLengths(info(vcfFiltered)$simpleRepeat) == 0 &
                                  elementLengths(info(vcfFiltered)$microsat) == 0)]

    novelVcfFile = gsub("nodups", "nodups.novel", vcfFileterdFile, fixed = T)
    cat(sprintf("Writing %s ...\n", novelVcfFile))
    writeVcf(novelVcf, novelVcfFile, index = T)

    # Common
    commonVcf = vcfFiltered[unlist((!is.na(info(vcfFiltered)$TG_ALL_AF) &
                                    (info(vcfFiltered)$TG_ALL_AF >= 0.01 |
                                     info(vcfFiltered)$TG_EAS_AF >= 0.01 |
                                     info(vcfFiltered)$TG_EUR_AF >= 0.01 |
                                     info(vcfFiltered)$TG_AFR_AF >= 0.01 |
                                     info(vcfFiltered)$TG_AMR_AF >= 0.01 |
                                     info(vcfFiltered)$TG_SAS_AF >= 0.01)) |
                                   ((!is.na(info(vcfFiltered)$TA_AF) & info(vcfFiltered)$TA_AF >= 0.01) |
                                    (!is.na(info(vcfFiltered)$EA_AF) & info(vcfFiltered)$EA_AF >= 0.01) |
                                    (!is.na(info(vcfFiltered)$AA_AF) & info(vcfFiltered)$AA_AF >= 0.01))
                                   )]

    commonVcfFile = gsub("nodups", "nodups.common", vcfFileterdFile, fixed = T)
    cat(sprintf("Writing %s ...\n", commonVcfFile))
    writeVcf(commonVcf, commonVcfFile, index = T)
}


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

germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
gdacDir = file.path(germDir, "GDAC")
tcgaDir = file.path(germDir, "TCGA")
eqtlDir = file.path(germDir, "EQTL")
tableDir = file.path(germDir, "Table")
rdataDir = file.path(germDir, "RData")

gdacStdDataDate = "2015_02_04"
gdacStdDataDateNoUb = gsub("_", "", gdacStdDataDate)
gdacStdDataDir = file.path(gdacDir, sprintf("stddata__%s", gdacStdDataDate))

gdacAnalysesDate = "2014_10_17"
gdacAnalysesDateNoUb = gsub("_", "", gdacAnalysesDate)
gdacAnalysesDir = file.path(gdacDir, sprintf("analyses__%s", gdacAnalysesDate))

args = commandArgs(TRUE)
canType = args[1]
#canType = "BLCA"

tcgaCanDir = file.path(tcgaDir, canType)
dir.create(tcgaCanDir, showWarnings = F, recursive = T)

eqtlCanDir = file.path(eqtlDir, canType)
dir.create(eqtlCanDir, showWarnings = F, recursive = T)

imageFile = file.path(rdataDir, sprintf("GCC-Germline-Run_MatrixEQTL-%s.RData", canType))
#load(imageFile)
#save.image(imageFile)


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
## Read TCGA Broad GDAC data
##
ncanType = ifelse(canType == "CRC", "COADREAD", canType)
gdacStdDataCanDir = file.path(gdacStdDataDir, ncanType, gdacStdDataDateNoUb)
gdacAnalysesCanDir = file.path(gdacAnalysesDir, ncanType, gdacAnalysesDateNoUb)
canVcfDir = file.path(germDir, "VCF/CANCERS", canType)


##
## Create gene location file for MatrixEQTL using GAF hg19 from NCI (https://tcga-data.nci.nih.gov/docs/GAF)
##
gafFile = file.path(germDir, "GAF/GAF.hg19.June2011.bundle/outputs/TCGA.hg19.June2011.gaf")
eqtlGeneLocFile = file.path(dirname(gafFile), gsub("gaf", "gene_locations.txt", basename(gafFile), fixed = T))


##
## Load residual log transformed gene expression data
##
#geneExpRsdEurFile = file.path(tcgaCanDir, sprintf("Residual_Gene_Expression-Tumor-EUR-%s.txt", canType))
#geneExpRsdEurDf = read.delim(geneExpRsdEurFile, header = T, as.is = T)
#sharedSnpExpScnaMethylEurPatients = colnames(geneExpRsdEurDf)[2:ncol(geneExpRsdEurDf)]
#rm(geneExpRsdEurDf)

geneLog2ExpRsdEurFile = file.path(tcgaCanDir, sprintf("Residual_Log2_Gene_Expression-Tumor-EUR-%s.new.txt", canType))
geneLog2ExpRsdEurDf = read.delim(geneLog2ExpRsdEurFile, header = T, as.is = T)
sharedSnpExpScnaMethylEurPatients = colnames(geneLog2ExpRsdEurDf)[2:ncol(geneLog2ExpRsdEurDf)]
rm(geneLog2ExpRsdEurDf)


##
## Create genotypes file for MatrixEQTL
##
genoDf = data.frame()
genoChrFiles = mixedsort(Sys.glob(file.path(germDir, sprintf("VCF/CANCERS/%s/*/chrs/*/*.genotypes.txt", canType, canType))))
for (i in 1:length(genoChrFiles)) {
    genoChrFile = genoChrFiles[i]
    print(genoChrFile)
    genoChrDf = read.delim(genoChrFile, header = T, as.is = F)
    genoDf = rbind(genoDf, genoChrDf)
}

genoEurDf = data.frame(id = rownames(genoDf), genoDf[, sharedSnpExpScnaMethylEurPatients])
genoEurFile = file.path(eqtlCanDir, "Genotypes.txt")
write.table(genoEurDf, genoEurFile, sep = "\t", row.names = F, col.names = T, quote = F)

rm(genoChrDf)
rm(genoDf)
rm(genoEurDf)


varLocDf = data.frame()
varLocChrFiles = mixedsort(Sys.glob(file.path(germDir, sprintf("VCF/CANCERS/%s/*/chrs/*/*.genotype_locations.txt", canType, canType))))
for (i in 1:length(varLocChrFiles)) {
    varLocChrFile = varLocChrFiles[i]
    print(varLocChrFile)
    varLocChrDf = read.delim(varLocChrFile, header = T, as.is = F)
    varLocDf = rbind(varLocDf, varLocChrDf)
}

varLocFile = file.path(eqtlCanDir, "Genotype_Locations.txt")
write.table(varLocDf, varLocFile, sep = "\t", row.names = F, col.names = T, quote = F)

rm(varLocChrDf)
rm(varLocDf)


##
## Read clinical data to create covariates file for MatrixEQTL
##
clinFields = c("Hybridization REF", "yearstobirth", "gender")
clinFile = file.path(gdacStdDataCanDir,
                     sprintf("gdac.broadinstitute.org_%s.Clinical_Pick_Tier1.Level_4.%s00.0.0", ncanType, gdacStdDataDateNoUb),
                     sprintf("%s.clin.merged.picked.txt", ncanType))
clinTmpDf = read.delim(clinFile, header = F, as.is = T)
clinTmpDf = clinTmpDf[clinTmpDf[,1] %in% clinFields,]
colnames(clinTmpDf) = gsub("-", "_", toupper(clinTmpDf[1,]), fixed = T)
colnames(clinTmpDf)[1] = "id"
clinTmpDf = clinTmpDf[-1,]
clinTmpDf[1, "id"] = "age"
clinTmpDf[2, "id"] = "gender"
clinTmpDf[2, 2:ncol(clinTmpDf)] = ifelse(clinTmpDf[2, 2:ncol(clinTmpDf)] == "female", 0, 1)
clinDf = clinTmpDf[, c("id", sharedSnpExpScnaMethylEurPatients)]

# If all the samples are the same gender then do not include gener information in the covariates.txt
genders = clinDf[2, 2:ncol(clinDf)]
genders = as.integer(genders[!is.na(genders)])
if (abs(max(genders) - min(genders)) < (.Machine$double.eps ^ 0.5)) {
    clinDf = clinDf[-2,]
}
covariatesFile = file.path(eqtlCanDir, "Covariates.txt")
write.table(clinDf, covariatesFile, sep = "\t", row.names = F, col.names = T, quote = F)


##
## eQTL analysis
##
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
output_file_name_cis = file.path(eqtlDir, canType, "cis-eQTLs-Residual_Log2_Gene_Expression.txt")
output_file_name_tra = file.path(eqtlDir, canType, "trans-eQTLs-Residual_Log2_Gene_Expression.txt")
#output_file_name_cis = file.path(eqtlDir, canType, "cis-eQTLs-Residual_Gene_Expression.txt")
#output_file_name_tra = file.path(eqtlDir, canType, "trans-eQTLs-Residual_Gene_Expression.txt")
pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 0.3e-2
#pvOutputThreshold_tra = 0 # for cis-eQTL only
cisDist = 1e6;

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(genoEurFile);

# Filtering SNPs based on MAF
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
    slice = snps[[sl]];
    maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
    maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

cat(sprintf('SNPs before filtering: %d\n',nrow(snps)))
snps$RowReorder(maf>=0.05);
cat(sprintf('SNPs before filtering: %d\n',nrow(snps)))

# Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
#gene$LoadFile(geneExpRsdEurFile);
gene$LoadFile(geneLog2ExpRsdEurFile);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariatesFile)
#cvrt$LoadFile(character()) # no covariates

## Run the analysis
snpsPosDf = read.table(varLocFile, header = T, as.is = T)
genePosDf = read.table(eqtlGeneLocFile, header = T, as.is = T)

matEqtl = Matrix_eQTL_main(snps = snps, 
                           gene = gene, 
                           cvrt = cvrt,
                           output_file_name = output_file_name_tra,
                           pvOutputThreshold = pvOutputThreshold_tra,
                           useModel = useModel, 
                           errorCovariance = numeric(), 
                           verbose = TRUE, 
                           output_file_name.cis = output_file_name_cis,
                           pvOutputThreshold.cis = pvOutputThreshold_cis,
                           snpspos = snpsPosDf, 
                           genepos = genePosDf,
                           cisDist = cisDist,
                           pvalue.hist = "qqplot",
                           min.pv.by.genesnp = TRUE,
                           noFDRsaveMemory = FALSE);

# eQTL results:
cat('Analysis done in: ', matEqtl$time.in.sec, ' seconds', '\n');
qqPlotFile = file.path(eqtlCanDir, sprintf("Q_Q_Plot-%s.pdf", canType))
pdf(qqPlotFile)
plot(matEqtl)
dev.off()


# Store R image data
save.image(imageFile)

rm(list=ls())

require(reshape2)
require(sqldf)
require(scales)
require(ggplot2)
require(RSQLite)

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

baseFontSize = 12
germDir = file.path(rootDir, "BiO/Research/GCC/Germline")
canTypes = c("CRC")
canType = canTypes[1]

sqldf()

for (canType in canTypes) {
    if (canType == "CRC") {
        ncanType = "COADREAD"
    } else {
        ncanType = canType
    }

    figDir = file.path(germDir, "Figure/Expression_Analysis", canType)
    dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
    canVcfDir = file.path(germDir, "VCF/CANCERS", canType)  
    sqlite3Db = file.path(canVcfDir, sprintf("%s.db", canType))
    sqldf(sprintf("ATTACH '%s' AS %s;", sqlite3Db, canType))

    ## Read gene expression data
    # GA data
    geneExpGaFile = file.path(germDir, sprintf("GDAC/gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014031600.0.0/%s.rnaseqv2__illuminaga_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", ncanType, ncanType))
    geneExpGaLines = readLines(geneExpGaFile)
    geneExpGaLines = geneExpGaLines[-2]
    geneExpGaDf = read.delim(textConnection(geneExpGaLines), header = TRUE, as.is = TRUE)
    colnames(geneExpGaDf) = c(colnames(geneExpGaDf)[1], as.vector(sapply(colnames(geneExpGaDf[,2:ncol(geneExpGaDf)]), function(x) {
                                                                         paste(strsplit(x, "\\.")[[1]][1:4], collapse="_")
})))
    gaSampleIds = grep("TCGA", colnames(geneExpGaDf), value = TRUE)

    # HiSeq data
    geneExpHiseqFile = file.path(germDir, sprintf("GDAC/gdac.broadinstitute.org_%s.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2014031600.0.0/%s.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", ncanType, ncanType))
    geneExpHiseqLines = readLines(geneExpHiseqFile)
    geneExpHiseqLines = geneExpHiseqLines[-2]
    geneExpHiseqDf = read.delim(textConnection(geneExpHiseqLines), header = TRUE, as.is = TRUE)
    colnames(geneExpHiseqDf) = c(colnames(geneExpHiseqDf)[1], as.vector(sapply(colnames(geneExpHiseqDf[,2:ncol(geneExpHiseqDf)]), function(x) {
                                                                               paste(strsplit(x, "\\.")[[1]][1:4], collapse="_")
})))
    hiseqSampleIds = grep("TCGA", colnames(geneExpHiseqDf), value = TRUE)
    sharedRnaSeqSampleIds = intersect(gaSampleIds, hiseqSampleIds)

    # Merge GA and HiSeq data
    geneExpDf = cbind(geneExpGaDf[,-c(1, which(colnames(geneExpGaDf) %in% sharedRnaSeqSampleIds))], geneExpHiseqDf[,-1])
    rownames(geneExpDf) = geneExpGaDf[,1]

    geneLog2ExpDf = log2(geneExpDf + 1)
    meanLog2Exps = rowMeans(geneLog2ExpDf)
    geneLog2FcExpDf = geneLog2ExpDf - meanLog2Exps

    ## Extract SNPs to be analyzed
    # Over-represented known SNPs
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

    novel05Snps1WithExpDf = data.frame()
  
    for (j in 1:nrow(novel05Snps1Df)) {
        print(j)
        r = novel05Snps1Df[j,]
        if (is.na(r$feature)) {
            print(sprintf("No genes nearby. Skipped"))
            next
        }
        mutSampleIds = colnames(r)[grepl("'\\S{1}/\\S{1}'", r, perl = TRUE) & 
                                   !grepl("0/0", r, perl = TRUE) & 
                                   !grepl("\\./\\.", r, perl = TRUE)]
        mutPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", mutSampleIds)
        wdtSampleIds = colnames(r)[grepl("0/0", r, perl = TRUE)]
        wdtPatientIds = gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", wdtSampleIds)

        #geneNames = sapply(strsplit(r$feature, ",")[[1]], function(x) { strsplit(x, "\\(")[[1]][1] })
        geneNames = strsplit(gsub("\\(.*?\\)", "", r$feature), ",")[[1]]

        for (geneName in geneNames) {
            if (geneName == "CAAP1") geneName = "C9orf82"
            varStem = sprintf("%s_%d_%s_%s_%s", r$chrom, r$pos, r$ref, r$alt, geneName)
            figStem = sprintf("chr%s:%d:%s>%s", r$chrom, r$pos, r$ref, r$alt)
            selGeneExps = geneExpDf[grep(paste(geneName, "|", sep=""), rownames(geneExpDf), fixed = T),]
            if (nrow(selGeneExps) != 1) {
                print(sprintf("No expression data for %s", geneName))
                novel05Snps1WithExpDf = rbind(novel05Snps1WithExpDf,
                                                       data.frame(r, 
                                                                  gene = geneName, 
                                                                  mean_mutant_gene_expression = NA,
                                                                  mean_wildtype_gene_expression = NA,
                                                                  mean_total_gene_expression = NA,
                                                                  median_mutant_gene_expression = NA,
                                                                  median_wildtype_gene_expression = NA,
                                                                  median_total_gene_expression = NA,
                                                                  mean_mutant_gene_expression_log2_fold_change = NA,
                                                                  mean_wildtype_gene_expression_log2_fold_change = NA,
                                                                  mean_total_gene_expression_log2_fold_change = NA,
                                                                  median_mutant_gene_expression_log2_fold_change = NA,
                                                                  median_wildtype_gene_expression_log2_fold_change = NA,
                                                                  median_total_gene_expression_log2_fold_change = NA,
                                                                  wt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  tt_pvalue_gene_expression_mutant_vs_wildtype = NA,
                                                                  wt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  tt_pvalue_gene_expression_mutant_vs_total = NA,
                                                                  wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = NA,
                                                                  wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA,
                                                                  tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = NA))
                next
            }
      
            mutGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% mutPatientIds)]
            wdtGeneExps = selGeneExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneExps)) %in% wdtPatientIds)]

            ## RSEM
            selGeneExpMutVsWdtDf = data.frame()
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), rsem = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsWdtDf = rbind(selGeneExpMutVsWdtDf, data.frame(type=rep("WT", ncol(wdtGeneExps)), rsem = as.numeric(wdtGeneExps[1,])))

            wtPvalueGeneExpMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))
            ttPvalueGeneExpMutVsWdt = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(wdtGeneExps[1,]))

            pvalueCondGeneExpMutVsWdt = (wtPvalueGeneExpMutVsWdt < 0.05 | ttPvalueGeneExpMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsWdt) & pvalueCondGeneExpMutVsWdt == TRUE) {
                mutVsWdtExpPlot = ggplot(selGeneExpMutVsWdtDf, aes(factor(type), rsem)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsWdt, ttPvalueGeneExpMutVsWdt)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsWdtExpPlotFile = file.path(figDir, sprintf("%s-MutVsWdt-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtExpPlotFile, plot = mutVsWdtExpPlot, width = 4, height = 5)
            }

            #
            selGeneExpMutVsTotDf = data.frame()
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneExps)), log2FoldChange = as.numeric(mutGeneExps[1,])))
            selGeneExpMutVsTotDf = rbind(selGeneExpMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneExps)), log2FoldChange = as.numeric(selGeneExps[1,])))

            wtPvalueGeneExpMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))
            ttPvalueGeneExpMutVsTot = my.t.test.p.value(as.numeric(mutGeneExps[1,]), as.numeric(selGeneExps[1,]))

            pvalueCondGeneExpMutVsTot = (wtPvalueGeneExpMutVsTot < 0.05 | ttPvalueGeneExpMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneExpMutVsTot) & pvalueCondGeneExpMutVsTot == TRUE) {
                mutVsTotExpPlot = ggplot(selGeneExpMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneExpMutVsTot, ttPvalueGeneExpMutVsTot)) +
                                   scale_y_continuous(name = "RSEM\n")
                                   mutVsTotExpPlotFile = file.path(figDir, sprintf("%s-MutVsTot-RSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotExpPlotFile, plot = mutVsTotExpPlot, width = 4, height = 5)
            }

            ## Log2 Fold Change of RSEM
            selGeneLog2FcExps = geneLog2FcExpDf[grep(paste(geneName, "|", sep=""), rownames(geneLog2FcExpDf), fixed = T),]
            mutGeneLog2FcExps = selGeneLog2FcExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneLog2FcExps)) %in% mutPatientIds)]
            wdtGeneLog2FcExps = selGeneLog2FcExps[, which(gsub("(TCGA_\\S{2}_\\S{4})_\\S{3}", "\\1", colnames(selGeneLog2FcExps)) %in% wdtPatientIds)]

            selGeneLog2FcExpMutVsWdtDf = data.frame()
            selGeneLog2FcExpMutVsWdtDf = rbind(selGeneLog2FcExpMutVsWdtDf, 
                                               data.frame(type=rep(figStem, ncol(mutGeneLog2FcExps)), log2FoldChange = as.numeric(mutGeneLog2FcExps[1,])))
            selGeneLog2FcExpMutVsWdtDf = rbind(selGeneLog2FcExpMutVsWdtDf, 
                                               data.frame(type=rep("WT", ncol(wdtGeneLog2FcExps)), log2FoldChange = as.numeric(wdtGeneLog2FcExps[1,])))

            wtPvalueGeneLog2FcExpMutVsWdt = my.wilcox.test.p.value(as.numeric(mutGeneLog2FcExps[1,]), as.numeric(wdtGeneLog2FcExps[1,]))
            ttPvalueGeneLog2FcExpMutVsWdt = my.t.test.p.value(as.numeric(mutGeneLog2FcExps[1,]), as.numeric(wdtGeneLog2FcExps[1,]))

            pvalueCondGeneLog2FcExpMutVsWdt = (wtPvalueGeneLog2FcExpMutVsWdt < 0.05 | ttPvalueGeneLog2FcExpMutVsWdt < 0.05)

            if (!is.na(pvalueCondGeneLog2FcExpMutVsWdt) & pvalueCondGeneLog2FcExpMutVsWdt == TRUE) {
                mutVsWdtLog2FcExpPlot = ggplot(selGeneLog2FcExpMutVsWdtDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", wtPvalueGeneLog2FcExpMutVsWdt, ttPvalueGeneLog2FcExpMutVsWdt)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsWdtLog2FcExpPlotFile = file.path(figDir, sprintf("%s-MutVsWdt-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsWdtLog2FcExpPlotFile, plot = mutVsWdtLog2FcExpPlot, width = 4, height = 5)
            }

            #
            selGeneLog2FcExpMutVsTotDf = data.frame()
            selGeneLog2FcExpMutVsTotDf = rbind(selGeneLog2FcExpMutVsTotDf, data.frame(type=rep(figStem, ncol(mutGeneLog2FcExps)), log2FoldChange = as.numeric(mutGeneLog2FcExps[1,])))
            selGeneLog2FcExpMutVsTotDf = rbind(selGeneLog2FcExpMutVsTotDf, data.frame(type=rep("Total", ncol(selGeneLog2FcExps)), log2FoldChange = as.numeric(selGeneLog2FcExps[1,])))

            wtPvalueGeneLog2FcExpMutVsTot = my.wilcox.test.p.value(as.numeric(mutGeneLog2FcExps[1,]), as.numeric(selGeneLog2FcExps[1,]))
            ttPvalueGeneLog2FcExpMutVsTot = my.t.test.p.value(as.numeric(mutGeneLog2FcExps[1,]), as.numeric(selGeneLog2FcExps[1,]))

            pvalueCondGeneLog2FcExpMutVsTot = (wtPvalueGeneLog2FcExpMutVsTot < 0.05 | ttPvalueGeneLog2FcExpMutVsTot < 0.05)

            if (!is.na(pvalueCondGeneLog2FcExpMutVsTot) & pvalueCondGeneLog2FcExpMutVsTot == TRUE) {
                mutVsTotLog2FcExpPlot = ggplot(selGeneLog2FcExpMutVsTotDf, aes(factor(type), log2FoldChange)) +
                ggtitle(sprintf("%s\n", geneName)) +
                geom_boxplot() +
                theme(plot.title   = element_text(size = baseFontSize + 2, face="bold"),
                      axis.title.y = element_text(size = baseFontSize, face="plain", family="sans", angle = 90),
                      axis.title.x = element_text(size = baseFontSize, face="plain", family="sans"),
                      axis.text.x  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black", angle = 0, hjust = NULL),
                      axis.text.y  = element_text(size = baseFontSize, face="plain", family="sans", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(color = "black", fill="white"),
                      legend.title = element_text(size = baseFontSize, face="plain", , family="sans", hjust = 0),
                      legend.text  = element_text(size = baseFontSize, family="sans"),
                      legend.direction = "horizontal",
                      legend.position = "bottom") +
                                   scale_x_discrete(name=sprintf("\nWilcoxon signed-rank test: %f\nStudent's t-test: %f", 
                                                                 wtPvalueGeneLog2FcExpMutVsTot, 
                                                                 ttPvalueGeneLog2FcExpMutVsTot)) +
                                   scale_y_continuous(name = "Log2 fold change of RSEM\n")
                                   mutVsTotLog2FcExpPlotFile = file.path(figDir, sprintf("%s-MutVsTot-Log2FoldChangeRSEM.pdf", varStem))
                                   ggsave(filename = mutVsTotLog2FcExpPlotFile, plot = mutVsTotLog2FcExpPlot, width = 4, height = 5)
            }

            novel05Snps1WithExpDf = rbind(novel05Snps1WithExpDf,
                                                   data.frame(r, 
                                                              gene = geneName, 
                                                              mean_mutant_gene_expression = mean(as.numeric(mutGeneExps)),
                                                              mean_wildtype_gene_expression = mean(as.numeric(wdtGeneExps)),
                                                              mean_total_gene_expression = mean(as.numeric(selGeneExps)),
                                                              median_mutant_gene_expression = median(as.numeric(mutGeneExps)),
                                                              median_wildtype_gene_expression = median(as.numeric(wdtGeneExps)),
                                                              median_total_gene_expression = median(as.numeric(selGeneExps)),
                                                              mean_mutant_gene_expression_log2_fold_change = mean(as.numeric(mutGeneLog2FcExps)),
                                                              mean_wildtype_gene_expression_log2_fold_change = mean(as.numeric(wdtGeneLog2FcExps)),
                                                              mean_total_gene_expression_log2_fold_change = mean(as.numeric(selGeneLog2FcExps)),
                                                              median_mutant_gene_expression_log2_fold_change = median(as.numeric(mutGeneLog2FcExps)),
                                                              median_wildtype_gene_expression_log2_fold_change = median(as.numeric(wdtGeneLog2FcExps)),
                                                              median_total_gene_expression_log2_fold_change = median(as.numeric(selGeneLog2FcExps)),
                                                              wt_pvalue_gene_expression_mutant_vs_wildtype = wtPvalueGeneExpMutVsWdt,
                                                              tt_pvalue_gene_expression_mutant_vs_wildtype = ttPvalueGeneExpMutVsWdt,
                                                              wt_pvalue_gene_expression_mutant_vs_total = wtPvalueGeneExpMutVsTot,
                                                              tt_pvalue_gene_expression_mutant_vs_total = ttPvalueGeneExpMutVsTot,
                                                              wt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = wtPvalueGeneLog2FcExpMutVsWdt,
                                                              tt_pvalue_gene_expression_log2_fold_change_mutant_vs_wildtype = ttPvalueGeneLog2FcExpMutVsWdt,
                                                              wt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = wtPvalueGeneLog2FcExpMutVsTot,
                                                              tt_pvalue_gene_expression_log2_fold_change_mutant_vs_total = ttPvalueGeneLog2FcExpMutVsTot))
        }
    }
    novel05Snps1WithExpFile = file.path(canVcfDir, "novel05Snps1WithExp.txt")
    write.table(novel05Snps1WithExpDf, novel05Snps1WithExpFile, sep = "\t", quote = FALSE, row.names = FALSE)
}

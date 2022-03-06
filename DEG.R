## Differentially expressed (DE) genes in alcohol-associated hepatocellular carcinoma
## The normalized gene expression values of the GSE59259 dataset was downloaded from NCBI's GEO database and log transformed. The dataset included 8 pairs of HCC tissue samples (i.e. hepatocellular carcinoma) and their corresponding normal liver tissue from eight donors (total of 16 samples). The `limma` package was used to a fit linear model for each genebank ID in the dataset, contrasting the gene expression in HCC versus liver samples with a paired design. The Bonferroni-Hochberg (BH) method was applied to control the false discovery rate and obtain the list of top DE genes ranked by their BH-adjusted pValues. The expression of top 5 DE genes across study samples was depicted using a heatmap plot.
silent <- suppressPackageStartupMessages

silent(library(limma))
silent(library(GEOquery))
silent(library(stringr))
silent(library(pheatmap))
silent(library(dplyr))

## retrieving and processing data
gse <- getGEO("GSE59259", GSEMatrix=TRUE)[[1]]
exprs(gse) <- log(exprs(gse)+0.001)
sampleTypes <- factor(pData(gse)[ , "tissue:ch1"]) ## samples: 'liver' or 'hcc'
levels(sampleTypes) <- c("hcc", "liver")
sampleTypes <- relevel(sampleTypes, ref="liver")

## DE analysis
experiment <- factor(c(1:8, 1:8)) ## samples are paired
design <- model.matrix(~experiment+sampleTypes+0)
fit <- lmFit(gse, design=design)
fit <- eBayes(fit) ##  smoothing the standard errors
deGenes <- topTable(fit, coef="sampleTypeshcc", adjust.method="BH", p.value=0.01, number = Inf)
deGeneNames <- gsub("ExemplarFor|\'", '', deGenes$DESCRIPTION) ## extracting gene_name
deGeneNames <- str_extract(deGeneNames, 'gene_name\\b(.+?);')
deGeneNames <- gsub("gene_name|;| ", '', deGeneNames)
deGenes[, "geneName"] <- deGeneNames

## saving deGenes to 'de.csv'
de <- deGenes[, c("geneName", "adj.P.Val")]
de <- distinct(de, geneName, .keep_all=TRUE) # multiple genebank IDs per gene_name
de <- de[de$geneName != "N/A",] # remove genebank IDs without gene_name
print("==================================================")
print("Writing DE genes to de.csv ...")
write.table(de, file="de.csv", quote=FALSE, sep=",", row.names=FALSE)

## Number of differentially expressed genes
print("==================================================")
print(paste("A total of ", nrow(de), " distinct genes were differentially expressed in hepatocellular carcinoma tissues as compared to the normal liver tissues"))
print("==================================================")
## Expression of top 5 DE genes
topFiveIds <- row.names(de)[1:5] # IDs of top five gene
topFiveEx <- as.matrix(exprs(gse)[topFiveIds, ]) # expression of top five genes
rownames(topFiveEx) <- de[ ,"geneName"][1:5]
annotations <- data.frame(Donor=paste0("Donor ", experiment), Tissue=c(rep("Liver", 8), rep("HCC", 8)))
rownames(annotations) <- colnames(topFiveEx)
heatmap <- pheatmap(topFiveEx, cluster_cols=FALSE, cluster_rows=FALSE, scale="row", annotation_col=annotations)

save.pheatmap <- function(object, file, width=1500, height=800, dpi=200) {
    png(file, width=width, height=height, res=dpi)
    grid::grid.newpage()
    grid::grid.draw(object$gtable)
    dev.off()
}
print("Saving heatmap of top 5 DE genes.....")
save.pheatmap(heatmap, "heatmap.png")
print("done!")

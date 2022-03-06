silent <- suppressPackageStartupMessages
source('./Functions.R')
silent(library(TCGAbiolinks))
silent(library(maftools))
silent(library(data.table))
silent(library(survival))

## Get maf and clinical data
print("1. Downloading maf file of the TCGA-LIC project (muse variant calling pipeline")
maf <- GDCquery_Maf("LIHC", pipelines="muse")

print("")
print("2. Downloading the clinical data of the TCGA-LIC project")
clinicalData <- get.clinicalData(project="TCGA-LIHC")

print("3. Parsing the MAF file....")
lihc <- read.maf(maf, clinicalData, removeDuplicatedVariants=TRUE, useAll=TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile=NULL, gisticDelGenesFile=NULL, gisticScoresFile=NULL, cnLevel="all", cnTable=NULL, isTCGA=TRUE, vc_nonSyn=NULL, verbose=TRUE)
## top 3 mutated genes
print("|======================================================================|")
top3 <- top.n.genes(n=3, maf=lihc)
print("top three most frequently mutated genes: ")
print(top3)

## Plots of top3 mutated genes
clinOs <- os(clinicalData)
sapply(top3, function(x) surv.km(gene=x, maf=lihc, clinicalData=clinOs))



# `utils::combn()`, which takes the vector of top n1
# mutated genes and generates all combinations of genes in the vector, taken k1 at a time.
# `combn` also applies the `impact` function to each combination. `impact` is a user
# defined function that takes genes in each combination, performs a log rank test and
# return genes names, the Pvalue and adjusted Pvalue of the log rank test.
print("|======================================================================|")
print(paste0("The line of code implementing combn(): ", "comb <- combn(topN, k1, impact)"))
print("|======================================================================|")
## Best k1 genes among top n1 most frequently mutated genes
print("The best combination k1=3, n1=3: ")
most.impact(n1=3, k1=3, maf=lihc, clinicalData=clinicalData)
print("|======================================================================|")
print("The best combination k1=3, n1=5: ")
most.impact(n1=10, k1=3, maf=lihc, clinicalData=clinicalData)
#print("The best combination k1=3, n1=100: ")
#most.impact(n1=100, k1=3, maf=lihc, clinicalData=clinicalData)

## Utility Functions

top.n.genes <- function(n, maf){
  #' get the top n most frequently mutated genes in a MAF file
  #'
  #' @param n Number of top most frequently mutated genes
  #' @param maf A MAF object (from maftools::read.maf())
  #' @return A vector of HUGO gene name
  stopifnot("MAF" %in% class(maf))
  genes = getGeneSummary(x = maf)[order(AlteredSamples, decreasing = TRUE)][1:n, Hugo_Symbol]
  genes <- as.list(genes)
  return(genes)
}

get.clinicalData <- function(project=NULL, dir=".", info = "patient"){
  #' Takes a TCGA project name and returns the clinical data of the project
  #'
  #' @param project name of the TCGA project. see a list of projects: TCGAbiolinks:::getGDCprojects()$project_id
  #' @dir directory for downloading clinical data
  #' @info The clinical.info that should be retrieved. The value will be passed to GDCprepare_clinic(clinical.info=).Options Clinical: drug, admin, follow_up,radiation, patient, stage_event or new_tumor_event Options Biospecimen: protocol, admin, aliquot, analyte, bio_patient, sample, portion, slide.
  if(is.null(project)){
    stop("please privide the name of a TCGA project (e.g. TCGA-LIHC) to get the clinical data. use TCGAbiolinks:::getGDCprojects()$project_id to see a list of all projects")
  }
  query <- TCGAbiolinks::GDCquery(project=project, data.category="Clinical", file.type="xml", legacy=FALSE)
  TCGAbiolinks::GDCdownload(query, directory=".")
  clinicalData <- TCGAbiolinks::GDCprepare_clinic(query, clinical.info = info, directory = ".")
    if(!"Tumor_Sample_Barcode" %in% colnames(clinicalData) & "bcr_patient_barcode" %in% colnames(clinicalData)){
    clinicalData[ , "Tumor_Sample_Barcode"] <- clinicalData[ , "bcr_patient_barcode"]
  }
  clinicalData <- data.table::setDT(clinicalData)
  return(clinicalData)
}

os <- function(clinicalData, statusCol="vital_status"){
  #' Calculates overall survival from clinical data of a TCGA project
  #'
  #' The function accepts as input clinicalData of a TCGA project and adds a column named Time to the clinical dataset and return the clinicaalData. Time (i.e. overall survival) = days_to_death (for patients with the "Dead" vital_status) or days_to_last_followup (for patients with the "alive" vital_status).
  #' Before calculating Time, samples without any follow up data and samples with a negative or zero value for the follow up data are filtered out (see https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/324).
  #'
  #' @param clinicalData clinical data of a TCGA project
  #' @param statusCol Name of the column containing event status (default: survival_status)
  #' @return the clinicalData with an extra column named Time, representing overall survival
  if(is.null(clinicalData)){
  stop("please provide a data.table containing clinicalData or use get.clinicalData() to download clinicalData for a TCGA project.")
  }
  if(!inherits(clinicalData, "data.table")){clinicalData <- data.table::setDT(clinicalData)}
  clinicalData$Time <- ifelse(clinicalData[ , ..statusCol] == "Dead", clinicalData[ , days_to_death], clinicalData[ , days_to_last_followup])
  clinicalData <- clinicalData[!duplicated(clinicalData$patient_id),]  ## ensuring no duplicates
  clinicalData <- clinicalData[Time > 0 & !is.na(Time),]
  return(clinicalData)
}

surv.km <- function(maf, gene = NULL, clinicalData=NULL, project=NULL, time = "os", statusCol = "vital_status", col = c("maroon", "royalblue")) {
  #' A function to plot km curves for a given gene name (a modification of maftools::mafSurvival())
  #'
  #' @param maf A MAF object (maftools::read.maf())
  #' @param gene A HUGO gene symbol.
  #' @param clinicalData the clinical data of a TCGA project. If missed, FUN uses the project parameter to retrieve the data using get.clinicalData()
  #' @param project Either provide clinicalData or the name of a TCGA project
  #' @param time If os, function calculates & uses overall survival. Alternatively, name of a column containing time to event.
  #' @param statCol Name of the column in clinicalData containing status of the event
  #' @param col A vector of length two. Color of km lines
  if (is.null(gene)) {
        stop("Provide a gene name for survival analysis")
  }
  if(is.null(clinicalData)){
      clinicalData <- get.clinicalData(project = project)
  } else {
      clinicalData <- data.table::setDT(clinicalData)
  }
  if (time == "os") {
        clinicalData <- os(clinicalData, statusCol)
  } else {
        if (!time %in% colnames(clinicalData)) {
            stop("time should either be os (for overall survival) or name of the column with time to event data in clinialData")
        }
        colnames(clinicalData)[colnames(clinicalData) %in% time] = "Time"
    }
  if (!statusCol %in% colnames(clinicalData)) {
        print(colnames(clinicalData))
        stop(paste0("Column", statusCol, "not found in clinical data. Provide column name containing vital status (Dead or Alive)"))
  }
  ## identifying samples with mutations in the gene
  geneTSB = genesToBarcodes(maf = maf, genes = gene, justNames = TRUE)[[1]]
  if(length(geneTSB) == 0){stop("No sample is mutated for this gene")}

  clinicalData$Time = suppressWarnings(as.numeric(clinicalData$Time))
  clinicalData$Status = ifelse(clinicalData[ , ..statusCol] == "Dead", 1, 0)
  clinicalData$Group = ifelse(clinicalData[, Tumor_Sample_Barcode] %in% geneTSB, "Mutant", "WT")
  ## survival analysis
  surv.km = survival::survfit(formula=survival::Surv(time=Time, event=Status) ~ Group, data=clinicalData, conf.type="log-log")
  res = summary(surv.km)
  surv.diff = survival::survdiff(formula=survival::Surv(time=Time, event=Status) ~ Group, data = clinicalData)
  surv.diff.pval = signif(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1), digits=3)
  surv.cox = survival::coxph(formula=survival::Surv(time=Time, event=Status) ~ Group, data = clinicalData)
  hr = signif(1/exp(stats::coef(surv.cox)), digits = 3)
  surv.dat=data.table::data.table(Group=res$strata, Time=res$time, survProb=res$surv, survUp=res$upper, survLower = res$lower)
  surv.dat$Group = gsub(pattern="Group=", replacement="", x=surv.dat$Group)
  ## ploting km
  png(paste0(gene, "_km.png"))
  par(mar = c(4, 4, 2, 2))
  x_lims = pretty(surv.km$time)
  y_lims = seq(0, 1, 0.2)
  plt <- plot(NA, xlim=c(min(x_lims), max(x_lims)), ylim = c(0, 1), frame.plot=FALSE, axes=FALSE, xlab=NA, ylab=NA)
  abline(h=y_lims, v=x_lims, lty=2, col=grDevices::adjustcolor(col="gray70", alpha.f=0.5), lwd=0.75)
  points(surv.dat[Group %in% "Mutant", Time], y=surv.dat[Group %in% "Mutant", survProb], pch=8, col = col[1])
  points(surv.dat[Group %in% "WT", Time], y=surv.dat[Group %in% "WT", survProb], pch=8, col=col[2])
  lines(surv.km[1], col=col[1], lty = 1, lwd=2, conf.int = FALSE)
  lines(surv.km[2], col=col[2], lty=1, lwd=2, conf.int = FALSE)
  axis(side=1, at=x_lims)
  axis(side =2, at=y_lims, las = 2)
  mtext(text="Survival probability", side = 2, line = 2.5)
  legend(x = "topright", legend = c("Mutant", "WT"), col=col, bty="n", lwd=2, pch = 8)
  title(main=NA, sub=paste0("P-value: ", surv.diff.pval, "; ", "HR: ", hr), cex.main=1, font.main = 4, col.main="black", cex.sub=1, font.sub=3, col.sub=ifelse(test=surv.diff.pval < 0.05, yes="red", no="black"), line=2.5, adj=0)
  title(main=paste0(gene," : Mutant", " v/s ", "WT"), adj=0, font.main=4)
  dev.off()
}

## `most.impact()` implemented using `utils::combn()`
most.impact <- function(n1, k1, maf, clinicalData){
  #'
  #' performs a log rank test comparing the survival of patients with or without mutations in a gene set
  #'
  #' @param n1 The number of top mutated genes in MAF object
  #' @param K1 The size of each gene set
  #' @param maf A MAF object (maftools::read.maf())
  #' @param clinicalData clinical data of the MAF file, containing vital_status, days_to_death, days_to_the-last_followup, Tumor_sample_Barcode columns
  #' @return Returns a list object containing name of the gene set with the lowest pvalue of the log rank test

  stopifnot("Provide a MAF object (maftools::read.maf())"= !missing(maf) & inherits(maf, "MAF"))
  stopifnot("Provide a clinicalata (get.clinicalData())"= !missing(clinicalData))
  topN <- unlist(top.n.genes(n=n1, maf=maf))
  clinicalData <- os(clinicalData) # time = overall survival
  Time = clinicalData$Time
  Status = ifelse(clinicalData[ , vital_status] == "Dead", 1, 0)
  ## given each gene set, identify mutated samples & calculate a pvalue
  impact <- function(x){
    tsb <- genesToBarcodes(genes=x, maf=maf,  justNames=TRUE, verbose=FALSE)
    ## samples with mutations in all members of a gene set
    mutSamples <- Reduce(intersect, tsb)
    Group <- ifelse(clinicalData[, Tumor_Sample_Barcode] %in% mutSamples, "Mutant", "WT")
    ## if any sample had all gene mutated?
    if(length(unique(Group)) == 2){
        surv.diff <- survdiff(formula=Surv(time=Time, event=Status) ~ Group)
        surv.diff.pval <- signif(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1), digits=3)
        return(list(x, surv.diff.pval))
    } else {
        return(list(x, NA))
    }
  }
  comb <- combn(topN, k1, impact)
  pval.adj <- p.adjust(p=comb[2,], method="BH")
  min.idx <- which.min(pval.adj)
  set <- comb[[1 , min.idx]]
  if(sum(!is.na(comb[2,])) == 0){
    message("no sample was mutated for all members of gene sets of length ", k, " in the list of ", n, "most mutated genes")
  } else {
    return(list(set=set, impact=comb[2, min.idx] , impact.adj=pval.adj[min.idx]))
  }
}

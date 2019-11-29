####### A GSEA-based metric for identifying differential gene set enrichment#######
####### By James Joly, USC (Advisor: Nicholas Graham, PhD)                  #######

rm(list = ls())

library(dplyr)
library(foreach)
library(doSNOW)
cl <- makeCluster(4, type="SOCK") # for 4 cores machine
registerDoSNOW (cl)

### Do you want to make mountain plots?
Make.Plots <- FALSE

###### Read in data

#Read in gene expression data
#Genes should be first column, named "Gene"
#Samples should be columsn 2 - N
setwd("/Users/JamesJoly/Documents/Research/Bioinformatics/Differential\ Pathway\ GSEA/GSEA\ on\ Corr\ bt\ DPEA\ Pathway\ to\ DepMap/")
data_in <- read.table(file = "DEPMAP Spearman Correlations v4 - all genes.csv", header=T, sep = ",")


#Read in gene sets
setwd("/Users/JamesJoly/Documents/Research/Bioinformatics/Differential\ Pathway\ GSEA/")
#Enter gene sets to be compared (e.g. glycolysis and oxphos)
metabolic_gene_sets <- read.table(file = "Human Core Glycolysis OxPhos gene sets.txt", header = T, sep = "\t")
# Enter background gene sets (e.g. all metabolic pathways)
all_pathways <- read.table(file = "KEGG_Metabolic_pathways.txt", sep = "\t", header = T)

# change directory for output files
setwd("/Users/JamesJoly/Documents/Research/Bioinformatics/Differential\ Pathway\ GSEA/GSEA\ on\ Corr\ bt\ DPEA\ Pathway\ to\ DepMap/")
output.file.name <- "Differential Core Glycolysis OxPhos - DepMap Correlations Glyc OxPhos DGSEA.csv"
nperm = 1000

### What GSEA metric do you want to use? Classic or Weighted? Set one of the below to true
Classic = FALSE
Weighted = TRUE

if (Classic == TRUE){
  score.weight = 0
}
if (Weighted == TRUE){
  score.weight = 1
}
if ( (Classic == TRUE)& (Weighted == TRUE)){
  stop("Set either Classic or Weighted to TRUE, not both") 
}
if ( (Classic == FALSE)& (Weighted == FALSE)){
  stop("Please specify which statistic to use, set one of Classic or Weighted to TRUE") 
}


GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL){
  tag.indicators <- sign(match(gene.list, gene.set, nomatch = 0))
  no.tag.indicator <- 1 - tag.indicators
  N <- length(gene.list)
  Nh <- numhits_pathway 
  Nm <- N - Nh
  if (weighted.score.type == 0){
    correl.vector <- rep(1,N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicators == 1])
  norm.tag <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicators * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicators))
} #for real ES

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = score.weight, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- numhits_pathway
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(ES)
  
} #for permutation ES

cells <- colnames(data_in)
cells <- cells[-1]

pathways <-colnames(metabolic_gene_sets)

### Annotate gene sets
data_in$Glycolysis.Pathway <- 0
data_in$OxPhos.Pathway <- 0
data_in$Hypoxia.Pathway <- 0

#OxPhos
for (i in 1:nrow(data_in)){
  if (data_in$Gene[i] %in% metabolic_gene_sets$OxPhos.Pathway){
    data_in$OxPhos.Pathway[i] = "X";
  }
}
#Glycolysis
for (i in 1:nrow(data_in)){
  if (data_in$Gene[i] %in% metabolic_gene_sets$Glycolysis.Pathway){
    data_in$Glycolysis.Pathway[i] = "X";
  }
}


all_pathways_unique <- colnames(all_pathways)

annotations <- matrix(data = 0, nrow = nrow(data_in), ncol = length(all_pathways_unique))
colnames(annotations) <- all_pathways_unique

annotations <- as.data.frame(annotations)

annotations <- cbind(data_in$Gene,annotations)
colnames(annotations) <- c("Gene", all_pathways_unique)
annotations <- as.matrix(annotations)
### Annotate gene sets
for (j in 1:length(all_pathways_unique)){
  temp.pathway <- all_pathways[,all_pathways_unique[j]]
  for (i in 1:nrow(annotations)){
    if (annotations[i,"Gene"] %in% temp.pathway){
      annotations[i,j+1] = "X";
    }
  }
}
annotations <- as.data.frame(annotations)
data_in <- merge(data_in, annotations, by = "Gene")

data_in <- na.omit(data_in)

Results.All.Cells <- matrix(data = NA, nrow = 0, ncol = 7)
colnames(Results.All.Cells) <- c("Cell","Pathway","KS","KS_Normalized",
                                 "Frequency of Random KS Scores > KS Score observed","Position at Max",
                                 "FDR q-value")

combined_pathways <- c(pathways,all_pathways_unique)

rm(annotations)
data_in2 <- array(data = NA)

for (u in 1:length(cells)){
  loop.time <- Sys.time()
  
  data_in2 <- cbind(subset(data_in, select = all_pathways_unique),
                    select(data_in, contains("Pathway")),
                    select(data_in, cells[u]))  #select one cell type and the genes and pathways
  data_in2 <- data_in2[order(-data_in2[,cells[u]]),] #sort by descending order for the rank metric
  rownames(data_in2) <- 1:nrow(data_in2) #reorder row indices for counting in for loop below


## Assuming first two columns in data table are Genes and Rank Metric (e.g. Foldchange, SNR)

Results <- matrix(data = NA, nrow = length(combined_pathways), ncol = 7)
colnames(Results) <- c("Cell","Pathway","KS","KS_Normalized",
                       "Frequency of Random KS Scores > KS Score observed","Position at Max",
                      "FDR q-value")
Results <- as.data.frame(Results)
Results$Pathway <- combined_pathways
Results$Cell <- cells[u]

ions <- nrow(data_in2)



##for plotting
if (Make.Plots == TRUE){
  ks_results_plot <- matrix(data = NA, nrow = nrow(data_in), ncol = length(pathways))
  colnames(ks_results_plot) <- pathways
  ks_results_plot <- as.data.frame(ks_results_plot)
}
gene.list <- 1:ions
rank_metric <- data_in2[,cells[u]] #Save the rank metric

pos_gene_set <- array(data = 0, dim = nrow(data_in2), dimnames = NULL);

## Calculate Real KS Statistic
for (i in 1:length(combined_pathways)){
  data_in3 <- data_in2[,combined_pathways[i]]
  numhits_pathway <- sum(data_in3 == "X"); #check to see if there is anything in the column (e.g. X)
  if (numhits_pathway > 1){
    pos_gene_set <- which(data_in2[,combined_pathways[i]] %in% c("X"))
    if (combined_pathways[i] == "Glycolysis.Pathway"){
      pos_gene_set.Glycolysis <- pos_gene_set
    }
    if (combined_pathways[i] == "OxPhos.Pathway"){
      pos_gene_set.OxPhos <- pos_gene_set
    }
    #Calculate Real KS statistic and save the maximum (pos or neg) as well as the position of the maximum to check
    #pos_gene_set <- as.integer(pos_gene_set)        
    KS_real <- GSEA.EnrichmentScore(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
    Results[Results$Pathway == combined_pathways[i],]$KS <- KS_real$ES;
    Results[Results$Pathway == combined_pathways[i],]$`Position at Max` <- KS_real$arg.ES;
    if (Make.Plots == TRUE){
      ks_results_plot[,combined_pathways[i]] = KS_real$RES
    }
  }
}


rand_order_gene_set <- cbind(select(data_in, ends_with("Pathway")), subset(data_in, select = all_pathways_unique))
rand_order_gene_set <- as.matrix(rand_order_gene_set)
## Permutations using the same random shuffled gene order for each pathway



KSRandomArray <- list()
print("Calculating permutations...")
KSRandomArray <- foreach(L = 1:nperm, .combine = "list") %dopar% {
  temp.KSRandomArray <- matrix(data = 0, nrow = 1, ncol = length(combined_pathways))
    #Shuffle Genes
    rand_order_gene_set <- rand_order_gene_set[sample(1:nrow(rand_order_gene_set)),] #Shuffle genes 
    for (i in 1:length(combined_pathways)) {
      pos_gene_set <- which(rand_order_gene_set[,combined_pathways[i]] %in% c("X"))
      numhits_pathway <- sum(rand_order_gene_set[,combined_pathways[i]] == "X")
      temp.KSRandomArray[,i] <- GSEA.EnrichmentScore2(gene.list, pos_gene_set, weighted.score.type = score.weight, correl.vector = rank_metric)
    }
    temp.KSRandomArray
}

KSRandomArray <- data.frame(matrix(unlist(KSRandomArray), nrow = nperm, byrow = T))
colnames(KSRandomArray) <- combined_pathways
KSRandomArray <- na.omit(KSRandomArray)

rm(rand_order_gene_set)
print("Normalizing Scores...")
KSRandomArray <- as.data.frame(KSRandomArray)
###normalize the GSEA distribution
KSRandomArray.Norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = ncol(KSRandomArray))
colnames(KSRandomArray.Norm) <- colnames(KSRandomArray)
avg <- 0
KSRandomArray.temp <- 0
for (i in 1:ncol(KSRandomArray.Norm)){
  avg <- 0
  KSRandomArray.temp <- KSRandomArray[,i]
  pos.temp <- KSRandomArray.temp[which(KSRandomArray.temp >= 0)]
  neg.temp <- KSRandomArray.temp[which(KSRandomArray.temp < 0)]
  
  avg.pos <- mean(pos.temp)
  avg.neg <- mean(neg.temp)
  
  norm.pos.temp <- pos.temp / avg.pos
  norm.neg.temp <- neg.temp / avg.neg * -1
  
  norm.perms <- c(norm.pos.temp,norm.neg.temp)
  
  KSRandomArray.Norm[,i] <- norm.perms 
  
}


KSRandomArray$DGSEA <- KSRandomArray$Glycolysis.Pathway - KSRandomArray$OxPhos.Pathway

KSRandomArray.Norm <- as.data.frame(KSRandomArray.Norm)
pos.DGSEA.perms <- KSRandomArray[which(KSRandomArray$DGSEA >= 0),]$DGSEA
neg.DGSEA.perms <- KSRandomArray[which(KSRandomArray$DGSEA < 0),]$DGSEA

norm.pos.DGSEA.perms <- pos.DGSEA.perms / mean(pos.DGSEA.perms)
norm.neg.DGSEA.perms <- neg.DGSEA.perms / mean(neg.DGSEA.perms) * -1

norm.DGSEA.perms <- c(norm.pos.DGSEA.perms,norm.neg.DGSEA.perms)

KSRandomArray.Norm$DGSEA <- norm.DGSEA.perms


###Make a DGSEA KSRandomArray for Glyc - All Pathways Unique
DGSEA_RandomArray_A <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(all_pathways_unique))
DGSEA_RandomArray_A <- as.data.frame(DGSEA_RandomArray_A)
DGSEA_RandomArray_B <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(all_pathways_unique))
DGSEA_RandomArray_B <- as.data.frame(DGSEA_RandomArray_B)

DGSEA_Real_A <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(all_pathways_unique)))

DGSEA_Real_B <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(all_pathways_unique)))


for (i in 1:length(all_pathways_unique)){
  identifier <- paste("DGSEA A", all_pathways_unique[i])
  colnames(DGSEA_RandomArray_A)[i] <- identifier
  DGSEA_RandomArray_A[,i] <- KSRandomArray[,"Glycolysis.Pathway"] - KSRandomArray[,all_pathways_unique[i]]
  colnames(DGSEA_Real_A)[i] <- identifier
  DGSEA_Real_A[,i] <- Results[Results$Pathway == "Glycolysis.Pathway",]$KS -Results[Results$Pathway == all_pathways_unique[i],]$KS
  
  identifier <- paste("DGSEA B", all_pathways_unique[i])
  colnames(DGSEA_RandomArray_B)[i] <- identifier
  DGSEA_RandomArray_B[,i] <- KSRandomArray[,"OxPhos.Pathway"] - KSRandomArray[,all_pathways_unique[i]]
  colnames(DGSEA_Real_B)[i] <- identifier
  DGSEA_Real_B[,i] <-  Results[Results$Pathway == all_pathways_unique[i],]$KS - Results[Results$Pathway == "OxPhos.Pathway",]$KS
}

DGSEA_Real <- cbind(DGSEA_Real_A,DGSEA_Real_B)
DGSEA_RandomArray <- cbind(DGSEA_RandomArray_A, DGSEA_RandomArray_B)

#normalize the DGSEA distribution
DEEPS <- colnames(DGSEA_RandomArray)
DEEPS.norm <- matrix(data = NA, nrow = nrow(KSRandomArray), ncol = length(DEEPS))
colnames(DEEPS.norm) <- DEEPS
Real.DEEPS.NES <- matrix(data = NA, nrow = 1, ncol = length(DEEPS))
colnames(Real.DEEPS.NES) <- DEEPS
for (i in 1:length(DEEPS)){
  avg <- 0
  
  DGSEA_RandomArray.temp <- DGSEA_RandomArray[,DEEPS[i]]
  pos.temp <- DGSEA_RandomArray.temp[which(DGSEA_RandomArray.temp >= 0)]
  neg.temp <- DGSEA_RandomArray.temp[which(DGSEA_RandomArray.temp < 0)]
  
  avg.pos <- mean(pos.temp)
  avg.neg <- mean(neg.temp)
  
  norm.pos.temp <- pos.temp / avg.pos
  norm.neg.temp <- neg.temp / avg.neg * -1
  
  norm.deeps <- c(norm.pos.temp,norm.neg.temp)
  if (is.finite(avg.pos)){
    DEEPS.norm[,DEEPS[i]] <- norm.deeps
  }
  real_dES <- DGSEA_Real[1,DEEPS[i]]
  if (real_dES >= 0){
    Real.DEEPS.NES[1,DEEPS[i]] <- real_dES / avg.pos
  } else if (real_dES < 0){
    Real.DEEPS.NES[1,DEEPS[i]] <- real_dES / avg.neg * -1
  }
}

#Remove NA
DEEPS.norm <- DEEPS.norm[,!is.na(colSums(DEEPS.norm))]
Real.DEEPS.NES <- Real.DEEPS.NES[,!is.na(colSums(Real.DEEPS.NES))]
#Now we have a distribution of DGSEA NES across all pathways and a distribution of Real DGSEA NES
#We can use this for FDR calcs.

DGSEA <- matrix(data = NA, nrow = 1, ncol = 7)
colnames(DGSEA) <- c("Cell","Pathway","KS","KS_Normalized",
                    "Frequency of Random KS Scores > KS Score observed","Position at Max",
                    "FDR q-value")
DGSEA <- as.data.frame(DGSEA)
DGSEA$Pathway <- c('DGSEA')
DGSEA$Cell <- cells[u]

DGSEA$KS<- Results[Results$Pathway == "Glycolysis.Pathway",]$KS - Results[Results$Pathway == "OxPhos.Pathway",]$KS
Results <- rbind(DGSEA,Results)

pathways2 <- unique(Results$Pathway)

#Calculate Normalized KS Statistic - sum all KS from permutations and average
for (i in 1:length(pathways2)){
  avg.perms <- 0
  #Normalize by taking the average permutation
  if (Results[Results$Pathway == pathways2[i],]$KS > 0){
    KSRandomArray.temp <- KSRandomArray[,pathways2[i]]
    KSRandomArray.temp <- KSRandomArray.temp[which(KSRandomArray.temp > 0)] #take only pos values
    avg.perms <- sum(KSRandomArray.temp) / length(KSRandomArray.temp)
    Results[Results$Pathway == pathways2[i],]$KS_Normalized = signif(Results[Results$Pathway == pathways2[i],]$KS / avg.perms, digits = 3);
  } else if (Results[Results$Pathway == pathways2[i],]$KS < 0){
    KSRandomArray.temp <- KSRandomArray[,pathways2[i]]
    KSRandomArray.temp <- KSRandomArray.temp[which(KSRandomArray.temp < 0)] #take only neg values
    avg.perms <- sum(KSRandomArray.temp) / length(KSRandomArray.temp) * -1 #multiply by -1 to turn the value positive, so we aren't dividing by neg
    Results[Results$Pathway == pathways2[i],]$KS_Normalized = signif(Results[Results$Pathway == pathways2[i],]$KS / avg.perms, digits = 3);
  }
  #Determine how many times Random Shuffle KS outperforms real KS
  count <- 0;
    if (Results[Results$Pathway == pathways2[i],]$KS > 0) {
      count <- sum(KSRandomArray[,pathways2[i]] > Results[Results$Pathway == pathways2[i],]$KS)
    } else if (Results[Results$Pathway == pathways2[i],]$KS < 0) {
      count <- sum(KSRandomArray[,pathways2[i]] < Results[Results$Pathway == pathways2[i],]$KS)
    }
    Results[Results$Pathway == pathways2[i],]$`Frequency of Random KS Scores > KS Score observed` = count / length(KSRandomArray.temp); 
  
}


Total.NES.DGSEA <- as.vector(DEEPS.norm)


percent.real.NES.pos <- sum(Real.DEEPS.NES > 0) / length(Real.DEEPS.NES)
percent.real.NES.neg <- sum(Real.DEEPS.NES < 0) / length(Real.DEEPS.NES)

print("Calculating FDR...")

### Change Total.NES.DGSEA to the permutation NES.DGSEA for A-B
if (Results[Results$Pathway == "DGSEA",]$KS > 0){
  temp <- sum(norm.pos.DGSEA.perms > Results[Results$Pathway == "DGSEA",]$KS_Normalized)
  percent.temp <- temp / length(norm.pos.DGSEA.perms)
  if (percent.temp / percent.real.NES.pos < 1){
    Results[Results$Pathway == "DGSEA",]$`FDR q-value` <- percent.temp / percent.real.NES.pos
  } else {Results[Results$Pathway == "DGSEA",]$`FDR q-value` = 1}
} else if (Results[Results$Pathway == "DGSEA",]$KS < 0){
  temp <- sum(norm.neg.DGSEA.perms < Results[Results$Pathway == "DGSEA",]$KS_Normalized)
  percent.temp <- temp / length(norm.neg.DGSEA.perms)
  if (percent.temp / percent.real.NES.neg < 1){
    Results[Results$Pathway == "DGSEA",]$`FDR q-value` <- percent.temp / percent.real.NES.neg
  } else {Results[Results$Pathway == "DGSEA",]$`FDR q-value` = 1}
}

Total.NES.GSEA <- as.vector(KSRandomArray.Norm)

#Start at 2 b/c DGSEA is pathways2[1]

num.pos <- Results[which(Results$KS_Normalized >= 0),]
num.neg <- Results[which(Results$KS_Normalized < 0),]
percent.real.NES.pos <- nrow(num.pos) / nrow(Results)
percent.real.NES.neg <- nrow(num.neg) / nrow(Results)
for(i in 2:length(pathways2)){
  KSRandomArray.Norm.temp <- KSRandomArray.Norm[,pathways2[i]]
  KSRandomArray.Norm.temp.pos <- KSRandomArray.Norm.temp[which(KSRandomArray.Norm.temp >= 0)]
  KSRandomArray.Norm.temp.neg <- KSRandomArray.Norm.temp[which(KSRandomArray.Norm.temp < 0)]
  if (Results[Results$Pathway == pathways2[i],]$KS > 0){
    temp <- sum(KSRandomArray.Norm.temp.pos > Results[Results$Pathway == pathways2[i],]$KS_Normalized)
    percent.temp <- temp / length(KSRandomArray.Norm.temp.pos)
    if (percent.temp / percent.real.NES.pos < 1){
      Results[Results$Pathway == pathways2[i],]$`FDR q-value` <- percent.temp / percent.real.NES.pos
      } else {Results[Results$Pathway == pathways2[i],]$`FDR q-value` = 1}
    } 
  else if (Results[Results$Pathway == pathways2[i],]$KS < 0){
    temp <- sum(KSRandomArray.Norm.temp.neg < Results[Results$Pathway == pathways2[i],]$KS_Normalized)
    percent.temp <- temp / length(KSRandomArray.Norm.temp.neg)
    if (percent.temp /percent.real.NES.neg < 1){
      Results[Results$Pathway == pathways2[i],]$`FDR q-value` <- percent.temp / percent.real.NES.neg
    } else {Results[Results$Pathway == pathways2[i],]$`FDR q-value` = 1}
  }
}

Results.All.Cells <- rbind(Results.All.Cells,Results)



####################################################
####### Below is only for data visualizaiton #######

  ## Create GSEA plot
  # Save default for resetting
if (Make.Plots == TRUE){
  p <- {
    def.par <- par(no.readonly = TRUE)

    # Create a new device of appropriate size
    dev.new(width = 3, height = 3)
   
    # Create a division of the device
    gsea.layout <- layout(matrix(c(1, 2, 3)), heights = c(2, 0.5, 2))
    layout.show(gsea.layout)
    

    # Create plots
    par(mar = c(0, 5, 2, 2))
    gsea.hit.indices <- pos_gene_set
    enrichment.score.range <- c(min(KS_real$RES), max(KS_real$RES))
    
    gsea.hit.indices.Glycolysis <- pos_gene_set.Glycolysis
    gsea.hit.indices.OxPhos <- pos_gene_set.OxPhos
    
    gsea.es.profile.Glycolysis <- ks_results_plot[gsea.hit.indices.Glycolysis,"Glycolysis.Pathway"]
    gsea.es.profile.OxPhos <- ks_results_plot[gsea.hit.indices.OxPhos,"OxPhos.Pathway"]
    
    enrichment.score.range <- c(min(ks_results_plot), max(ks_results_plot))
    
    
    plot(c(1, gsea.hit.indices.Glycolysis, length(data_in2[,cells[u]])),
         c(0, gsea.es.profile.Glycolysis, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
         xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
         ylim = enrichment.score.range,
         #main = list(gsea.gene.set, font = 1, cex = 1),
         panel.first = {
           abline(h = seq(round(enrichment.score.range[1], digits = 1),
                          enrichment.score.range[2], 0.1),
                  col = "gray95", lty = 2)
           abline(h = 0, col = "gray50", lty = 2)
          }
         )
    lines(c(1, gsea.hit.indices.OxPhos, length(data_in2[,cells[u]])),
         c(0, gsea.es.profile.OxPhos, 0), type = "l", col = "blue", lwd = 1.5, xaxt = "n",
         xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
         ylim = enrichment.score.range,
         #main = list(gsea.gene.set, font = 1, cex = 1),
         panel.first = {
           abline(h = seq(round(enrichment.score.range[1], digits = 1),
                          enrichment.score.range[2], 0.1),
                  col = "gray95", lty = 2)
           abline(h = 0, col = "gray50", lty = 2)
         }
    )
  
    plot.coordinates <- par("usr")
    gsea.enrichment.score <- Results[Results$Pathway == "DGSEA",]$KS
    gsea.normalized.enrichment.score <- Results[Results$Pathway == "DGSEA",]$KS_Normalized
    gsea.p.value <- signif(Results[Results$Pathway == "DGSEA",]$`Frequency of Random KS Scores > KS Score observed`, digits = 3)
    if(gsea.enrichment.score < 0) {
      text(length(data_in2[,cells[u]]) * 0.01, plot.coordinates[3] * 0.98,
           paste("DGSEA","\np-value:", gsea.p.value, "\nES:",
                 gsea.enrichment.score, "\nNES:",
                 gsea.normalized.enrichment.score), adj = c(0, 0))
      text(length(data_in2[,cells[u]]) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.02),
           paste("Glycolysis"), col = "red", adj = c(1,1))
      text(length(data_in2[,cells[u]]) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.10),
           paste("OxPhos"), col = "blue", adj = c(1,1))
      
    } else {
      text(length(data_in2[,cells[u]]) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
           paste("DGSEA","\np-value:", gsea.p.value, "\nES:",
                 gsea.enrichment.score, "\n NES:",
                 gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
      text(length(data_in2[,cells[u]]) * 0.01, plot.coordinates[3] * 0.90,
           paste("Glycolysis"), col = "red", adj = c(0,0))
      text(length(data_in2[,cells[u]]) * 0.01, plot.coordinates[3] * 0.6,
           paste("OxPhos"), col = "blue", adj = c(0,0))
    }
    
    par(mar = c(0, 5, 0, 2))
    plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
         ylab = "", xlim = c(1, length(data_in2[,cells[u]])))
    abline(v = gsea.hit.indices.Glycolysis, lwd = 0.75, col = "red")
    abline(v = gsea.hit.indices.OxPhos, lwd = 0.75, col = "blue")
    
    
    metric.range <- c(min(data_in2[,cells[u]]),max(data_in2[,cells[u]]))
    
    par(mar = c(5, 5, 0, 2))
    rank.metric <- rle(round(data_in2[,cells[u]], digits = 2))
    plot(data_in2[,cells[u]], type = "n", xaxs = "i",
         xlab = "Rank in ordered gene list", xlim = c(0, length(data_in2[,cells[u]])),
         ylim = metric.range, yaxs = "i",
         ylab = "Ranking Metric",
         panel.first = abline(h = seq(metric.range[1] / 2,
                                      metric.range[2] - metric.range[1] / 4,
                                      metric.range[2] / 2), col = "gray95", lty = 2))
    
    barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
            xlab = "Rank in ordered gene list", xlim = c(0, length(data_in2[,cells[u]])),
            ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
            #ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), 
            space = 0, add = TRUE)
    box()
    dev.print(pdf,file = paste(cells[u],"DGSEA.pdf"), height = 4, width = 4)
    dev.off()
  }
  par(def.par)
}

print(paste("Cell: ", u))

end.loop.time <- Sys.time()
total.loop.time <- end.loop.time - loop.time
print(paste("Time per cell:" , total.loop.time))
}

write.csv(Results.All.Cells, file = output.file.name, row.names = FALSE)

stopCluster(cl)

DGSEA.Results <- Results.All.Cells[Results.All.Cells$Pathway == "DGSEA",]

#uncomment these lines if you want GSEA results across all gene sets
#GSEA.Results <- Results.All.Cells[Results.All.Cells$Pathway != "DGSEA",]
#write.csv(GSEA.Results, file = paste("GSEA Only", output.file.name), row.names = FALSE)

### DGSEA Comparisons have Glyc, OxPhos, and DGSEA
DGSEA.Results <- Results.All.Cells[Results.All.Cells$Pathway == "DGSEA" |
                                    Results.All.Cells$Pathway == "Glycolysis.Pathway" | 
                                    Results.All.Cells$Pathway == "OxPhos.Pathway",]
write.csv(DGSEA.Results, file = paste("DGSEA Results", output.file.name), row.names = FALSE)


# This draft started by taking a chunk from Run_plastma_deconts_and_preds.R

# for each library, install it if needed
load.library <- function(lib){
    if (!require(lib, character.only = TRUE)) install.packages(lib, repos='http://cran.us.r-project.org')
    library(lib, character.only = TRUE)
}
load.library("splitstackshape")
load.library("tidyverse")
load.library("biomformat") 
load.library("vegan")
load.library("glmnet")
load.library("torch")
# library(microDecon)
# install_torch()

## Load packages ##
load.library("dplyr")
load.library("doMC")
load.library("tibble")
load.library("parallel")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

load.bioc.library <- function(lib){
    if (!require(lib, character.only = TRUE)) BiocManager::install(lib)
    library(lib, character.only = TRUE)
}
load.bioc.library("edgeR")
load.bioc.library("limma")
load.bioc.library("snm")
load.bioc.library("gbm")


sessionInfo()

#### read meta data ####

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

sampleSet = row.names(metadataPSMatchedDPQCFiltered)

message("Only include these ", length(sampleSet), " samples:")
message(paste(sampleSet, collapse = ", "))


#### find data files ####

# SCRuB files in ../results/data/decontaminate/SCRuB/trial_*/scrub_output.csv
# infiles = c(
#     dir("../results/decontamination/SCRuB", pattern="decontaminated.csv", full.names = T, recursive = T),
#     dir("../results/decontamination/SCRuB_noPlate", pattern="output.csv", full.names = T, recursive = T))

decontaminationFolder = "../results/decontamination"
infiles = dir(decontaminationFolder, pattern="decontaminated.csv", full.names = T, recursive = T)

#### methods ####

div_ind <- 'shannon'
message("Using diversity index: ", div_ind)

calcDiversity <- function(table, infile){
    message("Calculating diversity...")
    diver = diversity(table, index=div_ind)
    
    outfile = sub("decontaminated.csv", "decontaminated_shannon.csv", infile)
    message("Saving file: ", outfile)
    write.csv(diver, outfile)
    return(diver)
}

vsnm <- function(qcData, infile){
    imgfile = sub("decontaminated.csv", "decontaminated_voomMeanVarTrend.pdf", infile)
    message("Any images will be saved to: ", imgfile)
    pdf(imgfile)

    message("Running vsnm...")
    # above this line is stuff Ivory added.
    
    numCores <- parallel::detectCores()
    registerDoMC(cores=numCores)
    
    qcMetadata <- metadataPSMatchedDPQCFiltered # ADAPT THIS AS NEEDED
    # qcData <- qqq # ADAPT THIS AS NEEDED
    
    # Set up design matrix
    covDesignNorm <- model.matrix(~0 + disease_type_consol + #randomized_disease_type_consol + #disease_type_consol +
                                      host_age + # host_age should be numeric
                                      sex, # sex should be a factor
                                  data = qcMetadata)
    
    # Check row dimensions
    dim(covDesignNorm)[1] == dim(qcData)[1]
    
    # The following corrects for column names that are incompatible with downstream processing
    colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
    
    # Set up counts matrix
    counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
    
    # Quantile normalize and plug into voom
    dge <- edgeR::DGEList(counts = counts)
    vdge <<- limma::voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE, 
                  normalize.method="quantile")
    
    # List biological and normalization variables in model matrices
    bio.var <- model.matrix(~disease_type_consol, #randomized_disease_type_consol, #disease_type_consol,
                            data=qcMetadata)
    
    adj.var <- model.matrix(~host_age +
                                sex,
                            data=qcMetadata)
    
    colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
    colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
    print(dim(adj.var))
    print(dim(bio.var))
    print(dim(t(vdge$E)))
    print(dim(covDesignNorm))
    
    snmDataObjOnly <- snm::snm(raw.dat = vdge$E, 
                          bio.var = bio.var, 
                          adj.var = adj.var, 
                          rm.adj=TRUE,
                          verbose = TRUE,
                          diagnose = TRUE)
    snmData <<- t(snmDataObjOnly$norm.dat)
    
    ## above: vsnm function exactly as taken from Run_plasma_deconts_and_preds.R
    ## below: my addition.
    dev.off()
    outfile = sub("decontaminated.csv", "decontaminated_vsnm.csv", infile)
    message("Saving file: ", outfile)
    
    commentLines = readLines(infile, n=30) %>% grep(pattern="#METHODS", value = TRUE)
    writeLines(commentLines, outfile)
    suppressWarnings({
      write.table(snmData, file=outfile, append=TRUE, quote=F, sep=",")
    })
}

#### do for each file ####

shortNames = sub(decontaminationFolder, "", infiles)
shortNames = sub("trial_", "t", shortNames)
diversitySummary = data.frame(row.names=sampleSet)

for (infile in infiles){
    message("Processing file: ", infile)
    indata = read.csv(infile, row.names = 1, comment.char = "#")
    diver = calcDiversity(indata[sampleSet, ], infile)
    diversitySummary[,infile] = diver[sampleSet]
    normed <- vsnm(indata[sampleSet, ], infile)
}


#### diversity summary ####

suppressWarnings({
    dir.create("../results/shannon_diversity", recursive = TRUE) 
})

shortNames = sub(decontaminationFolder, "", names(diversitySummary))
shortNames = sub("trial_", "t", shortNames)
names(diversitySummary) = shortNames

summaryFile = "../results/shannon_diversity/shannon_summary.csv"
message("Saving diversity summary info to: ", summaryFile)
write.csv(diversitySummary, summaryFile)

heatFile = "../results/shannon_diversity/shannon_heatmap.png"
message("Saving diversity heatmap: ", heatFile)
png(heatFile)
heatmap(as.matrix(diversitySummary))
dev.off()


boxFile = "../results/shannon_diversity/shannon_heatmap.png"
message("Saving diversity boxplot: ", boxFile)
png(boxFile)
boxplot(diversitySummary, las=2, col="skyblue")
title("Shannon Diversity")
text(names(diversitySummary), x=c(1:4)-.2, y=2, srt=90)
dev.off()



# example
# Processing file: ../results/decontamination/SCRuB_noPlate/trial_2/scrub_decontaminated.csv
# Saving file: ../results/decontamination/SCRuB_noPlate/trial_2/scrub_decontaminated_diversity.csv
# Saving file: ../results/decontamination/SCRuB_noPlate/trial_2/scrub_decontaminated_vsnm.csv


# scrubbed_normalized <- vsnm(scrub_df[row.names(metadataPSMatchedDPQCFiltered), ])
# 
# microdecon_normalized <- vsnm(microdec_df[row.names(metadataPSMatchedDPQCFiltered), ])
# 
# # scrubbed_normalized <- vsnm(scrub_df[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ), ] )
# # 
# raw_inp <- unique_samps[row.names(metadataPSMatchedDPQCFiltered),]
# raw_inp <- raw_inp[, colSums(raw_inp)>500]
# # raw_inp <- group_to_genus(raw_inp)
# raw_normalized <- vsnm(raw_inp)
# 
# dec_normalized <- vsnm(decontammed_data[row.names(metadataPSMatchedDPQCFiltered), ])#
# # dec_normalized <- vsnm(decontammed_data[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# 
# dec_standard_normalized <- vsnm(decontammed_data_standard[row.names(metadataPSMatchedDPQCFiltered), ])#
# 
# # dec_standard_normalized <- vsnm(decontammed_data_standard[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# dec_lb_normalized <- vsnm(decontammed_data_low_bm[row.names(metadataPSMatchedDPQCFiltered), ])#
# # dec_lb_normalized <- vsnm(decontammed_data_low_bm[paste0( 'X' , row.names(metadataPSMatchedDPQCFiltered) ),])
# 
# 
# 
# restrictive_normalized <-  vsnm( restrictive[row.names(metadataPSMatchedDPQCFiltered) %>% as.character(), ])#

message("Done!")

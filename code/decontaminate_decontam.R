
# This draft started by taking a chunk of Run_plastma_deconts_and_preds.R

if (!require(vegan)) install.packages(vegan, repos='http://cran.us.r-project.org')
library(vegan)
if (!require(splitstackshape)) install.packages("splitstackshape", repos='http://cran.us.r-project.org')
library(splitstackshape)
# required for SCRuB
if (!require(tidyverse)) install.packages("tidyverse", repos='http://cran.us.r-project.org')
library(tidyverse)
if (!require(decontam)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("decontam")
}
library(decontam)
# if (!require(glmnet)) install.packages("glmnet", repos='http://cran.us.r-project.org')
# library(glmnet)
# if (!require(torch)) install.packages("torch", repos='http://cran.us.r-project.org')
# library(torch)
# install_torch()
# if (!require(SCRuB)) devtools::install_github("shenhav-and-korem-labs/SCRuB")
# library(SCRuB)

sessionInfo()

scriptResultsDir = "../results/decontamination/decontam"
suppressWarnings({
    dir.create(scriptResultsDir, recursive = TRUE)
})

#### read data ####

dec_ind <- read.csv('../data/Fig4_plasma/Metadata-Plasma-For-Decontam-With-Negative-And-Positive-Controls.csv')

full_df = read.csv("../data/inGit/dataset/Nejman-et-al-2020/full_df.csv", row.names = 1)
message("full_df.csv has dimentions:")
dim(full_df)

metadata <- read.csv('../data/inGit/dataset/Nejman-et-al-2020/metadata.csv', comment.char = "#")
row.names(metadata) = metadata$X.SampleID
message("metadata has dimentions:")
dim(metadata)


# These lines were not needed when using R version 4.2.2 (2022-10-31), 
# but without them I hit errors when using R version 3.6.3 (2020-02-29)
metadata$sample_plate = as.character(metadata$sample_plate)
metadata$sample_type = as.character(metadata$sample_type)
metadata$sample_well = as.character(metadata$sample_well)

sample_intersect = intersect(row.names(metadata), 
                             row.names(full_df))

full_df = full_df[sample_intersect, ]
metadata = metadata[sample_intersect, ]

#### run through decontam ####

numTrials = 10
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) numTrials = args[1]

for (trialNumber in 1:numTrials){
    message("Trial ", trialNumber, " of ", numTrials)
    
    seed = trialNumber * 29
    message(paste('Using seed:', seed))
    set.seed(seed)


    control_types = c('control blank library prep', 'control blank DNA extraction')
    metadata$is_control = metadata$sample_type %in% control_types
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # 
    # plate_meta_split = split(metadata, f=metadata$sample_plate)
    # message("SCRuB the ", length(plate_meta_split), " plates separately.")
    
    decontam_out_list = list()
    
    # for (plate in names(plate_meta_split)){
        # message("SCRuB-ing plate ", plate)
    #     plate_meta = plate_meta_split[[plate]]
    #     plate_data = full_df[row.names(plate_meta),]
    #     
    #     # limit metadata to the EXACT three columns permitted.
    #     plate_meta = plate_meta[,c("is_control", "sample_type", "sample_well")]
    #     
    #     
    #     # when I passed a table with no controls, got error: "Error in rowSums(.) : 'x' must be numeric"
    #     # so skip any plate with 0 controls
    #     numControls = sum(plate_meta$is_control)
    #     message("Plate has controls: ", numControls)
    #     message(paste(row.names(plate_meta)[plate_meta$is_control], collapse=", "))
    #     
    #     if (numControls > 0 ){
    #         scrub_output = SCRuB::SCRuB(plate_data, metadata = plate_meta)
    #         
    #         # decontaminated data
    #         scrub_out_list[[plate]] <- scrub_output$decontaminated_samples
    #         
    #     }else{
    #         message("Skipping plate")
    #         scrub_out_list[[plate]] <- plate_data
    #     }
    # }
    # scrub_df = do.call("rbind", scrub_out_list)
    # 
    # message("scrub_df has dimensions:")
    # dim(scrub_df)
    # 
    # message("Filtering down to columns with sums > 500")
    # scrub_df <- scrub_df[, (colSums(scrub_df) > 500) %>% which]
    # 
    # message("scrub_df has dimensions:")
    # dim(scrub_df)
    
    
    
    
    
    
    
    
    
    
    
    #### write methods and results ####
    METHODS=c(decontaminationTool="SCRuB-byPlate", 
              blankType=paste(control_types, collapse=","), 
              numberBlanks=sum(metadata$is_control), 
              trialNumber=trialNumber, 
              seed=seed)
    
    resDir = file.path(scriptResultsDir, paste0("trial_", trialNumber))
    suppressWarnings({dir.create(resDir)})
    file = file.path(resDir, "decontam_decontaminated.csv")
    message("Saving cleaned data as: ", file)
    
    # write methods and results to the same file
    commentLines = paste0("#METHODS ", paste(names(METHODS), METHODS, sep="="))
    writeLines(commentLines, file)
    suppressWarnings({
      write.table(scrub_df, file=file, append=TRUE, quote=F, sep=",")
      })
    
}

message("Done!")

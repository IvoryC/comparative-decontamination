
# This draft started by taking the first 334 lines of Run_plastma_deconts_and_preds.R

library(tidyverse)

sessionInfo()

resDir = "../results/decontamination/raw/trial_0"
suppressWarnings({
    dir.create(resDir, recursive = TRUE)
})

#### read data ####

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

counts = full_df[sample_intersect, ]
meta = metadata[sample_intersect, ]

#### decontamination step ####

# nope

#### general process ####
# just to match any filtering and/or formating that was done in the decontamination scripts.

message("counts has dimensions:")
dim(counts)

message("Filtering down to columns with sums > 500")
counts_df <- counts[, (colSums(counts) > 500) %>% which]

message("counts_df has dimensions:")
dim(counts_df)



#### write methods and results ####
METHODS=c(decontaminationTool="no-decontamination", 
          blankType=NA, 
          numberBlanks=0, 
          trialNumber=0, 
          seed=NA)

file=file.path(resDir, "not_decontaminated.csv")
message("Saving non-scrubbed (raw) data as: ", file)

# write methods and results to the same file
commentLines = paste0("#METHODS ", paste(names(METHODS), METHODS, sep="="))
writeLines(commentLines, file)
suppressWarnings({
  write.table(counts_df, file=file, append=TRUE, quote=F, sep=",")
})

message("Done!")

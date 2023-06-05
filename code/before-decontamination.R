
# This draft started by taking the first several lines of Run_plastma_deconts_and_preds.R
# Assumes the working dir is capsule<capsule_id>/code/

library(tidyverse)
library(biomformat) 

sessionInfo()

#### a place to save data ####
suppressWarnings({
    dir.create("../results/before-decontamination", recursive = TRUE) 
})

#### as is files ####

# dec_ind <- read.csv('../data/Fig4_plasma/Metadata-Plasma-For-Decontam-With-Negative-And-Positive-Controls.csv')
# metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

# not sure if this gets used at all
# snmDataKrakenCFDecontamDPQC <-read.csv('../data/Fig4_plasma/Kraken-Plasma-Voom-SNM-Age-And-Sex-Data.csv', row.names=1)


#### full_df ####

data <- read_biom(biom_file = '../data/Fig4_plasma/136205_47212_analysis_Metagenomic_Woltkav011Databaseoptgenomeqiitadbsqpwoltkarep200Rep200BIOMnonebiom.biom')

taxa_names <- c()

for(i in 1:length(data$rows)) taxa_names <- c(taxa_names, data$rows[[i]]$id)

samp_names <- c()

for(i in 1:length(data$columns)) samp_names <- c(samp_names, data$columns[[i]]$id)

full_df <- matrix( unlist(data$data), byrow=TRUE, nrow=length(data$data) ) 
colnames(full_df) <- samp_names
row.names(full_df) <- taxa_names
full_df <- full_df %>% t()

write.csv(full_df,
          "../results/before-decontamination/full_df.csv")

#### metadata ####

metadata <- read.csv('../data/Fig4_plasma/47212_47212_analysis_mapping.txt', sep='\t')

remov_left <- function(x, n){
    substr(x, n, nchar(x))
}

unique_metadata <- metadata %>%
    filter(as.character(X.SampleID) %in% row.names(full_df)) %>%
    mutate(general_id = remov_left(as.character(X.SampleID), 7)) %>%
    group_by(general_id) %>% sample_n(1)

#### unique_samps ####

unique_samps <- full_df[as.character(unique_metadata$X.SampleID), ]

write.csv(unique_samps,
          "../results/before-decontamination/unique_samps.csv")

#### not-used--final_df ####
final_df <- full_df[which(row.names(full_df) %>% str_detect('Control') == F ), ] %>%
    rbind( full_df[metadata %>% filter(str_detect(X.SampleID, 'Control'), 
                                       #X.SampleID %>% str_detect('Control'),
                                       X.SampleID %>% str_detect('12691') ) %>% 
                       pull(X.SampleID) %>% unique() %>% as.character, ] ) 

#### removed_copies ####

removed_copies <- metadata %>% filter( ( str_detect(X.SampleID, 'Control') == F ) ) %>% #&(str_detect(description, 'HIV') == F ) ) %>%
    rbind( metadata %>% filter(str_detect(X.SampleID, 'Control'), 
                               X.SampleID %>% str_detect('12691') ) )

write.csv(removed_copies,
          "../results/before-decontamination/removed_copies.csv")

#### unique_metadata ####

unique_metadata <- as.data.frame(unique_metadata)
row.names(unique_metadata) <- unique_metadata[, 'X.SampleID'] %>% as.character()

write.csv(unique_metadata,
          "../results/before-decontamination/unique_metadata.csv")

#### well_dists ####

# well_dists <- unique_metadata %>%
#     mutate(well_loc= sample_well %>% substr(1,1) %>% sapply( function(x) which(LETTERS==x)[1]), 
#            indices_loc = sample_well%>% substr(2,3) %>% as.integer ) %>%
#     select(well_loc, indices_loc) %>%
#     dist(method = 'euclidean') %>% as.matrix()
# 
# well_dists = round(well_dists, digits = 2)
# 
# write.csv(well_dists,
#           "../results/before-decontamination/well_dists.csv")

message("Done!")
date()
sessionInfo()


#----delete these notes sometime soon------#
#### objects we need
# unique_samps
# unique_metadata
# X.SampleID #not sure if this is an obj...actually I don't think it is.
# well_dists

## additonal for microdecon
# metadataPSMatchedDPQCFiltered

## additional for decontam
# dec_ind 
# dec_ind$decontam_prevalence_use
# removed_copies
# removed_copies$well_conc
# full_df

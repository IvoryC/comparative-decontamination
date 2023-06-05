
#### libraries ####

library(vegan)
if (!require(ape)) install.packages("ape")
library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)

# install.packages("remotes")
# remotes::install_github("nchlis/pca.utils")
# library("pca.utils")
#
# not sure why I have problems loading this library. This is the function we use from it:
string = "Plesae cite github('nchlis/pca.utils') see: https://rdrr.io/github/nchlis/pca.utils/man/project_pca.html"
message(string); print(string)
project_pca <- function (Xnew = NULL, pc = NULL) 
{
  return(scale(Xnew, pc$center, pc$scale) %*% pc$rotation)
}


sessionInfo()

#### categories of interest ####

# default
categories = c("Control", "SKCM") 

# command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  args2 = args %>% strsplit(split="=")
  key = sapply(args2, function(x) x[[1]])
  argVals = sapply(args2, function(x) x[[2]])
  names(argVals) = key
  if ("categories" %in% key) {
    categories = strsplit(argVals["categories"], split=",")[[1]]
    message("Using command line arg value for 'categories':", args)
  }else{
    message("Using default value for 'categories'.")
  }
}
message("In this run, we will aim to distinguish between 'categories': ", paste(categories, collapse = " vs "))

#### read meta data ####

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

#### find input data files ####

decontaminationFolder = "../results/decontamination"
inputPattern = "_vsnm.csv"
infiles = dir(decontaminationFolder, pattern=inputPattern, full.names = T, recursive = T)
message("Found ", length(infiles), " input files.")

for (infile in infiles){
  
  message("Reading data from file: ", infile)
  data = read.csv(infile, row.names = 1, comment.char = "#")
  
  # save output with matching sub-directory structure
  predictionFolderBase = "../results/prediction"
  indir = dirname(infile)
  predictionFolderBase = sub(decontaminationFolder, predictionFolderBase, indir)
  predictionFolder = file.path(predictionFolderBase, paste(categories, collapse = "_vs_"))
  suppressWarnings(dir.create(predictionFolder, recursive = T))
  message("Using output dir: ", predictionFolder)
  
  # Save all images for this input file to one pdf
  plotFile = file.path(predictionFolder, sub(inputPattern, "_leave-1-out-plots.pdf", basename(infile)))
  pdf(plotFile)
  
  fileSummary = data.frame()
  
  #### select test data ####
  
  hasCat = metadataPSMatchedDPQCFiltered$disease_type_consol %in% categories
  nameHasCat = row.names(metadataPSMatchedDPQCFiltered)[which(hasCat)]
  data = data[nameHasCat,]
  
  for (test1 in nameHasCat){
    
    message("Leave out sample ", test1, " as the test.")
    
    train = filter(data, rownames(data) != test1)
    
    #### train ####
    
    ##### prcomp #####
    
    res = stats::prcomp(train)
    
    axes = res$x
    
    pcDF = merge(axes, metadataPSMatchedDPQCFiltered, by=0)
    row.names(pcDF) = pcDF$Row.names
    pcDF = pcDF %>% select(-Row.names)
    
    #### classic biplots
    #
    # plot1 = ggplot2::ggplot(data=pcDF[nameHasCat,]) +
    #     geom_point(mapping = aes(x=PC1, y=PC2, color = disease_type)) +
    #     ggtitle("prcomp")
    # plot1
    # 
    # plot2 = ggplot2::ggplot(data=pcDF[nameHasCat,]) +
    #     geom_point(mapping = aes(x=PC3, y=PC4, color = disease_type)) +
    #     ggtitle("prcomp")
    # plot2
    
    axisSet = colnames(axes)
    if (length(axisSet) > 12) axisSet = axisSet[1:12]
    
    d4 = pcDF %>% 
      filter(disease_type_consol %in% categories) %>% 
      select(axisSet, disease_type_consol) %>% 
      gather("value", key="PC", -disease_type_consol) %>%
      mutate(axisNum = as.numeric(gsub("PC", "", PC))) %>%
      arrange(axisNum) 
    
    # reset the factor levels so that the plots appear in order
    d5 = d4 %>%
      arrange(axisNum) %>%
      mutate(PC = factor(PC, levels=unique(PC))) %>%
      select(-axisNum)
    
    plot4 = ggplot2::ggplot(data=d5) +
      geom_boxplot(mapping = aes(y=value, x=disease_type_consol, color = disease_type_consol)) +
      ggtitle("prcomp values split by sample type", subtitle = paste0("made without sample: ", test1)) +
      facet_wrap(~PC, nrow=3)
    #plot4
    
    ##### prcomp choose best #####
    # Choose the single best axis for separating the data.
    axisPs = sapply(axisSet, function(ax) {
      t.test(data=d4 %>% filter(PC == ax), 
             value ~ disease_type_consol)$p.value
    })
    bestAxis = names(axisPs)[which(axisPs == min(axisPs))]
    message("Axis ", bestAxis, " is the best at separating ", paste(categories, collapse=" from "))
    
    # highlight the best axis in the plot
    plot6 = plot4 +
      geom_rect(data = subset(d5, PC == bestAxis), 
                fill = NA, colour = "red", 
                xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)
    
    # save plot
    #plotFile = sub(inputPattern, "_prcomp-boxplot.png", infile)
    #ggsave(plotFile, plot=plot6)
    # printed images are directed to pdf
    print(plot6)
    
    
    #### test sample ####
    
    ###### sanity check ######
    
    # from https://rdrr.io/github/nchlis/pca.utils/man/project_pca.html
    # install.packages("remotes")
    # remotes::install_github("nchlis/pca.utils")
    #
    # example: pca.utils::project_pca(Xnew = NULL, pc = NULL)
    train1 = row.names(train)[1]
    pro.check = project_pca(Xnew = data[train1,], pc = res)
    sanity = pro.check[,"PC1"] == res$x[,"PC1"][train1]
    # TRUE #--yay!
    if (!sanity) warning("Oh no! we might not be sane! Check the project_pca method.")
    
    ### the project_pca function is very straight forward:
    # function (Xnew = NULL, pc = NULL) 
    # {
    #     return(scale(Xnew, pc$center, pc$scale) %*% pc$rotation)
    # }
    
    ###### actual test ######
    
    pro.test = project_pca(Xnew = data[test1,], pc = res)
    pro.test.val = pro.test[,bestAxis]
    
    
    claims = sapply(categories, function(cat){
      catVals = pcDF %>% 
        filter(disease_type_consol == cat) %>%
        select(bestAxis) %>% 
        unlist()
      t.test(catVals, pro.test.val, var.equal = T)$p.value
    })
    claims = claims[order(claims, decreasing = T)]
    
    # which category had the highest p? that's our prediction.
    predict = names(claims)[1]
    confidenceP = claims[2]
    
    realCategory = metadataPSMatchedDPQCFiltered[test1,"disease_type_consol"]
    isCorrect = predict == realCategory
    correctOrNot = ifelse(isCorrect, "CORRECT", paste0("WRONG (is actually ", realCategory, ")"))
    
    message("We predict (with confidence p=", round(confidenceP,3),") that sample [", test1, "] belongs in category [", predict, "], which is ", correctOrNot, "." )
    
    
    
    fileSummary = rbind(fileSummary, 
                        c(sampleID=test1, actual=realCategory, prediction=predict, confidenceP=confidenceP))
    
  } #done test with sample
  
  names(fileSummary) = c("sampleID", "actual", "prediction", "confidenceP")
  levels(fileSummary$prediction) = levels(fileSummary$actual)
  fileSummary$isCorrect = fileSummary$actual == fileSummary$prediction
  summr = paste("Predictions were correct for", sum(fileSummary$isCorrect), "of", nrow(fileSummary), "leave-1-out samples.")
  message(summr)
  
  tableFile = file.path(predictionFolder, sub(inputPattern, "_leave-1-out-prediction-summary.txt", basename(infile)))
  message("Saving file: ", tableFile)
  
  # append METHODS comments about this process to the comments that were in the infile
  METHODS_KEY="#METHODS "
  commentLinesIN = readLines(infile, n=30) %>% grep(pattern=METHODS_KEY, value = TRUE)
  commentLinesOUT = c(commentLinesIN, 
                      paste0(METHODS_KEY, "predictionCategories=", paste0(categories, collapse=",")))
  generalComment = paste("#", summr)
  writeLines(tableFile, text = c(generalComment, commentLinesOUT))
  
  # write the main table to the same file
  write.table(fileSummary, file=tableFile, quote=F, row.names = F, sep="\t", append=T)
  rm(fileSummary, commentLines, additionalComment)
  
  dev.off()
  
  message("Done with file ", infile, ".")
  
}

message("Done!")

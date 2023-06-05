
# summarize performance

#### libraries ####

if (!require(ggalluvial)) install.packages("ggalluvial")
library(ggalluvial)

#used earlier in pipeline
library(ggplot2)
library(dplyr)
library(tidyr)

#### find data ####

predictionFolder = "../results/prediction"
inputPattern = "prediction-summary.txt"
infiles = dir(predictionFolder, pattern=inputPattern, full.names = T, recursive = T)
message("Found ", length(infiles), " prediction files.")

resultsDir = "../results/accuracySummary"
dir.create(resultsDir, recursive = TRUE)

#### read and record accuracy ####

acc = data.frame()
preds = list()
predAll = read.delim2(infiles[1], comment.char = "#")[c("sampleID", "actual")]
allMethods = c()
shortFileId = paste0("file", 1:length(infiles))
names(infiles) = shortFileId
for (i in 1:length(infiles)){
  infile = infiles[i]
  commentLines = readLines(infile, n=30)  %>%
    grep(pattern="#METHODS ", value = TRUE) %>%
    gsub(pattern="#METHODS ", replacement="") %>%
    strsplit(split="=")
  methods = sapply(commentLines, "[[", 2)
  names(methods) = sapply(commentLines, "[[", 1)
  allMethods = rbind(allMethods, methods)
  nBlanks = methods["numberBlanks"] %>% as.numeric()
  df = read.delim2(infile, comment.char = "#")
  newRow = c(accuracy=sum(df$isCorrect)/nrow(df), 
             numCorrect=sum(df$isCorrect),
             methods)
  keys = union(names(newRow), names(acc))
  if (nrow(acc)==0) acc = data.frame(t(newRow))
  else {
    newcols = setdiff(keys, names(acc))
    for (nc in newcols) acc[nc] = NA
    acc = rbind(acc[,keys], newRow[keys])
  }
  preds[[infile]] = df
  df = df[,c("sampleID",  "prediction")]
  names(df)[2] = shortFileId[i]
  predAll = merge(predAll, df, by="sampleID")
  # names(df)[2] = paste0("prediction_", nBlanks, "_t", methods["trialNumber"])
  # mtemp = merge(predAll, df, by="sampleID")
  # predAll = mtemp
}
names(preds) = infiles

acc = acc %>% mutate_at(c("accuracy", "numCorrect", "numberBlanks"), as.numeric)

allMethods=data.frame(allMethods)
row.names(allMethods) = shortFileId
allMethods$numberBlanks = as.numeric(allMethods$numberBlanks)

#### scatter plot  and boxplot ####
s0 = ggplot(acc) +
  geom_hline(yintercept=c(0,1), color="white", linewidth=2) +
  #geom_point(aes(x=numberBlanks, y=accuracy, color=decontaminationTool, shape=decontaminationTool)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1)) +
  scale_x_continuous(labels=unique(acc$numberBlanks), breaks = unique(acc$numberBlanks)) +
  xlab("number of blanks for decontamination") +
  ylab("accuracy") +
  ggtitle("Accuracy by classifaction task") +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 25)) +
  theme(axis.text = element_text(size = 15))+
  facet_wrap("predictionCategories", nrow=2)

s1 = s0 + geom_point(aes(x=numberBlanks, y=accuracy, color=decontaminationTool, shape=decontaminationTool))

imgFile1 = file.path(resultsDir, "accuracy-scatterplot.png")
message("Saving image to: ", imgFile1)
ggsave(imgFile1, s1, device = png, units="in", width=6, height = 5)

s2 = s0 + geom_boxplot(aes(x=numberBlanks, y=accuracy, group=numberBlanks, color=numberBlanks))

imgFile2 = file.path(resultsDir, "accuracy-boxplot.png")
message("Saving image to: ", imgFile2)
ggsave(imgFile2, s2, device = png, units="in", width=6, height = 5)


#### heatmap ####

if (FALSE){
  
  predStyles = gather(predAll, key="file", value="prediction", -sampleID)
  predStyles2 = merge(predStyles, allMethods, by.x="file", by.y=0, all.x = T)
  predStyles2 = predStyles2 %>% group_by(numberBlanks)
  
  tile = ggplot(predStyles2 %>% group_by(decontaminationTool, blankType, numberBlanks), 
                aes(x=file, y=sampleID, fill=prediction)) +
    geom_tile() +
    facet_wrap(vars(numberBlanks),nrow = 1, scales="free_x") +
    ggtitle(paste(unique(preds[[1]]$prediction), collapse=" vs "))
  tile
  
  imgFileTile = file.path(resultsDir, "sample-prediction-collorgrid.png")
  message("Saving image to: ", imgFileTile)
  ggsave(imgFileTile, tile, device = png, units="in", width=10, height = 10)
  
}else{
  message("Ignoring the heatmap section.")
}

#### alluvial ####

if(FALSE){
  rawFile = "../results/prediction/raw/trial_0/not_decontaminated_leave-1-out-prediction-summary.txt"
  predNone = read.delim2(rawFile, comment.char = "#")
  predNone = predNone[,c("sampleID",  "prediction")]
  names(predNone)[2] = "predictionRaw"
  predNone$numberBlanks = 0
  
  predAll = preds[[3]][,c("sampleID",  "actual", "prediction")]
  names(predAll)[3] = "prediction90"
  predHalf = preds[[2]][,c("sampleID",  "prediction")]
  names(predHalf)[2] = "prediction45"
  predFew = preds[[1]][,c("sampleID",  "prediction")]
  names(predFew)[2] = "prediction23"
  
  m0 = merge(predNone, predAll, by="sampleID")
  m1 = merge(m0, predHalf, by="sampleID")
  m2 = merge(m1, predFew, by="sampleID")
  
  
  p2 = m2 %>% 
    group_by(predictionRaw, actual, prediction90, prediction45, prediction23) %>% 
    summarise(Freq = n())
  
  ## helpful pages
  # https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
  # https://felixfan.github.io/ggplot2-remove-grid-background-margin/
  
  a1 = ggplot(p2,
              aes(y = Freq, 
                  axis1=predictionRaw, axis2 = actual, axis3 = prediction90, axis4=prediction45, axis5=prediction23,
                  #axis1 = actual, axis2 = prediction90, axis3=prediction45, axis4=prediction23,
                  fill = actual)) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum))) + 
    theme(axis.text=element_text(size=12)) +
    scale_x_continuous(breaks = 1:5, labels = c("raw", "truth", "90 blanks", "45 blanks", "23 blanks")) + 
    ggtitle("prediction of leave-1-out samples")
  a1    
  
  imgFile2 = file.path(resultsDir, "accuracy-alluvial.png")
  message("Saving image to: ", imgFile2)
  ggsave(imgFile2, a1, device = png, units="in", width=9, height = 7)
  
}else{
  message("Ignoring alluvial plot section.")
}

message("Done!")

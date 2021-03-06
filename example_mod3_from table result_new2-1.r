
library("ggplot2")
library("dplyr")
library("parallel")
library("NBPSeq")
library("qvalue")
library("RNAsense")
library("SummarizedExperiment")

#important! the table for reading should have the following col.names: gene/genename(geneID and symble)/samples
#sample name following genotype(big letter)_rep(small letter)_time(1:8)
normalizedReads <- read.table("final_reads_for GEO.txt",header=TRUE,sep = "\t",fill = TRUE)#read normalized reads
hpf<-c(2.5,3,3.5,4,4.5,5,5.5,6)
genotype<-c("WT","MZsox","MZspg","MZsoxspg")
startTime <- Sys.time() # get start time
analyzeConditions <- c("MZspg","WT") #rename groups names
ini <- TRUE

# mydata is data in standard format
thCount <- 100  
nrcores <- 6
FC = 2 #fold change
pValueSwitch = 0.15 # p value for switchlists
pValueFC=0.01 #p value for the fold change
experimentStepDetection = "WT" #switch list from which time series
#times <- list(c(2.5,3,3.5,4,4.5,5,5.5,6))
times <- c(2.5,3,3.5,4,4.5,5,5.5,6)

genotype1<-paste("^", analyzeConditions[1],"_",sep = "")
genotype2<-paste("^", analyzeConditions[2],"_",sep = "")

reads1<-normalizedReads[grepl(pattern = genotype1, names(normalizedReads))]
reads2<-normalizedReads[grepl(pattern = genotype2, names(normalizedReads))]
genenames<-normalizedReads[grepl(pattern = "name", names(normalizedReads))]
mydata <- cbind(genenames,reads1,reads2)

#cut with thCount
vec2Keep <- which(vapply(1:dim(mydata)[1],function(i) !Reduce("&",mydata[i,3:dim(mydata)[2]]<thCount), c(TRUE)))
mydata <- mydata[vec2Keep,]
genenames <- genenames[vec2Keep,]

#sort data
mydata<-arrange(mydata,name)

#transfer data to summarizedexperiment class  
colnames <- names(mydata)[-c(1,2)]
mydataSE <- SummarizedExperiment(assays=SimpleList(mydata=as.matrix(mydata[,-c(1,2)])),
                                 rowData = DataFrame(name=mydata$name, genename=mydata$genename),
                                 colData = DataFrame(condition = do.call(rbind, strsplit(colnames, "_"))[,1],
                                                     time = do.call(rbind, strsplit(colnames, "_"))[,2],
                                                     replicate = do.call(rbind, strsplit(colnames, "_"))[,3],
                                                     row.names=colnames))

#getFC
resultFC <- RNAsense:::getFC(dataset = mydataSE, myanalyzeConditions = analyzeConditions, cores = nrcores, mytimes=times)
head(resultFC)

#write.table(resultFC, "Foldchange.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "", dec = ".", row.names = FALSE, qmethod = c("escape", "double"))


# generate and output vulcano plot
VULCANO<-ggplot(subset(resultFC, FCdetect!="none" & time==3 & logFoldChange**2 < 12**2), 
                aes(x=logFoldChange, y=-log10(pValue), color=FCdetect)) + 
  xlab("log2(Fold Change)") + geom_point(shape=20) + 
  scale_color_dMod()
png(file = paste0(Sys.Date(), "_VULCANO_colorPlot.png"), width = 600, height = 500)
print(VULCANO)
dev.off()

# get Switch
#resultSwitch <- RNAsense:::getSwitch(dataset = mydataSE, experimentStepDetection = experimentStepDetection, 
#pValueSwitch = pValueSwitch, cores = nrcores, mytimes = times)
resultSwitch <- RNAsense:::getSwitchCorrect(dataset = mydataSE, experimentStepDetection = experimentStepDetection, 
                                            pValueSwitch = pValueSwitch, cores = nrcores, mytimes = times,chooseFirst = TRUE)
  

head(resultSwitch)

resultCombined <- RNAsense:::combineResults(resultSwitch, resultFC)
head(resultCombined)

# generate and output SSGS plot
pdf(file = paste0(Sys.Date(), "_SSGS_colorPlot.pdf"), width = 7, height = 6)
SSGS<-RNAsense:::plotSSGS(resultCombined, c(2.5,3,3.5,4,4.5,5,5.5,6), myanalyzeConditions = analyzeConditions, experimentStepDetection = experimentStepDetection)
print(SSGS)
dev.off()

RNAsense:::outputGeneTables(resultCombined)





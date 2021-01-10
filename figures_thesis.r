
library("ggplot2")
library("dplyr")
library("parallel")
library("NBPSeq")
library("qvalue")
library("RNAsense")
library("SummarizedExperiment")
library("dMod")
library("cowplot")

path <- "~/Dokumente/Doktorarbeit/Thesis/Figures/Zebrafish/"

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

# get Switch
#resultSwitch <- RNAsense:::getSwitch(dataset = mydataSE, experimentStepDetection = experimentStepDetection, 
#pValueSwitch = pValueSwitch, cores = nrcores, mytimes = times)
resultSwitch <- RNAsense:::getSwitchCorrect(dataset = mydataSE, experimentStepDetection = experimentStepDetection, 
                                            pValueSwitch = pValueSwitch, cores = nrcores, mytimes = times,chooseFirst = TRUE)
  

head(resultSwitch)

################################################

resultCombined <- RNAsense:::combineResults(resultSwitch, resultFC)
head(resultCombined)


# generate and output vulcano plot
VULCANO<-ggplot(subset(resultFC[seq(1,dim(resultFC)[1], by=1),], FCdetect!="none" & time==3 & logFoldChange**2 < 12**2), 
                aes(x=logFoldChange, y=-log10(pValue), color=FCdetect)) + 
  xlab("log2 (fold change)") + geom_point(shape=20) + theme_bw(base_size = 12) +
  ylab("-log10 (pvalue)") +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  scale_color_dMod(name="Fold change detection", breaks=c("MZspg<WT", "MZspg=WT", "MZspg>WT"), 
                   labels=c("MZspg < WT (downregulated)", "MZspg = WT (no effect)", "MZspg > WT (upregulated)"))
pdf(file = paste0(path, "VULCANO_colorPlot.pdf"), width = 7, height = 4)
print(VULCANO)
dev.off()


data_switchplot <- do.call(rbind, lapply(c("slc35a5", "sox3", "sox2", "dap", "thraa", "triob", "nanog", "sox19b"), function(mygene){
  cbind(subset(resultSwitch, genename==mygene),
        data.frame(time = c(2.5,3,3.5,4,4.5,5,5.5,6),
                   value = as.numeric(assays(mydataSE[which(rowData(mydataSE)$genename==mygene),
                                                      which(colData(mydataSE)$condition=="WT" &
                                                              colData(mydataSE)$replicate=="B1")])[[1]][1,])))
}))

data_switchplot <- do.call(rbind, lapply(unique(data_switchplot$genename), function(gene){
  sub <- subset(data_switchplot, genename==gene)
  sub$value <- sub$value/max(sub$value)
  sub$timepoint <- sub$timepoint+0.5
  ystep1 <- do.call(c, lapply(1:7, function(i){mean(sub$value[1:i])}))
  ystep2 <- do.call(c, lapply(1:7, function(i){mean(sub$value[(i+1):8])}))
  ymean <- mean(sub$value)
  cbind(sub, ystep1 = c(ystep1, NA), ystep2 = c(ystep2, NA), ymean = c(ymean, rep(NA, 7)))
}))

# data_switchplot1 <- subset(data_switchplot, genename%in%c("slc35a5", "sox3"))
# P1 <- ggplot(data_switchplot1, aes(x=time, y=value)) + #geom_smooth(color="grey", linetype="dashed", fill=NA) +
#   geom_segment(data = data_switchplot1, aes(x=2, xend=time+0.25, y=ystep1, yend=ystep1, color=as.factor(time)), linetype="dashed") +
#   geom_segment(data = data_switchplot1, aes(x=time+0.25, xend=time+0.25, y=ystep1, yend=ystep2, color=as.factor(time)), linetype="dashed") +
#   geom_segment(data = data_switchplot1, aes(x=time+0.25, xend=6.5, y=ystep2, yend=ystep2, color=as.factor(time)), linetype="dashed") +
#   geom_point() + facet_wrap(~genename, scales="fixed") + theme_dMod(base_size = 10) + xlab("Time [hpf]") + 
#   ylab("Expr. level [normalized]") +
#   scale_y_continuous(breaks=c(0,0.5,1)) +
#   scale_x_continuous(breaks=c(3,4,5,6)) +
#   theme(panel.grid.minor = element_blank(), legend.position = "bottom") +
#   scale_color_dMod(name="Step at time point")

data_switchplot2 <- subset(data_switchplot, genename%in%c("slc35a5", "sox3", "sox19b", "nanog"))
P2 <- ggplot(data_switchplot2, aes(x=time, y=value)) + #geom_smooth(color="grey", linetype="dashed", fill=NA) +
  geom_segment(data = data_switchplot2, aes(x=2.2, xend=time+0.25, y=ystep1, yend=ystep1, color=time==timepoint, linetype=time==timepoint)) +
  geom_segment(data = data_switchplot2, aes(x=time+0.25, xend=time+0.25, y=ystep1, yend=ystep2, color=time==timepoint, linetype=time==timepoint)) +
  geom_segment(data = data_switchplot2, aes(x=time+0.25, xend=6.5, y=ystep2, yend=ystep2, color=time==timepoint, linetype=time==timepoint)) +
  geom_point() + facet_wrap(~genename, scales="fixed", nrow=1) + theme_dMod(base_size = 10) + xlab("Time [hpf]") + 
  ylab("Expr. level [normalized]") +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  scale_x_continuous(breaks=c(3,4,5,6)) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.key.width = unit(3, "line")) +
  scale_color_manual(name="One-step models", values= c("red", "grey"), breaks=c(TRUE, FALSE), labels=c("Best", "Candidate")) + 
  scale_linetype_manual(name="One-step models", values=c("solid", "longdash"), breaks=c(TRUE, FALSE), labels=c("Best", "Candidate"))

P3 <- ggplot(data_switchplot, aes(x=time, y=value)) + #geom_smooth(color="grey", linetype="dashed", fill=NA) +
  #geom_vline(data = data_switchplot, aes(xintercept=timepoint, color=switch)) +
  geom_segment(data = data_switchplot, aes(x=2.2, xend=6.5, y=ymean, yend=ymean), color="blue") +
  geom_point() + facet_wrap(~genename, scales="fixed") + theme_dMod(base_size = 10) + xlab("Time [hpf]") + 
  geom_segment(data = subset(data_switchplot, time==timepoint), aes(x=2.2, xend=time+0.25, y=ystep1, yend=ystep1), color="red") +
  geom_segment(data = subset(data_switchplot, time==timepoint), aes(x=time+0.25, xend=time+0.25, y=ystep1, yend=ystep2), color="red") +
  geom_segment(data = subset(data_switchplot, time==timepoint), aes(x=time+0.25, xend=6.5, y=ystep2, yend=ystep2), color="red") +
  ylab("Expr. level [normalized]") +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  scale_x_continuous(breaks=c(3,4,5,6)) +
  theme(panel.grid.minor = element_blank(), legend.position = "right") +
  scale_color_dMod(name="", breaks=c("down", "none", "up"), labels=c("Switch down", "No switch", "Switch up"))

pdf(file = paste0(path, "Switchplot.pdf"), width = 7, height = 7)
ggdraw() +
  #draw_plot(P1, x=0, y=.6, width=.5, height=.4) +
  draw_plot(P2, x=0, y=.6, width=1, height=.4) +
  draw_plot(P3, x=0, y=0, width=1, height=.6) +
  draw_plot_label(c("A", "B"), x=c(0,0), y=c(1,.6), size=12)
dev.off()


# generate and output SSGS plot
source("myplotSSGS.R")
pdf(file = paste0(path, "SSGS_colorPlot.pdf"), width = 7, height = 8)
SSGS<- myplotSSGS(resultCombined, c(2.5,3,3.5,4,4.5,5,5.5,6), 
                          myanalyzeConditions = analyzeConditions,
                          experimentStepDetection = experimentStepDetection) + 
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2)) +
  xlab("Time point of fold change") + ylab("Time point of switch")
  
print(SSGS)
dev.off()

# RNAsense:::outputGeneTables(resultCombined)





myplotSSGS <- function(myresultCombined = resultCombined, mytimes = times, myanalyzeConditions = analyzeConditions, experimentStepDetection = experimentStepDetection){
  # auxiliary function for application of fisher test
  getFT <- function(myresult=result, myswitch="up",switchtime=3, xaxis="2.5hpf", identifier="FCdown"){
    if(identifier == "FCdown"){
      a=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                     grepl(xaxis, FCdown)))[1]
      b=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                     grepl(xaxis, FCdown)))[1]
      c=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                     !grepl(xaxis, FCdown)))[1]
      d=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                     !grepl(xaxis, FCdown)))[1]
    } else {
      a=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                     grepl(xaxis, FCup)))[1]
      b=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                     grepl(xaxis, FCup)))[1]
      c=dim(subset(myresult, switch==myswitch & timepoint==switchtime &
                     !grepl(xaxis, FCup)))[1]
      d=dim(subset(myresult, !(switch==myswitch & timepoint==switchtime) &
                     !grepl(xaxis, FCup)))[1]
    }
    
    pvalue_suppress <- fisher.test(rbind(c(a,b), c(c,d)), alternative="less")$p.value
    pvalue_enhance <- fisher.test(rbind(c(a,b), c(c,d)), alternative="greater")$p.value
    if(pvalue_suppress < 1e-20) cluster <- "1E-20 Suppression" else
      if(pvalue_suppress < 1e-10) cluster <- "1E-10 Suppression" else
        if(pvalue_suppress < 1e-5) cluster <- "1E-05 Suppression" else
          if(pvalue_suppress < 1e-2) cluster <- "1E-02 Suppression" else
            if(pvalue_enhance < 1e-20) cluster <- "1E-20 Enhancement" else
              if(pvalue_enhance < 1e-10) cluster <- "1E-10 Enhancement" else
                if(pvalue_enhance < 1e-5) cluster <- "1E-05 Enhancement" else
                  if(pvalue_enhance < 1e-2) cluster <- "1E-02 Enhancement" else cluster <- "none"
    data.frame(pvalue_enhance=pvalue_enhance, pvalue_suppress=pvalue_suppress, cluster=cluster, myswitch=myswitch)
  }
  
  out <- do.call(rbind, lapply(mytimes, function(t){
    do.call(rbind, lapply(paste0(format(mytimes, nsmall = 1),"hpf"), function(x){
      do.call(rbind, lapply(c("FCdown", "FCup"), function(ident){
        do.call(rbind, lapply(c("up", "down"), function(myswitch){
          cbind(getFT(myresult=subset(myresultCombined, experiment==experimentStepDetection), myswitch=myswitch, switchtime = t, xaxis = x, identifier = ident),
                time=t, xaxis=x, identifier=ident, experiment=experimentStepDetection)
        }))
      }))
    }))
  }))
  out$cluster <- factor(out$cluster, levels = c("none", "1E-20 Enhancement","1E-10 Enhancement", "1E-05 Enhancement","1E-02 Enhancement",
                                                "1E-20 Suppression","1E-10 Suppression", "1E-05 Suppression", "1E-02 Suppression"))
  
  mylabeller <- c("FCdown" = paste0(myanalyzeConditions[1]," > ",myanalyzeConditions[2]),
                  "FCup" = paste0(myanalyzeConditions[2]," > ",myanalyzeConditions[1]),
                  "up" = "up",
                  "down" = "down")
  
  P <- ggplot(out, aes(x=xaxis, y=time, fill=cluster)) + geom_tile(color="black") +
    xlab(paste0("Stage-specific gene sets ",experimentStepDetection)) + ylab("Switch Time") +theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
    scale_x_discrete(breaks=paste0(format(seq(2.5,6, by=0.5),nsmall=1), "hpf"),labels=format(seq(2.5,6, by=0.5),nsmall=1)) +
    scale_y_reverse(breaks=seq(2.5,6, by=0.5)) +
    facet_grid(myswitch~identifier, labeller = as_labeller(mylabeller)) +
    scale_fill_manual(name="Fisher's exact test",
                      values=c("none" = "grey", "1E-20 Enhancement" = "darkred", "1E-10 Enhancement" = "red", "1E-05 Enhancement"= "orange",
                               "1E-02 Enhancement" = "yellow", "1E-20 Suppression"= "violet", "1E-10 Suppression" = "darkblue",
                               "1E-05 Suppression" = "blue", "1E-02 Suppression" = "lightblue" ),
                      breaks = c("none", "1E-20 Enhancement", "1E-10 Enhancement", "1E-05 Enhancement", "1E-02 Enhancement",
                                 "1E-20 Suppression", "1E-10 Suppression", "1E-05 Suppression", "1E-02 Suppression"),
                      labels = c("none", "p<1e-20 Enhancement", "p<1e-10 Enhancement", "p<1e-05 Enhancement", "p<1e-02 Enhancement",
                                 "p<1e-20 Suppression", "p<1e-10 Suppression", "p<1e-05 Suppression", "p<1e-02 Suppression"))
  return(P)
}
#' Students T-test and the analysis of variance (ANOVA) on Omics data sets
#'
#' MassOmics R-package also includes standard statistical analyses such as the Students T-test and ANOVA, coupled with false discovery rate (FDR), Bonferroni correction, and Tukey HSD test to adjust p-values when conducting multiple comparisons. To open the user-interface for statistical analysis, click <Statistics> and <Run ANVOA & T-test>. The implementation of Students T-test or ANOVA is determined automatically by number of input conditions. After selecting the type of statistical corrections, p value filtering, and log transformation, click <submit> and a pop-out window will appear to allow the user to classify samples into each group. The final statistical result is saved in the chosen working directory.
#' 
#' @return None
#' @examples
#' Omics_htest()
#'
#' @export
Omics_htest <- function(){ htest<-function (conditions=1:2, signif.level = 0.05, save = TRUE, log = TRUE, stat_type="False disover rate")
{
  
  report_folder <- tk_choose.dir(caption = "Select the Working directory")
  setwd(report_folder)
  
  data <- tk_choose.files(caption = "Select metabolite or pathways profile in .csv")
  final.df <- read.csv(data)
  
  final.stat <- final.df
  abc <- c(LETTERS)
  
  selected.samples <- numeric()
  for (q in 1:length(conditions)) {
    cat(paste("Select samples from condition:", q, "\n"))
    samples <- select.list(names(final.stat), multiple = TRUE,
                           title = paste("samples in condition", q))
    columns <- names(final.stat)
    equal <- 0
    equal_list <- 0
    for (t in 1:length(samples)) {
      equal <- which(columns == samples[t])
      if (length(equal) != 0) {
        equal_list <- rbind(equal_list, equal)
      }
    }
    equal_list <- (equal_list[-1])
    names(final.stat)[equal_list] <- abc[q]
    if (length(selected.samples) == 0) {
      selected.samples <- equal_list
    }
    else {
      selected.samples <- c(selected.samples, equal_list)
    }
  }
  facNames <- names(final.stat)[selected.samples]
  replicates <- factor(facNames)
  if (length(conditions) > 2) {
    cat(paste("ANOVA & TukeyHSD is in progress..", "\n"))
    yMat2 <- final.stat[selected.samples]
    
    
    if(log==TRUE) {
      yMat2[is.na(yMat2)]<-0
      yMat2<-log(yMat2)
      yMat2[yMat2=="-Inf"]<-0
    }
    
    aov.list <- apply(yMat2, 1, function(x) aov(as.numeric(x) ~ replicates))
    aov.summary <-  lapply(aov.list, anova)
    
    
    ############ p value correction #####################
    
    if(stat_type=="Non"){
      pvalues <- data.frame(aov.summary)
      t <- 5
      for (i in 1:nrow(yMat2)) {
        yMat2$pvalues[i] <- pvalues[1, t]
        t <- t + 5
      }
      pvalues_list<-as.numeric(yMat2$pvalues)
      ANOVA.df<-cbind(final.df, ANOVA=pvalues_list)
      ANOVA.df<-subset(ANOVA.df, ANOVA.df$ANOVA < signif.level)
      ANOVA.df <- ANOVA.df[order(ANOVA.df$ANOVA, decreasing = F),]
      
      if (log==TRUE) { sheet <- "ANOVA (log).csv"} else {sheet <- "ANOVA.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv(  ANOVA.df , file =  FileName , row.names = FALSE)
        cat(paste("ANOVA",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
    }
    
    if(stat_type=="TukeyHSD"){
      
      TukeyHSDlist <- lapply(aov.list,TukeyHSD)
      pvalues <- data.frame(aov.summary)
      
      TukeyHSD.df <-NULL
      remove(TukeyHSD.df )
      
      t <- 5
      for (i in 1:nrow(yMat2)) {
        yMat2$pvalues[i] <- pvalues[1, t]
        t <- t + 5
        if (!exists("TukeyHSD.df")){
          TukeyHSD.df <- data.frame(TukeyHSDlist[[i]]$replicates[,4])
          rownames(TukeyHSD.df)<-rownames(TukeyHSDlist[[i]]$replicates)
          colnames(TukeyHSD.df)<-i
        } else {
          TukeyHSD_L <- data.frame(TukeyHSDlist[[i]]$replicates[,4])
          rownames(TukeyHSD_L)<-rownames(TukeyHSDlist[[i]]$replicates)
          TukeyHSD.df <- merge(TukeyHSD.df, TukeyHSD_L,by="row.names", all=TRUE)
          rownames(TukeyHSD.df) <- TukeyHSD.df$Row.names
          TukeyHSD.df$Row.names <- NULL
          colnames(TukeyHSD.df)[i]<-i
        }
      }
      TukeyHSD.df1<-t(TukeyHSD.df)
      rownames(TukeyHSD.df1)<-NULL
      TukeyHSD.df1<-data.frame(TukeyHSD.df1)
      TukeyHSD.true<-TukeyHSD.df1<signif.level
      TukeyHSD.Sig<-data.frame(apply(TukeyHSD.true, 1, function(x) table(x)["TRUE"]))
      TukeyHSD.Sig[is.na(TukeyHSD.Sig)]<-0
      names(TukeyHSD.Sig)<-"#.TukeyHSD.sig"
      TukeyHSD.data<-cbind.data.frame(TukeyHSD.df1, TukeyHSD.Sig, yMat2[ncol(yMat2)])
      rownames(TukeyHSD.data)<- rownames(yMat2)
      
      finalResults.df<-merge(final.df,TukeyHSD.data, by =0)
      finalResults.df[1]<-NULL
      
      finalResults.df <- subset( finalResults.df, as.numeric(finalResults.df$pvalues) <
                                   signif.level)
      finalResults.df <-  finalResults.df[order( finalResults.df$"#.TukeyHSD.sig", decreasing = T),]
      
      if (log==TRUE) { sheet <- "ANOVA & TukeyHSD Result (log).csv"} else {sheet <- "ANOVA & TukeyHSD Result.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv(finalResults.df, file =  FileName , row.names = FALSE)
        cat(paste("ANOVA & TukeyHSD complete!",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
    }
    
    if(stat_type=="False discovery rate"){
      
      require(qvalue)
      
      pvalues <- data.frame(aov.summary)
      t <- 5
      for (i in 1:nrow(yMat2)) {
        yMat2$pvalues[i] <- pvalues[1, t]
        t <- t + 5
      }
      pvalues_list<-as.numeric(yMat2$pvalues)
      FDR <- qvalue(pvalues_list)$qvalue ###### mistake need to fix......
      FDR.df<-cbind(final.df, ANOVA=pvalues_list, False_discovery_rate = FDR)
      FDR.df<-subset(FDR.df, FDR.df$False_discovery_rate < signif.level)
      FDR.df<-subset(FDR.df, FDR.df$ANOVA < signif.level)
      FDR.df <- FDR.df[order(FDR.df$False_discovery_rate, decreasing = F),]
      
      if (log==TRUE) { sheet <- "ANOVA & False discovery rate (log).csv"} else {sheet <- "ANOVA & False discovery rate.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv( FDR.df, file =  FileName , row.names = FALSE)
        cat(paste("ANOVA & False discovery rate",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
      
    }
    
    if(stat_type=="Bonferroni correction"){
      pvalues <- data.frame(aov.summary)
      t <- 5
      for (i in 1:nrow(yMat2)) {
        yMat2$pvalues[i] <- pvalues[1, t]
        t <- t + 5
      }
      pvalues_list <-as.numeric(yMat2$pvalues)
      number <- length(pvalues_list)
      Bonferroni_correction <-signif.level/number
      BC.df <-cbind(final.df, ANOVA_Bonferroni=pvalues_list)
      BC.df <-subset(BC.df, BC.df$ANOVA_Bonferroni < Bonferroni_correction)
      BC.df <- BC.df[order(BC.df$ANOVA_Bonferroni, decreasing = F),]
      
      if (log==TRUE) { sheet <- "ANOVA & Bonferroni correction (log).csv"} else {sheet <- "ANOVA & Bonferroni correction.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv( BC.df , file =  FileName , row.names = FALSE)
        cat(paste("ANOVA & Bonferroni correction",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
      
    }
    
  }
  else {
    cat(paste("t-test is in progress..", "\n"))
    yMat2 <- final.stat[selected.samples]
    if(log==TRUE) {
      yMat2[is.na(yMat2)]<-0
      yMat2<-log(yMat2)
      yMat2[yMat2=="-Inf"]<-0
    }
    
    list.p <- apply(yMat2, 1, function(x) t.test(as.numeric(x) ~
                                                   replicates))
    list.p <- unlist(list.p)
    t <- 3
    for (i in 1:nrow(yMat2)) {
      yMat2$pvalues[i] <- list.p[t]
      t <- t + 11
    }
    
    if(stat_type=="Non"){
      finalResult.df <- data.frame(final.df, Ttest <-as.numeric(yMat2$pvalues))
      finalResult.df <- subset(finalResult.df, finalResult.df$Ttest < signif.level)
      finalResult.df <-  finalResult.df[order(finalResult.df$Ttest, decreasing = F),]
      finalResult.df[1] <- NULL
      
      if (log==TRUE) {sheet <- "TTest result (log).csv"} else {sheet <- "TTest result.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv(finalResult.df, file =  FileName, row.names = FALSE)
        cat(paste("t-test complete! The file ", sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
    }
    
    if(stat_type=="False discovery rate"){
      require(qvalue)
      pvalues_list<-as.numeric(yMat2$pvalues)
      FDR <- qvalue(pvalues_list)$qvalue
      FDR.df<-cbind(final.df, TTest=pvalues_list, False_discovery_rate = FDR)
      FDR.df<-subset(FDR.df, FDR.df$False_discovery_rate < signif.level)
      FDR.df<-subset(FDR.df, FDR.df$TTest < signif.level)
      FDR.df <- FDR.df[order(FDR.df$False_discovery_rate, decreasing = F),]
      
      if (log==TRUE) { sheet <- "TTest & False discovery rate (log).csv"} else {sheet <- "TTest & False discovery rate.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv( FDR.df, file =  FileName , row.names = FALSE)
        cat(paste("TTest & False discovery rate",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
    }
    
    if(stat_type=="Bonferroni correction"){
      
      pvalues_list <-as.numeric(yMat2$pvalues)
      number <- length(pvalues_list)
      Bonferroni_correction <-signif.level/number
      BC.df <-cbind(final.df, ANOVA_Bonferroni=pvalues_list)
      BC.df <-subset(BC.df, BC.df$ANOVA_Bonferroni < Bonferroni_correction)
      BC.df <- BC.df[order(BC.df$ANOVA_Bonferroni, decreasing = F),]
      
      if (log==TRUE) { sheet <- "TTtest & Bonferroni correction (log).csv"} else {sheet <- "TTest & Bonferroni correction.csv"}
      if (save == TRUE) {
        store <- paste(report_folder, "\\", sheet, sep = "")
        FileName <- tclvalue(tkgetSaveFile(initialfile= store))
        write.csv( BC.df , file =  FileName , row.names = FALSE)
        cat(paste("TTest & Bonferroni correction",sheet ,"was saved in the folder",
                  report_folder, "\n"))
      }
    }
    
    if(stat_type=="TukeyHSD") {stop("TukeyHSD only works on 3 or more conditions")}
    
  }
}
##### user interphase

require(tcltk)
tclRequire("BWidget")
ConditionN <- tclVar("")
P_value <- tclVar("0.05")
gg <- tktoplevel()
tkwm.title(gg,"ANOVA or T-TEST")
c.entry <- tkentry(background="white",gg, textvariable=ConditionN)
p.entry <- tkentry(background="white",gg, textvariable=P_value)
reset <- function() {
  tclvalue(ConditionN)<-""
  tclvalue(P_value)<-""
}


##### Choice stat correction #####

stat <- c("Non","False discovery rate","Bonferroni correction","TukeyHSD")
comboBox <- tkwidget(gg,"ComboBox",editable=FALSE,values=stat,textvariable=tclVar("False discovery rate"),width=18)

reset.but <- tkbutton(gg, text="Reset", command=reset)

cb <- tkcheckbutton(gg)
cbValue <- tclVar("1")
tkconfigure(cb,variable=cbValue)

submit <- function() {
  c <- as.numeric(tclvalue(ConditionN))
  conditions<-1:c
  signif.level<-as.numeric(tclvalue(P_value))
  cbVal <- as.character(tclvalue(cbValue))
  if (cbVal=="1")
    log<-TRUE
  if (cbVal=="0")
    log<-FALSE
  stat_choose <-stat[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
  htest(conditions=conditions,signif.level=signif.level, log = log, stat_type = stat_choose)
}


submit.but <- tkbutton(gg, text="submit", command=submit)
quit.but <- tkbutton(gg, text = "Close Session",
                     command = function() {
                       tkdestroy(gg)
                     }
)

tkgrid(tklabel(gg,text="Run T-TEST (conditions=2)  |  Run ANOVA (conditions>2)"),columnspan=3, pady = 10)
tkgrid(tklabel(gg,text="Number of conditions"), c.entry, pady= 10, padx= 10)
tkgrid(tklabel(gg,text="Statistical correction"),comboBox, pady= 10, padx= 10)
tkgrid(tklabel(gg,text="P-value"), p.entry, tklabel(gg,text="Log transform"),cb, pady= 10, padx= 10)
tkgrid(submit.but, reset.but, quit.but, pady= 10, padx= 10)
}
#' Partial least squares-discriminant analysis (PLS-DA) analysis on Omics data sets
#'
#' Partial least squares-discriminant analysis (PLS-DA) is used to compress the dimension of multivariate data, and thereby display outliers and relationships between samples or classification of conditions. Under <Statistics> click <Run PLSDA> and a user interface window will appear to adjust a number of conditions, sphere radius, text size, log transfor-mation, and model validations. To plot 3 dimensional PLSDA, click <submit> and a pop-out window will allow the user to classify sam-ples into each group. All parameters can be changed on the user inter-face of PLSDA and replotted instantly via the <Re-Plot> function.
#'
#' @return None
#' @examples
#' Omics_PCA()
#'
#' @export
Omics_PCA <- function(){
  
  PCA<-function(conditions=1:3, ball.size=0.05, text.size = 1, log=TRUE){
    
    
    
    require(rgl)
    
    report_folder <<- tk_choose.dir(caption = "Select the working directory")
    setwd(report_folder)
    original.df<-  tk_choose.files(caption = "Select metabolite or pathways profile in .csv")
    final.df <- read.csv(original.df)
    final.df[is.na(final.df)]<-0
    final.stat <- final.df
    abc <- 1:100
    
    options("menu.graphics"=TRUE)
    selected.samples <- numeric()
    for (q in 1:length(conditions)) {
      cat(paste("Select samples from condition:", q, "\n"))
      samples <- select.list(names(final.stat), multiple = TRUE,
                             title = paste("Samples in condition", q))
      columns <- names(final.stat)
      equal <- 0
      equal_list <- 0
      for (t in 1:length(samples)) {
        equal <- which(columns == samples[t])
        if (length(equal) != 0) {
          equal_list <- rbind(equal_list, equal)
        }
      }
      equal_list <- (equal_list[-1])
      names(final.stat)[equal_list] <- abc[q]
      if (length(selected.samples) == 0) {
        selected.samples <- equal_list
      }
      else {
        selected.samples <- c(selected.samples, equal_list)
      }
    }
    group <- names(final.stat)[selected.samples]
    
    data.df<- final.stat[,selected.samples]
    names(data.df)<-names(final.df[,selected.samples])
    
    group<-as.numeric(group)
    replicates <- group[!duplicated(group)]
    colur.panel<-rainbow(length((replicates)))
    colour<-NULL
    for(i in 1:length(group)){
      colour<- c(colour,colur.panel[group[i]])
    }
    colour_replot<<-colour
    cat(paste("PCA analysis is in progress..", "\n"))
    
    
    data.replot<<-data.df
    
    if(log==TRUE){
      
      data.df[is.na(data.df)]<-0
      data.df<-log(data.df)
      data.df[data.df=="-Inf"]<-0
    }
    
    pc1<-prcomp(data.df)
    PC<-pc1$r
    #biplot(pc1)
    
    ### 3D PCA plot ##
    
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d(PC[,1],PC[,2],PC[,3], radius=ball.size, color=colour, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    mtext3d(paste("PC1 ","(",round(summary(pc1)$importance[2,"PC1"]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(summary(pc1)$importance[2,"PC2"]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(summary(pc1)$importance[2,"PC3"]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    
    ifelse(text.size==0, print("No sample name"),text3d(PC[,1],PC[,2],PC[,3], text=rownames(PC), col=colour,adj=c(0.5,2),cex=text.size))
    
    title3d("", "", "", "", "", col='black')
    
    # setup coordinate
    #snap <- par3d( c("userMatrix") )
    #par3d("windowRect")  # To see window position
    #dput(snap)
    
    
    
    ##### loading plots
    
    
    
    
    
    pc <- princomp(data.df[1:nrow(data.df),], cor=TRUE, scores=TRUE)
    open3d()
    
    par3d("windowRect"= c(161, 134, 944, 759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    loading.name<<-final.df$Name
    text3d(pc$score[,1:3],texts=loading.name,cex=text.size*0.7)
    text3d(pc$loadings[,1:3], texts=rownames(pc$loadings), col=colour,cex=text.size*0.7)
    spheres3d(pc$loadings[,1:3], radius=ball.size*3, color=colour, alpha=1, shininess=20)
    coords <- NULL
    colour2<- NULL
    for (i in 1:nrow(pc$loadings)) {
      coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
    }
    lines3d(coords, col=c("red"), lwd=1)
    aspect3d(1, 1, 1)
    
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    mtext3d(paste("PC1 ","(",round(summary(pc1)$importance[2,"PC1"]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(summary(pc1)$importance[2,"PC2"]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(summary(pc1)$importance[2,"PC3"]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    
    
    
    
  }
  
  
  replot<-function(data.replot=data.replot,ball.size=0.01, text.size = 1,colour1=colour_replot, log=TRUE){
    
    if(log==TRUE){
      data.replot<-log(data.replot)
      data.replot[data.replot=="-Inf"]<-0
    }
    
    
    pca<-prcomp(data.replot)
    pca1<-pca$r
    
    require(rgl)
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d(pca1[,1],pca1[,2],pca1[,3], radius=ball.size, color=colour1, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    mtext3d(paste("PC1 ","(",round(summary(pca)$importance[2,"PC1"]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(summary(pca)$importance[2,"PC2"]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(summary(pca)$importance[2,"PC3"]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    ifelse(text.size==0, print("No sample name"),text3d(pca1[,1],pca1[,2],pca1[,3], text=rownames(pca1), col=colour1,adj=c(0.5,2),cex=text.size))
    title3d("", "", "", "", "", col='black')
    
    
    #loading ploting
    pc <- princomp(data.replot[1:nrow(data.replot),], cor=TRUE, scores=TRUE)
    
    open3d()
    par3d("windowRect"= c(161, 134, 944, 759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    ifelse(text.size==0, print("No sample name"), text3d(pc$score[,1:3],texts= loading.name,cex=text.size*0.7))
    ifelse(text.size==0, print("No sample name"),text3d(pc$loadings[,1:3], texts=rownames(pc$loadings), col=colour1,cex=text.size*0.7))
    spheres3d(pc$loadings[,1:3], radius=ball.size*3, color=colour1, alpha=1, shininess=20)
    coords <- NULL
    for (i in 1:nrow(pc$loadings)) {
      coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
    }
    lines3d(coords, col=c("red"), lwd=1)
    
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    mtext3d(paste("PC1 ","(",round(summary(pca)$importance[2,"PC1"]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(summary(pca)$importance[2,"PC2"]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(summary(pca)$importance[2,"PC3"]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
  }
  
  
  
  require(tcltk)
  ConditionN <- tclVar("")
  ball.size <- tclVar("0.05")
  text.size <- tclVar("1")
  PL <- tktoplevel()
  tkwm.title(PL,"PCA analysis")
  c.entry <- tkentry(background="white",PL, textvariable=ConditionN)
  p.entry <- tkentry(background="white",PL, textvariable=ball.size)
  t.entry <- tkentry(background="white",PL, textvariable=text.size)
  
  reset <- function() {
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    replot(ball.size=ball.size,text.size=text.size,data.replot=data.replot,colour1=colour_replot,log=log)
  }
  
  reset.but <- tkbutton(PL, text="Re-Plot", command=reset)
  
  cb <- tkcheckbutton(PL)
  cbValue <- tclVar("1")
  tkconfigure(cb,variable=cbValue)
  
  submit <- function() {
    c <- as.numeric(tclvalue(ConditionN))
    conditions<-1:c
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    PCA(conditions=conditions,ball.size=ball.size,text.size=text.size,log=log)
  }
  submit.but <- tkbutton(PL, text="submit", command=submit)
  
  
  SaveImage<-function(){
    
    time<-format(Sys.time(), "%H:%M:%S")
    time2<-gsub(":",".",time)
    pca.name<-paste("PCA(",time2,").pdf",sep="")
    rgl.postscript(pca.name,"pdf")
    print(paste(pca.name," is save under ", report_folder,sep=""))
  }
  saveImage.but <-tkbutton(PL, text="Save PDF", command=SaveImage)
  
  quit.but <- tkbutton(PL, text = "Close Session",
                       command = function() {
                         tkdestroy(PL)
                         rgl.close()
                       }
  )
  tkgrid(tklabel(PL,text="Enter the Parameters"),columnspan=4, pady = 10)
  tkgrid(tklabel(PL,text="Number of conditions"), c.entry, pady= 10, padx= 10)
  tkgrid(tklabel(PL,text="Sphere Radius"), p.entry, pady= 10, padx= 10)
  tkgrid(tklabel(PL,text="Text size"), t.entry, tklabel(PL,text="Log transform"),cb, pady= 10, padx= 10)
  tkgrid(submit.but, reset.but, saveImage.but, quit.but, pady= 10, padx= 10)
}

#' Principal component analysis (PCA) analysis on Omics data sets
#'
#' Principal component analysis (PCA) is used to compress the dimension of multivariate data, and thereby display outliers and relationships between samples or classification of conditions. Under <Statistics> click <Run PCA> and a user interface window will appear to adjust a number of conditions, sphere radius, text size, log transfor-mation, and model validations. To plot 3 dimensional PCA, click <submit> and a pop-out window will allow the user to classify sam-ples into each group. All parameters can be changed on the user inter-face of PCA and replotted instantly via the <Re-Plot> function.
#' 
#' @return None
#' @examples
#' Omics_PLSDA()
#'
#' @export
#' 
Omics_PLSDA <- function(){
  
  PLSDA<-function(conditions=1:2, ball.size=0.7, text.size = 1, log=TRUE, valid=TRUE){
    
    
    require(mixOmics)
    require(rgl)
    require(tcltk)
    
    report_folder <<- tk_choose.dir(caption = "Select the working directory")
    setwd(report_folder)
    original.df<-  tk_choose.files(caption = "Select metabolite or pathways profile in .csv")
    final.df <- read.csv(original.df)
    final.df[is.na(final.df)]<-0
    final.stat <- final.df
    abc <- 1:100
    
    options("menu.graphics"=TRUE)
    selected.samples <- numeric()
    for (q in 1:length(conditions)) {
      cat(paste("Select samples from condition:", q, "\n"))
      samples <- select.list(names(final.stat), multiple = TRUE,
                             title = paste("Samples in condition", q))
      columns <- names(final.stat)
      equal <- 0
      equal_list <- 0
      for (t in 1:length(samples)) {
        equal <- which(columns == samples[t])
        if (length(equal) != 0) {
          equal_list <- rbind(equal_list, equal)
        }
      }
      equal_list <- (equal_list[-1])
      names(final.stat)[equal_list] <- abc[q]
      if (length(selected.samples) == 0) {
        selected.samples <- equal_list
      }
      else {
        selected.samples <- c(selected.samples, equal_list)
      }
    }
    group <- names(final.stat)[selected.samples]
    
    data.df <- final.stat[,selected.samples]
    names(data.df)<-names(final.df[,selected.samples])
    
    group<-as.numeric(group)
    group1<<-group
    replicates <- group[!duplicated(group)]
    colur.panel<-rainbow(length((replicates)))
    colour<-NULL
    for(i in 1:length(group)){
      colour<- c(colour,colur.panel[group[i]])
    }
    colour_replot<<-colour
    
    
    
    plsda.replot<<-data.df
    
    if(log==TRUE){
      data.df[is.na(data.df)]<-0
      data.df<-log(data.df)
      data.df[data.df=="-Inf"]<-0
    }
    
    plsr<-pls(t(data.df),mode="regression", group, ncomp=4)
    
    
    #plsr<-plsda(t(data.df), group, ncomp=3)
    
    
    #barplot(pls.val$Q2.total,main= "Leave-one-out cross validation (Q2)")
    
    
    
    
    if(valid==TRUE){
      nSample <- length(data.df[1,])
      pls.val<-perf(plsr, validation = "Mfold", folds = nSample/2)
      win.graph(width = 4.5, height = 4 ,pointsize = 10)
      counts<-rbind(pls.val$R2[1:4],pls.val$Q2.total[1:4])
      rownames(counts)<-c("R2","Q2")
      colnames(counts)<-c("PC1","PC2","PC3","PC4")
      val.plot<-barplot(counts, main="Leave-one-out cross validation", ylim=c(0,1.1), beside=T, legend=rownames(counts),
                        args.legend = list(x = "topright", bty = "n",inset=c(-0.08,-0.05)))
      text(val.plot, counts, labels = round(counts,2), pos = 3,cex=0.9)
    }
    
    
    
    ### 3D PCA plot ##
    
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d(plsr$variates$X[,1],plsr$variates$X[,2],plsr$variates$X[,3], radius=ball.size, color=colour, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    ifelse(text.size == 0, print("No sample name"), text3d(plsr$variates$X[,1],plsr$variates$X[,2],plsr$variates$X[,3], text=rownames(plsr$variates$X), col=colour,adj=c(0.5,2),cex=text.size)) # Add names #
    
    mtext3d(paste("PC1 ","(",round(plsr$explained_variance$X[1]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(plsr$explained_variance$X[2]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(plsr$explained_variance$X[3]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    title3d("", "", "", "", "", col='black')
    cat(paste("PLSDA analysis is Completed", "\n"))
    # setup coordinate
    #snap <- par3d( c("userMatrix") )
    #par3d("windowRect") see position
    #dput(snap)
    return(plsr)
  }
  
  
  replot<-function(plsda.replot1=plsda.replot,ball.size=0.7, text.size = 1,color1=colour_replot, group=group1, log = TRUE,valid=TRUE){
    
    
    
    if(log==TRUE){
      plsda.replot1[is.na(plsda.replot1)]<-0
      plsda.replot1<-log(plsda.replot1)
      plsda.replot1[plsda.replot1=="-Inf"]<-0
    }
    
    
    plsr1<-pls(t(plsda.replot1), mode="regression",group, ncomp=4)
    
    
    
    
    if(valid==TRUE){
      nSample <- length(plsda.replot1[1,])
      pls.val<-perf(plsr, validation = "Mfold", folds = nSample/2)
      win.graph(width = 4.5, height = 4 ,pointsize = 10)
      counts<-rbind(pls.val$R2[1:4],pls.val$Q2.total[1:4])
      rownames(counts)<-c("R2","Q2")
      colnames(counts)<-c("PC1","PC2","PC3","PC4")
      val.plot<-barplot(counts, main="Leave-one-out cross validation", ylim=c(0,1.1), beside=T, legend=rownames(counts),
                        args.legend = list(x = "topright", bty = "n",inset=c(-0.08,-0.05)))
      text(val.plot, counts, labels = round(counts,2), pos = 3,cex=0.9)
    }
    
    
    
    
    require(rgl)
    rgl.close()
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d(plsr1$variates$X[,1],plsr1$variates$X[,2],plsr1$variates$X[,3], radius=ball.size, color=color1, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=F,col="black")
    axis3d("X+",nticks = 5)
    axis3d("Y",nticks = 5)
    axis3d("Z++",nticks = 5)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    ifelse(text.size == 0, print("No sample name"), text3d(plsr1$variates$X[,1],plsr1$variates$X[,2],plsr1$variates$X[,3], text=rownames(plsr1$variates$X), col=color1,adj=c(0.5,2),cex=text.size)) # Add names #
    mtext3d(paste("PC1 ","(",round(plsr$explained_variance$X[1]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(plsr$explained_variance$X[2]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(plsr$explained_variance$X[3]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    
    title3d("", "", "", "", "", col='black')
    cat(paste("PLSDA analysis is Completed", "\n"))
    return(plsr1)
  }
  
  require(tcltk)
  ConditionN <- tclVar("")
  ball.size <- tclVar("0.5")
  text.size <- tclVar("1")
  PL <- tktoplevel()
  tkwm.title(PL,"PLSDA analysis")
  c.entry <- tkentry(background="white",PL, textvariable=ConditionN)
  p.entry <- tkentry(background="white",PL, textvariable=ball.size)
  t.entry <- tkentry(background="white",PL, textvariable=text.size)
  
  reset <- function() {
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    
    ValVal <- as.character(tclvalue(ValValue))
    if (ValVal=="1")
      valid<-TRUE
    if (ValVal=="0")
      valid<-FALSE
    
    replot(ball.size=ball.size,text.size=text.size,plsda.replot=plsda.replot,color1=colour_replot,log=log,valid=valid)
  }
  reset.but <- tkbutton(PL, text="Re-Plot", command=reset)
  
  
  cb <- tkcheckbutton(PL)
  cbValue <- tclVar("1")
  tkconfigure(cb,variable=cbValue)
  
  
  Val <- tkcheckbutton(PL)
  ValValue <- tclVar("1")
  tkconfigure(Val,variable=ValValue)
  
  
  
  submit <- function() {
    c <- as.numeric(tclvalue(ConditionN))
    conditions<-1:c
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    
    
    ValVal <- as.character(tclvalue(ValValue))
    if (ValVal=="1")
      valid<-TRUE
    if (ValVal=="0")
      valid<-FALSE
    
    
    PLSDA(conditions=conditions,ball.size=ball.size,text.size=text.size,log=log,valid=valid)
  }
  submit.but <- tkbutton(PL, text="submit", command=submit)
  
  
  SaveImage<-function(){
    
    time<-format(Sys.time(), "%H:%M:%S")
    time2<-gsub(":",".",time)
    pca.name<-paste("PLSDA(",time2,").pdf",sep="")
    rgl.postscript(pca.name,"pdf")
    print(paste(pca.name," is save under ", report_folder,sep=""))
  }
  saveImage.but <-tkbutton(PL, text="Save PDF", command=SaveImage)
  
  quit.but <- tkbutton(PL, text = "Close Session",
                       command = function() {
                         tkdestroy(PL)
                         rgl.close()
                       }
  )
  tkgrid(tklabel(PL,text="Enter the Parameters"),columnspan=3, pady = 10)
  tkgrid(tklabel(PL,text="Number of conditions"), c.entry, pady= 10, padx= 10)
  tkgrid(tklabel(PL,text="Sphere Radius"),p.entry, tklabel(PL,text="Validation"),Val , pady= 10, padx= 10)
  tkgrid(tklabel(PL,text="Text size"), t.entry, tklabel(PL,text="Log transform"),cb, pady= 10, padx= 10)
  tkgrid(submit.but, reset.but, saveImage.but, quit.but, pady= 10, padx= 10)
}

######################################### Omics_sPCA #########################################
Omics_sPCA <- function(){
  require(rgl)
  require (tcltk )
  require("mixOmics")
  
  
  PCA<-function(conditions=1:3, ball.size=0.01, text.size = 0.5, log=TRUE, ncomp1=10, ncomp2=10, ncomp3=10,  cen=TRUE, rescale=TRUE){
    
    
    report_folder <<- tk_choose.dir(caption = "Select the working directory")
    setwd(report_folder)
    original.df<-  tk_choose.files(caption = "Select metabolite or pathways profile in .csv")
    final.df <- read.csv(original.df)
    final.df[is.na(final.df)]<-0
    final.stat <- final.df
    abc <- 1:100
    
    
    selected.samples <- numeric()
    for (q in 1:length(conditions)) {
      cat(paste("Select samples from condition:", q, "\n"))
      samples <- select.list(names(final.stat), multiple = TRUE,
                             title = paste("Samples in condition", q))
      columns <- names(final.stat)
      equal <- 0
      equal_list <- 0
      for (t in 1:length(samples)) {
        equal <- which(columns == samples[t])
        if (length(equal) != 0) {
          equal_list <- rbind(equal_list, equal)
        }
      }
      equal_list <- (equal_list[-1])
      names(final.stat)[equal_list] <- abc[q]
      if (length(selected.samples) == 0) {
        selected.samples <- equal_list
      }
      else {
        selected.samples <- c(selected.samples, equal_list)
      }
    }
    group <- names(final.stat)[selected.samples]
    
    data.df<- final.stat[,selected.samples]
    names(data.df)<-names(final.df[,selected.samples])
    
    group<-as.numeric(group)
    replicates <- group[!duplicated(group)]
    colur.panel<-rainbow(length((replicates)))
    colour<-NULL
    for(i in 1:length(group)){
      colour<- c(colour,colur.panel[group[i]])
    }
    colour_replot<<-colour
    cat(paste("sPCA analysis is in progress..", "\n"))
    
    
    data.replot<<-data.df
    
    if(log==TRUE){
      data.df[is.na(data.df)]<-0
      data.df<-log(data.df)
      data.df[data.df=="-Inf"]<-0
    }
    
    #  spca.results <- spca(t(data.df), ncomp=3, center = TRUE, scale =TRUE, keepX =c(ncomp1,ncomp2,ncomp3))
    spca.results <- spca(t(data.df), ncomp=3, center = cen, scale =rescale, keepX =c(ncomp1,ncomp2,ncomp3))
    
    
    # pc1<-prcomp(data.df)
    # PC<-pc1$r
    #biplot(pc1)
    
    ### 3D PCA plot ##
    
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d( spca.results$x[,1], spca.results$x[,2], spca.results$x[,3], radius=ball.size, color=colour, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=FALSE,col="black")
    axis3d("X+",nticks = 6)
    axis3d("Y",nticks = 6)
    axis3d("Z++",nticks = 6)
    mtext3d(paste("PC1 ","(",round(spca.results$varX[1]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(spca.results$varX[2]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(spca.results$varX[3]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    
    ifelse(text.size==0, print("No sample name"),text3d(spca.results$x[,1],spca.results$x[,2],spca.results$x[,3], text= spca.results$names$indiv, col= colour,adj=c(0.5,2),cex=text.size))
    
    title3d("", "", "", "", "", col='black')
    
    cat(paste("sPCA analysis is Complected", "\n"))
    
  }
  
  
  replot<-function(data.replot=data.replot,ball.size=0.01, text.size = 1,colour1=colour_replot, log=TRUE, ncomp1=10, ncomp2=10, ncomp3=10,cen=TRUE, rescale=TRUE){
    
    if(log==TRUE){
      data.replot[is.na(data.replot)]<-0
      data.replot<-log(data.replot)
      data.replot[data.replot=="-Inf"]<-0
    }
    
    spca.results_R <- spca(t(data.replot), ncomp=3, center = cen, scale =rescale, keepX =c(ncomp1,ncomp2,ncomp3))
    
    
    require(rgl)
    open3d()
    par3d("windowRect"= c(1018,134,1801,759))
    myframe <-structure(list(userMatrix = structure(c(0.686889827251434, -0.209069654345512,
                                                      0.696040630340576, 0, 0.0137451644986868, -0.953820705413818,
                                                      -0.300063371658325, 0, 0.726631700992584, 0.21567764878273, -0.652295589447021,
                                                      0, 0, 0, 0, 1), .Dim = c(4L, 4L)), FOV = 30, zoom = 1.20000004768372), .Names = c("userMatrix",
                                                                                                                                        "FOV", "zoom"))
    par3d(myframe)
    bg3d("white")
    spheres3d( spca.results_R$x[,1], spca.results_R$x[,2], spca.results_R$x[,3], radius=ball.size, color=colour1, alpha=1, shininess=20)
    aspect3d(1, 1, 1)
    axes3d(edges=c('x++', 'y--', 'z-+'),labels=FALSE,col="black")
    axis3d("X+",nticks = 6)
    axis3d("Y",nticks = 6)
    axis3d("Z++",nticks = 6)
    mtext3d(paste("PC1 ","(",round(spca.results_R$varX[1]*100,1),"%",")",sep=""),edge="x+",line=3,las=2)
    mtext3d(paste("PC2 ","(",round(spca.results_R$varX[2]*100,1),"%",")",sep=""),edge="Y",line=3,las=2)
    mtext3d(paste("PC3 ","(",round(spca.results_R$varX[3]*100,1),"%",")",sep=""),edge="z++",line=3,las=2)
    
    grid3d(c("x", "y+", "z+"),lty="dash",col="grey")
    ifelse(text.size==0, print("No sample name"),text3d(spca.results_R$x[,1],spca.results_R$x[,2],spca.results_R$x[,3], text= spca.results_R$names$indiv, col=colour1,adj=c(0.5,2),cex=text.size))
    title3d("", "", "", "", "", col='black')
  }
  
  require(tcltk)
  ConditionN <- tclVar("3")
  ball.size <- tclVar("0.01")
  text.size <- tclVar("0.5")
  ncomp1.var <- tclVar("10")
  ncomp2.var <- tclVar("10")
  ncomp3.var <- tclVar("10")
  
  PL <- tktoplevel()
  tkwm.title(PL,"Sparse PCA analysis")
  c.entry <- tkentry(background="white",PL, textvariable=ConditionN, width=17)
  p.entry <- tkentry(background="white",PL, textvariable=ball.size, width=17)
  t.entry <- tkentry(background="white",PL, textvariable=text.size, width=17)
  ncomp1.entry <- tkentry(background="white",PL, textvariable= ncomp1.var, width=17)
  ncomp2.entry <- tkentry(background="white",PL, textvariable= ncomp2.var, width=17)
  ncomp3.entry <- tkentry(background="white",PL, textvariable= ncomp3.var, width=17)
  
  reset <- function() {
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    ncomp1.var<-as.numeric(tclvalue(ncomp1.var))
    ncomp2.var<-as.numeric(tclvalue(ncomp2.var))
    ncomp3.var<-as.numeric(tclvalue(ncomp3.var))
    
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="1")
      rescaleb<-TRUE
    if (rbVal=="0")
      rescaleb<-FALSE
    
    cenbVal <- as.character(tclvalue(cenbValue))
    if (cenbVal=="1")
      cenb<-TRUE
    if (cenbVal=="0")
      cenb<-FALSE
    
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    replot(ball.size=ball.size,text.size=text.size,data.replot=data.replot,colour1=colour_replot,log=log, ncomp1=ncomp1.var, ncomp2=ncomp2.var, ncomp3=ncomp3.var, rescale=rescaleb,  cen= cenb)
  }
  
  reset.but <- tkbutton(PL, text="Re-Plot", command=reset)
  
  #log
  cb <- tkcheckbutton(PL)
  cbValue <- tclVar("1")
  tkconfigure(cb,variable=cbValue)
  #rescale
  rb <- tkcheckbutton(PL)
  rbValue <- tclVar("1")
  tkconfigure(rb,variable=rbValue)
  #centre
  cenb <- tkcheckbutton(PL)
  cenbValue <- tclVar("1")
  tkconfigure(cenb,variable=cenbValue)
  
  
  submit <- function() {
    c <- as.numeric(tclvalue(ConditionN))
    conditions<-1:c
    ball.size<-as.numeric(tclvalue(ball.size))
    text.size<-as.numeric(tclvalue(text.size))
    ncomp1.var<-as.numeric(tclvalue(ncomp1.var))
    ncomp2.var<-as.numeric(tclvalue(ncomp2.var))
    ncomp3.var<-as.numeric(tclvalue(ncomp3.var))
    
    
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="1")
      rescaleb<-TRUE
    if (rbVal=="0")
      rescaleb<-FALSE
    
    cenbVal <- as.character(tclvalue(cenbValue))
    if (cenbVal=="1")
      cenb<-TRUE
    if (cenbVal=="0")
      cenb<-FALSE
    
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="1")
      log<-TRUE
    if (cbVal=="0")
      log<-FALSE
    PCA(conditions=conditions,ball.size=ball.size,text.size=text.size,log=log, ncomp1=ncomp1.var, ncomp2=ncomp2.var, ncomp3=ncomp3.var, rescale=rescaleb,  cen= cenb)
  }
  submit.but <- tkbutton(PL, text="submit", command=submit)
  
  
  SaveImage<-function(){
    
    time<-format(Sys.time(), "%H:%M:%S")
    time2<-gsub(":",".",time)
    pca.name<-paste("PCA(",time2,").pdf",sep="")
    rgl.postscript(pca.name,"pdf")
    print(paste(pca.name," is save under ", report_folder,sep=""))
  }
  saveImage.but <-tkbutton(PL, text="Save PDF", command=SaveImage)
  
  quit.but <- tkbutton(PL, text = "Close Session",
                       command = function() {
                         tkdestroy(PL)
                         rgl.close()
                       }
  )
  tkgrid(tklabel(PL,text="Enter the Parameters"),columnspan=4, pady = 10)
  tkgrid(tklabel(PL,text="Number of conditions"), c.entry, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(PL,text="Number of variables for PC1"),  ncomp1.entry , pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(PL,text="Number of variables for PC2"),  ncomp2.entry , pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(PL,text="Number of variables for PC3"),  ncomp3.entry , tklabel(PL,text="Centralise data"),cenb, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(PL,text="Sphere Radius"), p.entry, tklabel(PL,text="Scale data"),rb, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(PL,text="Text size"), t.entry, tklabel(PL,text="Log transform"),cb, pady= 10, padx= 10, sticky="w")
  tkgrid(submit.but, reset.but, saveImage.but, quit.but, pady= 10, padx= 10)
}
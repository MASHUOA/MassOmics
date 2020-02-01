
######################################### DataCorrection #########################################
DataCorrection <- function(){
  ## normalize data by mutiple internal stadnards
  
  NormMultIS <-function(cutoff= 0.3, IS.index = c(1,2,3,4)){
    
    library(tcltk)
    print(cutoff)
    print(IS.index)
    
    folder <- tk_choose.dir(caption = "Select working directory")
    setwd(folder)
    dir()
    dataR <- tk_choose.files(caption = "Select metabolite profile in .csv format (e.g. GC-MS Result.csv)")
    info <- tk_choose.files(caption = "Select subject informatiom in .csv format (e.g. info.csv)")
    
    data.df <- read.csv(dataR)
    data.df1 <- t(data.df[,-c(1:which(names(data.df)=="Name"))])
    
    data.df1[is.na(data.df1)] <- 0
    info.df <- read.csv(info)
    info.df$Name <-make.names(info.df$Name)
    
    
    merge.df <- merge(info.df, data.df1, by.x = "Name", by.y = "row.names"  , all= TRUE)
    nTrue <- length(merge.df[is.na(merge.df[,ncol(merge.df)])==TRUE])
    
    if(nTrue> 0) { tkmessageBox(message = paste("ERROR:",nTrue, "Sample names differ between spreadsheets"), icon = "warning", type = "ok"); stop("ERROR: Sample names differ")}
    
    #Batch <- merge.df$Batch
    Type <- merge.df$Type
    Data <-merge.df[,-c(1:length(names(info.df)))]
    IS.names<-as.character(data.df$Name[IS.index])
    AA<-gsub(",", "",toString(paste(IS.names,'\n',sep="")))
    
    res <- tkmessageBox(title = "Your select internal standards",
                        message = AA, icon = "question", type = "yesno", default = "yes")
    
    if(as.character(res)=="no"){stop("select internal standards again")}
    
    
    HairMetab <-Data
    IS<-Data[,IS.index]
    nMetab<-dim(HairMetab)[2]
    cormat<-matrix(NA,nMetab,length(IS.index))
    for(i in 1:nMetab){
      for(j in 1:length(IS.index)){
        cormat[i,j]<-cor(HairMetab[Type=='QC',i],IS[Type=='QC',j], use='pairwise.complete.obs')
      }
    }
    maxcors<-apply(cormat,1,max)
    
    #abline(v=0.7, col=2)
    norm.index<-apply(cormat,1,which.max)
    
    # show which internal standard are used
    
    under<-maxcors< cutoff
    colnames(cormat) <- IS.names
    MetabolteName<-as.character(data.df$Name)
    ChosenIS <- norm.index
    ChosenIS[under]<-NA
    
    
    hist(maxcors, main ="Correlation")
    
    corMatrix <- cbind(MetabolteName, cormat, ChosenIS)
    write.csv(corMatrix, file= "Chosen IS for each metabolites.csv", row.names = FALSE)
    dev.new.OS()
    boxplot(data.df$Ret.Time~norm.index, main='Chosen Internal Standards', ylab='Retention time')
    
    
    HairMetabReplace<- HairMetab
    for(i in 1:nMetab){
      if(under[i]==TRUE) next()
      HairMetabReplace[,i]<-HairMetab[,i]/IS[,norm.index[i]]*median(IS[,norm.index[i]])
    }
    
    final.df <- cbind(data.df[,c(1:which(names(data.df)=="Name"))],t(HairMetabReplace))
    names(final.df) <-c( names(data.df[,c(1:which(names(data.df)=="Name"))]), merge.df$Name)
    
    
    
    
    do.call(data.frame,lapply(final.df, function(x) replace(x, is.infinite(x),NA)))
    final.df[is.na(final.df)]<-0
    
    store<- "2_GC-MS Result_NormIS.csv"
    FileName <- tclvalue(tkgetSaveFile(initialfile=store))
    write.csv(final.df, file=FileName , row.names = FALSE)
    
  }
  
  ## interphase
  
  library(tcltk2)
  tclRequire("BWidget")
  # parameters
  
  
  
  #function
  submitA <- function() {
    a <- as.numeric(strsplit(tclvalue(avar), ",")[[1]])
    b <- as.numeric(tclvalue(bvar))
    NormMultIS(IS.index=a, cutoff=b)
  }
  
  #Button
  
  
  
  win1 <- tktoplevel()
  tkwm.title(win1,"Data correction")
  win1$env$nb <- tk2notebook(win1, tabs = c("   Intro   ", "  Scaling  ","   Contaminants  ", "  Multiple IS   ", "  Batch effect  "))
  tkpack(win1$env$nb, fill = "both", expand = TRUE)
  win1$env$tb1 <- tk2notetab(win1$env$nb, "  Multiple IS   ")
  
  submitA.but <- tkbutton(win1, text="Run Analysis", command=submitA, width= 10)
  
  quit.but <- tkbutton(win1 , text = "Close Session",
                       command = function() {
                         tkdestroy(kk)
                       }
  )
  
  open.but <- tkbutton(win1 ,text="Open Folder",
                       width= 10,
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )
  
  
  avar <- tclVar(c("31, 50, 98, 159"))
  bvar <- tclVar("0")
  
  a.entry <- tkentry(background="white",win1, textvariable=avar, width=20)
  b.entry <- tkentry(background="white",win1, textvariable=bvar, width=20)
  
  
  
  frameNM<- tkframe(win1$env$tb1,borderwidth=2, relief="groove",padx=5,pady=5)
  tkgrid(tklabel(frameNM,text="Positions for IS"), a.entry, pady= 8, padx= 10, sticky="w")
  tkgrid(tklabel(frameNM,text="Correlation cut-off"), b.entry, pady= 8, padx= 10, sticky="w")
  
  frameButton<- tkframe(win1$env$tb1)
  tkgrid(tklabel(frameButton,text=""), submitA.but, open.but , pady= 10, padx= 10, sticky="w")
  
  
  #tkgrid(tklabel(kk,image=logo_1))
  tkgrid(tklabel(win1$env$tb1,text=""), pady= 5, padx= 20, sticky="w")
  tkgrid(tklabel(win1$env$tb1,text="Settings for multiple internal standard (IS) normalization"), pady = 10, padx= 20, sticky="w")
  tkgrid(frameNM, pady= 0, padx= 30)
  tkgrid(tklabel(win1$env$tb1,text=""), pady= 5, padx= 20, sticky="w")
  tkgrid(frameButton, pady= 10, padx= 30)
  
  ######################################### Remove contamint #########################
  
  RemoveContaminant <- function(A=0.4,B=2){
    #Setup workding directory
    folder <- tk_choose.dir(caption = "Select working directory")
    setwd(folder)
    dir()
    dataR <- tk_choose.files(caption = "Select metabolite profile in .csv format (e.g. GC-MS Result.csv)")
    info <- tk_choose.files(caption = "Select subject informatiom in .csv format (e.g. info.csv)")
    
    data.df <<- read.csv(dataR)
    data.df1 <- t(data.df[,-c(1:which(names(data.df)=="Name"))])
    data.df1[is.na(data.df1)] <- 0
    info.df <- read.csv(info)
    info.df$Name <-make.names(info.df$Name)
    
    merge.df <- merge(info.df, data.df1, by.x = "Name", by.y = "row.names"  , all= TRUE)
    
    
    
    nTrue <- length(merge.df[is.na(merge.df[,ncol(merge.df)])==TRUE])
    merge.df[,1:5]
    if(nTrue> 0) { tkmessageBox(message = paste("ERROR:",nTrue, "Sample names differ between spreadsheets"), icon = "warning", type = "ok"); stop("ERROR: Sample names differ")}
    
    Batch <<- merge.df$Batch
    Type <<- merge.df$Type
    Data <<-merge.df[,-c(1:length(names(info.df)))]
    
    # Beatrix codes
    IdentifyContaminants<-function(Data, Batch, Type, A,  B){
      nMetab<-dim(Data)[[2]]
      Batch<-factor(Batch)
      nBatch<-length(unique(Batch))
      outcome<-data.frame(Crit1=rep(NA,nMetab),Crit2=rep(NA, nMetab), PropBelow=rep(NA, nMetab), RatioAbove=rep(NA, nMetab))
      # print(c(A, B))
      for(i in 1:nMetab){
        blankmeans<-sapply(split(Data[Type=='Blank',i],Batch[Type=='Blank']),mean)
        samples<-split(Data[Type=='Sample',i], Batch[Type=='Sample'])
        BatchLength<-sapply(samples, length)
        TotBelow<-rep(NA, nBatch)
        RatioAbove<-rep(NA, nBatch)
        for(j in 1:nBatch){
          TotBelow[j]<-sum(samples[[j]]< blankmeans[[j]])
          RatioAbove[j]<-mean((samples[[j]][samples[[j]]>blankmeans[[j]]])/blankmeans[[j]])
        }
        outcome$PropBelow[i]<-sum(TotBelow)/sum(Type=='Sample')
        outcome$RatioAbove[i]<-mean(RatioAbove[blankmeans>0], na.rm=T)
        if(outcome$PropBelow[i]<A & mean(samples$`1`)!=0){outcome$Crit1[i]<-'PASS' } else{outcome$Crit1[i]<-'FAIL'}
        if(outcome$PropBelow[i]>0 |mean(samples$`1`)==0){
          if(!is.na(outcome$RatioAbove[i]) & outcome$RatioAbove[i]> B & (1-outcome$PropBelow[i])> (1-A)/2 & mean(samples$`1`)!=0){outcome$Crit2[i]<-'PASS'}
          else{outcome$Crit2[i]<-'FAIL'}
        }}
      outcome[,2][is.na(outcome[,2])]<-'PASS'
      return(outcome)
    }
    Finaloutcome<-IdentifyContaminants(Data, Batch, Type, A,  B)
    
    x<-(1-A)*100
    names(Finaloutcome) <- c(paste( x ,"PercenAboveBlank"), paste( B, "RatioAboveBlank"),"PercenAboveBlank","RatioAboveBlank" )
    
    
    
    subtract<-sapply(Finaloutcome[3],function(x) (1-x)*100)
    
    Finaloutcome[,3] <-c(round(subtract,1))
    Finaloutcome[,4] <-c(round(Finaloutcome[4],1))
    findResult.df <- cbind(Finaloutcome, data.df)
    
    store<-"1_GC-MS Result_RemoveContaminants.csv"
    FileName <- tclvalue(tkgetSaveFile(initialfile=store))
    write.csv(findResult.df, file=FileName )
    
  }
  
  plotfigure <-function(metaNo=30){
    Batch1 <- Batch[Type=="Blank"|Type=="Sample"]
    Type1 <- Type[Type=="Blank"|Type=="Sample"]
    MetaName<-data.df$Name[metaNo]
    Data1 <- Data[Type=="Blank"|Type=="Sample",]
    plot(Batch1,Data1[,metaNo], col=factor(Type1), ylab="level", main=   MetaName ,xlab ="Batch")
    points(Batch1[Type=='Blank'], Data1[,metaNo][Type1=='Blank'],pch=19)
    legend('topright',pch=c(19,1), col=c(1:2),legend=c('Blank','Sample'),cex=0.8)
  }
  
  
  
  
  
  ## Interphase for Remove contaminat
  
  
  
  xvar <- tclVar("60")
  vvar <- tclVar("2")
  yvar <- tclVar("1")
  
  
  
  # Run function
  submit <- function() {
    x <- (100-as.numeric(tclvalue(xvar)))/100
    v <- as.numeric(tclvalue(vvar))
    RemoveContaminant(A=x,B=v)
  }
  
  # Plot function
  
  plotBatchfigure <- function() {
    yfigure <- as.numeric(tclvalue(yvar))
    plotfigure(metaNo= yfigure)
  }
  
  
  #Tip
  Tips<-function () {
    require(tcltk)
    kk <- tktoplevel()
    tktitle(kk) <- "Excel format for the subject information"
    image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/SubjectInfo.gif", sep=""))
    FrontImage<-tcl("image",  "create", "photo", "MaxC2.image", file=image.path)
    tkgrid(tklabel(kk,image=FrontImage))
  }
  
  
  # Reset function
  reset <- function() {
    tclvalue(xvar)<-"60"
    tclvalue(vvar)<-"2"
    tclvalue(yvar)<-"1"
  }
  
  reset.but <- tkbutton(win1, text="Reset", command=reset, width=9)
  
  submit.but <- tkbutton(win1, text="Run Analysis", command=submit, width= 10)
  
  
  quit.but <- tkbutton(win1, text = "Close Session",
                       command = function() {
                         tkdestroy(kk)
                       }
  )
  
  open.but <- tkbutton(win1,text="Open Folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )
  
  plot.but <- tkbutton(win1, text = "Plot Figure",
                       command = plotBatchfigure, width=9)
  
  Tips.but<-tkbutton(win1, text="Example", command=Tips,width=9)
  
  ## Userinterphase structure
  
  win1$env$tb2 <- tk2notetab(win1$env$nb, "   Contaminants  ")
  
  frameMiddle<- tkframe(win1$env$tb2 ,borderwidth=2, relief="groove",padx=5,pady=5)
  
  x.entry <- tkentry(background="white",frameMiddle, textvariable=xvar, width=15)
  v.entry <- tkentry(background="white",frameMiddle, textvariable=vvar, width=15)
  y.entry <- tkentry(background="white",frameMiddle, textvariable=yvar, width=15)
  
  
  
  tkgrid(tklabel(frameMiddle,text="% samples above blanks"), x.entry,pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameMiddle,text="Fold changes"),v.entry, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameMiddle,text="Plot metabolite using ID"), y.entry, pady= 10, padx= 10, sticky="w")
  
  
  #tkgrid(tklabel(frameMiddle,text="Retention Time (min,max)"), y.entry, z.entry, tklabel(frameMiddle,text="Sec"), pady= 8, padx= 10, sticky="w")
  
  # Button
  frameButton<- tkframe(win1$env$tb2)
  tkgrid(tklabel(frameButton,text=""), open.but, submit.but, plot.but, pady= 10, padx= 10, sticky="w")
  
  
  #reset.but, open.but, quit.but,
  tkgrid(tklabel(win1$env$tb2,text=""), pady= 5, padx= 20, sticky="w")
  # Interphase structure
  #tkgrid(tklabel(kk,image=logo_1))
  tkgrid(tklabel(win1$env$tb2 ,text="Fitler Criteria to remove contaminants"), pady = 10, padx= 20, sticky="w")
  #tkgrid(frameAll)
  tkgrid(frameMiddle, pady= 0, padx= 30)
  tkgrid(frameButton, pady= 10, padx= 30)
  
  
  ################################################Batch correction######################
  BatchCorrection <- function(funMeanMedian= "median", SampleType = "QC", datatype="original", Log_trans=F, vis_norm_result=F){
    
    ## Batch normalization ##
    
    library(tcltk)
    library("rfviz")
    message(paste("Applying",paste(SampleType,collapse = "+"),"'s",funMeanMedian,"on samples"))
    
    folder <- tk_choose.dir(caption = "Select working directory")
    setwd(folder)
    #dir()
    if (datatype=="original"){
      dataR <- tk_choose.files(caption = "Select metabolite profile in .csv format (e.g. GC-MS Result.csv)")
      info <- tk_choose.files(caption = "Select subject informatiom in .csv format (e.g. info.csv)")
      
      data.df <- read.csv(dataR)
      data.df1 <- t(data.df[,-c(1:which(names(data.df)=="Name"))])
      data.df1[is.na(data.df1)] <- 0
      data.df1[(data.df1==-Inf)] <- 0
      if (Log_trans) {
        data.df1<-log(data.df1)
        data.df1[(data.df1==-Inf)] <- 0
      }
      
      
      info.df <- read.csv(info)
      info.df <- info.df[,c("Name", "Type", "Batch")]
      info.df$Name <-make.names(info.df$Name)
      
      merge.df <- merge(info.df[,c("Name", "Type", "Batch")], data.df1, by.x = "Name", by.y = "row.names"  , all.y= TRUE)
      merge.df <- merge.df[!is.na(merge.df$Batch),]
      
      nTrue <- nrow(merge.df[is.na(merge.df[,ncol(merge.df)]),])
      
      if(nTrue> 0) { tkmessageBox(message = paste("ERROR:",nTrue, "Sample names differ between spreadsheets"), icon = "warning", type = "ok"); stop("ERROR: Sample names differ")}
      
      Type <- merge.df$Type
      Data <- merge.df[,-c(1:length(names(info.df)))]
      Batch <-merge.df$Batch
      SampleNorm <- Data
      
      nMetab<-dim(SampleNorm)[[2]]
      nBatch <- length(table(Batch))
    }
    
    if (datatype=="Customized"){
      dataR <- tk_choose.files(caption = "Select metabolite profile in .csv format (e.g. GC-MS Result.csv)")
      
    }
    
    if(funMeanMedian=="median"){funM <-median} 
    
    if(funMeanMedian=="mean"){funM <-mean} 
    
    if(funMeanMedian=="SERRF"){funM <-SERRF}
    
    if((funMeanMedian=="median")|(funMeanMedian=="mean")){
      
      for( i in 1:nMetab){
        M<-median(Data[Type==SampleType,i], na.rm=T)
        for(j in 1:nBatch){
          m<-funM(Data[Type==SampleType & Batch==j,i], na.rm=T) # change here for mean or median
          if(m==0) m<-mean(SampleNorm[Type==SampleType & Batch==j,i], na.rm=T)
          if(m>0){SampleNorm[Batch==j,i] <- Data[Batch==j,i]*M/m } else {print(c(round(i,digits = 0),round(j,digits = 0),m, M/m))}
        }
      }
    }else if(funMeanMedian=="SERRF"){
      
      SERRF(input=input,Predict_level=SampleT)
      
      
      
    }
    
    
    
    # Write results
    final.df <- cbind(data.df[,c(1:which(names(data.df)=="Name"))],t(SampleNorm))
    names(final.df) <-c( names(data.df[,c(1:which(names(data.df)=="Name"))]), merge.df$Name)
    
    
    
    final.df[is.na(final.df)]<-0
    
    
    store<- "3_GC-MS Result_NormBatch.csv"
    FileName <- tclvalue(tkgetSaveFile(initialfile=store))
    write.csv(final.df, file=FileName , row.names = FALSE)
    
    
    if(vis_norm_result){
      
    if(SampleType=="QC"){
      ## Plot QC
      means.before<-apply(Data[Type=="QC",],2,mean, na.rm=T)
      sds.before<-apply(Data[Type=="QC",],2,sd, na.rm=T)
      cv.before<-sds.before/means.before
      
      means.post<-apply(SampleNorm[Type=="QC",],2,mean, na.rm=T)
      sds.post<-apply(SampleNorm[Type=="QC",],2,sd, na.rm=T)
      cv.post<-sds.post/means.post
      dev.new.OS()
      reg1 <- lm(cv.before~cv.post)
      plot(cv.post, cv.before, ylab='Coefficent of variation (CV) before batch correction',
           xlab='Coefficent of variation (CV) after batch correction' ,main='QC',col="red")
      abline(reg1)
      abline(a=0,b=1, col = "lightgray", lty = 3)
      mtext(paste("Number of imrpoved metabolites:", sum(cv.s.before > cv.s.post,na.rm = T))) 
      
    }
    ## Plot sample
    means.s.before<-apply(Data[Type=="Sample",],2,mean, na.rm=T)
    sds.s.before<-apply(Data[Type=="Sample",],2,sd, na.rm=T)
    cv.s.before<-sds.s.before/means.s.before
    
    means.s.post<-apply(SampleNorm[Type=="Sample",],2,mean, na.rm=T)
    sds.s.post<-apply(SampleNorm[Type=="Sample",],2,sd, na.rm=T)
    cv.s.post<-sds.s.post/means.s.post
    dev.new.OS()
    reg1 <- lm(cv.s.before~cv.s.post)
    plot(cv.s.post, cv.s.before, ylab='Coefficent of variation (CV) before batch correction',
         xlab='Coefficent of variation (CV) after batch correction' ,main='Samples',col="red")
    abline(reg1)
    abline(a=0,b=1, col = "lightgray", lty = 3)
    mtext(paste("Number of imrpoved metabolites:", sum(cv.s.before > cv.s.post,na.rm = T))) 
    }

  }
  
  ## interphase for BatchCorrection
  
  #library(tcltk2)
  #tclRequire("BWidget")
  
  #win1 <- tktoplevel()
  #tkwm.title(win1,"Data correction")
  #win1$env$nb <- tk2notebook(win1, tabs = c("   Contaminants  ", "  Multiple IS   ", "  Batch effect  "))
  #tkpack(win1$env$nb, fill = "both", expand = TRUE)
  
  
  win1$env$tb3 <- tk2notetab(win1$env$nb, "  Batch effect  ")
  
  ## Setup
  SampleType <- c("QC","Sample","QC+Sample")
  #StatType <- c("median","mean","polynomial")
  StatType <- c("SERRF","median","mean")
  comboBox1 <- tkwidget(background="white",win1 ,"ComboBox",editable=FALSE,values=SampleType ,textvariable=tclVar("QC"),width=9)
  comboBox2 <- tkwidget(background="white",win1 ,"ComboBox",editable=FALSE,values=StatType ,textvariable=tclVar("SERRF"),width=9)
  
  zero_imputaion_check <- tkcheckbutton(win1)
  zero_imputaion<- tclVar("1")
  tkconfigure(zero_imputaion_check,variable=zero_imputaion)
  
  vis_norm_result_check <- tkcheckbutton(win1)
  vis_norm_result<- tclVar("1")
  tkconfigure(vis_norm_result_check,variable=vis_norm_result)
  ## Function
  submitC <- function() {
    zero_imputaion=tclvalue(zero_imputaion)
    if (zero_imputaion=="0"){zero_imputaion<-FALSE}
    if (zero_imputaion=="1"){zero_imputaion<-TRUE}
    vis_norm_result=tclvalue(vis_norm_result)
    if (vis_norm_result=="0"){vis_norm_result<-FALSE}
    if (vis_norm_result=="1"){vis_norm_result<-TRUE}
    SampleT <- strsplit(SampleType[[as.numeric(tclvalue(tcl(comboBox1,"getvalue")))+1]],"\\+")[[1]]
    StatT <- StatType[[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1]]
    if (StatT %in% c("median","mean")){
      BatchCorrection(funMeanMedian= StatT , SampleType =  SampleT,vis_norm_result=vis_norm_result)
    }else if (StatT %in% c("SERRF")){
      SERRF(Predict_level = SampleT,datatype="MASSOMICS",zero_imputaion=zero_imputaion,vis_norm_result=vis_norm_result)
    }
    
  }
  
  submitd <- function() {
    zero_imputaion=tclvalue(zero_imputaion)
    if (zero_imputaion=="0"){zero_imputaion<-FALSE}
    if (zero_imputaion=="1"){zero_imputaion<-TRUE}
    SampleT <- base::strsplit(SampleType[[as.numeric(tclvalue(tcl(comboBox1,"getvalue")))+1]],"\\+")[[1]]
    StatT <- StatType[[as.numeric(tclvalue(tcl(comboBox2,"getvalue")))+1]]
    input=Data_prep_area_SERRF()
    SERRF(input=input,Predict_level=SampleT,zero_imputaion=zero_imputaion)
  }
  
  submitC.but <- tkbutton(win1, text="MassOmics", command=submitC , width= 12)
  submitd.but <- tkbutton(win1, text="MSDIAL", command=submitd , width= 12)
  open.but <- tkbutton(win1,text="Open Folder", width= 10,
                       command= function(){dir = getwd()
                       shell.exec(dir)}
  )
  
  frameCorrection<- tkframe(win1$env$tb3  ,borderwidth=2, relief="groove",padx=5,pady=5)
  tkgrid(tklabel(frameCorrection,text="Centralised batches based on"), comboBox1,pady= 10, padx= 10, sticky="W")
  tkgrid(tklabel(frameCorrection,text="Normalization method"),comboBox2, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameCorrection,text="Perform Zero value imputation (1/2*min)"), zero_imputaion_check, sticky="w")
  tkgrid(tklabel(frameCorrection,text="Check normalization result"), vis_norm_result_check, sticky="w")
  frameButtonC<- tkframe(win1$env$tb3)
  tkgrid(tklabel(frameButtonC,text=""),  submitC.but,submitd.but, open.but , pady= 10, padx= 10, sticky="w")
  
  tkgrid(tklabel(win1$env$tb3,text=""), pady= 5, padx= 20, sticky="w")
  tkgrid(tklabel(win1$env$tb3 ,text="Methods for minimise the batch effect"), pady = 10, padx= 20, sticky="w")
  tkgrid(frameCorrection, pady= 0, padx= 30)
  tkgrid(tklabel(win1$env$tb3,text="Select the source of file to be normalized"), pady= 5, padx= 20, sticky="w")
  tkgrid(frameButtonC, pady= 10, padx= 40)
  
  #################### Intro ###########################
  win1$env$tb0 <- tk2notetab(win1$env$nb, "   Intro   ")
  
  frameButtonE<- tkframe(win1$env$tb0)
  tkgrid(tklabel(frameButtonE,text=""),  Tips.but , pady= 5, padx= 5, sticky="w")
  
  tkgrid(tklabel(win1$env$tb0 ,text="Suggest GC-MS analysis setup (each batch)"), pady= 5, padx= 20, sticky="w")
  image.path1c<-c(paste(file.path(path.package(package="MassOmics")),"/R/Intro.gif", sep=""))
  logo_3<-tcl("image",  "create", "photo", "MaxC34.image", file=image.path1c)
  tkgrid(tklabel(win1$env$tb0,image=logo_3))
  tkgrid(tklabel(win1$env$tb0 ,text="Requirements:"), pady= 5, padx= 5, sticky="w")
  tkgrid(tklabel(win1$env$tb0,text="QCs, negative controls, and samples are organised in the\n subject information spreadsheet as shown by the example"), pady= 5, padx= 5, sticky="w")
  tkgrid( frameButtonE, pady= 5, padx= 5)
  
  #################### Scaling ###########################
  win1$env$tb5 <- tk2notetab(win1$env$nb, "  Scaling  ")
  
  Scaling<-function(){
    ScalingT <- ScaleType[[as.numeric(tclvalue(tcl(comboBox_scaling,"getvalue")))+1]]
    ImputeT <- ImputeType[[as.numeric(tclvalue(tcl(comboBox_Impute,"getvalue")))+1]]
    Data_scaling_impute(ScalingT=ScalingT,Impute=ImputeT)
  }
    
  ImputeType <- c("NULL","KNN","Minimum", "Zero")
  
  ScaleType <- c("NULL","Log","Mean","Pareto")
  comboBox_Impute <- tkwidget(background="white",win1 ,"ComboBox",editable=FALSE,values=ImputeType ,textvariable=tclVar("KNN"),width=9)
  comboBox_scaling <- tkwidget(background="white",win1 ,"ComboBox",editable=FALSE,values=ScaleType ,textvariable=tclVar("Log"),width=9)
  
  frameButtonF<- tkframe(win1$env$tb5)
  
  Scaling.but<-tkbutton(win1, text="Choose file", command=Scaling,width=9)
  tkgrid(tklabel(frameButtonF,text="Imputation method"),comboBox_Impute, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameButtonF,text="scaling method"),comboBox_scaling, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameButtonF,text="Run"),  Scaling.but , pady= 10, padx= 10, sticky="w")
  
  tkgrid(tklabel(win1$env$tb5 ,text="Select Massomics integration result to have a data scaling"), pady= 5, padx= 20, sticky="w")
  
  tkgrid(frameButtonF, pady= 5, padx= 5)
  
}

Data_scaling_impute<-function(path=NULL,ScalingT=NULL,Impute=NULL){
  library(tcltk)
  library(tcltk2)
  library(data.table)
  library(impute)
  if (is.null(path)){path=tk_choose.files(filter = matrix(c( "CSV", ".csv","excel file", ".xls", "excel file", ".xlsx"),3, 2, byrow = TRUE),
                                          caption  = "Choose integration result to have a data scaling")}
  
  data=fread(path)
  datacol=colnames(data)[!(colnames(data) %in% c("CAS","Ref.Ion" , "Total.ID(n=3)","Library.match","Ret.Time" ,"RT.shfit.Lower","Warnings" ,"Name"))]
  datafcol=colnames(data)[(colnames(data) %in% c("CAS","Ref.Ion" , "Total.ID(n=3)","Library.match","Ret.Time" ,"RT.shfit.Lower","Warnings" ,"Name"))]
  
  data=as.data.frame(data)
  for (col in datacol){
    data[,col]=as.numeric(as.character(data[,col]))  
  }
  
  # impute missing value.
 
  e=data[,datacol]
  e[e==0]<-NA
  e[e=="Inf"]<-NA
  e[e=="-Inf"]<-NA
  e[e==""]<-NA
  e[e=="NaN"]<-NA
  e=as.matrix(e)
  if((sum(is.na(e)))>0 && !is.null(Impute)){
    
    switch(Impute,
           Minimum= {e[is.na(e)] = min(e[!is.na(e)])},
           
           Zero= {e[is.na(e)] = 0},
           
           KNN= {e=try(impute.knn(e)$data,silent = T)}
           
           )
  }
  
  if(!is.null(ScalingT)){
    
    switch(ScalingT,
           Log= {
                e<-log(e)
                e[(e==-Inf)] <- 0
                },
           
           Mean= {e=e/mean(na.omit(e))},
           
           Pareto = {
                     for (row in nrow(e)){
                     e[row,]=e[row,]/sqrt(sd(e[row,]))
                     }
                    }
           
    )
  }
  
  data=cbind(data[,datafcol],e)
  
  write.csv(data,file = paste0(tools::file_path_sans_ext(path)," Impute_",ifelse(is.null(Impute),"None",Impute)," Scaling_",ifelse(is.null(ScalingT),"None",ScalingT),".csv"),row.names = F)
  
  message("Done")
  
  }
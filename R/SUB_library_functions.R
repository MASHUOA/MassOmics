######################################### Liz_analysis #########################################
Liz_analysis <- function(workdir=tk_choose.dir(caption = "Select working directory")){
  
  #Require library
  require(xcms)
  require(lattice)
  require(plyr)
  require(flux)
  library(tcltk)
  
  WorkingDir <-function(workdir = tk_choose.dir(caption = "Select working directory")){
    workdir <<- workdir
    setwd(workdir)
  }
  
  
  Id_report<-function(      workdir=workdir,
                            MS.L= tk_choose.files(caption="Select MS library (e.g. SVB_total or NIST) in .msl"),
                            amdis.report = read.csv(tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                                    multi = FALSE, caption = "Select the AMDIS report in .TXT"),
                                                    sep = "\t", stringsAsFactors=FALSE),
                            File.name = "Identification Report (Liz version).csv",
                            MsLibrary="InHouse",
                            Ret.Time.Filter=2){
    
    workdir <- workdir
    setwd(workdir)
    
    ## Extract Reference ion, CAS # from MS library
    lib.txt<-readLines(MS.L)
    lib.txt1<-lib.txt[-grep("^[ ]*$",lib.txt)]
    att.nms<-unique(sapply(strsplit(lib.txt1[grep(":",lib.txt1)],":"),function(x) x[1]))
    entry.nms<-sapply(strsplit(lib.txt1[grep("NAME:",lib.txt1)],"ME:"), function(x) x[2])
    rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))
    
    starts<-grep("NAME:", lib.txt)
    stops<-c(starts[1:(length(starts)-1)]+diff(starts)-1,length(lib.txt))
    for (i in 1:length(starts)){
      tmp<-lib.txt[starts[i]:stops[i]]
      sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
    }
    
    ## for Nist library
    
    if(MsLibrary=="NIST"){
      substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
      }
      reference_ion<-NULL
      
      for (i in 1:length(starts)){
        tmp<-lib.txt[starts[i]:stops[i]]
        reference_ion<-rbind(reference_ion, (as.numeric(substrRight(lapply(strsplit(grep("1000)",tmp,value=TRUE)," 1000)"), function(x) x[1]),3)))[1])
        
        sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
      }
      
      rez<-data.frame(rez,stringsAsFactors=FALSE)
      rez$NAME<-entry.nms
      
      
      libr<-cbind(rez,reference_ion)## returns list of peaks
      libr$NAME <- gsub("?", "", libr$NAME, fixed = TRUE)
      libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
      
      
      ############ Step 2. Generate the list of metabolites and clean up
      AmRep<-amdis.report
      # AmRep<-AmRep[!AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)+Ret.Time.Filter & AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)-Ret.Time.Filter,]## check this one ##
      AmRep$Name <- gsub("?", "", AmRep$Name, fixed = TRUE)
      AmRep$Name <- gsub("^ ", "", AmRep$Name, perl=T)
      
      ## RT statistics from the AReport
      RT.stats<-t(sapply(split(AmRep$RT, AmRep$Name),function(x) c(RT.median=round(median(x,na.rm=T),3),
                                                                   RT.mean=round(mean(x,na.rm=TRUE),3),
                                                                   RT.sd=sd(x,na.rm=T))))
      
      
      ID.stats<-t(sapply(split(AmRep$Weighted, AmRep$Name),function(x) c(AMDIS.percentage.match=round(median(x,na.rm=T)),
                                                                         Number.of.ID=length(x[!is.na(x)]))))
      
      
      #AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,ID.stats, Exp.RT=0,Diff.RT=0, ref_ion=0))##data frame for list of metabolites
      AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,ID.stats, ref_ion=0, CAS=0))##data frame for list of metabolites
      AmRep.RT.stats$ref_ion <-  libr$reference_ion[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
      #AmRep.RT.stats$Exp.RT <- as.numeric(t(sapply(split(AmRep$Expec..RT, AmRep$Name), median)))
      #AmRep.RT.stats$Diff.RT <- AmRep.RT.stats$Exp.RT- AmRep.RT.stats$RT.median
      AmRep.RT.stats$CAS <-  libr$CAS[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
      AmRep.RT.stats$Name <- rownames(RT.stats)
      AmRep.RT.stats <- AmRep.RT.stats[order(AmRep.RT.stats$RT.median, decreasing = F),]
      File.nameNIST <- "Identification Report_NIST (Liz version).csv"
      File.nameNIST <- tclvalue(tkgetSaveFile(initialfile=File.nameNIST))
      write.csv(AmRep.RT.stats, file = File.nameNIST,row.names = FALSE)
      print(paste("Identificaiton report was created and save in",workdir, sep=" "))
      
    } else {
      
      rez<-data.frame(rez,stringsAsFactors=FALSE)
      rez$NAME<-entry.nms
      libr<-rez## returns list of peaks
      
      
      ############ Step 2. Generate the list of metabolites and clean up
      AmRep<-amdis.report
      AmRep<-AmRep[!AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)+Ret.Time.Filter & AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)-Ret.Time.Filter,]## check this one ##
      AmRep$Name <- gsub("?", "", AmRep$Name, fixed = TRUE)
      AmRep$Name <- gsub("^ ", "", AmRep$Name, perl=T)
      
      ## RT statistics from the AReport
      RT.stats<-t(sapply(split(AmRep$RT, AmRep$Name),function(x) c(RT.median=round(median(x,na.rm=T),3),
                                                                   RT.mean=round(mean(x,na.rm=TRUE),3),
                                                                   RT.sd=sd(x,na.rm=T))))
      
      ID.stats<-t(sapply(split(AmRep$Weighted, AmRep$Name),function(x) c(AMDIS.percentage.match=round(median(x,na.rm=T)),
                                                                         Number.of.ID=length(x[!is.na(x)]))))
      
      AmRep.RT.stats<-NULL
      AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,ID.stats, Exp.RT=0,Diff.RT=0,ref_ion=0,CAS=0))##data frame for list of metabolites
      AmRep.RT.stats$ref_ion <- libr$RSN[match(rownames(RT.stats),libr$NAME)]
      AmRep.RT.stats$Exp.RT <- libr$RT[match(rownames(RT.stats),libr$NAME)]
      AmRep.RT.stats$CAS <- libr$CASNO[match(rownames(RT.stats),libr$NAME)]
      AmRep.RT.stats$Diff.RT <- as.numeric(AmRep.RT.stats$Exp.RT)- AmRep.RT.stats$RT.median
      AmRep.RT.stats$Name <- rownames(RT.stats)
      AmRep.RT.stats <- AmRep.RT.stats[order(AmRep.RT.stats$RT.median, decreasing = F),]
      File.name <- tclvalue(tkgetSaveFile(initialfile=File.name))
      write.csv(AmRep.RT.stats, file = File.name,row.names = FALSE)
      print(paste("Identificaiton report was created and save in",workdir, sep=" "))
    }
  }
  
  
  
  ### MERAGE TWO FILE TOGETHER ## remove
  
  merage_RT_standard<-function(workdir=workdir,
                               File.name = "Identification Report (TMS) with StandardMix_RT.csv"
  ){
    workdir <- workdir
    setwd(workdir)
    data.df=read.csv(tk_choose.files(caption = "Load Samples Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    standard.df=read.csv(tk_choose.files(caption = "Load Standard Mixture Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    RT_standard_names <- tk_select.list(names(standard.df), multiple = FALSE, title = paste("Select Rentention time"))
    SelectStandard<-standard.df[,RT_standard_names]
    standard.df.final<-data.frame(cbind(Name=standard.df$Name,StandardMixture_RT=SelectStandard))
    Merge.data <- merge(data.df, standard.df.final, by.x = "Name",
                        by.y = "Name", all = TRUE)
    
    Merge.data  <- Merge.data [order(Merge.data $RT.median, decreasing = F),]
    File.name <- tclvalue(tkgetSaveFile(initialfile=File.name))
    write.csv(Merge.data, file = File.name,row.names = FALSE)
    print(paste("Identificaiton report with StandardMix RT was created and save in",workdir))
    
  }
  
  
  
  Tips<-function () {
    require(tcltk)
    kk <- tktoplevel()
    tktitle(kk) <- "GC-Chromatogram of TMS standard mixtures"
    image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/TMSstandardMix.gif", sep=""))
    FrontImage<-tcl("image",  "create", "photo", "MaxC2.image", file=image.path)
    tkgrid(tklabel(kk,image=FrontImage))
  }
  
  
  RunMCF<-function(save = TRUE,
                   output = "Final Data Process Result",
                   Ion.bin= 0.5,
                   workdir=workdir,
                   RT.bin = 5,
                   peak="Peak Height",
                   FileName="Final Result (Backup)",
                   GGReport="Slow"
  ) {
    
    require(xcms)
    require(lattice)
    require(plyr)
    setwd(workdir)
    final.check.data <- read.csv(tk_choose.files(caption = "Load a corrected Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    
    check.ref_ion<-table(is.na(final.check.data$ref_ion))["TRUE"]
    check.ref_ion[is.na(check.ref_ion)]<-0
    if(check.ref_ion>=1){tkmessageBox(message = "ERROR: Missing reference ions in the Identification Report!", icon = "warning", type = "ok"); stop("Missing reference ion")}
    
    
    
    ion.lib<-final.check.data
    files <- c(dir())
    info <- file.info(files)
    isDir <- info$isdir
    conditions_pre1 <- c(files[isDir == TRUE])
    
    conditions_pre<-c(conditions_pre1 ,"Manually choose folder")
    
    conditions <- tk_select.list(conditions_pre, multiple = FALSE,
                                 title = paste("Folder contains all cdf files"))
    
    if (conditions=="Manually choose folder"){ conditions<-tk_choose.dir(default = "", caption = "A folder contains all cdf files")}
    
    ###
    library_file <- ion.lib[,c("Name","ref_ion","RT.median","CAS")]
    final.df<-NULL
    remove(final.df)
    Graphic.df<-data.frame()
    
    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]
    for (q in 1:length(filenames)) {
      file1 <- filenames[q]
      name.file <- file1
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      surefinal <- data.frame(Name=ion.lib[,"Name"],Base.Peak=0)
      raw_data <- xcmsRaw(filename = paste(conditions,
                                           filenames[q], sep = "\\"))
      
      # plotChrom(raw_data) # plot TIC_remove because Liz don't like it
      
      for (h in 1:nrow(library_file)) {
        metabolite <-  library_file[h, ]
        R.Ion<-metabolite$ref_ion
        R.Time<-metabolite$RT.median*60
        IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(R.Time-RT.bin,R.Time+RT.bin))
        abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
        ifelse(peak=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-auc(abundance$rt, abundance$intensity*10))
        surefinal$Base.Peak[h] <-  Peakvalue
        
        if (GGReport=="Slow"){
          Graphic.df<-rbind.fill(Graphic.df,data.frame(Metabolite.Names=paste(metabolite$Name," (m/z:",R.Ion,")",sep="") ,Retention.Time=(abundance$rt)/60,
                                                       Intensity=abundance$intensity ,Sample.Names=name.file ))
        }
        
      }
      surefinal <- surefinal[!duplicated(surefinal$Name),]
      names(surefinal)[2] <- name.file
      
      if (!exists("final.df")){
        final.df <- surefinal
        confirmation <- paste("File", q, "(", name.file,")", "done!", sep = " ")
        print(confirmation)
      } else {
        final.df <- cbind(final.df,surefinal[-1])
        confirmation <- paste("File", q, "(", name.file,
                              ")", "done!", sep = " ")
        print(confirmation)
      }
    }
    
    
    # False discover rate
    
    finalX.df<-cbind(False.Positive.ID=NA ,AMDIS.percentage.match =ion.lib$AMDIS.percentage.match , Ref.Ion=ion.lib$ref_ion, Ret.Time=ion.lib$RT.median,final.df)
    
    # Remove the repeat aboundace
    dulpicated_aboundace<-unique(finalX.df[duplicated(finalX.df[,6]),6])
    Repeate.aboundance<-NULL
    if (length(dulpicated_aboundace)!=0){
      for(n in 1: length(dulpicated_aboundace)){
        Repeat.All<- finalX.df[finalX.df[,6]==dulpicated_aboundace[n],]
        Repeate.aboundance <-rbind(Repeate.aboundance, Repeat.All[-which.max(Repeat.All$AMDIS.percentage.match),])
      }
      finalX.df[row.names(finalX.df) %in% row.names(Repeate.aboundance),"False.Positive.ID"]<-"Repeated Aboundance"
    }
    
    # Remove the repated RT
    dulpicated_RT<-unique(finalX.df$Ret.Time[duplicated(finalX.df$Ret.Time)])
    
    Repeate.RT<-NULL
    if (length(dulpicated_RT)!=0){
      for(n in 1: length(dulpicated_RT)){
        Repeat.All<- finalX.df[finalX.df$Ret.Time==dulpicated_RT[n],]
        Repeate.RT <-rbind(Repeate.RT, Repeat.All[-which.max(Repeat.All$AMDIS.percentage.match),])
      }
      if (is.na(match("Repeated Aboundance",Repeate.RT$False.Positive.ID))==FALSE){
        R.Abd.RT<-Repeate.RT[Repeate.RT$False.Positive.ID %in% "Repeated Aboundance",]
        finalX.df[row.names(finalX.df) %in% row.names(R.Abd.RT),"False.Positive.ID"]<-"Repeated Aboundance & RT"
      }
      if (is.na(match(NA,Repeate.RT$False.Positive.ID))==FALSE){
        R.RT<-Repeate.RT[Repeate.RT$False.Positive.ID %in% NA,]
        finalX.df[row.names(finalX.df) %in% row.names(R.RT),"False.Positive.ID"]<-"Repeated RT"
      }
    }
    
    
    finalF.df<-cbind(CAS=ion.lib$CAS, Num.ID=ion.lib$Number.of.ID,  finalX.df)
    finalF.df <- finalF.df[order(finalF.df$Ret.Time , decreasing = F),]
    row.names(finalF.df) <- 1:nrow(finalF.df)
    finalF.df$False.Positive.ID[is.na(finalF.df$False.Positive.ID)]<-""
    
    
    ## Save Files
    ifelse(peak=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
    sheet <- paste(output," (RT bin=", RT.bin,"s_", ValueType ,")", sep="")
    store <- paste(workdir, "\\", sheet, ".csv", sep = "")
    FileName <- tclvalue(tkgetSaveFile(initialfile=store))
    write.csv(finalF.df, file = FileName, row.names = FALSE)
    if (GGReport=="Slow"){write.csv(Graphic.df,file="GC-Peak Diagnosis Report.csv")}
    print(paste("Data process is completed and saved as: ",sheet,".csv", sep=""))
    
    
  }
  
  
  
  
  ## Help page
  
  
  
  TMS.hlep.page<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")
    
    tkinsert(txt,"end","Step 1\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end"," -Prepare a working folder containing:\n")
    tkinsert(txt,"end","      AMDIS Batch report in .txt formate\n")
    tkinsert(txt,"end","      A folder contains all GC-MS raw data in .cdf formate\n")
    tkinsert(txt,"end","      A TMS MS library in .msl formate\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end"," -Suggested AMIDS setting\n")
    tkinsert(txt,"end","      Minimum match factor : 80\n")
    tkinsert(txt,"end","      Deconv. : one, High, Very High, Low\n")
    tkinsert(txt,"end","      Batch Job : include only first 1 hit\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","Step 2\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","Step 3\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }
  
  
  
  ## Building an interfance for lIZ mcf anlaysis
  
  
  library(tcltk)
  require(rgl)
  tclRequire("BWidget")
  
  
  #Setup
  TMS <- tktoplevel()
  frameOverall <- tkframe(TMS)
  tkwm.title(TMS,"GC-MS data analysis (Liz version)")
  
  #All the functions
  workingdirctory<-function(){WorkingDir()}
  Id_report1<-function(){
    Ret.Time.Filter <- as.numeric(tclvalue(RT.filter.var))
    NbVal <- as.character(tclvalue(NbValue))
    if (NbVal=="0"){MsLibrary_s<-"InHouse"}
    if (NbVal=="1"){MsLibrary_s<-"NIST"}
    Id_report(workdir=workdir,MsLibrary=MsLibrary_s)
  }
  Tips1<-function(){Tips()}
  
  RunMCF1<-function(){
    RT.bin_s <-as.numeric(tclvalue(SliderValue))
    Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="0"){GGReport_s<-"Slow"}
    if (cbVal=="1"){GGReport_s<-"Fast"}
    Mzbin<-as.numeric(tclvalue(MzBinValue))
    RunMCF(RT.bin = RT.bin_s,workdir=workdir, peak=Peakvar,GGReport=GGReport_s,Ion.bin=Mzbin)
  }
  merage_RT_standard1<-function(){merage_RT_standard(workdir=workdir)}
  TMS.hlep.page1<-function(){TMS.hlep.page()}
  
  
  
  ### all the button function
  Wdir.but <- tkbutton(TMS, text="Choose", command=workingdirctory,width=6)
  ID_report.but <- tkbutton(TMS, text="Run", command=Id_report1,width=6)
  #Tips.but<-tkbutton(TMS, text="Tips", command=Tips1,width=6)
  integrate.but<-tkbutton(TMS, text="Run", command= RunMCF1,width=6)
  #mergeSM.but<-tkbutton(TMS, text="Run", command= merage_RT_standard1,width=6)
  help.but<-tkbutton(TMS, text="Help", command= TMS.hlep.page1,width=6)
  
  ## fast or slow mode
  cb <- tkcheckbutton(TMS)
  cbValue <- tclVar("0")
  tkconfigure(cb,variable=cbValue)
  
  ## NIST library mode
  Nb <- tkcheckbutton(TMS)
  NbValue <- tclVar("0")
  tkconfigure(Nb,variable=NbValue)
  
  
  
  ## RT filter
  RT.filter.var <- tclVar("2")
  RT.entry <- tkentry(TMS, textvariable=RT.filter.var,width=8)
  
  ## MZ bin
  MzBinValue <- tclVar("0.5")
  MzBin.entry <- tkentry(TMS, textvariable=MzBinValue ,width=8)
  
  ## Upper
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper,text="Setup Working Directory"),Wdir.but, sticky="w")
  tkgrid(tklabel(frameUpper,text="AutoSelect Reference ions (NIST)"),Nb,sticky="w")
  #tkgrid(tklabel(frameUpper,text="Add StandardMix RT (option)"),mergeSM.but,sticky="w ")
  tkgrid(tklabel(frameUpper,text="Retention time filter"), RT.entry,tklabel(frameUpper,text="min  "),sticky="w")
  tkgrid(tklabel(frameUpper,text="Create an Identification Report                     "),ID_report.but,sticky="w")
  
  ## Middle
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  
  
  ## Slider
  
  SliderValue <- tclVar("5")
  SliderValueLabel <- tklabel(TMS,text=as.character(tclvalue(SliderValue)),widt=2)
  tkconfigure(SliderValueLabel,textvariable=SliderValue)
  slider <- tkscale(frameMid, from=1, to=10,
                    showvalue=F, variable=SliderValue,
                    resolution=1, orient="horizontal")
  
  
  ## choose between peakheight or area
  
  Peak <- c("Peak Height","Peak Area")
  comboBox <- tkwidget(TMS,"ComboBox",editable=FALSE,values=Peak,textvariable=tclVar("Peak Height"),width=14)
  
  
  # tkgrid(tklabel(frameMid,text="Manually correct the Identification Report"),Tips.but,sticky="w")
  tkgrid(tklabel(frameMid,text="Fast Mode (No Chromatograms)"),cb,sticky="w")
  tkgrid(tklabel(frameMid,text="Peak height or Peak area"),comboBox,sticky="w")
  tkgrid(tklabel(frameMid,text="Retention time bin (Default = 5 s) "),slider,SliderValueLabel,tklabel(frameMid,text="s"),sticky="w")
  tkgrid(tklabel(frameMid,text="Ion Mass Bin (m/z)"),  MzBin.entry,sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic GC-peak integration"),integrate.but,sticky="w")
  
  ## lower
  
  font.size.var <- tclVar("0.5")
  font.size.entry <- tkentry(TMS, textvariable=font.size.var,width=4)
  
  
  rb1 <- tkradiobutton(TMS)
  rb2 <- tkradiobutton(TMS)
  rbValue <- tclVar("PeakReport")
  tkconfigure(rb1,variable=rbValue,value="PeakReport")
  tkconfigure(rb2,variable=rbValue,value="IonExtractor_Window")
  
  OnOK <- function(){
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="PeakReport"){
      font.size<- as.numeric(tclvalue(font.size.var))
      ifelse(sum(dir()=="GC-Peak Diagnosis Report.csv")==1,PeakReport(font.size=font.size),stop("Please select correct working directory"))
    }
    if (rbVal=="IonExtractor_Window"){
      IonExtractor_Window()}
  }
  OK.but <- tkbutton(TMS,text="Run",command=OnOK,width=6)
  
  frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  #tkgrid(tklabel(frameLower,text="Check TMS results"))
  tkgrid(tklabel(frameLower,text="Overlay all chromatograms                     "),rb1,tklabel(frameLower,text="Font Size "),font.size.entry, sticky="w")
  tkgrid(tklabel(frameLower,text="IonExtractor Single Mode "), rb2,OK.but,sticky="w")
  
  
  quit.but <- tkbutton(TMS, text = "Close Session",
                       command = function() {
                         tkdestroy(TMS)
                         rgl.close()
                       }
  )
  
  ## Open working directory
  
  open.but <- tkbutton(TMS,text="Open folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )
  
  ## load images
  image.path1<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_SIlas_logo.gif", sep=""))
  logo<-tcl("image",  "create", "photo", "MaxC3.image", file=image.path1)
  
  
  
  
  button <- tkframe(frameOverall,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text="                              "),quit.but,open.but, tklabel(button ,text="                    "),help.but)
  
  
  ## Userinterphase
  
  tkgrid(tklabel(frameOverall,image=logo))
  tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="                           STEP 1. Setup"),sticky="w")
  tkgrid(frameUpper,pady= 10, padx= 10)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="                           STEP 2. Integration"),sticky="w")
  tkgrid(frameMid,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="                           STEP 3. Check GC-MS results"),sticky="w")
  tkgrid(frameLower,pady= 10, padx= 10)
  tkgrid(frameOverall)
  #tkgrid(quit.but,help.but, pady= 10, padx= 10)
  tkgrid(button,pady= 10, padx= 10)
  tkgrid(tklabel(frameOverall,text="    "))
}

#' create_sub_library
#'
#' This is a function for peak integration according to the IDs in identification summary
#'
#' @param output specify the output file name initial 
#' @param Ion.bin the mz tolerance (Da) for peak integration
#' @param workdir locate the working Dir
#' @param intensity_type the type of signal to be reported \code{"Peak Height"} or \code{"Peak Area"}
#' @param GGReport the reporting mode \code{"Slow"} plot every integrated peak, or \code{"Fast"} generate an overall report 
#'
#' @return None
#'
#' @examples
#' GCMS_integration()
#'
#' @export
######################################### create_sub_library #########################################
create_sub_library <- function(workdir=tk_choose.dir(caption = "Select working directory")){
  
  require(xcms)
  require(lattice)
  require(plyr)
  require(flux)
  
  WorkingDir <-function(workdir = tk_choose.dir(caption = "Select working directory")){
    workdir <<- workdir
    setwd(workdir)
  }
  
  
  Id_report<-function(      workdir=workdir,
                            MS.L= tk_choose.files(caption="Select SUb-NIST library in .MSL"),
                            amdis.report = read.csv(tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                                    multi = FALSE, caption = "Select the AMDIS report in .TXT"),
                                                    sep = "\t", stringsAsFactors=FALSE),
                            File.name = "Identification Report (Liz version).csv",
                            Ret.Time.Filter=2){
    
    workdir <- workdir
    setwd(workdir)
    lib.txt<-readLines(MS.L)
    lib.txt<-lib.txt[-grep("^[ ]*$",lib.txt)]
    att.nms<-unique(sapply(strsplit(lib.txt[grep(":",lib.txt)],":"),function(x) x[1]))
    entry.nms<-sapply(strsplit(lib.txt[grep("NAME:",lib.txt)],"ME:"), function(x) x[2])
    rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))
    
    starts<-grep("NAME:", lib.txt)
    stops<-c(starts[1:(length(starts)-1)]+diff(starts)-1,length(lib.txt))
    
    
    # function for select right digit
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    reference_ion<-NULL
    
    for (i in 1:length(starts)){
      tmp<-lib.txt[starts[i]:stops[i]]
      reference_ion<-rbind(reference_ion, (as.numeric(substrRight(lapply(strsplit(grep("1000)",tmp,value=TRUE)," 1000)"), function(x) x[1]),3)))[1])
      
      sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
    }
    
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez$NAME<-entry.nms
    
    
    libr<-cbind(rez,reference_ion)## returns list of peaks
    libr$NAME <- gsub("?", "", libr$NAME, fixed = TRUE)
    libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
    
    
    ############ Step 2. Generate the list of metabolites and clean up
    AmRep<-amdis.report
    # AmRep<-AmRep[!AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)+Ret.Time.Filter & AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)-Ret.Time.Filter,]## check this one ##
    AmRep$Name <- gsub("?", "", AmRep$Name, fixed = TRUE)
    AmRep$Name <- gsub("^ ", "", AmRep$Name, perl=T)
    
    ## RT statistics from the AReport
    RT.stats<-t(sapply(split(AmRep$RT, AmRep$Name),function(x) c(RT.median=round(median(x,na.rm=T),3),
                                                                 RT.mean=round(mean(x,na.rm=TRUE),3),
                                                                 RT.sd=sd(x,na.rm=T))))
    
    
    ID.stats<-t(sapply(split(AmRep$Net, AmRep$Name),function(x) c(AMDIS.percentage.match=round(median(x,na.rm=T)),
                                                                  Number.of.ID=length(x[!is.na(x)]))))
    
    
    #AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,ID.stats, Exp.RT=0,Diff.RT=0, ref_ion=0))##data frame for list of metabolites
    
    AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,ID.stats, ref_ion=0, CAS=0))##data frame for list of metabolites
    AmRep.RT.stats$ref_ion <-  libr$reference_ion[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    #AmRep.RT.stats$Exp.RT <- as.numeric(t(sapply(split(AmRep$Expec..RT, AmRep$Name), median)))
    #AmRep.RT.stats$Diff.RT <- AmRep.RT.stats$Exp.RT- AmRep.RT.stats$RT.median
    AmRep.RT.stats$CAS <-  libr$CAS[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    AmRep.RT.stats$Name <- rownames(RT.stats)
    AmRep.RT.stats <- AmRep.RT.stats[order(AmRep.RT.stats$RT.median, decreasing = F),]
    write.csv(AmRep.RT.stats, file = File.name,row.names = FALSE)
    print(paste("Identificaiton report was created and save in",workdir, sep=" "))
  }
  
  
  
  
  ### MERAGE TWO FILE TOGETHER ## remove
  
  
  
  
  
  
  Tips<-function () {
    require(tcltk)
    kk <- tktoplevel()
    tktitle(kk) <- "GC-Chromatogram of TMS standard mixtures"
    image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/TMSstandardMix.gif", sep=""))
    FrontImage<-tcl("image",  "create", "photo", "MaxC2.image", file=image.path)
    tkgrid(tklabel(kk,image=FrontImage))
  }
  
  
  RunTMS<-function(save = TRUE,
                   output = "Final Data Process Result",
                   Ion.bin= 0.5,
                   workdir=workdir,
                   RT.bin = 5,
                   peak="Peak Height",
                   GGReport="Slow"
  ) {
    
    require(xcms)
    require(lattice)
    require(plyr)
    setwd(workdir)
    final.check.data <- read.csv(tk_choose.files(caption = "Load a corrected Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    
    check.ref_ion<-table(is.na(final.check.data$ref_ion))["TRUE"]
    check.ref_ion[is.na(check.ref_ion)]<-0
    if(check.ref_ion>=1){print("ERROR: Missing reference ions in the Identification Report")}
    
    
    ion.lib<-final.check.data
    files <- c(dir())
    info <- file.info(files)
    isDir <- info$isdir
    conditions_pre <- c(files[isDir == TRUE])
    
    
    conditions <- tk_select.list(conditions_pre, multiple = FALSE,
                                 title = paste("Folder contains all cdf files"))
    
    
    
    ###
    library_file <- ion.lib[,c("Name","ref_ion","RT.median")]
    final.df<-NULL
    remove(final.df)
    Graphic.df<-data.frame()
    
    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]
    for (q in 1:length(filenames)) {
      file1 <- filenames[q]
      name.file <- file1
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      surefinal <- data.frame(Name=ion.lib[,"Name"],Base.Peak=0)
      raw_data <- xcmsRaw(filename = paste(conditions,
                                           filenames[q], sep = "\\"))
      
      # plotChrom(raw_data) # plot TIC_remove because Liz don't like it
      
      for (h in 1:nrow(library_file)) {
        metabolite <-  library_file[h, ]
        R.Ion<-metabolite$ref_ion
        R.Time<-metabolite$RT.median*60
        IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(R.Time-RT.bin,R.Time+RT.bin))
        abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
        ifelse(peak=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-auc(abundance$rt, abundance$intensity*10))
        surefinal$Base.Peak[h] <-  Peakvalue
        
        ## create overlay chromatograms
        
        if (GGReport=="Slow"){
          Graphic.df<-rbind.fill(Graphic.df,data.frame(Metabolite.Names=paste(metabolite$Name," (m/z:",R.Ion,")",sep="") ,Retention.Time=(abundance$rt)/60,
                                                       Intensity=abundance$intensity ,Sample.Names=name.file ))
        }
      }
      surefinal <- surefinal[!duplicated(surefinal$Name),]
      names(surefinal)[2] <- name.file
      
      if (!exists("final.df")){
        final.df <- surefinal
        confirmation <- paste("File", q, "(", name.file,")", "done!", sep = " ")
        print(confirmation)
      } else {
        final.df <- cbind(final.df,surefinal[-1])
        confirmation <- paste("File", q, "(", name.file,
                              ")", "done!", sep = " ")
        print(confirmation)
      }
    }
    final.df<-cbind( CAS=ion.lib$CAS, Num.ID=ion.lib$Number.of.ID, AMDIS.percentage.match =ion.lib$AMDIS.percentage.match , Ref.Ion=ion.lib$ref_ion, Ret.Time=ion.lib$RT.median,final.df)
    final.df <- final.df[order(final.df$Ret.Time , decreasing = F),]
    row.names(final.df) <- 1:nrow(final.df)
    
    
    
    if (save == TRUE) {
      ifelse(peak=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
      sheet <- paste(output," (RT bin=", RT.bin,"s_", ValueType ,")", sep="")
      store <- paste(workdir, "\\", sheet, ".csv", sep = "")
      write.csv(final.df, file = store, row.names = FALSE)
      
      if (GGReport=="Slow"){write.csv(Graphic.df,file="GC-Peak Diagnosis Report.csv")}
      print(paste("Data process is completed and saved as: ",sheet,".csv", sep=""))
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  ## Help page
  
  
  
  TMS.hlep.page<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")
    
    tkinsert(txt,"end","Step 1\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end"," -Prepare a working folder containing:\n")
    tkinsert(txt,"end","      AMDIS Batch report in .txt formate\n")
    tkinsert(txt,"end","      A folder contains all GC-MS raw data in .cdf formate\n")
    tkinsert(txt,"end","      A TMS MS library in .msl formate\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end"," -Suggested AMIDS setting\n")
    tkinsert(txt,"end","      Minimum match factor : 80\n")
    tkinsert(txt,"end","      Deconv. : one, High, Very High, Low\n")
    tkinsert(txt,"end","      Batch Job : include only first 1 hit\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","Step 2\n")
    tkinsert(txt,"end","  \n")
    tkinsert(txt,"end","Step 3\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }
  
  
  ## HOw to make MSP library
  Help1<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")
    
    tkinsert(txt,"end"," -Generate MSP library\n")
    tkinsert(txt,"end","      Open: Lib2NIST Converter\n")
    tkinsert(txt,"end","           Select: NIST05~9a.L library\n")
    tkinsert(txt,"end","            Click: Define subsite\n")
    tkinsert(txt,"end","                   under subset Source choose entry number from file in .txt\n")
    tkinsert(txt,"end","                   Make sure SUbset Type is ID numbers\n")
    tkinsert(txt,"end","            Click: OK\n")
    tkinsert(txt,"end","      Choose a correct NIST library under Input Libraries window\n")
    tkinsert(txt,"end","      Choose an appropriate Output location\n")
    tkinsert(txt,"end","      Click: Convert\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }
  
  
  ## HOw to convert MSP to MSL library
  Help2<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")
    
    
    tkinsert(txt,"end"," -Convert MSP library to MSL library\n")
    tkinsert(txt,"end","    Open: AMDIS software\n")
    tkinsert(txt,"end","          click: library\n")
    tkinsert(txt,"end","          Click: library Transfer...\n")
    tkinsert(txt,"end","          Click: Files...\n")
    tkinsert(txt,"end","          Click: Load Library...(under source) and select MSP library\n")
    tkinsert(txt,"end","          Click: Creat New Library...(under Destination)\n")
    tkinsert(txt,"end","                 and creat a new MS library in MSL formate \n")
    tkinsert(txt,"end","          Select all compounds from MSP library and transfer to MSP library\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }
  
  
  ## HOw to generate AMDIS report
  Help3<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")
    
    tkinsert(txt,"end"," -Use MSL library generated previously for compound identificaitons\n")
    tkinsert(txt,"end"," -Suggested AMIDS setting\n")
    tkinsert(txt,"end","      Minimum match factor : ~45\n")
    tkinsert(txt,"end","      Deconv. : one, low, Very Low, High\n")
    tkinsert(txt,"end","      Batch Job : include only first 1 hit\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }
  
  ## combine library
  
  combine.ms.library <- function(){
    require(tcltk)
    workdir <- tk_choose.dir(caption = "Select folder where MS library is located")
    setwd(workdir)
    NIST.library <- readLines(tk_choose.files(multi = FALSE, caption = "Select the NIST library in .msl"))
    Inhouse.library <- readLines(tk_choose.files(multi = FALSE, caption = "Select the In-house library in .msl"))
    lib_md <- gsub("NAME:","NAME:(NIST)",NIST.library)
    combin.library <- cbind(c(Inhouse.library, lib_md ))
    File.nameNIST <- "Inhouse_NIST_library.msl"
    File.nameNIST <- tclvalue(tkgetSaveFile(initialfile=File.nameNIST))
    write.table(combin.library , file= File.nameNIST,row.names=FALSE, col.names=FALSE, quote =FALSE)
  }
  
  
  ## Building an interfance for lIZ mcf anlaysis
  
  library(tcltk)
  library(tcltk2)
  require(rgl)
  tclRequire("BWidget")
  
  TMS <- tktoplevel()
  frameOverall <- tkframe(TMS)
  tkwm.title(TMS,"Create a sub-NIST library")
  
  
  
  ## fast or slow mode
  
  cb <- tkcheckbutton(TMS)
  cbValue <- tclVar("0")
  tkconfigure(cb,variable=cbValue)
  
  
  workingdirctory<-function(){WorkingDir()}
  Id_report1<-function(){
    Ret.Time.Filter <- as.numeric(tclvalue(RT.filter.var))
    Id_report(workdir=workdir,Ret.Time.Filter=Ret.Time.Filter)
  }
  Tips1<-function(){Tips()}
  RunTMS1<-function(){
    RT.bin_s <-as.numeric(tclvalue(SliderValue))
    Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="0"){GGReport_s<-"Slow"}
    if (cbVal=="1"){GGReport_s<-"Fast"}
    RunTMS(RT.bin = RT.bin_s,workdir=workdir, peak=Peakvar,GGReport=GGReport_s)
  }
  merage_RT_standard1<-function(){merage_RT_standard(workdir=workdir)}
  TMS.hlep.page1<-function(){TMS.hlep.page()}
  ChemstationLibraryEntry1<-function(){
    
    
    ChemstationLibraryEntry(workdir = getwd())
  }
  
  
  
  MSP_generate_lib2nist1<-function(){
    cbVal <- as.character(tclvalue(RI_retrievecbValue))
    if (cbVal=="0"){RI_retrieve<-FALSE}
    if (cbVal=="1"){RI_retrieve<-TRUE}
    RI_column_type=RI_column_typecbValue[[as.numeric(tclvalue(tcl(RI_column_typecb,"getvalue")))+1]]
    message(paste("Retrieve RI from NIST website",RI_retrieve,", Column: ",RI_column_type))
    cbVal <- as.character(tclvalue(convert_MSPtoMSLcbcbValue))
    if (cbVal=="0"){convert<-FALSE}
    if (cbVal=="1"){convert<-TRUE}
    entryfile=tcltk::tk_choose.files(default = "merged_library_entries_File.csv",
                                     multi = F, caption = "Select the merged library entries file in csv format")
    message(paste("Generate msp for",entryfile,"\nusing lib2nist.exe in",as.character(tclGetValue("lib2nistpath")),"\nGenerate MSL library as well:",convert))
    lib2nistpath=as.character(tclGetValue("lib2nistpath"))
    MSP_generate_lib2nist(unique_Entry_no=entryfile,lib2nistpath = lib2nistpath,convert_MSPtoMSL=convert,RI_retrieve=RI_retrieve,RI_column_type=RI_column_type)
  }
  
  RI_retrievecb <- tkcheckbutton(TMS)
  RI_retrievecbValue <- tclVar("0")
  tkconfigure(RI_retrievecb,variable=RI_retrievecbValue)
  
  RI_column_typecbValue <- (c("5ms","wax"))
  RI_column_typecb <- tkwidget(background="white",TMS ,"ComboBox",editable=FALSE,values=RI_column_typecbValue ,textvariable=tclVar("5ms"),width=9)
  
  msp_to_msl1<-function(){
    cbVal <- as.character(tclvalue(RI_retrievecbValue))
    if (cbVal=="0"){RI_retrieve<-FALSE}
    if (cbVal=="1"){RI_retrieve<-TRUE}
    RI_column_type=RI_column_typecbValue[[as.numeric(tclvalue(tcl(RI_column_typecb,"getvalue")))+1]]
    message(paste("Retrieve RI from NIST website",RI_retrieve,", Column: ",RI_column_type))
    msp_to_msl(RI_retrieve=RI_retrieve,RI_column_type=RI_column_type)
  }
  
  #tclSetValue("convert_MSPtoMSL",0)
  convert_MSPtoMSLcb <- tkcheckbutton(TMS)
  convert_MSPtoMSLcbcbValue <- tclVar("0")
  tkconfigure(convert_MSPtoMSLcb,variable=convert_MSPtoMSLcbcbValue)
  
  convert_MSPtoMSLfun.but <- tkbutton(TMS, text="MSP->MSL", command=msp_to_msl1,width=10)
  
  lib2nistpath="C:\\NIST17\\MSSEARCH"
  tclSetValue("lib2nistpath",lib2nistpath)
  lib2nistpathentry <- tkentry(background="white",TMS,width="12",textvariable="lib2nistpath",validate="key")
  ### all the button function
  Wdir.but <- tkbutton(TMS, text="Choose", command=workingdirctory,width=6)
  ID_report.but <- tkbutton(TMS, text="Run", command=Id_report1,width=6)
  #Tips.but<-tkbutton(TMS, text="Tips", command=Tips1,width=6)
  integrate.but<-tkbutton(TMS, text="Run", command= RunTMS1,width=6)
  #mergeSM.but<-tkbutton(TMS, text="Run", command= merage_RT_standard1,width=6)
  help.but<-tkbutton(TMS, text="Help", command= TMS.hlep.page1,width=6)
  
  combin.but <- tkbutton(TMS, text="Run", command= combine.ms.library, width=6)
  
  
  
  LibraryEntry.but<-tkbutton(TMS, text="Run", command=ChemstationLibraryEntry1 ,width=6)
  
  GeneratMSP.but1<-tkbutton(TMS, text="RUN", command=MSP_generate_lib2nist1 ,width=6)
  help.but2<-tkbutton(TMS, text="Tips", command=Help2 ,width=6)
  help.but3<-tkbutton(TMS, text="Tips", command=Help3 ,width=6)
  
  RT.filter.var <- tclVar("2")
  RT.entry <- tkentry(TMS, textvariable=RT.filter.var,width=6)
  
  
  ## Upper 1
  frameUpper1 <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper1,text="Setup Working Directory        "),Wdir.but,sticky="w")
  tkgrid(tklabel(frameUpper1,text="Generate entry #"), LibraryEntry.but, sticky="w")
  
  
  ## Upper 2
  frameUpper2 <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper2,text="Retrieve RI from NIST"), RI_retrievecb,sticky="W")
  tkgrid(tklabel(frameUpper2,text="Select COlumn typr for RI retrieving"), RI_column_typecb,sticky="W")
  tkgrid(tklabel(frameUpper2,text="Specify NIST folder"),lib2nistpathentry,sticky="w")
  tkgrid(tklabel(frameUpper2,text="Generate MSP library"),GeneratMSP.but1, sticky="w")
  tkgrid(tklabel(frameUpper2,text="Convert MSP to MSL library"),convert_MSPtoMSLcb, sticky="w")
  tkgrid(tklabel(frameUpper2,text="Already have the MSP library?"),convert_MSPtoMSLfun.but, sticky="w")
  
  ## Upper 2b
  
  combine.ms.library
  frameUpper2b <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper2b,text="Combine two libraries              "), combin.but, sticky="w")
  
  
  
  ## Upper 3
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  
  #tkgrid(tklabel(frameUpper,text="Add StandardMix RT (option)"),mergeSM.but,sticky="w")
  #tkgrid(tklabel(frameUpper,text="Retention time filter"), RT.entry,tklabel(frameUpper,text="min"),sticky="w")
  tkgrid(tklabel(frameUpper,text="Generate AMDIS report"),help.but3, sticky="w")
  tkgrid(tklabel(frameUpper,text="Create an Identification Report "),ID_report.but,sticky="w")
  
  
  ## Middle
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  
  
  ## Slider
  
  SliderValue <- tclVar("5")
  SliderValueLabel <- tklabel(TMS,text=as.character(tclvalue(SliderValue)),width=2)
  tkconfigure(SliderValueLabel,textvariable=SliderValue)
  slider <- tkscale(frameMid, from=1, to=10,
                    showvalue=F, variable=SliderValue,
                    resolution=1, orient="horizontal")
  
  
  
  
  
  
  ## choose between peakheight or area
  
  Peak <- c("Peak Height","Peak Area")
  comboBox <- tkwidget(TMS,"ComboBox",editable=FALSE,values=Peak,textvariable=tclVar("Peak Height"),width=10)
  
  
  # tkgrid(tklabel(frameMid,text="Manually correct the Identification Report"),Tips.but,sticky="w")
  tkgrid(tklabel(frameMid,text="Fast Mode (No Overlay chromatograms)"),cb,sticky="w")
  tkgrid(tklabel(frameMid,text="Peak height or Peak area"),comboBox,sticky="w")
  tkgrid(tklabel(frameMid,text="Retention time bin (Default = 5 s) "),slider,SliderValueLabel,tklabel(frameMid,text="s"),sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic GC-peak integration"),integrate.but,sticky="w")
  
  ## lower
  
  font.size.var <- tclVar("0.5")
  font.size.entry <- tkentry(TMS, textvariable=font.size.var,width=4)
  
  
  rb1 <- tkradiobutton(TMS)
  rb2 <- tkradiobutton(TMS)
  rbValue <- tclVar("PeakReport")
  tkconfigure(rb1,variable=rbValue,value="PeakReport")
  tkconfigure(rb2,variable=rbValue,value="IonExtractor_Window")
  
  OnOK <- function(){
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="PeakReport"){
      font.size<- as.numeric(tclvalue(font.size.var))
      ifelse(sum(dir()=="GC-Peak Diagnosis Report.csv")==1,PeakReport(font.size=font.size),stop("Please select correct working directory"))
    }
    if (rbVal=="IonExtractor_Window"){
      IonExtractor_Window()}
  }
  OK.but <- tkbutton(TMS,text="RUN",command=OnOK,width=6)
  
  frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  #tkgrid(tklabel(frameLower,text="Check TMS results"))
  tkgrid(tklabel(frameLower,text="Overlay all chromatograms "),rb1,tklabel(frameLower,text="Font Size "),font.size.entry, sticky="w")
  tkgrid(tklabel(frameLower,text="IonExtractor Single Mode "), rb2,OK.but,sticky="w")
  
  
  quit.but <- tkbutton(TMS, text = "Close Session",
                       command = function() {
                         tkdestroy(TMS)
                         rgl.close()
                       }
  )
  
  ## Open working directory
  
  open.but <- tkbutton(TMS,text="Open folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )
  
  ## load images
  image.path1<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_SIlas_logo.gif", sep=""))
  logo<-tcl("image",  "create", "photo", "MaxC3.image", file=image.path1)
  
  
  
  
  button <- tkframe(frameOverall,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text="                              "),quit.but,open.but, tklabel(button ,text="                    "),help.but)
  
  
  ## Userinterphase
  
  tkgrid(tklabel(frameOverall,image=logo))
  tkgrid(tklabel(frameOverall,text="                          STEP 1. Setup and generate entry #"),sticky="w")
  tkgrid(frameUpper1,pady= 10, padx= 10)
  tkgrid(tklabel(frameOverall,text="                          STEP 2. Make a MS library"),sticky="w")
  tkgrid(frameUpper2,pady= 10, padx= 10)
  tkgrid(tklabel(frameOverall,text="                          Optional: Combine NIST and in-house libraries"),sticky="w")
  tkgrid(frameUpper2b,pady= 10, padx= 10)
  #tkgrid(tklabel(frameOverall,text="                                            STEP 3. Generate ID report"),sticky="w")
  #tkgrid(frameUpper,pady= 10, padx= 10)
  #tkgrid(tklabel(frameOverall,text="    "))
  #tkgrid(tklabel(frameOverall,text="                                            STEP 4. Integration"),sticky="w")
  #tkgrid(frameMid,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall,text="    "))
  #tkgrid(tklabel(frameOverall,text="                                            STEP 5. Check GC-MS results"),sticky="w")
  #tkgrid(frameLower,pady= 10, padx= 10)
  tkgrid(frameOverall)
  #tkgrid(quit.but,help.but, pady= 10, padx= 10)
  tkgrid(button,pady= 10, padx= 10)
}

ChemstationLibraryEntry<-function(workdir=NULL,rootdir=dirname(workdir)){
  # author: vfan
  # date: 24/07/2014
  
  # merge multiple ChemStation .xls output files into a combined file containg
  # FileName, CAS.Number, Name, Retention.Time.min, Width, Peak.Area, Baseline.Height,
  
  
  
  # Peak.Height, Library.Match.Quality, Mol.Weight.amu
  
  # Note: you need the xlsx library
  # IMPORTANT: this code ASSUMES that the relevant info starts at row 9 in the
  #            ChemStation outfiles if the version of ChemStation changes, you will
  #            NEED to check if this still correct
  
  # Questions? email Vicky: v.fan@auckland.ac.nz
  
  library(tcltk)
  if (!dir.exists(workdir)){
  workdir=tcltk::tk_choose.dir(caption = "Select Working Directory that contains .D file or lib report csv")  
  }
  
  message(workdir)
  setwd(workdir)
  
  # Modified by Morgan
  
  files <- dir()
  info <- file.info(files)
  isDir <- info$isdir
  conditions_pre <- c(files[isDir == TRUE])
  
  #chemStationFileLoc<- tk_select.list(conditions_pre, multiple = FALSE,title = paste("Folder contains all xls files"))
  chemStationFileLoc=workdir
  
  
  
  # --- METHODS ---
  
  # need to copy down "Retention.Time.min", "Width", "Peak.Area", "Baseline.Height", "Peak.Height"
  # as only the 1st line of the series contains it.
  # function to copy down the above 5 fields for the identified compounds for each compound
  # assume that eat Hit.No ALWAYS starts at 1
  
  # arguments: infile.df - the merged and renamed infile
  
  processDF <- function(infile.df){
    toProcess.df <- infile.df
    
    compound.data = c()
    for(i in 1:nrow(toProcess.df)){
      # if first compound of set, copy data
      if( is.na(toProcess.df[i,"Hit.No"])==T){}else{
        if(toProcess.df[i,"Hit.No"] == 1){
          compound.data <- toProcess.df[i,c("Retention.Time.min", "Width", "Peak.Area", "Baseline.Height", "Peak.Height")]
        } else{
          toProcess.df[i,c("Retention.Time.min", "Width", "Peak.Area", "Baseline.Height", "Peak.Height")] <- compound.data
        }
      }} # end for
    return(toProcess.df)
  } # end processDF
  
  
  # --- END METHODS ---
  
  
  
  library(rJava)
  library(data.table)
  library(xlsx)
  
  # Set "chekStationFileLoc" to directory location to where ChemStation outfiles are.
  # There should ONLY be the ChemStation xls files in this directory
  
  #setwd(chemStationFileLoc)
  message("Data analysis in progress...")
  allFiles <- dir(chemStationFileLoc)
  
  excelFiles <- allFiles[grepl(".xls", allFiles)]
  csvFilesinfolder <- allFiles[grepl(".csv", allFiles)]
  csvFiles <-paste0(conditions_pre,"/MSRep.csv")
  potentialexcelFiles<-paste0(conditions_pre,"/MSRep.xls")
  potentialexcelFiles=potentialexcelFiles[file.exists(potentialexcelFiles)]
  excelFiles <- c(excelFiles ,potentialexcelFiles)
  csvFiles <-csvFiles[file.exists(csvFiles)]
  csvFiles<-c(csvFiles,csvFilesinfolder)
  
  newMergedFile <- c()
  libraryname=""
  message(paste(length(excelFiles),"Xls files found, first 5 files:"))
  message(paste(excelFiles[1:5],collaps=" "))
  #message(paste(length(excelFiles),"files in total"))
  message(paste(length(csvFiles),"csv files found, first 5 files:"))
  message(paste(csvFiles[1:5],collaps=" "))
  #message(paste(length(csvFiles),"files in total"))
  message("")
  TB <- txtProgressBar(min = 0, max = length(excelFiles)+length(csvFiles), style = 3, width = 50)
  if (length(excelFiles)>=1){
      for(i in 1:length(excelFiles)){
    
    infile <- xlsx::read.xlsx(file=paste(chemStationFileLoc,excelFiles[i], sep = "/"), sheetName="LibRes", startRow=9)
    if (!is.null(infile)){
      if(libraryname==""){libraryname=as.character(infile[1,"Library"])}
      libraryname=unique(c(libraryname,infile[,"Library"]))
      fileName <- gsub("", replacement="", x=excelFiles[i])
      
      #infile.reduced <- infile[,-c(1,3,8,13,14)]
      infile.reduced <- infile[,-c(3)]
      colnames(infile.reduced) <- c("Compound", "Retention.Time.min", "Peak.Area", "Baseline.Height",
                                    "Peak.Height", "Width", "Hit.No", "Name", "Library.Match.Quality",
                                    "Mol.Weight.amu", "CAS.Number", "Library", "Entry.No.Lib")
      infile.reduced=infile.reduced[infile.reduced[["Entry.No.Lib"]]!="",]
      infile.reduced$FileName <- fileName
      
      infile.reduced <- infile.reduced[c("Compound", "Hit.No", "FileName", "CAS.Number", "Name",
                                         "Retention.Time.min", "Width", "Peak.Area",
                                         "Baseline.Height", "Peak.Height", "Library.Match.Quality",
                                         "Mol.Weight.amu", "Entry.No.Lib", "Library")]
      
      newMergedFile <- rbind(newMergedFile, infile.reduced)
    }
    setTxtProgressBar(TB, i)
  } # end for(i)
  newMergedFile <- processDF(newMergedFile)
  }
  if (length(csvFiles)>=1){
  for(i in 1:length(csvFiles)){
    
    #infile <- read.csv(file=paste(chemStationFileLoc,csvFiles[i], sep = "\\"), skip=8,header=T)
    #infile <- fread(file=paste(chemStationFileLoc,csvFiles[i], sep = "/"), skip=8,header=T,sep=",",select = c("CPD" , "RT", "Scan number" , "Area" , "Baseline Heigth"  ,"Absolute Heigth" , "Peak Width min","Hit Number" , "Quality"  ,  "Mol Weight" ,"CAS Number" , "Library", "Entry Number Library" ))
    infile <- read.delim(file=paste(chemStationFileLoc,csvFiles[i], sep = "/"), skip=8,quote = "",stringsAsFactors = F,header = F)
    infile <- as.list(t(infile))
    infile1 <- sapply(infile, function(x){

      ((strsplit(x,split ='(,)(?=(?:[^"]|"[^"]*")*$)',perl = TRUE)[[1]][1:14]))
      })
    infile1<-t(infile1)
    colnames(infile1)<-infile1[1,]
    infile1<-infile1[-1,]
    infile1[,1:10]=as.numeric(infile1[,1:10])
    infile1[,11:14]=as.character(infile1[,11:14])
    infile=na.omit(infile1)
    if('&'(ncol(infile)==14,colnames(infile)[1]=="CPD")){
      
      if(libraryname==""){libraryname=as.character(infile[1,"Library"])}
      libraryname=unique(c(libraryname,infile[,"Library"]))
      fileName <- gsub("", replacement="", x=csvFiles[i])
      
      #infile.reduced <- infile[,-c(1,3,8,13,14)]
      infile.reduced <- data.frame(infile[,-c(3)],stringsAsFactors = F)
      colnames(infile.reduced) <- c("Compound", "Retention.Time.min", "Peak.Area", "Baseline.Height",
                                    "Peak.Height", "Width", "Hit.No", "Library.Match.Quality",
                                    "Mol.Weight.amu", "CAS.Number", "Library","Entry.No.Lib", "Name")
      
      infile.reduced=infile.reduced[infile.reduced[,"Entry.No.Lib"]!="",]
      infile.reduced$FileName <- fileName
      
      infile.reduced <- infile.reduced[c("Compound", "Hit.No", "FileName", "CAS.Number", "Name",
                                         "Retention.Time.min", "Width", "Peak.Area",
                                         "Baseline.Height", "Peak.Height", "Library.Match.Quality",
                                         "Mol.Weight.amu", "Entry.No.Lib", "Library")]
      
      newMergedFile <- rbind(newMergedFile, infile.reduced)
    }
    setTxtProgressBar(TB, i+length(excelFiles))
  }
  }
  #newData <- processDF(infile.reduced)
  newMergedFile$Name=paste0("\"",as.character(newMergedFile$Name),"\"")
  newMergedFile$Name=as.character(gsub("\"","",newMergedFile$Name))
  newMergedFile$Name=as.character(gsub(",","",newMergedFile$Name))
  #newMergedFile$Name=as.character(windows_filename2(newMergedFile$Name))
  newData <- newMergedFile

  mergedOutFile = paste0(rootdir,"/merged_library_entries_File.csv")
  
  # remove "Compound" and "Hit.No"
  write.csv(newData[,-c(1,2)], mergedOutFile, row.names=FALSE, na="", quote=F)
  
  # write out "Entry.No.Lib" with duplicates removed in a separate file
  #unique.Entry.No <- unique(newData[,"Entry.No.Lib"])
  libraryname<<-libraryname
  library_list<-unique(newData[,c("Entry.No.Lib", "Library")])
  
  for (libraryname in unique(library_list$Library)){
  unique.Entry.No <- library_list[library_list$Library==libraryname,"Entry.No.Lib"]
  #unique.Entry.No<-c(libraryname,unique.Entry.No)
  libraryname<-gsub("\\\\","/",libraryname)
  libraryname<-gsub(base::dirname(libraryname),"",libraryname)
  libraryname<-gsub("/","",libraryname)
  write.table(unique.Entry.No, paste0(rootdir,"/",libraryname,"_unique_Entry_no.txt"), row.names=FALSE, na="", quote=FALSE, sep="\t", col.names=FALSE)
  }

  # mergedOutFile is the name of the new outfile
  #mergedOutFile = "mergedFile.txt"
  
  # write out mergedOutFile to the current working directory (chemStationFileLoc - set above)
  #write.table(newMergedFile, mergedOutFile, sep="\t", row.names=FALSE, quote=FALSE)
  message("")
  message(paste("The Chemstation results are created and save in",rootdir, sep=" "))
  message(paste("Library(ies) found in the ID results: ",unique(library_list$Library)))
  message(paste("Done."))
}



#' msp_to_msl
#'
#' This is a function will convert the MSP library to MSL library. The default export folder will be the same as of source msp file.
#'
#' @param mspfile specify the mspfile file name with its location 
#' 
#' @return MSL library file Path
#'
#' @examples
#' msp_to_msl()
#'
#' @export
######################################### create_sub_library #########################################



msp_to_msl<-function(mspfile=NULL,RI_retrieve=T,RI_column_type=c("5ms","wax"),save_new_msp=RI_retrieve){
  library(tcltk)
  standard_msp_col=c("Name","MW","RI","Comment","CASNO","Formula","Num peaks")
  standard_MSL_col=c("NAME","MW","RI","COMMENT","CASNO","FORM","NUM PEAKS")
  if (is.null(mspfile)){mspfile=tk_choose.files(caption = "Select the .msp library file",multi = F)}
  mspcontent=readLines(mspfile)
  mspcontent1<-mspcontent[-grep("^[ ]*$",mspcontent)]
  att.nms<-unique(sapply(strsplit(mspcontent1[grep(":",mspcontent1)],":"),function(x) x[1]))
  entry.nms<-sapply(strsplit(mspcontent1[grep("Name:",mspcontent1)],"Name:"), function(x) x[2])
  cas.nms<-sapply(strsplit(mspcontent1[grep("CASNO:",mspcontent1)],"CASNO:"), function(x) x[2])
  cas.nms<-casno_reformating_string(cas.nms,hypon = F)
  comments.nms<-sapply(strsplit(mspcontent1[grep("Comment:",mspcontent1)],"Comment:"), function(x) x[2])
  rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))
  starts<-grep("Name:", mspcontent)
  stops<-c(starts[1:(length(starts)-1)]+diff(starts)-1,length(mspcontent))
  for (i in 1:length(starts)){
    tmp<-mspcontent[starts[i]:stops[i]]
    sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
  }
  rez[,"Comment"]=comments.nms
  rez<-data.frame(rez,stringsAsFactors =F )
  rez[,"SOURCE"]=mspfile
  test_col=testcolume(rez,c("Name","Comment","MW","RI","Formula","CASNO","SOURCE","Num.peaks"))
  if (length(test_col$failcol)>0){
    for (col in test_col$failcol){
      
      rez[[col]]=NA
      
    }
  }
  
  rez<-rez[,c("Name","Comment","MW","RI","Formula","CASNO","SOURCE","Num.peaks")]
  rez<-casno_reformating(rez,hypon = F)
  if (RI_retrieve){
    message(paste("Retrieving RI information... Selected colomn type:",RI_column_type))
    library(reticulate)
    library(pbapply)
    source_python(paste0(file.path(path.package(package="MassOmics")),"/src/lookup_ri.py"))
    retrieve_CASNO<-unique(rez$CASNO)
    retrieve_RI<-pblapply(retrieve_CASNO,function(x, RI_column_type){
      library(reticulate)
      library(MassOmics)
      source_python(paste0(file.path(path.package(package="MassOmics")),"/src/lookup_ri.py"))
      get_ri(x, RI_column_type)
    }, RI_column_type,cl=autoStopCluster(detectCores()))
    retrieve_df<-data.frame(stringsAsFactors = F,CASNO=retrieve_CASNO,Retrieved_RI=unlist(retrieve_RI))
    rez_merge<-merge(rez,retrieve_df,by="CASNO")
    rez_merge$RI<-rez_merge$Retrieved_RI
    org_rez<-rez
    rez<-rez_merge[,c("Name","Comment","MW","RI","Formula","CASNO","SOURCE","Num.peaks")]
  }
  if(save_new_msp){
    mspcontent_new<-character(0)
    for (i in 1:length(starts)){
      tmp<-mspcontent[starts[i]:stops[i]]
      rez$RI[rez$CASNO==cas.nms[i]]
      strsplit(tmp,":")
      find_RI<-sapply(strsplit(tmp,":"), function(x){
        if (length(x)>0)
        {if (x[1]=="RI"){T}else{F}
        }else{F}
      })
      rt_RI<-unique(rez$RI[rez$CASNO==cas.nms[i]])
      if(rt_RI!=-1){
      if (sum(find_RI)==1){
        tmp[which(find_RI)]=paste0("RI:",rt_RI)
      }else{
        #message(paste("find",i))
        tmp=c(tmp[1],paste0("RI:",rt_RI),tmp[2:length(tmp)])
      }
      }else if(sum(find_RI)==1){
        tmp<-tmp[-which(find_RI)]
      }
      
      mspcontent_new<-c(mspcontent_new,tmp)
    }
    newmspfile<-paste0(gsub(".MSP","",mspfile,ignore.case = T),"_RI_",RI_column_type,".MSP")
    writeLines(mspcontent_new,newmspfile)
    message("New MSP file have been created!")
  }
  
  
  colnames(rez)<-c("NAME","COMMENT","MW","RI","FORM","CASNO","SOURCE","NUM PEAKS")
  
      
  SumMSL=parSapply(cl=autoStopCluster(makeCluster(detectCores())),1:nrow(rez),generateMSL_par,rez,starts,stops,mspcontent)
  SumMSL=paste0(SumMSL,collapse = "\n")
  sink(paste0(gsub(".MSP","",mspfile),".MSL"))
  cat(SumMSL)
  sink()
  message("MSL file have been converted!")
  
}

generateMSL_par<- function(i,rez,starts,stops,mspcontent){
  entry=data.frame(rez[i,],stringsAsFactors = F)
  #rownames(entry)=toupper(rownames(entry))
  #entry=entry[!is.na(entry[1,]),]
  colnames(entry)=colnames(rez)
  line=1
  writeLinetmp=list()
  
  entry <- entry[,colSums(is.na(entry))<nrow(entry)]
  for (col in colnames(entry)){
    if (entry[[col]]!=""){
      writeLinetmp[[line]]=paste0(col,": ",entry[[col]])
      line=line+1
    } 
  }
  
  tmp<-mspcontent[starts[i]:stops[i]]
  tmp<-tmp[(which(tmp==paste0("Num peaks: ",as.numeric(entry[1,"NUM PEAKS"])))+1):length(tmp)]
  tmp<-tmp[tmp!=""]
  ionslist=data.frame(t(sapply(tmp, function(x) strsplit(x," ")[[1]])),stringsAsFactors = F)
  ionslist[,2]=as.numeric(ionslist[,2])
  ionslist[,2]=as.character(round(ionslist[,2]/max(ionslist[,2])*1000,digits=0))
  
  ionnumperline=1 
  for (ions in 1:nrow(ionslist)){
    if(ionnumperline==1){
      linetemp=""
      linetemp=paste0("(",paste0(sprintf("% 4s", ionslist[ions,]),collapse = " ")  ,")") 
      ionnumperline=ionnumperline+1
    }else if(ionnumperline==6){
      writeLinetmp[[line]]=linetemp
      line=line+1
      linetemp=""
      linetemp=paste0("(",paste0(sprintf("% 4s", ionslist[ions,]),collapse = " ")  ,")") 
      ionnumperline=1+1
    } else{
      linetemp=paste0(linetemp," ",paste0("(",paste0(sprintf("% 4s", ionslist[ions,]),collapse = " ")  ,")")) 
      ionnumperline=ionnumperline+1
      
    }
  }
  
  writeLinetmp[[line]]=linetemp
  line=line+1
  writeLinetmp[[line]]=""
  writeLinetmp=unlist(writeLinetmp)
  writeLinetmp=paste0(writeLinetmp,collapse = "\n")
  return(writeLinetmp)
}

MSP_generate_lib2nist<-function(unique_Entry_no="",lib2nistpath="C:\\NIST17\\MSSEARCH",convert_MSPtoMSL=T,RI_retrieve=F,RI_column_type="5ms"){
  library("tcltk")
  if(lib2nistpath==""){ lib2nistpath=  tk_choose.dir(caption = "Select lib2nist directory") }
  
  if(unique_Entry_no==""){
    unique_Entry_no=tcltk::tk_choose.files(default = "merged_library_entries_File.csv",
                                           multi = F, caption = "Select the merged library entries file in csv format")
  }
  unique_Entry_no=gsub("\\\\","/",unique_Entry_no)
  wdpath=unique(base::dirname(unique_Entry_no))
  output_filename=gsub(pattern = paste0(wdpath,"/"),replacement = "",x=tools::file_path_sans_ext(unique_Entry_no))
  
  if (wdpath=="."){wdpath=getwd()}
  
  unique_Entry_no_name=gsub(paste0(wdpath,"/"),"",unique_Entry_no)
  
  newData<-read.csv(unique_Entry_no,stringsAsFactors = F)
  library_list<-unique(newData[,c("Entry.No.Lib", "Library")])
  
  for (libraryname in unique(library_list$Library)){
    message(libraryname)
    unique.Entry.No <- library_list[library_list$Library==libraryname,"Entry.No.Lib"]
    #unique.Entry.No<-c(libraryname,unique.Entry.No)
    libraryname_org=libraryname
    libraryname<-gsub("\\\\","/",libraryname)
    libraryname<-gsub(base::dirname(libraryname),"",libraryname)
    libraryname<-gsub("/","",libraryname)
    writeLines(as.character(unique.Entry.No), paste0(wdpath,"/",libraryname,"_unique_Entry_no.txt"))
    
  #message(as.character(unique.Entry.No))
  copyresult=file.copy(paste0(wdpath,"/",libraryname,"_unique_Entry_no.txt"), paste0(lib2nistpath,"/",libraryname,"_unique_Entry_no.txt"), overwrite = TRUE)
  #unique_entries_lines<-readLines(unique_Entry_no[i])
  
  
  writeini=NULL
  writeini[1]="[Directory]"
  writeini[2]="Input=" 
  writeini[3]=paste0("Output=",gsub("/","\\\\",wdpath))
  writeini[4]=paste0("NIST=",lib2nistpath)
  writeini[5]="Type=0" 
  writeini[6]="[Output]"  
  writeini[7]="Text=1"
  writeini[8]="TextFileType=0"  
  writeini[9]="DB=0" 
  writeini[10]="CalcMW=0"  
  writeini[11]="IncludeSynonyms=0"
  writeini[12]="KeepIDs=0" 
  writeini[13]="LinkMOLfile=0" 
  writeini[14]="MzAdd=0" 
  writeini[15]="MzMpy=1" 
  writeini[16]="NeedSubset=1"  
  writeini[17]="MsmsOnly=0"  
  writeini[18]="Msms2008-Compat=0"  
  writeini[19]="[Subset]" 
  writeini[20]="Type=1"  
  writeini[21]="ForAll=1"  
  writeini[22]="FromFile=1"  
  writeini[23]=paste0("IdFile=",paste0(lib2nistpath,"/",libraryname,"_unique_Entry_no.txt"))
  
  writelineresult=writeLines(writeini,paste0(lib2nistpath,"/","Mylib.ini"))
  
  #writelineresult=writeLines(writeini,paste0(lib2nistpath,"/","Mylibb.ini"))
  
  
  
  
  
  #cmdline=paste0("cd /d ",gsub("/","\\\\",lib2nistpath),"&&",paste0("lib2nist /log3 Mylib.log Mylib.ini C:\\Database\\W11N17MAIN.L \"",gsub("/","\\\\",wdpath),"\" =Mylib"))
  cmdline=paste0("cd /d \"",gsub("/","\\\\",lib2nistpath),"\"&&","lib2nist /log3 Mylib.log Mylib.ini ",libraryname_org, " \"",gsub("/","\\\\",wdpath),"\" =Mylib_",libraryname,".MSP")
  message(cmdline)
  shell(cmdline)
  message("done")
  #readini=readLines(paste0(lib2nistpath,"/","Mylib.ini"))
  
  if(convert_MSPtoMSL){
    mspfile=paste0(gsub("/","\\\\",wdpath),"\\Mylib_",libraryname,".MSP")
    msp_to_msl(mspfile=mspfile,RI_retrieve=RI_retrieve,RI_column_type=RI_column_type)
  }
    
  }
  
  
}

MSP_generate_lib2nist_a<-function(unique_Entry_no="",lib2nistpath="C:/NIST17/MSSEARCH",convert_MSPtoMSL=T,RI_retrieve=F,RI_column_type="5ms"){
  library("tcltk")
  wdpath=base::dirname(unique_Entry_no)
  if(lib2nistpath==""){ lib2nistpath=  tk_choose.dir(caption = "Select lib2nist directory") }
  
  if(unique_Entry_no==""){unique_Entry_no=tk_choose.files(multi = F,caption = "Select lib2nist directory",filter =  matrix(c( "Text", ".txt"),1, 2, byrow = TRUE))
  wdpath=base::dirname(unique_Entry_no)
  }
  if (wdpath=="."){wdpath=getwd()}
  
  unique_Entry_no_name=gsub(paste0(wdpath,"/"),"",unique_Entry_no)
  copyresult=file.copy(unique_Entry_no, paste0(lib2nistpath,"/",unique_Entry_no_name), overwrite = TRUE)
  
  writeini=NULL
  writeini[1]="[Directory]"
  writeini[2]="Input=" 
  writeini[3]=paste0("Output=",gsub("/","\\\\",wdpath,"\\\\Mylibs",))
  writeini[4]=paste0("NIST=",lib2nistpath)
  writeini[5]="Type=0" 
  writeini[6]="[Output]"  
  writeini[7]="Text=1"
  writeini[8]="TextFileType=0"  
  writeini[9]="DB=0" 
  writeini[10]="CalcMW=0"  
  writeini[11]="IncludeSynonyms=0"
  writeini[12]="KeepIDs=0" 
  writeini[13]="LinkMOLfile=0" 
  writeini[14]="MzAdd=0" 
  writeini[15]="MzMpy=1" 
  writeini[16]="NeedSubset=1"  
  writeini[17]="MsmsOnly=0"  
  writeini[18]="Msms2008-Compat=0"  
  writeini[19]="[Subset]" 
  writeini[20]="Type=1"  
  writeini[21]="ForAll=1"  
  writeini[22]="FromFile=1"  
  writeini[23]=paste0("IdFile=",paste0(lib2nistpath,"/",unique_Entry_no_name))
  
  writelineresult=writeLines(writeini,paste0(lib2nistpath,"/","Mylib.ini"))
  
  #writelineresult=writeLines(writeini,paste0(lib2nistpath,"/","Mylibb.ini"))
  
  
  
  
  
  cmdline=paste0("cd /d ",gsub("/","\\\\",lib2nistpath),"&&",paste0("lib2nist /log3 Mylib.log Mylib.ini C:\\Database\\W11N17MAIN.L \"",gsub("/","\\\\",wdpath),"\" =Mylib"))
  
  shell(cmdline)
  
  readini=readLines(paste0(lib2nistpath,"/","Mylib.ini"))
  
  if(convert_MSPtoMSL){
    mspfile=paste0(wdpath,"/Mylib.MSP")
    msp_to_msl(mspfile=mspfile,RI_retrieve=RI_retrieve,RI_column_type=RI_column_type)
    
  }
  
}
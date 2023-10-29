#' GCMS_integration
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


GCMS_integration<-function(output = "GC-MS Result",
                           Ion.bin= 0.5,
                           workdir=getwd(),
                           intensity_type=c("Peak Area","Peak Height"),
                           GGReport=c("Fast","Slow")
) {
  require(flux)
  require(parallel)
  require(xcms)
  require(lattice)
  require(plyr)
  require(tcltk)
  require(mzR)
  require(Rcpp)
  require(data.table)
  setwd(workdir)
  final.check.data <- read.csv(paste0(tcltk::tk_choose.files(caption = "Load a Summary Report.csv",  multi = FALSE),collapse = " "),stringsAsFactors=FALSE)
  
  final.check.data <- final.check.data[final.check.data$Total.ID!=0,]
  check.Ref.ion <-table(is.na(final.check.data$Ref.ion))["TRUE"]
  check.Ref.ion[is.na(check.Ref.ion)]<-0
  if(check.Ref.ion>=1){tkmessageBox(message = "ERROR: Missing reference ions in the Identification Report!", icon = "warning", type = "ok"); stop("Missing reference ion")}
  
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
  
  library_file_org=ion.lib[,c("Name","Ref.ion", "RT.median", "RT.shfit.Lower", "RT.shfit.upper","Peak.width")]
  library_file_org<-library_file_org[!is.na(library_file_org$RT.median),]
  final.df<-NULL
  remove(final.df)
  Graphic.df<-data.frame()
  filenames <- dir(conditions)
  filenames <- filenames[grep(c(".cdf$|.mzXML$"),filenames, ignore.case=TRUE)]
  thread <- round(nrow(library_file_org)/200,digits = 0)
  if (thread>=6) {thread=6}
  if (thread<1) {thread=1}
  
  for (q in 1:length(filenames)) {
    findScanTime=NULL
    res<-try(findScanTime <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "/"),))
    if(class(res) == "try-error"){
     #scanextract=stringr::str_extract_all(res[1],"scan.....")[[1]]
     #scanextract=gsub("scan","",scanextract)
     #scanextract=gsub(" ","",scanextract)
     #scanextract=as.numeric(scanextract)
     
     filepath <- xcmsSource(paste(conditions, filenames[q], sep = "/"))
     rawdata <- loadRaw(filepath, includeMSn = F)
     mzrange_file<-range(rawdata$mz,na.rm = T)
     scanrange_file<-max(which((min(which(is.na(rawdata$mz)))>=rawdata$scanindex)))-1
     res<-try(findScanTime <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "/"),scanrange = c(1,scanrange_file)))
     #res<-try(findScanTime <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "/"),scanrange = c(1,length(rawdata$scanindex))))
     #res<-try(findScanTime <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "/"),scanrange = 1:10))
    } 
    
    
    if(class(res) != "try-error"){
    ScanPerSec <- round(length(findScanTime@ scantime)/(findScanTime@ scantime[length(findScanTime@ scantime)]-findScanTime@ scantime[1]),2)
    scanrangeL <- min(findScanTime@scantime)
    scanrangeU <- max(findScanTime@scantime)
    message(paste("Scan rate:", ScanPerSec,"scan/second"))
    
    library_file <- library_file_org
    library_file$RT.shfit.Lower[is.na(library_file_org$RT.shfit.Lower)]<-0 ## testing
    library_file$RT.shfit.Lower<-library_file_org$RT.shfit.Lower + library_file_org$Peak.width/ScanPerSec
    
    library_file$RT.shfit.upper[is.na(library_file_org$RT.shfit.upper)]<-0 ## testing
    library_file$RT.shfit.upper<-library_file_org$RT.shfit.upper +library_file_org$Peak.width/ScanPerSec
    library_file$Ref.ion=as.numeric(as.character(library_file$Ref.ion))
    library_file<-library_file[!is.na(library_file$RT.median),]
    
    file1 <- filenames[q]
    name.file <- file1
    name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
    surefinal <- data.frame(Name=ion.lib[,"Name"],Base.Peak=0)
    raw_data <- findScanTime
    
    if (GGReport=="Slow"){
    
      
        integrationresult=parLapply(cl=autoStopCluster(makeCluster(thread)),
                                             1:nrow(library_file),
                                             Par_peakintegrate_slow,
                                             library_file,
                                             findScanTime,
                                             raw_data,
                                             intensity_type,
                                             Ion.bin,
                                             name.file)
        
        
        Graphic.dflist<-list()
        for (h in 1:nrow(library_file)){
          surefinal$Base.Peak[h]<-integrationresult[[h]]["Peakvalue"]
          Graphic.dflist[[h]]<-as.data.frame(integrationresult[[h]]["Graphic.df"])
          #message(h)
        }
        surefinal$Base.Peak=unlist(surefinal$Base.Peak)
        Graphic.df<-rbind(Graphic.df,data.table::rbindlist(Graphic.dflist))
        #Graphic.df<-dplyr::bind_rows(Graphic.dflist)
  
    }
    
    
    if (GGReport=="Fast"){
    
    
    surefinal$Base.Peak=unlist(parLapply(cl=autoStopCluster(makeCluster(thread)),
                                         1:nrow(library_file),
                                         Par_peakintegrate,
                                         library_file,
                                         findScanTime,
                                         raw_data,
                                         intensity_type,
                                         Ion.bin))
    }
    
    
    
    surefinal <- surefinal[!duplicated(surefinal$Name),]
    names(surefinal)[2] <- name.file
    
    if (!exists("final.df")){
      final.df <- as.data.frame(surefinal)
      colnames(final.df)=c("Name",filenames[q])
      confirmation <- paste("File", q, "(", name.file,")", "done!", sep = " ")
      message(confirmation)
    } else {
      colnames=colnames(final.df)
      final.df <- cbind(final.df,surefinal[-1])
      
      colnames(final.df)=c(colnames,filenames[q])
      confirmation <- paste("File", q, "(", name.file,
                            ")", "done!", sep = " ")
      message(confirmation)
    }
    
    write.csv(final.df,file = "Integration report tempfile.csv",row.names = F)
    
    }
  }
  
  message("Estimating errors")
  
  
  
  ### setup warnign message
  
  if(class(ion.lib$CAS)!="NULL"){
    finalX.df <- cbind(CAS=ion.lib$CAS, Ref.Ion=ion.lib$Ref.ion, Total.ID=ion.lib$Total.ID, Library.match=ion.lib$Library.match, Ret.Time=ion.lib$RT.median, RT.shfit.Lower=ion.lib$RT.shfit.Lower, Warnings=NA, final.df)
    
    finalX.df <- finalX.df[order(finalX.df$Ret.Time , decreasing = FALSE),]
    parameterLlibraryMatch <- 40
    parameterRTShif <- 60
    parameterSimilarRT <- 0.01
    parameterN <-0.05 # Percentage
    
    # warning matrix
    warning.df <- data.frame(matrix(vector(), length(finalX.df[,1]), 6, dimnames=list(c(), c("LowID", "Library.match.d", "RT.shift.d", "Repeated.aboudnace", "Repeated.RT", "Close.RT"))), stringsAsFactors=F)
    
    # filter low match filter
    low.match<- which(finalX.df$Library.match <  parameterLlibraryMatch)
    warning.df[row.names(warning.df) %in% low.match,"Library.match.d"]<-"Library.Match<40%/"
    
    
    # Filter low RT shift
    high.RT.shift <- which(finalX.df$RT.shfit.Lower >  parameterRTShif )
    warning.df[row.names(warning.df) %in%  high.RT.shift,"RT.shift.d"]<-"RT.shift>60s/"
    
    # Filter low Number of ID
    
    low.nID <- which(finalX.df$Total.ID  <  parameterN*length(filenames))
    
    warning.df[row.names(warning.df) %in%   low.nID,"LowID"]<-"Total.ID<5%/"
    low.nID <- which(finalX.df$Total.ID  ==1 )
    warning.df[row.names(warning.df) %in%   low.nID,"LowID"]<-"Total.ID=1/"
    
    total.n<-length(filenames)
    ID_name<-paste("Total.ID(n=",total.n,")",sep="")
    names(finalX.df)[3]<-ID_name
    
    
    
    # Remove the repeat aboundace
    last.postion <- 8 + length(filenames)
    dulpicated_aboundace<-unique(finalX.df[duplicated(finalX.df[, last.postion]), last.postion])
    Repeate.aboundance<-NULL
    if (length(dulpicated_aboundace)!=0){
      for(n in 1: length(dulpicated_aboundace)){
        Repeat.All<- finalX.df[finalX.df[,last.postion]==dulpicated_aboundace[n],]
        Repeate.aboundance <-rbind(Repeate.aboundance, Repeat.All[-which.max(Repeat.All$Library.match),])
      }
      warning.df[row.names(warning.df) %in% row.names(Repeate.aboundance),"Repeated.aboudnace"]<-"Repeated.Level/"
    }
    
    # Remove the repated RT
    dulpicated_RT<-unique(finalX.df$Ret.Time[duplicated(finalX.df$Ret.Time)])
    
    Repeate.RT<-NULL
    if (length(dulpicated_RT)!=0){
      for(n in 1: length(dulpicated_RT)){
        Repeat.All<- finalX.df[finalX.df$Ret.Time==dulpicated_RT[n],]
        Repeate.RT <-rbind(Repeate.RT, Repeat.All[-which.max(Repeat.All$Library.match),])
      }
      warning.df[row.names(warning.df) %in% row.names(Repeate.RT),"Repeated.RT"]<-"Repeated.RT/"
    }
    
    # Identfiy similar RT
    close.RT<-NULL
    for(i in 1:(length(finalX.df[,1])-1)){
      n=i+1
      diff <-finalX.df$Ret.Time[n] - finalX.df$Ret.Time[i]
      if(diff <=  parameterSimilarRT & diff!=0 ){
        close.df <- finalX.df[i:n,]
        close.RT <- rbind(close.RT,  close.df[-which.max(close.df$Library.match),])
        warning.df[row.names(warning.df) %in% row.names(close.RT),"Close.RT"]<-"Similar.RT/"
      }
    }
    
    warning.list <- paste(warning.df[,1], warning.df[,2], warning.df[,3], warning.df[,4], warning.df[,5], warning.df[,6], sep = ' ')
    warning.list <- gsub("NA","", warning.list)
    warningZ.list <- gsub("  ","", warning.list)
    warningY.list <- gsub(" ","", warning.list)
    warningX.list <- gsub(".$", "",warningY.list)
    finalX.df$Warnings <-warningX.list
    finalF.df <- finalX.df
    finalZ.df<-finalX.df
    finalZ.df$Warnings <-warningZ.list
    
    erroreport<- finalZ.df
    finalF_errorplot <-erroreport[order(erroreport$Name , decreasing = FALSE),]
    greystrip<-gsub(" ", "grey95",finalF_errorplot$Warnings)
    greystrip[greystrip!="grey95"]<-"firebrick1"
    redwarning<<-greystrip
    
  } else {finalF.df=cbind( Ref.Ion=ion.lib$Ref.ion, Ret.Time=ion.lib$RT.mean,RT.shift=ion.lib$RT.shift, final.df)
  finalF.df <- finalF.df[order(finalF.df$Ret.Time , decreasing = FALSE),]
  }
  
  ## Save Files
  message("Save file")
  if (GGReport=="Slow"){
    colnames(Graphic.df)=stringr::str_replace(colnames(Graphic.df),"Graphic.df.","")
    write.csv(Graphic.df,file="GC-Peak Diagnosis Report.csv")
    }
  ifelse(intensity_type=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
  sheet <- paste(output,"(", ValueType ,")", sep="")
  store <- paste(workdir, "\\", sheet, ".csv", sep = "")
  FileName <- tclvalue(tkgetSaveFile(initialfile=store))
  write.csv(finalF.df, file = FileName, row.names = FALSE)
  
  message(paste("Data process is completed and saved as: ",sheet,".csv", sep=""))
}

#' mz range detection
#'
#' This is a function that help user to determine the sharing mz range acrossing data files. 
#'
#' @param QCfiles specify the data file(s) to find the mz range
#' 
#' @return a vector contains mz lower/upper limit 
#'
#' @examples
#' MZrange_autodetect()
#'
#' @export


MZrange_autodetect<-function(){
  mz_L=0
  mz_U=1000
  require(xcms)
  require(lattice)
  require(plyr)
  require(tcltk)
  require(mzR)
  require(Rcpp)
  QCfiles = tk_choose.files(multi = T, caption = "Select a data file to detect the mz range")
  matchtype=c("mzxml$","cdf$")
  #datafiles=gsub(paste0(dirname(datafiles)[1],"/"),"",datafiles)
  
  QCfiles=QCfiles[grep(pattern =paste(matchtype,collapse="|"),x=QCfiles,ignore.case = T)]
  openCDFpar<-function(h,QCfiles){
    require(xcms)
    try(findScanTime <- xcmsRaw(filename = QCfiles[h])) 
    findScanTime@mzrange
    
  }
  
  #mzranges<-parLapply(cl=autoStopCluster(makeCluster(detectCores())),1:length(QCfiles),openCDFpar,QCfiles)
  for (QCfile in QCfiles){
    findScanTime <- xcmsRaw(filename = QCfile)
    if (mz_L<=min(findScanTime@mzrange)) mz_L=min(findScanTime@mzrange)
    if (mz_U>=max(findScanTime@mzrange)) mz_U=max(findScanTime@mzrange)  
  }
  
  
  
  #mz_L <- tclVar(mz_L)
  #mz_U <- tclVar(mz_U)
  #tkconfigure(MZrange_L,variable=mz_L)
  #tkconfigure(MZrange_U,variable=mz_U)
  #MZrange_L <- tkentry(TMS, textvariable=mz_L ,width=6)
  #MZrange_U <- tkentry(TMS, textvariable=mz_U ,width=6)
  #tkgrid(tklabel(frameOverall,text="    "))
  #try(tclSetValue("mz_L", mz_L))
  #try(tclSetValue("mz_U", mz_U))
  return(c(mz_L,mz_U))
}


#' AMDIS report file parsing and summary
#'
#' This is a function to process the AMDIS_report.txt file. specify a in-house library or a subset of NIST library to find the most intense quantification ion within the GC-MS experiment m/z range. The script will invoke selection windows to let the user specify the files location.  
#'
#' @param workdir locate the working Dir
#' @param MS.L specify the library file
#' @param MsLibrary specify the library file origin, could be "NIST" or "InHouse"
#' @param amdis.report specify the identification report file
#' @param Ret.Time.Filter set the retetion time filter window around the expected retention time (in +/-min).
#' @param RT.shift.limt set a threshole in seconds to define a molecule have 
#' @param mz_L lower mz limit for quantification ion selection
#' @param mz_U upper mz limit for quantification ion selection
#' @param generate_rt_shift_graph Set TRUE to generate the retetion time shift plots for each molecules that have retention time discrepencies greater than the \code{"RT.shift.limt"}
#' @param RTcorrection Set TRUE to enable a retention time correction before multiple peaks resolving
#' @return None
#'
#' @examples
#' amdis_id_Summary()
#'
#' @export
amdis_id_Summary<-function(workdir= NULL,
                           MS.L= NULL,
                           amdis.report = NULL,
                           File.name = "Summary Report.csv",
                           MsLibrary=c("NIST", "InHouse"),
                           Ret.Time.Filter=2.5,
                           RT.shift.limt = 60,
                           mz_L=38,
                           mz_U=550,
                           generate_rt_shift_graph=F,
                           generate_rt_shift_graphs=F,
                           RTcorrection=F
){
  library("tcltk")
  library("tcltk2")
  # require("Metab")
  library("pander")
  library("parallel")
  #library("OneR")
  library("dplyr")
  library("ggplot2")
  library("ggpubr")
  #workdir <- workdir
  if (is.null(workdir)){workdir= tk_choose.dir(caption = "Select working directory")}

  if (is.null(amdis.report)){amdis.report = tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                            multi = FALSE, caption = "Select the AMDIS report in .TXT",
                                                            filter=matrix(c("Text", ".txt", "All files", "*"),2, 2, byrow = TRUE))}
  
    if (is.null(MS.L)){MS.L= tk_choose.files(caption="Select MS library (e.g. SVB_total or NIST) in .msl",
                                           filter=matrix(c("MSL", ".msl", "All files", "*"),2, 2, byrow = TRUE))}
  setwd(workdir)
  amdis.report = data.table::fread(amdis.report,sep = "\t", stringsAsFactors=FALSE)
  message("amdis_id_Summary parameters")
  message(paste("RT filter: ", Ret.Time.Filter, " min", sep=""))
  message(paste("Detect mutiple peaks when RT range is greater than: ", RT.shift.limt, " s", sep=""))
  message(paste("Quant mass will be selected within: ", mz_L, "-", mz_U,sep=""))
  message(paste("Perform RT correction before multi-peak resolving: ",RTcorrection ,sep=""))
  # Extract Reference ion, CAS # from MS library
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
  
  # ExtractIon for NIST library
  if(MsLibrary=="NIST"){
    message("Running in NIST mode, will generate Ref ions for each metabolites within mz detection range")
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    reference_ion<-NULL
    
    #test<-buildLib(AmdisLib=MS.L, folder=workdir, save = F, output = "ion_lib.csv", verbose = F,mz_L = mz_L,mz_U = mz_U)
    
    
    #reference_ion=parLapply(cl=autoStopCluster(makeCluster(detectCores()-1)), 1:length(starts),parse_msl_par,lib.txt,starts,stops,mz_L,mz_U)
    
    #reference_ion=unlist(reference_ion)
    
    #sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
    reference_ion=0
    
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez$NAME<-entry.nms
    libr<-cbind(rez,reference_ion)## returns list of peaks
    libr$NAME <- gsub("?", "", libr$NAME, fixed = TRUE)
    libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
    libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
    libr$NAME <- gsub("[*:%^!;$]", "", libr$NAME, perl=T)
    
    # Generate Summary report
    AMDIS.report <-amdis.report
    AMDIS.report$Name <- gsub("?", "", AMDIS.report$Name, fixed = TRUE)
    AMDIS.report$Name <- gsub("^ ", "", AMDIS.report$Name, perl=T)
    AMDIS.report$Name <- gsub("[*:%^!*?&;$]", "", AMDIS.report$Name, perl=T)
    AMDIS.report$Width <- gsub(">", "",  AMDIS.report$Width, perl=T)
    AMDIS.report$Width <- as.numeric(gsub("scans", "",  AMDIS.report$Width, perl=T))
    #AMDIS.reportcas<-AMDIS.report[AMDIS.report$CAS %in% libr$CASNO,]
    AMDIS.report<-AMDIS.report[AMDIS.report$Name %in% gsub("^ ","",libr$NAME),]
    #AMDIS.report1<-AMDIS.report[AMDIS.report$Name %in% libr$NAME,]
    #
    
    if (RTcorrection){
      library(ggpubr)
      RT.correction.grid<-RT_correction_xcms()
      RT.correction.result<-RT_correction(df=AMDIS.report)
      AMDIS.report<-RT.correction.result[[1]]
      ggpubr::ggarrange(RT.correction.result[[2]],RT.correction.result[[3]],
                        labels = c("Before","After"),common.legend = T,
                        ncol = 2, nrow =1,vjust=1.5) %>%  ggexport(filename = paste(getwd(),"/retention time correction.png",sep=""),width = 1800, height = 900,verbose = NULL,res = 100)
      
    }
    Metabolite.list <-split(AMDIS.report$RT, AMDIS.report$Name)
    RT.stats<-t(sapply(Metabolite.list,function(x) c(RT.median=round(median(x,na.rm=TRUE),3),
                                                     RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5))*60,2),
                                                     RT.shift=round(((quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)-median(x,na.rm=TRUE))*60,2)
    )))
    colnames(RT.stats)<-c("RT.median","RT.shfit.Lower", "RT.shfit.upper")
    RT.stats=RT.stats[is.na(RT.stats[,"RT.median"])!=T,]
    Peak.width <- sapply(split(AMDIS.report$Width, AMDIS.report$Name), function(x) median(x))
    
    
    ID.stats<-t(sapply(split(AMDIS.report$Weighted, AMDIS.report$Name),function(x) c(Library.match=round(mean(x,na.rm=T)),
                                                                                     Total.ID=length(x[!is.na(x)]))))
    #message("Running in NIST mode, will generate Ref ions for each metabolites within mz detection range")
    AMDIS.report.RT.stats <- data.frame(cbind(Name=0,Ref.ion=0, ID.stats, RT.stats, Peak.width, CAS=0))
    #message(paste(AMDIS.report.RT.stats$Ref.ion[1:40],sep="\n"))
    #AMDIS.report.RT.stats$Ref.ion <-  libr$reference_ion[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    #matchID=match(rownames(RT.stats),libr$NAME)
    matchID=match(strtrim(rownames(RT.stats),170),strtrim(libr$NAME,170))
    #AMDIS.report.RT.stats=AMDIS.report.RT.stats[which(!is.na(matchID)),]
    #matchID70=match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))
    #intersect(matchID,matchID70)
    AMDIS.report.RT.stats$Ref.ion <-  unlist(parallel::parLapply(cl=autoStopCluster(makeCluster(detectCores()-1)), 
                                                                 matchID,
                                                                 parse_msl_par,lib.txt,starts,stops,mz_L,mz_U))
    
    #for (ids in matchID){
    # Ref.iontest= parse_msl_par(ids,lib.txt,starts,stops,mz_L,mz_U)
      
    #}
    #match("2-Cyclopropylethynylcyclopropane",libr$NAME)
    #message(AMDIS.report.RT.stats$Ref.ion)
    #which(match(rownames(RT.stats),libr$NAME)==NA)
    AMDIS.report.RT.stats$CAS <-  libr$CAS[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    AMDIS.report.RT.stats$Name <- rownames(RT.stats)
    
    #total.n<-length(unique(AMDIS.report$FileName))
    #ID_name<-paste("Total.ID(n=",total.n,")",sep="" )
    #names(AMDIS.report.RT.stats)[3]<-ID_name
    
    
    ## Remove point outside the range
    RT.list<-lapply(Metabolite.list,function(x) c(RT.shift=(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5),
                                                  RT.shift=(quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)))
    Remove.outlier<-function(x,y){
      x[x >= y[1]]
      x[x <= y[2]]
    }
    
    Remove.metabolite <-mapply(Remove.outlier, Metabolite.list, RT.list )
    
    # Plot all retention time of all metabolites
    
    n.total <- length(Metabolite.list)
    png(paste(getwd(),"/retention time of all metabolites.png",sep=""),width = 30,height = 30, bg = "white",units = "in",res = 150)
    par(mar=c(5,10,5,5))
    boxplot(Metabolite.list, xlab= "Retention time (min)", cex.axis=0.2,  horizontal = TRUE, las=2, xaxt="n", main = paste("Retention time distribution of all metabolites ","(n=",n.total,")",sep=""))
    axis(1,cex.axis=1,las=1)
    dev.off()
    
    
    ## solve multimodel issue
    
    
    
    AMDIS.report.split<-NULL
    deletenamelist<-NULL
    save.folder<-paste("RTshift_warnings_",RT.shift.limt,"s_","shift",sep="")
    if (dir.exists(save.folder)!=T) try(dir.create (save.folder))
    
    for (component in 1:length(Remove.metabolite)){
      if(`&`((RT.stats[component,2] + RT.stats[component,3]) > RT.shift.limt, T)){
        #names(Remove.metabolite)[component]
        RT_libsimilarity=NULL
        #plot figure
        names.var<- windows_filename(names(Remove.metabolite)[component])
        RT_libsimilarity<-AMDIS.report[AMDIS.report$Name==names(Remove.metabolite)[component],]
        #nbins=6
        #res <- try(RT_libsimilarity$peakgroup<-OneR::bin(RT_libsimilarity$RT,nbins = nbins,labels = c(paste("Peak",1:nbins)), method = "clusters",na.omit = T),silent = TRUE)
        #if (class(res) != "try-error"){
        #while(`&`(sum(grep("Peak",RT_libsimilarity$peakgroup))==0,nbins>=1)){
        #  nbins=nbins-1
        #  res <- try(RT_libsimilarity$peakgroup<-OneR::bin(RT_libsimilarity$RT,nbins = nbins,labels = c(paste("Peak",1:nbins)), method = "clusters",na.omit = T),silent = TRUE)
        #}
        #RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT
        #if (sum(grep("Peak",RT_libsimilarity$peakgroup))!=0){
        #  RT_libsimilarity$peakgroup_RT_median=NULL
        #  RT_libsimilarity.stats<-RT_libsimilarity %>% group_by(peakgroup) %>% summarise(peakgroup_RT_median = median(RT))
         # RT_libsimilarity=merge(RT_libsimilarity,RT_libsimilarity.stats)
          
         # }
        #RT_libsimilarity<-AMDIS.report[AMDIS.report$Name==names(Remove.metabolite)[component],]
        #RT_libsimilarity$peakgroupRT=median(RT_libsimilarity$peakgroup)
        #RT_libsimilarity$series=1:nrow(RT_libsimilarity)
        #try(RT_libsimilarity<-RT_libsimilarity[order(RT_libsimilarity$FileName),])

        #estimate number of model
        m.var <- data.frame(Remove.metabolite[component])[,1]
        metabolite.name<-names(Remove.metabolite[component])
        res=NULL
        res <- try(density.all <- density(m.var, bw = "SJ"))
        if(class(res)=="try-error"){density.all <- density.default(m.var)}
        density.var <- density.all$y
        density.x <-  density.all$x[density.var>0.01]
        density.var <-  density.var[density.var>0.01]
        
        modes <- NULL
        centre.p <- NULL
        for ( i in 2:(length(density.var)-1) ){
          if ( (density.var[i] > density.var[i-1]) & (density.var[i] > density.var[i+1]) ) {
            modes <- c(modes,i)
          }
          if ( (density.var[i] < density.var[i-1]) & (density.var[i] < density.var[i+1]) ) {
            centre.p <- c(centre.p,i)
          }
        }
        
        new.median<- density.x[modes]
        new.boundary<- density.x [centre.p]
        npeaks <-length(new.median)
        
        if(npeaks>1){
          
          vComposti <- vector("list", npeaks)
          m.var.sub <- m.var
          
          for (i in 1:(npeaks-1)){
            vComposti[i] <-list(m.var.sub [m.var.sub <new.boundary[i]])
            m.var.sub[m.var.sub<new.boundary[i]]<-NA
          }
          
          vComposti[npeaks] <- list(m.var.sub)
          
          RT.add.stats <- t(sapply(vComposti,function(x) c(Total.ID=length(x[!is.na(x)]),
                                                           RT.median=round(median(x,na.rm=TRUE),3),
                                                           RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*3))*60,2),
                                                           RT.shift=round(((quantile(x, prob = 0.75, na.rm=TRUE)+IQR(x,na.rm=TRUE)*3)-median(x,na.rm=TRUE))*60,2))))
          
          colnames(RT.add.stats)<-c("Total.ID","RT.median","RT.shfit.Lower", "RT.shfit.upper")
          AMDIS.report.RT.stats.add <- AMDIS.report.RT.stats[rownames(AMDIS.report.RT.stats)==metabolite.name,]
          AMDIS.report.RT.stats.add <- cbind(AMDIS.report.RT.stats.add[,1:3] ,RT.add.stats, AMDIS.report.RT.stats.add[,8:9], row.names = NULL)
          
          # re-name
          for (i in 1:npeaks){
            AMDIS.report.RT.stats.add$Name[i] <- paste(metabolite.name, " (split peak ", i, ")", sep="")
          }
          AMDIS.report.split <- rbind(AMDIS.report.split, AMDIS.report.RT.stats.add)
          deletenamelist<-c(deletenamelist,metabolite.name)

        }          
        if (generate_rt_shift_graph){          
          #png(filename=paste(save.folder,"/",names.var,".png",sep=""),width = 480,height = 720,)
          #par(mfrow = c(2, 1))
          #png(filename=paste(save.folder,"/",names.var,"lib_score.png",sep=""),width = 900,height = 1200)
          #par(mfrow = c(2, 3))
            #RT_libsimilarity<-AMDIS.report[AMDIS.report$Name==names(Remove.metabolite)[i],]
            #RT_libsimilarity$peakgroup="Outliner"
            #RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT

          names.var<- windows_filename(names(Remove.metabolite)[component])          

          
          #boxplotlibsim=as_ggplot(boxplot(Weighted~peakgroup,data=RT_libsimilarity, main=names(Remove.metabolite)[i], xlab="Retention time (min)", ylab="Library Score"))

          #rtdensityplot=plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",xlab= "Retention time (min)", main = names(Remove.metabolite)[i], lwd=1.5)
          #dev.off()
          #plot figure
          
            png(filename=paste(save.folder,"/",names.var,".png",sep=""))
            #png(filename=paste(save.folder,"/","aaa",".png",sep=""))
            
            plot(RT_libsimilarity$RT,RT_libsimilarity$Net,  ylab = "Net score",
                 xlab= "Retention time (min)", main = names(Remove.metabolite)[component], lwd=1.5)
            dev.off()
          
          
          #png(filename=paste(save.folder,"/",names.var,"_lib_score.png",sep=""),width = 400, height = 400)

          
          }
        if (generate_rt_shift_graphs){          
            #png(filename=paste(save.folder,"/",names.var,".png",sep=""),width = 480,height = 720,)
            #par(mfrow = c(2, 1))
            #png(filename=paste(save.folder,"/",names.var,"lib_score.png",sep=""),width = 900,height = 1200)
            #par(mfrow = c(2, 3))
            #RT_libsimilarity<-AMDIS.report[AMDIS.report$Name==names(Remove.metabolite)[i],]
            RT_libsimilarity$peakgroup="Outliner"
            RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT
            for (peak in 1: nrow(AMDIS.report.RT.stats.add)){
              RT_libsimilarity$peakgroup=ifelse(`&`(RT_libsimilarity$RT>=(AMDIS.report.RT.stats.add$RT.median[peak]-AMDIS.report.RT.stats.add$RT.shfit.Lower[peak]/60),
                                                    RT_libsimilarity$RT<=(AMDIS.report.RT.stats.add$RT.median[peak]+AMDIS.report.RT.stats.add$RT.shfit.upper[peak]/60)),paste("Peak",peak),RT_libsimilarity$peakgroup)
              RT_libsimilarity$peakgroup_RT_median=ifelse(`&`(RT_libsimilarity$RT>=(AMDIS.report.RT.stats.add$RT.median[peak]-AMDIS.report.RT.stats.add$RT.shfit.Lower[peak]/60),
                                                              RT_libsimilarity$RT<=(AMDIS.report.RT.stats.add$RT.median[peak]+AMDIS.report.RT.stats.add$RT.shfit.upper[peak]/60)),AMDIS.report.RT.stats.add$RT.median[peak],"")
            }
            names.var<- windows_filename(names(Remove.metabolite)[component])          
            boxplotlibsim <- ggplot(RT_libsimilarity, aes(x=peakgroup, y=Net, fill=peakgroup)) + 
              geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2)) + theme(
                panel.background = element_rect(fill = "white"),
                plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
                plot.background = element_rect(
                  fill = "grey90",
                  colour = "black"
                ))
            
            #boxplotlibsim=as_ggplot(boxplot(Weighted~peakgroup,data=RT_libsimilarity, main=names(Remove.metabolite)[i], xlab="Retention time (min)", ylab="Library Score"))
            
            #rtdensityplot=plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",xlab= "Retention time (min)", main = names(Remove.metabolite)[i], lwd=1.5)
            #dev.off()
            
            
            #png(filename=paste(save.folder,"/","aaa",".png",sep=""))
            
            Net.plot=ggplot(RT_libsimilarity, aes(x=RT, y=Net, color=peakgroup,size=2)) + geom_point() + theme(
              panel.background = element_rect(fill = "white"),
              plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
              plot.background = element_rect(
                fill = "grey90",
                colour = "black"
              ))
            
            ggpubr::ggarrange(boxplotlibsim,Net.plot,
                              labels = c("boxplot","Net score"),
                              ncol = 2, nrow =1,vjust=1.5) %>%  ggexport(filename = paste(save.folder,"/",names.var,"_lib_score.png",sep=""),width = 900, height = 400,verbose = NULL)
            #dev.off()
          }
      }
    }
    
    if(is.null(AMDIS.report.split)!=TRUE){
      
      deletenamelist=unique(deletenamelist)
      
      AMDIS.report.RT.statscombine<-AMDIS.report.RT.stats[!(AMDIS.report.RT.stats$Name %in% (deletenamelist)),]
      
      AMDIS.report.split=unique(AMDIS.report.split)
      
      AMDIS.report.RT.stats <- rbind(AMDIS.report.RT.statscombine, AMDIS.report.split)
    }
    
    AMDIS.report.RT.stats <- AMDIS.report.RT.stats[order(AMDIS.report.RT.stats$RT.median, decreasing = FALSE),]
    File.nameNIST <- "Summary report_NIST.csv"
    write.csv(AMDIS.report.RT.stats, file = File.nameNIST,row.names = FALSE)
    write.csv(libr, file = "Library summary.csv",row.names = FALSE)
    
    
    
    
    message(paste("Summary report_NIST.csv was created and save in",workdir, sep=" "))
  } else if (MsLibrary=="InHouse") {
    
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez$NAME<-entry.nms
    libr<-rez## returns list of peaks
    
    # Generate summary report
    AMDIS.report<-amdis.report
    AMDIS.report<-AMDIS.report[AMDIS.report$RT-AMDIS.report$`Expec. RT` < Ret.Time.Filter,]
    AMDIS.report<-AMDIS.report[AMDIS.report$RT-AMDIS.report$`Expec. RT`> -Ret.Time.Filter,]
    AMDIS.report$Name <- gsub("?", "", AMDIS.report$Name, fixed = TRUE)
    AMDIS.report$Name <- gsub("^ ", "", AMDIS.report$Name, perl=T)
    AMDIS.report$Name <- gsub("[*:%^!*?&;$]", "", AMDIS.report$Name, perl=T)
    AMDIS.report$Width <- gsub(">", "",  AMDIS.report$Width, perl=T)
    AMDIS.report$Width <- as.numeric(gsub("scans", "",  AMDIS.report$Width, perl=T))
    
    
    
    Metabolite.list <-split(AMDIS.report$RT, AMDIS.report$Name)
    
    RT.stats<-t(sapply(Metabolite.list,function(x) c(RT.median=round(median(x,na.rm=TRUE),3),
                                                     RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25)-IQR(x,na.rm=TRUE)*1.5))*60,2),
                                                     RT.shift=round(((quantile(x, prob = 0.75)+IQR(x,na.rm=TRUE)*1.5)-median(x,na.rm=TRUE))*60,2))))
    colnames(RT.stats)<-c("RT.median","RT.shfit.Lower", "RT.shfit.upper")
    
    Peak.width <- sapply(split(AMDIS.report$Width, AMDIS.report$Name), function(x) median(x))
    
    
    RT.stats[is.na(RT.stats)]<-0
    ID.stats<-t(sapply(split(AMDIS.report$Weighted, AMDIS.report$Name),function(x) c(Library.match=round(mean(x,na.rm=T)),
                                                                                     Total.ID=length(x[!is.na(x)]))))
    AMDIS.report.RT.stats<-NULL
    AMDIS.report.RT.stats <- as.data.frame(cbind(Name=0,Ref.ion=0,ID.stats,RT.stats,Expec.RT=0,Diff.RT=0, Peak.width, CAS=0))
    AMDIS.report.RT.stats$Ref.ion <- libr$RSN[match(rownames(RT.stats),libr$NAME)]
    AMDIS.report.RT.stats$Expec.RT <- round(as.numeric(libr$RT[match(rownames(RT.stats),libr$NAME)]),2)
    AMDIS.report.RT.stats$CAS <- libr$CASNO[match(rownames(RT.stats),libr$NAME)]
    AMDIS.report.RT.stats$Diff.RT <- AMDIS.report.RT.stats$Expec.RT- as.numeric(AMDIS.report.RT.stats$RT.median)
    AMDIS.report.RT.stats$Name <- rownames(RT.stats)
    
    
    
    
    
    ## Remove point outside the range
    RT.list<-lapply(Metabolite.list,function(x) c(RT.shift=(quantile(x, prob = 0.25)-IQR(x,na.rm=TRUE)*1.5),
                                                  RT.shift=(quantile(x, prob = 0.75)+IQR(x,na.rm=TRUE)*1.5)))
    Remove.outlier<-function(x,y){
      x[x >= y[1]]
      x[x <= y[2]]
    }
    
    Remove.metabolite <-mapply(Remove.outlier, Metabolite.list, RT.list )
    
    # Plot all retention time of all metabolites
    n.total <- length(Metabolite.list)
    par(mar=c(5,10,5,5))
    boxplot(Metabolite.list, xlab= "Retention time (min)", cex.axis=0.5,  horizontal = TRUE, las=2, xaxt="n", main = paste("Retention time distribution of all metabolites ","(n=",n.total,")",sep=""))
    axis(1,cex.axis=1,las=1)
    
    
    ## solve multimodel issue
    AMDIS.report.split<-NULL
    
    save.folder<-paste("RTshift_warnings_",RT.shift.limt,"s_","shift",sep="")
    if (dir.exists(save.folder)!=T) try(dir.create (save.folder))
    for (i in 1:length(Remove.metabolite)){
      if(RT.stats[i,2] + RT.stats[i,3] > RT.shift.limt ){
        
        # plot figure
        names.var<- gsub("/", "", names(Remove.metabolite)[i], perl=T)
        png(filename=paste(save.folder,"/",names.var,".png",sep=""))
        plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",
             xlab= "Retention time (min)", main = names(Remove.metabolite)[i],lwd=1.5)
        dev.off()
        
        #estimate number of model
        m.var <- data.frame(Remove.metabolite[i])[,1]
        metabolite.name<-names(Remove.metabolite[i])
        
        density.all <- density(m.var, bw = "SJ")
        density.var <- density(m.var, bw = "SJ")$y
        density.x <-  density.all$x[density.var>0.01]
        density.var <-  density.var[density.var>0.01]
        
        modes <- NULL
        centre.p <- NULL
        for ( i in 2:(length(density.var)-1) ){
          if ( (density.var[i] > density.var[i-1]) & (density.var[i] > density.var[i+1]) ) {
            modes <- c(modes,i)
          }
          if ( (density.var[i] < density.var[i-1]) & (density.var[i] < density.var[i+1]) ) {
            centre.p <- c(centre.p,i)
          }
        }
        
        new.median<- density.x[modes]
        new.boundary<- density.x [centre.p]
        npeaks <-length(new.median)
        
        if(npeaks>1){
          
          vComposti <- vector("list", npeaks)
          m.var.sub <- m.var
          
          for (i in 1:(npeaks-1)){
            vComposti[i] <-list(m.var.sub [m.var.sub <new.boundary[i]])
            m.var.sub[m.var.sub<new.boundary[i]]<-NA
          }
          
          vComposti[npeaks] <- list(m.var.sub)
          
          RT.add.stats <- t(sapply(vComposti,function(x) c(Total.ID=length(x[!is.na(x)]),
                                                           RT.median=round(median(x,na.rm=TRUE),3),
                                                           RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*3))*60,2),
                                                           RT.shift=round(((quantile(x, prob = 0.75, na.rm=TRUE)+IQR(x,na.rm=TRUE)*3)-median(x,na.rm=TRUE))*60,2))))
          
          colnames(RT.add.stats)<-c("Total.ID","RT.median","RT.shfit.Lower", "RT.shfit.upper")
          AMDIS.report.RT.stats.add <- AMDIS.report.RT.stats[rownames(AMDIS.report.RT.stats)==metabolite.name,]
          AMDIS.report.RT.stats.add <- cbind(AMDIS.report.RT.stats.add[,1:3] ,RT.add.stats, AMDIS.report.RT.stats.add[,8:11], row.names = NULL)
          
          # re-name
          for (i in 1:npeaks){
            AMDIS.report.RT.stats.add$Name[i] <- paste(metabolite.name, " (split peak ", i, ")", sep="")
          }
          AMDIS.report.split <- rbind(AMDIS.report.split, AMDIS.report.RT.stats.add)
        }
      }
    }
    
    if(is.null(AMDIS.report.split)!=TRUE){
      AMDIS.report.RT.stats <- rbind(AMDIS.report.RT.stats, AMDIS.report.split)
    }
    if ("CAS" %in% names(AMDIS.report.RT.stats)) {
      library(stringr)
      AMDIS.report.RT.stats$CAS<-gsub(" ","",AMDIS.report.RT.stats$CAS)
      AMDIS.report.RT.stats$CAS<-gsub("-","",AMDIS.report.RT.stats$CAS)
      AMDIS.report.RT.stats_CAS_width<-str_length(AMDIS.report.RT.stats$CAS)
      cas_last<-str_sub(AMDIS.report.RT.stats$CAS,AMDIS.report.RT.stats_CAS_width,AMDIS.report.RT.stats_CAS_width)
      cas_last_23<-str_sub(AMDIS.report.RT.stats$CAS,AMDIS.report.RT.stats_CAS_width-2,AMDIS.report.RT.stats_CAS_width-1)
      cas_last_rest<-str_sub(AMDIS.report.RT.stats$CAS,1,AMDIS.report.RT.stats_CAS_width-3)
      casfinal<-data.frame(a=cas_last_rest,b=cas_last_23,c=cas_last,stringsAsFactors = F)
      cas_final<-unlist(lapply(1:nrow(casfinal),function(x,casfinal){
        str_glue(casfinal$a[x],"-",casfinal$b[x],"-",casfinal$c[x])
      },casfinal))
      AMDIS.report.RT.stats$CAS<-cas_final
    }
    AMDIS.report.RT.stats <- AMDIS.report.RT.stats[order(AMDIS.report.RT.stats$RT.median, decreasing = FALSE),]
    File.name <- tclvalue(tkgetSaveFile(initialfile=File.name))
    write.csv(AMDIS.report.RT.stats, file = File.name,row.names = FALSE)
    print(paste("AMDIS summary report.csv is generated and save in",workdir, sep=" "))
  }
}

#' mgf report file parsing and summary
#'
#' This is a function to process the AMDIS_report.txt file. specify a in-house library or a subset of NIST library to find the most intense quantification ion within the GC-MS experiment m/z range. The script will invoke selection windows to let the user specify the files location.  
#'
#' @param workdir locate the working Dir
#' @param MS.L specify the library file
#' @param MsLibrary specify the library file origin, could be "NIST" or "InHouse"
#' @param amdis.report specify the identification report file
#' @param Ret.Time.Filter set the retetion time filter window around the expected retention time (in +/-min).
#' @param RT.shift.limt set a threshole in seconds to define a molecule have 
#' @param mz_L lower mz limit for quantification ion selection
#' @param mz_U upper mz limit for quantification ion selection
#' @param generate_rt_shift_graph Set TRUE to generate the retetion time shift plots for each molecules that have retention time discrepencies greater than the \code{"RT.shift.limt"}
#' @param RTcorrection Set TRUE to enable a retention time correction before multiple peaks resolving
#' @return None
#'
#' @examples
#' mgf_id_Summary()
#'
#' @export
mgf_id_Summary<-function(workdir= getwd(),
                         mgf.report = NULL,mgf.peak.rds = NULL,
                           Ret.Time.Filter=2.5,
                           RT.shift.limt = 30,
                           generate_rt_shift_graph=T,
                           generate_rt_shift_graphs=F
){
  library("tcltk")
  library("tcltk2")
  # require("Metab")
  library("pander")
  library("parallel")
  #library("OneR")
  library("dplyr")
  library("ggplot2")
  library("ggpubr")
  library(stringr)
  #workdir <- workdir
  if (is.null(workdir)){workdir= tk_choose.dir(caption = "Select working directory")}
  setwd(workdir) 
  if (is.null(mgf.report) ){
    mgf.report = tk_choose.files(default = "Select the mgf report",
                                                            multi = T, caption = "Select the mgf report in .mgf",
                                                            filter=matrix(c("MGF", ".mgf", "All files", "*"),2, 2, byrow = TRUE))
    peaks_list<-NULL
    for (mgf_file in mgf.report){
      lib.txt<-readLines(mgf_file)
      lib.txt1<-lib.txt[-grep("\\t",lib.txt)]
      att.nms<-unique(sapply(strsplit(lib.txt1[grep("=",lib.txt1)],"="),function(x) x[1]))
      entry.nms<-sapply(strsplit(lib.txt1[grep("TITLE=",lib.txt1)],"TITLE="), function(x) x[2])
      rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))
      
      starts<-grep("TITLE=", lib.txt)
      stops<-c(starts[1:(length(starts)-1)]+diff(starts)-1,length(lib.txt))
      for (i in 1:length(starts)){
        tmp<-lib.txt[starts[i]:stops[i]]
        sapply(strsplit(tmp[grep("=",tmp)],"="), function(x) rez[i,x[1]]<<-x[2])
      }
      rez<-as.data.frame(rez)
      rez$file<-stringr::str_remove(basename(mgf_file),".mgf")
      rez->peaks_list[[stringr::str_remove(basename(mgf_file),".mgf")]]
    }
    peaks_list_bind<-as.data.frame(do.call(rbind,peaks_list))
    peaks_list_bind$RT<-as.numeric(peaks_list_bind$RTINMINUTES)
    peaks_list_bind$Name<-peaks_list_bind$TITLE<-str_to_title(peaks_list_bind$TITLE)
  } else if(!is.null(mgf.peak.rds)) {
    peaks_list_raw<-readRDS(mgf.peak.rds)
    lapply(names(peaks_list_raw),function(x,peaks_list_raw){
      peaks_list_raw[[x]]->rez
      rez<-as.data.frame(rez)
      rez$file<-stringr::str_remove(basename(mgf_file),".mgf")
    })->peaks_list
      
    peaks_list_bind<-as.data.frame(do.call(rbind,peaks_list))
    peaks_list_bind$RT<-as.numeric(peaks_list_bind$RTINMINUTES)
    peaks_list_bind$Name<-peaks_list_bind$TITLE<-str_to_title(peaks_list_bind$TITLE)
    }else{
    stop("no mgf input.")
      }
  
  message("mgf_id_Summary parameters")
  
    Metabolite.list <-split(peaks_list_bind$RT, peaks_list_bind$Name)
    RT.stats<-t(sapply(Metabolite.list,function(x) {
      x<-as.numeric(x)
      c(RT.median=round(median(x,na.rm=TRUE),3),
        RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5)),2),
        RT.shift=round(((quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)-median(x,na.rm=TRUE)),2)
    )}
    ))
    colnames(RT.stats)<-c("RT.median","RT.shfit.Lower", "RT.shfit.upper")
    RT.stats=RT.stats[is.na(RT.stats[,"RT.median"])!=T,]
    
    
    ID.stats<-t(sapply(split(peaks_list_bind$RT, peaks_list_bind$Name),function(x) c(RT.first=x[1],
                                                                                     Total.ID=length(x[!is.na(x)]))))
    Ref.ion.stats<-t(sapply(split(peaks_list_bind$MODELION, peaks_list_bind$Name),function(x) c(Ref.ion=paste(unique(x),collapse = ";"),
                                                                                                Ref.ion.unique=length(unique(x)))))
    report.RT.stats <- data.frame(cbind(Name=0,Ref.ion.stats, ID.stats, RT.stats, CAS=0))
    
    report.RT.stats$Name <- rownames(RT.stats)
    
    RT.list<-lapply(Metabolite.list,function(x) c(RT.shift=(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5),
                                                  RT.shift=(quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)))
    Remove.outlier<-function(x,y){
      x[x >= y[1]]
      x[x <= y[2]]
    }
    
    Remove.metabolite <-mapply(Remove.outlier, Metabolite.list, RT.list )
    
    # Plot all retention time of all metabolites
    
    n.total <- length(Metabolite.list)
    png(paste(getwd(),"/retention time of all metabolites.png",sep=""),width = 30,height = 30, bg = "white",units = "in",res = 600)
    par(mar=c(5,10,5,5))
    boxplot(Metabolite.list, xlab= "Retention time (min)", cex.axis=0.2,  horizontal = TRUE, las=2, xaxt="n", main = paste("Retention time distribution of all metabolites ","(n=",n.total,")",sep=""))
    axis(1,cex.axis=1,las=1)
    dev.off()
    
    
    ## solve multimodel issue
    
    
    
    report.split<-NULL
    deletenamelist<-NULL
    save.folder<-paste("RTshift_warnings_",RT.shift.limt,"s_","shift",sep="")
    if (dir.exists(save.folder)!=T) try(dir.create (save.folder))
    
    for (component in 1:length(Remove.metabolite)){
      windows_filename<- function(stringX){
        stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\|]")
        stringX<-gsub("\"", "", stringX)
        if (nchar(stringX)>=100){stringX=strtrim(stringX,100)}
        return(stringX)
        
      }
      
      if(`&`((RT.stats[component,2] + RT.stats[component,3]) > RT.shift.limt, T)){
        RT_libsimilarity=NULL
        
        names.var<- windows_filename(names(Remove.metabolite)[component])
        RT_libsimilarity<-peaks_list_bind[peaks_list_bind$Name==names(Remove.metabolite)[component],]
        m.var <- data.frame(Remove.metabolite[component])[,1]
        metabolite.name<-names(Remove.metabolite[component])
        res=NULL
        res <- try(density.all <- density(m.var, bw = "SJ"))
        if(class(res)=="try-error"){density.all <- density.default(m.var)}
        density.var <- density.all$y
        density.x <-  density.all$x[density.var>0.01]
        density.var <-  density.var[density.var>0.01]
        
        modes <- NULL
        centre.p <- NULL
        for ( i in 2:(length(density.var)-1) ){
          if ( (density.var[i] > density.var[i-1]) & (density.var[i] > density.var[i+1]) ) {
            modes <- c(modes,i)
          }
          if ( (density.var[i] < density.var[i-1]) & (density.var[i] < density.var[i+1]) ) {
            centre.p <- c(centre.p,i)
          }
        }
        
        new.median<- density.x[modes]
        new.boundary<- density.x [centre.p]
        npeaks <-length(new.median)
        
        if(npeaks>1){
          
          vComposti <- vector("list", npeaks)
          m.var.sub <- m.var
          
          for (i in 1:(npeaks-1)){
            vComposti[i] <-list(m.var.sub [m.var.sub <new.boundary[i]])
            m.var.sub[m.var.sub<new.boundary[i]]<-NA
          }
          
          vComposti[npeaks] <- list(m.var.sub)
          
          RT.add.stats <- t(sapply(vComposti,function(x) c(Total.ID=length(x[!is.na(x)]),
                                                           RT.median=round(median(x,na.rm=TRUE),3),
                                                           RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*3)),2),
                                                           RT.shift=round(((quantile(x, prob = 0.75, na.rm=TRUE)+IQR(x,na.rm=TRUE)*3)-median(x,na.rm=TRUE)),2))))
          
          colnames(RT.add.stats)<-c("Total.ID","RT.median","RT.shfit.Lower", "RT.shfit.upper")
          report.RT.stats.add <- report.RT.stats[rownames(report.RT.stats)==metabolite.name,]
          
          report.RT.stats.add <- cbind(report.RT.stats.add[,1:4] ,RT.add.stats, CAS=report.RT.stats.add[,"CAS"], row.names = NULL)
          report.RT.stats.add$Name->report.RT.stats.add$Name_org
          # re-name
          for (i in 1:npeaks){
            report.RT.stats.add$Name[i] <- paste(metabolite.name, " (split peak ", i, ")", sep="")
          }
          report.split <- rbind(report.split, report.RT.stats.add)
          deletenamelist<-c(deletenamelist,metabolite.name)
          
        }          
        if (generate_rt_shift_graph){          
          names.var<- windows_filename(names(Remove.metabolite)[component])     
          
          png(filename=paste(save.folder,"/",names.var,".png",sep=""))
          plot(RT_libsimilarity$RT,log(as.numeric(RT_libsimilarity$INTEGRATEDAREA)),  ylab = "AREA (log)",
               xlab= "Retention time (min)", main = names(Remove.metabolite)[component], lwd=1.5)
          dev.off()
        }
        
      }
    }
    
    if(is.null(report.split)!=TRUE){
      
      deletenamelist=unique(deletenamelist)
      
      report.RT.statscombine<-report.RT.stats[!(report.RT.stats$Name %in% (deletenamelist)),]
      
      report.RT.resolved<-report.RT.stats[(report.RT.stats$Name %in% (deletenamelist)),]
      
      
      report.split=unique(report.split)
      report.split<-report.split[report.split$Total.ID>1,]
      report.split<-report.split[report.split$Name_org!="Unknown",]
      for ( i in 1:nrow(report.split)){
        peaks_list_bind_cpd<-peaks_list_bind[peaks_list_bind$Name==report.split$Name_org[i],]
        rtrow<-between(peaks_list_bind_cpd$RT,report.split$RT.median[i]-report.split$RT.shfit.Lower[i],report.split$RT.shfit.upper[i]+report.split$RT.median[i])
        report.split$Ref.ion[i]<-paste0(unique(peaks_list_bind_cpd$MODELION[rtrow]),collapse=";")
      }
      #report.split$Name_org<-NULL
      report.RT.statscombine$Name->report.RT.statscombine$Name_org
      report.RT.stats <- rbind(report.RT.statscombine, report.split)
    }
    
    report.RT.stats <- report.RT.stats[order(report.RT.stats$RT.median, decreasing = FALSE),]
    File.nameNIST <- "Summary_ID_report.csv"
    write.csv(report.RT.stats, file = File.nameNIST,row.names = FALSE)
    
    message(paste("Summary report was created and save in",workdir, sep=" "))
    message(paste("Summarizing area report from mgf report",workdir, sep=" "))
    
    summary_temp<-NULL
      report.RT.stats$RT.median<-as.numeric(report.RT.stats$RT.median)
      report.RT.stats$RT.shfit.upper<-as.numeric(report.RT.stats$RT.shfit.upper)
      report.RT.stats$RT.shfit.Lower<-as.numeric(report.RT.stats$RT.shfit.Lower)
    for (i in 1:nrow(report.RT.stats)){
      
      df<-NULL
      if(report.RT.stats$Total.ID[i]==1){
        peaks_list_bind[peaks_list_bind$RT == report.RT.stats$RT.first[i],]->df
        df[ df$MODELION==report.RT.stats$Ref.ion[i],]->df
        df$Name<-report.RT.stats$Name[i]
        summary_temp[[report.RT.stats$Name[i]]]<-df
      }else{
        peaks_list_bind[between(peaks_list_bind$RT ,
                        report.RT.stats$RT.median[i]-report.RT.stats$RT.shfit.Lower[i]-0.001,
                        report.RT.stats$RT.shfit.upper[i]+report.RT.stats$RT.median[i]+0.001),]->df
        df[ df$MODELION %in% str_split(report.RT.stats$Ref.ion[i],";")[[1]],]->df
        df$Name<-report.RT.stats$Name[i] 
        summary_temp[[report.RT.stats$Name[i]]]<-df
       }
        
      
        
    }
      summary_temp_bind<-do.call(rbind,summary_temp)
      summary_temp_bind$INTEGRATEDAREA<-as.numeric(summary_temp_bind$INTEGRATEDAREA)
      library(reshape2)
      summary_temp_bind_cast_ID_count<-dcast(summary_temp_bind,Name ~ file,value.var = "INTEGRATEDAREA",fun.aggregate = length)
      summary_temp_bind_cast_ID_sum<-dcast(summary_temp_bind,Name ~ file,value.var = "INTEGRATEDAREA",fun.aggregate = sum)
      
      File.nameNIST="Summary_ID_count_report.csv"
      write.csv(summary_temp_bind_cast_ID_count, file = File.nameNIST,row.names = FALSE)
      File.nameNIST="Summary_ID_area_report.csv"
      write.csv(summary_temp_bind_cast_ID_sum, file = File.nameNIST,row.names = FALSE)
      message(paste("Summarizing area report saved in",workdir, sep=" "))
      
}
#' MassHunter report file (.cefs) parsing and summary
#'
#' This is a function to process the AMDIS_report.txt file. specify a in-house library or a subset of NIST library to find the most intense quantification ion within the GC-MS experiment m/z range. The script will invoke selection windows to let the user specify the files location.  
#'
#' @param workdir locate the working Dir
#' @param MS.L specify the library file
#' @param MsLibrary specify the library file origin, could be "NIST" or "InHouse"
#' @param amdis.report specify the identification report file
#' @param Ret.Time.Filter set the retetion time filter window around the expected retention time (in +/-min).
#' @param RT.shift.limt set a threshole in seconds to define a molecule have 
#' @param mz_L lower mz limit for quantification ion selection
#' @param mz_U upper mz limit for quantification ion selection
#' @param generate_rt_shift_graph Set TRUE to generate the retetion time shift plots for each molecules that have retention time discrepencies greater than the \code{"RT.shift.limt"}
#' @param RTcorrection Set TRUE to enable a retention time correction before multiple peaks resolving
#' @return None
#'
#' @examples
#' MassHunter_id_Summary()
#'
#' @export
MassHunter_id_Summary<-function(workdir= NULL,
                           MS.L= NULL,
                           MassHunter.report = NULL,
                           File.name = "Summary Report.csv",
                           MsLibrary=c("NIST", "InHouse"),
                           Ret.Time.Filter=2.5,
                           RT.shift.limt = 60,
                           mz_L=38,
                           mz_U=550,
                           generate_rt_shift_graph=F,
                           generate_rt_shift_graphs=F,
                           RTcorrection=F
){
  library("tcltk")
  library("tcltk2")
  # require("Metab")
  library("pander")
  library("parallel")
  #library("OneR")
  library("dplyr")
  library("ggplot2")
  library("ggpubr")
  library(reticulate)
  #os <- import("os")
  #os$listdir(".")
  #py_install("fuzzywuzzy")
  #virtualenv_create("r-reticulate")
  #virtualenv_install("r-reticulate", "fuzzywuzzy")
  #system("pip install tk")
  if (is.null(workdir)){workdir= tk_choose.dir(caption = "Select working directory")}
  if (is.null(MassHunter.report)){MassHunter.report = tk_choose.dir(default = "",caption = "Select the MassHunter report folder contains .cef files")}
  if (is.null(MS.L)){MS.L= tk_choose.files(caption="Select MS library (e.g. SVB_total or NIST) in .msl",
                                           filter=matrix(c("MSL", ".msl", "All files", "*"),2, 2, byrow = TRUE))}
  
 #source_python(paste0(file.path(path.package(package="MassOmics")),"/src/masshunter-massOmics-summary-report.py"))
  #rm(parse_cef)
  #rm(return_df)
  setwd(workdir)
  source_python(paste0(file.path(path.package(package="MassOmics")),"/src/parse_Agilent_CEF.py"))
  return_df<-parse_cef(Folder = MassHunter.report)
  #workdir <- workdir
  
  
  

  
  
  setwd(workdir)
  MassHunter.report =as.data.frame (return_df)
  colnames(MassHunter.report)<-c("Name","n" ,"RT" , "ref","area","height","CAS" )
  message("MassHunter ID Summary parameters")
  message(paste("RT filter: ", Ret.Time.Filter, " min", sep=""))
  message(paste("Detect mutiple peaks when RT range is greater than: ", RT.shift.limt, " s", sep=""))
  message(paste("Quant mass will be selected within: ", mz_L, "-", mz_U,sep=""))
  message(paste("Perform RT correction before multi-peak resolving: ",RTcorrection ,sep=""))
  # Extract Reference ion, CAS # from MS library
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
  
  # ExtractIon for NIST library
  if(MsLibrary=="NIST"){
    message("Running in NIST mode, will generate Ref ions for each metabolites within mz detection range")
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    reference_ion<-NULL
    
    #test<-buildLib(AmdisLib=MS.L, folder=workdir, save = F, output = "ion_lib.csv", verbose = F,mz_L = mz_L,mz_U = mz_U)
    
    
    #reference_ion=parLapply(cl=autoStopCluster(makeCluster(detectCores()-1)), 1:length(starts),parse_msl_par,lib.txt,starts,stops,mz_L,mz_U)
    
    #reference_ion=unlist(reference_ion)
    
    #sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
    reference_ion=0
    
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez$NAME<-entry.nms
    libr<-cbind(rez,reference_ion)## returns list of peaks
    libr$NAME <- gsub("?", "", libr$NAME, fixed = TRUE)
    libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
    libr$NAME <- gsub("^ ", "", libr$NAME, perl=T)
    libr$NAME <- gsub("[*:%^!;$]", "", libr$NAME, perl=T)
    
    # Generate Summary report
    MassHunter.report <-MassHunter.report
    MassHunter.report$Name <- gsub("?", "", MassHunter.report$Name, fixed = TRUE)
    MassHunter.report$Name <- gsub("^ ", "", MassHunter.report$Name, perl=T)
    MassHunter.report$Name <- gsub("[*:%^!*?&;$]", "", MassHunter.report$Name, perl=T)
    MassHunter.report$Width <- 10
    #MassHunter.report$Width <- as.numeric(gsub("scans", "",  MassHunter.report$Width, perl=T))
    #MassHunter.reportcas<-MassHunter.report[MassHunter.report$CAS %in% libr$CASNO,]
    MassHunter.report<-MassHunter.report[MassHunter.report$Name %in% gsub("^ ","",libr$NAME),]
    #MassHunter.report1<-MassHunter.report[MassHunter.report$Name %in% libr$NAME,]
    #
    
    if (RTcorrection){
      library(ggpubr)
      RT.correction.grid<-RT_correction_xcms()
      RT.correction.result<-RT_correction(df=MassHunter.report)
      MassHunter.report<-RT.correction.result[[1]]
      ggpubr::ggarrange(RT.correction.result[[2]],RT.correction.result[[3]],
                        labels = c("Before","After"),common.legend = T,
                        ncol = 2, nrow =1,vjust=1.5) %>%  ggexport(filename = paste(getwd(),"/retention time correction.png",sep=""),width = 1800, height = 900,verbose = NULL,res = 100)
      
    }
    Metabolite.list <-split(MassHunter.report$RT, MassHunter.report$Name)
    RT.stats<-t(sapply(Metabolite.list,function(x) c(RT.median=round(median(x,na.rm=TRUE),3),
                                                     RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5))*60,2),
                                                     RT.shift=round(((quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)-median(x,na.rm=TRUE))*60,2)
    )))
    colnames(RT.stats)<-c("RT.median","RT.shfit.Lower", "RT.shfit.upper")
    RT.stats=RT.stats[is.na(RT.stats[,"RT.median"])!=T,]
    Peak.width <- sapply(split(MassHunter.report$Width, MassHunter.report$Name), function(x) median(x))
    
    
    ID.stats<-t(sapply(split(MassHunter.report$height, MassHunter.report$Name),function(x) c(Library.match=round(mean(x,na.rm=T)),
                                                                                     Total.ID=length(x[!is.na(x)]))))
    #message("Running in NIST mode, will generate Ref ions for each metabolites within mz detection range")
    MassHunter.report.RT.stats <- data.frame(cbind(Name=0,Ref.ion=0, ID.stats, RT.stats, Peak.width, CAS=0))
    #message(paste(MassHunter.report.RT.stats$Ref.ion[1:40],sep="\n"))
    #MassHunter.report.RT.stats$Ref.ion <-  libr$reference_ion[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    #matchID=match(rownames(RT.stats),libr$NAME)
    matchID=match(strtrim(rownames(RT.stats),170),strtrim(libr$NAME,170))
    #MassHunter.report.RT.stats=MassHunter.report.RT.stats[which(!is.na(matchID)),]
    #matchID70=match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))
    #intersect(matchID,matchID70)
    MassHunter.report.RT.stats$Ref.ion <-  unlist(parallel::parLapply(cl=autoStopCluster(makeCluster(detectCores()-1)), 
                                                                 matchID,
                                                                 parse_msl_par,lib.txt,starts,stops,mz_L,mz_U))
    
    #for (ids in matchID){
    # Ref.iontest= parse_msl_par(ids,lib.txt,starts,stops,mz_L,mz_U)
    
    #}
    #match("2-Cyclopropylethynylcyclopropane",libr$NAME)
    #message(MassHunter.report.RT.stats$Ref.ion)
    #which(match(rownames(RT.stats),libr$NAME)==NA)
    MassHunter.report.RT.stats$CAS <-  libr$CAS[match(strtrim(rownames(RT.stats),70),strtrim(libr$NAME,70))]
    MassHunter.report.RT.stats$Name <- rownames(RT.stats)
    
    #total.n<-length(unique(MassHunter.report$FileName))
    #ID_name<-paste("Total.ID(n=",total.n,")",sep="" )
    #names(MassHunter.report.RT.stats)[3]<-ID_name
    
    
    ## Remove point outside the range
    RT.list<-lapply(Metabolite.list,function(x) c(RT.shift=(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*1.5),
                                                  RT.shift=(quantile(x, prob = 0.75,na.rm=TRUE)+IQR(x,na.rm=TRUE)*1.5)))
    Remove.outlier<-function(x,y){
      x[x >= y[1]]
      x[x <= y[2]]
    }
    
    Remove.metabolite <-mapply(Remove.outlier, Metabolite.list, RT.list )
    
    # Plot all retention time of all metabolites
    
    n.total <- length(Metabolite.list)
    png(paste(getwd(),"/retention time of all metabolites.png",sep=""),width = 30,height = 30, bg = "white",units = "in",res = 150)
    par(mar=c(5,10,5,5))
    boxplot(Metabolite.list, xlab= "Retention time (min)", cex.axis=0.2,  horizontal = TRUE, las=2, xaxt="n", main = paste("Retention time distribution of all metabolites ","(n=",n.total,")",sep=""))
    axis(1,cex.axis=1,las=1)
    dev.off()
    
    
    ## solve multimodel issue
    
    
    
    MassHunter.report.split<-NULL
    deletenamelist<-NULL
    save.folder<-paste("RTshift_warnings_",RT.shift.limt,"s_","shift",sep="")
    if (dir.exists(save.folder)!=T) try(dir.create (save.folder))
    
    for (component in 1:length(Remove.metabolite)){
      if(`&`((RT.stats[component,2] + RT.stats[component,3]) > RT.shift.limt, T)){
        #names(Remove.metabolite)[component]
        RT_libsimilarity=NULL
        #plot figure
        names.var<- windows_filename(names(Remove.metabolite)[component])
        RT_libsimilarity<-MassHunter.report[MassHunter.report$Name==names(Remove.metabolite)[component],]
        #nbins=6
        #res <- try(RT_libsimilarity$peakgroup<-OneR::bin(RT_libsimilarity$RT,nbins = nbins,labels = c(paste("Peak",1:nbins)), method = "clusters",na.omit = T),silent = TRUE)
        #if (class(res) != "try-error"){
        #while(`&`(sum(grep("Peak",RT_libsimilarity$peakgroup))==0,nbins>=1)){
        #  nbins=nbins-1
        #  res <- try(RT_libsimilarity$peakgroup<-OneR::bin(RT_libsimilarity$RT,nbins = nbins,labels = c(paste("Peak",1:nbins)), method = "clusters",na.omit = T),silent = TRUE)
        #}
        #RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT
        #if (sum(grep("Peak",RT_libsimilarity$peakgroup))!=0){
        #  RT_libsimilarity$peakgroup_RT_median=NULL
        #  RT_libsimilarity.stats<-RT_libsimilarity %>% group_by(peakgroup) %>% summarise(peakgroup_RT_median = median(RT))
        # RT_libsimilarity=merge(RT_libsimilarity,RT_libsimilarity.stats)
        
        # }
        #RT_libsimilarity<-MassHunter.report[MassHunter.report$Name==names(Remove.metabolite)[component],]
        #RT_libsimilarity$peakgroupRT=median(RT_libsimilarity$peakgroup)
        #RT_libsimilarity$series=1:nrow(RT_libsimilarity)
        #try(RT_libsimilarity<-RT_libsimilarity[order(RT_libsimilarity$FileName),])
        
        #estimate number of model
        m.var <- data.frame(Remove.metabolite[component])[,1]
        metabolite.name<-names(Remove.metabolite[component])
        res=NULL
        res <- try(density.all <- density(m.var, bw = "SJ"))
        if(class(res)=="try-error"){density.all <- density.default(m.var)}
        density.var <- density.all$y
        density.x <-  density.all$x[density.var>0.01]
        density.var <-  density.var[density.var>0.01]
        
        modes <- NULL
        centre.p <- NULL
        for ( i in 2:(length(density.var)-1) ){
          if ( (density.var[i] > density.var[i-1]) & (density.var[i] > density.var[i+1]) ) {
            modes <- c(modes,i)
          }
          if ( (density.var[i] < density.var[i-1]) & (density.var[i] < density.var[i+1]) ) {
            centre.p <- c(centre.p,i)
          }
        }
        
        new.median<- density.x[modes]
        new.boundary<- density.x [centre.p]
        npeaks <-length(new.median)
        
        if(npeaks>1){
          
          vComposti <- vector("list", npeaks)
          m.var.sub <- m.var
          
          for (i in 1:(npeaks-1)){
            vComposti[i] <-list(m.var.sub [m.var.sub <new.boundary[i]])
            m.var.sub[m.var.sub<new.boundary[i]]<-NA
          }
          
          vComposti[npeaks] <- list(m.var.sub)
          
          RT.add.stats <- t(sapply(vComposti,function(x) c(Total.ID=length(x[!is.na(x)]),
                                                           RT.median=round(median(x,na.rm=TRUE),3),
                                                           RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*3))*60,2),
                                                           RT.shift=round(((quantile(x, prob = 0.75, na.rm=TRUE)+IQR(x,na.rm=TRUE)*3)-median(x,na.rm=TRUE))*60,2))))
          
          colnames(RT.add.stats)<-c("Total.ID","RT.median","RT.shfit.Lower", "RT.shfit.upper")
          MassHunter.report.RT.stats.add <- MassHunter.report.RT.stats[rownames(MassHunter.report.RT.stats)==metabolite.name,]
          MassHunter.report.RT.stats.add <- cbind(MassHunter.report.RT.stats.add[,1:3] ,RT.add.stats, MassHunter.report.RT.stats.add[,8:9], row.names = NULL)
          
          # re-name
          for (i in 1:npeaks){
            MassHunter.report.RT.stats.add$Name[i] <- paste(metabolite.name, " (split peak ", i, ")", sep="")
          }
          MassHunter.report.split <- rbind(MassHunter.report.split, MassHunter.report.RT.stats.add)
          deletenamelist<-c(deletenamelist,metabolite.name)
          
        }          
        if (generate_rt_shift_graph){          
          #png(filename=paste(save.folder,"/",names.var,".png",sep=""),width = 480,height = 720,)
          #par(mfrow = c(2, 1))
          #png(filename=paste(save.folder,"/",names.var,"lib_score.png",sep=""),width = 900,height = 1200)
          #par(mfrow = c(2, 3))
          #RT_libsimilarity<-MassHunter.report[MassHunter.report$Name==names(Remove.metabolite)[i],]
          #RT_libsimilarity$peakgroup="Outliner"
          #RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT
          
          names.var<- windows_filename(names(Remove.metabolite)[component])          
          
          
          #boxplotlibsim=as_ggplot(boxplot(Weighted~peakgroup,data=RT_libsimilarity, main=names(Remove.metabolite)[i], xlab="Retention time (min)", ylab="Library Score"))
          
          #rtdensityplot=plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",xlab= "Retention time (min)", main = names(Remove.metabolite)[i], lwd=1.5)
          #dev.off()
          #plot figure
          
          png(filename=paste(save.folder,"/",names.var,".png",sep=""))
          #png(filename=paste(save.folder,"/","aaa",".png",sep=""))
          
          plot(RT_libsimilarity$RT,RT_libsimilarity$Net,  ylab = "Net score",
               xlab= "Retention time (min)", main = names(Remove.metabolite)[component], lwd=1.5)
          dev.off()
          
          
          #png(filename=paste(save.folder,"/",names.var,"_lib_score.png",sep=""),width = 400, height = 400)
          
          
        }
        if (generate_rt_shift_graphs){          
          #png(filename=paste(save.folder,"/",names.var,".png",sep=""),width = 480,height = 720,)
          #par(mfrow = c(2, 1))
          #png(filename=paste(save.folder,"/",names.var,"lib_score.png",sep=""),width = 900,height = 1200)
          #par(mfrow = c(2, 3))
          #RT_libsimilarity<-MassHunter.report[MassHunter.report$Name==names(Remove.metabolite)[i],]
          RT_libsimilarity$peakgroup="Outliner"
          RT_libsimilarity$peakgroup_RT_median=RT_libsimilarity$RT
          for (peak in 1: nrow(MassHunter.report.RT.stats.add)){
            RT_libsimilarity$peakgroup=ifelse(`&`(RT_libsimilarity$RT>=(MassHunter.report.RT.stats.add$RT.median[peak]-MassHunter.report.RT.stats.add$RT.shfit.Lower[peak]/60),
                                                  RT_libsimilarity$RT<=(MassHunter.report.RT.stats.add$RT.median[peak]+MassHunter.report.RT.stats.add$RT.shfit.upper[peak]/60)),paste("Peak",peak),RT_libsimilarity$peakgroup)
            RT_libsimilarity$peakgroup_RT_median=ifelse(`&`(RT_libsimilarity$RT>=(MassHunter.report.RT.stats.add$RT.median[peak]-MassHunter.report.RT.stats.add$RT.shfit.Lower[peak]/60),
                                                            RT_libsimilarity$RT<=(MassHunter.report.RT.stats.add$RT.median[peak]+MassHunter.report.RT.stats.add$RT.shfit.upper[peak]/60)),MassHunter.report.RT.stats.add$RT.median[peak],"")
          }
          names.var<- windows_filename(names(Remove.metabolite)[component])          
          boxplotlibsim <- ggplot(RT_libsimilarity, aes(x=peakgroup, y=Net, fill=peakgroup)) + 
            geom_boxplot() +  geom_jitter(shape=16, position=position_jitter(0.2)) + theme(
              panel.background = element_rect(fill = "white"),
              plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
              plot.background = element_rect(
                fill = "grey90",
                colour = "black"
              ))
          
          #boxplotlibsim=as_ggplot(boxplot(Weighted~peakgroup,data=RT_libsimilarity, main=names(Remove.metabolite)[i], xlab="Retention time (min)", ylab="Library Score"))
          
          #rtdensityplot=plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",xlab= "Retention time (min)", main = names(Remove.metabolite)[i], lwd=1.5)
          #dev.off()
          
          
          #png(filename=paste(save.folder,"/","aaa",".png",sep=""))
          
          Net.plot=ggplot(RT_libsimilarity, aes(x=RT, y=Net, color=peakgroup,size=2)) + geom_point() + theme(
            panel.background = element_rect(fill = "white"),
            plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
            plot.background = element_rect(
              fill = "grey90",
              colour = "black"
            ))
          
          ggpubr::ggarrange(boxplotlibsim,Net.plot,
                            labels = c("boxplot","Net score"),
                            ncol = 2, nrow =1,vjust=1.5) %>%  ggexport(filename = paste(save.folder,"/",names.var,"_lib_score.png",sep=""),width = 900, height = 400,verbose = NULL)
          #dev.off()
        }
      }
    }
    
    if(is.null(MassHunter.report.split)!=TRUE){
      
      deletenamelist=unique(deletenamelist)
      
      MassHunter.report.RT.statscombine<-MassHunter.report.RT.stats[!(MassHunter.report.RT.stats$Name %in% (deletenamelist)),]
      
      MassHunter.report.split=unique(MassHunter.report.split)
      
      MassHunter.report.RT.stats <- rbind(MassHunter.report.RT.statscombine, MassHunter.report.split)
    }
    
    MassHunter.report.RT.stats <- MassHunter.report.RT.stats[order(MassHunter.report.RT.stats$RT.median, decreasing = FALSE),]
    if ("CAS" %in% names(MassHunter.report.RT.stats)) {
      library(stringr)
      MassHunter.report.RT.stats$CAS<-gsub(" ","",MassHunter.report.RT.stats$CAS)
      MassHunter.report.RT.stats$CAS<-gsub("-","",MassHunter.report.RT.stats$CAS)
      MassHunter.report.RT.stats_CAS_width<-str_length(MassHunter.report.RT.stats$CAS)
      cas_last<-str_sub(MassHunter.report.RT.stats$CAS,MassHunter.report.RT.stats_CAS_width,MassHunter.report.RT.stats_CAS_width)
      cas_last_23<-str_sub(MassHunter.report.RT.stats$CAS,MassHunter.report.RT.stats_CAS_width-2,MassHunter.report.RT.stats_CAS_width-1)
      cas_last_rest<-str_sub(MassHunter.report.RT.stats$CAS,1,MassHunter.report.RT.stats_CAS_width-3)
      casfinal<-data.frame(a=cas_last_rest,b=cas_last_23,c=cas_last,stringsAsFactors = F)
      cas_final<-unlist(lapply(1:nrow(casfinal),function(x,casfinal){
        str_glue(casfinal$a[x],"-",casfinal$b[x],"-",casfinal$c[x])
      },casfinal))
      MassHunter.report.RT.stats$CAS<-cas_final
    }
    File.nameNIST <- "Summary report_NIST.csv"
    write.csv(MassHunter.report.RT.stats, file = File.nameNIST,row.names = FALSE)
    write.csv(libr, file = "Library summary.csv",row.names = FALSE)
    
    
    
    
    message(paste("Summary report_NIST.csv was created and save in",workdir, sep=" "))
  } else if (MsLibrary=="InHouse") {
    
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez$NAME<-entry.nms
    libr<-rez## returns list of peaks
    
    # Generate summary report
    MassHunter.report<-MassHunter.report
    MassHunter.report<-MassHunter.report[MassHunter.report$RT-MassHunter.report$Expec..RT < Ret.Time.Filter,]
    MassHunter.report<-MassHunter.report[MassHunter.report$RT-MassHunter.report$Expec..RT > -Ret.Time.Filter,]
    MassHunter.report$Name <- gsub("?", "", MassHunter.report$Name, fixed = TRUE)
    MassHunter.report$Name <- gsub("^ ", "", MassHunter.report$Name, perl=T)
    MassHunter.report$Name <- gsub("[*:%^!*?&;$]", "", MassHunter.report$Name, perl=T)
    MassHunter.report$Width <- gsub(">", "",  MassHunter.report$Width, perl=T)
    MassHunter.report$Width <- as.numeric(gsub("scans", "",  MassHunter.report$Width, perl=T))
    
    
    
    Metabolite.list <-split(MassHunter.report$RT, MassHunter.report$Name)
    
    RT.stats<-t(sapply(Metabolite.list,function(x) c(RT.median=round(median(x,na.rm=TRUE),3),
                                                     RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25)-IQR(x,na.rm=TRUE)*1.5))*60,2),
                                                     RT.shift=round(((quantile(x, prob = 0.75)+IQR(x,na.rm=TRUE)*1.5)-median(x,na.rm=TRUE))*60,2))))
    colnames(RT.stats)<-c("RT.median","RT.shfit.Lower", "RT.shfit.upper")
    
    Peak.width <- sapply(split(MassHunter.report$Width, MassHunter.report$Name), function(x) median(x))
    
    
    RT.stats[is.na(RT.stats)]<-0
    ID.stats<-t(sapply(split(MassHunter.report$Weighted, MassHunter.report$Name),function(x) c(Library.match=round(mean(x,na.rm=T)),
                                                                                     Total.ID=length(x[!is.na(x)]))))
    MassHunter.report.RT.stats<-NULL
    MassHunter.report.RT.stats <- as.data.frame(cbind(Name=0,Ref.ion=0,ID.stats,RT.stats,Expec.RT=0,Diff.RT=0, Peak.width, CAS=0))
    MassHunter.report.RT.stats$Ref.ion <- libr$RSN[match(rownames(RT.stats),libr$NAME)]
    MassHunter.report.RT.stats$Expec.RT <- round(as.numeric(libr$RT[match(rownames(RT.stats),libr$NAME)]),2)
    MassHunter.report.RT.stats$CAS <- libr$CASNO[match(rownames(RT.stats),libr$NAME)]
    MassHunter.report.RT.stats$Diff.RT <- MassHunter.report.RT.stats$Expec.RT- as.numeric(MassHunter.report.RT.stats$RT.median)
    MassHunter.report.RT.stats$Name <- rownames(RT.stats)
    
    
    
    
    
    ## Remove point outside the range
    RT.list<-lapply(Metabolite.list,function(x) c(RT.shift=(quantile(x, prob = 0.25)-IQR(x,na.rm=TRUE)*1.5),
                                                  RT.shift=(quantile(x, prob = 0.75)+IQR(x,na.rm=TRUE)*1.5)))
    Remove.outlier<-function(x,y){
      x[x >= y[1]]
      x[x <= y[2]]
    }
    
    Remove.metabolite <-mapply(Remove.outlier, Metabolite.list, RT.list )
    
    # Plot all retention time of all metabolites
    n.total <- length(Metabolite.list)
    par(mar=c(5,10,5,5))
    boxplot(Metabolite.list, xlab= "Retention time (min)", cex.axis=0.5,  horizontal = TRUE, las=2, xaxt="n", main = paste("Retention time distribution of all metabolites ","(n=",n.total,")",sep=""))
    axis(1,cex.axis=1,las=1)
    
    
    ## solve multimodel issue
    MassHunter.report.split<-NULL
    
    save.folder<-paste("RTshift_warnings_",RT.shift.limt,"s_","shift",sep="")
    if (dir.exists(save.folder)!=T) try(dir.create (save.folder))
    for (i in 1:length(Remove.metabolite)){
      if(RT.stats[i,2] + RT.stats[i,3] > RT.shift.limt ){
        
        # plot figure
        names.var<- gsub("/", "", names(Remove.metabolite)[i], perl=T)
        png(filename=paste(save.folder,"/",names.var,".png",sep=""))
        plot(data.frame(Remove.metabolite[i])[,1],1:length(data.frame(Remove.metabolite[i])[,1]),  ylab = "Density",
             xlab= "Retention time (min)", main = names(Remove.metabolite)[i],lwd=1.5)
        dev.off()
        
        #estimate number of model
        m.var <- data.frame(Remove.metabolite[i])[,1]
        metabolite.name<-names(Remove.metabolite[i])
        
        density.all <- density(m.var, bw = "SJ")
        density.var <- density(m.var, bw = "SJ")$y
        density.x <-  density.all$x[density.var>0.01]
        density.var <-  density.var[density.var>0.01]
        
        modes <- NULL
        centre.p <- NULL
        for ( i in 2:(length(density.var)-1) ){
          if ( (density.var[i] > density.var[i-1]) & (density.var[i] > density.var[i+1]) ) {
            modes <- c(modes,i)
          }
          if ( (density.var[i] < density.var[i-1]) & (density.var[i] < density.var[i+1]) ) {
            centre.p <- c(centre.p,i)
          }
        }
        
        new.median<- density.x[modes]
        new.boundary<- density.x [centre.p]
        npeaks <-length(new.median)
        
        if(npeaks>1){
          
          vComposti <- vector("list", npeaks)
          m.var.sub <- m.var
          
          for (i in 1:(npeaks-1)){
            vComposti[i] <-list(m.var.sub [m.var.sub <new.boundary[i]])
            m.var.sub[m.var.sub<new.boundary[i]]<-NA
          }
          
          vComposti[npeaks] <- list(m.var.sub)
          
          RT.add.stats <- t(sapply(vComposti,function(x) c(Total.ID=length(x[!is.na(x)]),
                                                           RT.median=round(median(x,na.rm=TRUE),3),
                                                           RT.shift=round((median(x,na.rm=TRUE)-(quantile(x, prob = 0.25,na.rm=TRUE)-IQR(x,na.rm=TRUE)*3))*60,2),
                                                           RT.shift=round(((quantile(x, prob = 0.75, na.rm=TRUE)+IQR(x,na.rm=TRUE)*3)-median(x,na.rm=TRUE))*60,2))))
          
          colnames(RT.add.stats)<-c("Total.ID","RT.median","RT.shfit.Lower", "RT.shfit.upper")
          MassHunter.report.RT.stats.add <- MassHunter.report.RT.stats[rownames(MassHunter.report.RT.stats)==metabolite.name,]
          MassHunter.report.RT.stats.add <- cbind(MassHunter.report.RT.stats.add[,1:3] ,RT.add.stats, MassHunter.report.RT.stats.add[,8:11], row.names = NULL)
          
          # re-name
          for (i in 1:npeaks){
            MassHunter.report.RT.stats.add$Name[i] <- paste(metabolite.name, " (split peak ", i, ")", sep="")
          }
          MassHunter.report.split <- rbind(MassHunter.report.split, MassHunter.report.RT.stats.add)
        }
      }
    }
    
    if(is.null(MassHunter.report.split)!=TRUE){
      MassHunter.report.RT.stats <- rbind(MassHunter.report.RT.stats, MassHunter.report.split)
    }
    
    MassHunter.report.RT.stats <- MassHunter.report.RT.stats[order(MassHunter.report.RT.stats$RT.median, decreasing = FALSE),]
    File.name <- tclvalue(tkgetSaveFile(initialfile=File.name))
    write.csv(MassHunter.report.RT.stats, file = File.name,row.names = FALSE)
    print(paste("MassHunter summary report.csv is generated and save in",workdir, sep=" "))
  }
}


#' Peak diagnosis by overlay multiple chromatograms
#'
#'The user could visually confirm all peaks, and determine whether the peak is due to background noise by selecting <Overlay All Chromatograms> and clicking the <RUN> button in the box of step 3. A popout window will show a figure that plots chromatograms for all identified metabolites. You can also adjust the font size of metabolite names if the display is inappropriate. You must disable the <Fast Mode> to plot overlay of all chromatograms.
#'
#' @param font.size Size of report font
#' 
#' @return None
#'
#' @examples
#' PeakDiagnosis()
#'
#' @export
PeakDiagnosis <- function(path="GC-Peak Diagnosis Report.csv",font.size=0.5){
  if (!file.exists(path)){path=tcltk::tk_choose.files(default ="GC-Peak Diagnosis Report.csv", caption="choose the Peak Diagnosis Report.csv",multi = F)}
  IntegrationReport<-read.csv(path,stringsAsFactors = F)
  require(lattice)
  dev.new.OS()
  if(exists("redwarning")){
    xyp<-xyplot(Intensity~Retention.Time| factor(Metabolite.Names), groups=Sample.Names, type="l", data=IntegrationReport,
                plot.points=FALSE, scales = list(x ="free",y="free"), as.table=TRUE,
                main="GC-Peak overlay chromatogram",
                xlab="Rentention Time (min)",
                ylab="Peak Intensity",
                par.strip.text=list(cex=font.size),
                par.settings = list(axis.text=list(cex=0.5)),
                
                
                strip = function(..., which.panel, bg) {
                  bg.col = redwarning
                  strip.default(..., which.panel = which.panel,
                                bg = rep(bg.col, length = which.panel)[which.panel])}
                
    )
    print(xyp)
  } else {
    IntegrationReport<-read.csv(path,stringsAsFactors = F)
    dev.new.OS()
    
    xyp<-xyplot(Intensity~Retention.Time| factor(Metabolite.Names), groups=Sample.Names, type="l", data=IntegrationReport,
                plot.points=FALSE, scales = list(x ="free",y="free"), as.table=TRUE,
                main="GC-Peak Integration Report",
                par.strip.text=list(cex=font.size),
                par.settings = list(strip.background = list(col = "gray95"),axis.text=list(cex=0.5))
    )
    print(xyp)
  }
}

multimodel<-function(data,nbins=5,labels=1:nbins,method=c("clusters","length", "content"),na.omit=T){
  library("OneR")
  library("dplyr")
  
  
  
  
}

RT_correction<-function(df,topN=0.90,RT_std_limit=0.002,RT_drift_limit=0.2){
  
  library(dplyr)
  library(dplyr)
  library(ggplot2)
  df = data.table::data.table(df,stringsAsFactors = F)
  
  dfsum <- df %>% dplyr::group_by(Name) %>% dplyr::summarise(unqiuefilenum=length(unique(FileName)),medianRT=median(RT),stdRT=std.outliner.rm(RT),sdRT=sd(RT))
  
  dfsum$stdRT=ifelse(is.na(dfsum$stdRT),0,dfsum$stdRT)
  
  dfsum<-dfsum[dfsum$unqiuefilenum>=topN*max(dfsum$unqiuefilenum),]
  
  dfsum<-dfsum[dfsum$stdRT<=RT_std_limit,]
  
  dfdatafilesum<-df %>% dplyr::group_by(FileName,Name) %>% dplyr::summarise(medianRT=median(as.numeric(RT)))
  
  dfdatafilesum<-dfdatafilesum[dfdatafilesum$Name %in% dfsum$Name,]
  
  
  
  dfdatafilesum$Name=as.factor(dfdatafilesum$Name)
  
  beforecorrection<-ggplot2::ggplot(dfdatafilesum,aes(x=medianRT,y=FileName,color=Name)) + geom_point(size=1) +
    theme(axis.text.y=element_blank())
  
  
  
  df$newRT<- parSapply(cl=autoStopCluster(makeCluster(try(detectCores()))),1:nrow(df),RT_correction_LORESS,df,dfsum)
  df$newRT<-as.numeric(unlist(df$newRT))
  df$rtdiff=abs(df$RT-df$newRT)
  df$finalRT=ifelse('&'(!is.na(df$newRT),df$rtdiff<RT_drift_limit),df$newRT,df$RT)
  df$RT=df$finalRT
  
  dfdatafilesum<-df %>% dplyr::group_by(FileName,Name) %>% dplyr::summarise(medianRT=median(as.numeric(RT)))
  
  dfdatafilesum<-dfdatafilesum[dfdatafilesum$Name %in% dfsum$Name,]
  
  dfdatafilesum$Name=as.factor(dfdatafilesum$Name)
  aftercorrection<-ggplot2::ggplot(dfdatafilesum,aes(x=medianRT,y=FileName,color=Name)) + geom_point(size=1) +
    theme(axis.text.y=element_blank())
  
  
  return(list(df,beforecorrection,aftercorrection))
  }

RT_correction_LORESS<-function(i,df,dfsum){
  
  #df$FileName[i]
  #df$Name[i]
  refdf=df['&'(df$FileName==df$FileName[i],df$Name %in% dfsum$Name),]
  refcpd=unique(dfsum$Name[which(abs(dfsum$medianRT-df$RT[i])==min(abs(dfsum$medianRT-df$RT[i])))])
  refRTmedian=dfsum$medianRT[dfsum$Name %in% refcpd]
  refRT=unique(refdf$RT[refdf$Name %in% refcpd])
  RTdiff=(refRTmedian-refRT)[abs(refRTmedian-refRT)==min(abs(refRTmedian-refRT))]
  if (length(df$RT[i]+(RTdiff))==1) {
  return(df$RT[i]+(RTdiff))  
  }else{return("")}
  
  
}

RT_correction_xcms<-function(datafiles=NULL){
  library(tcltk)
  library(stringr)
  library(xcms)
  library(RColorBrewer)
  if(is.null(datafiles)){datafiles=tk_choose.files(caption = "select the MS data files for RT adjustment")}
  matchtype=c("mzxml$","cdf$")
  #datafiles=gsub(paste0(dirname(datafiles)[1],"/"),"",datafiles)

  datafiles=datafiles[grep(pattern =paste(matchtype,collapse="|"),x=datafiles,ignore.case = T)]
  message("Reading raw data...")
  raw_data <- readMSData(datafiles, mode = "onDisk", msLevel = 1,centroided =T, verbose = F )
  message("Reading raw data finished")
  mfp <- MatchedFilterParam(snthresh = 30, binSize = .5)
  message("Finding features...")
  closeAllConnections()
  cl=makeCluster(12)
  res <- findChromPeaks(raw_data, param = mfp)
  #head(chromPeaks(res))
  #table(chromPeaks(res)[, "sample"])
  
  ## Performing the peak grouping using the "peak density" method.
  p <- PeakDensityParam(sampleGroups = rep(1, length(datafiles)),maxFeatures=1000,minFraction = 0.5,minSamples = 1,binSize = 0.25)
  res <- groupChromPeaks(res, param = p)
  
  foverlap=overlappingFeatures(res, expandRt = 60)
  ## Perform the retention time adjustment using peak groups found in both
  ## files.minFraction equal to the fraction required as being a principle compound
  fgp <- xcms::PeakGroupsParam(minFraction = 0.95)
  ## Before running the alignment we can evaluate which features (peak groups)
  ## would be used based on the specified parameters.
  pkGrps <- xcms::adjustRtimePeakGroups(res, param = fgp)
  #pkGrps <- xcms::adjustRtime(res, param = fgp)
  ## We can also plot these to evaluate if the peak groups span a large portion
  ## of the retention time range.
  pkGrps<-pkGrps/60
  pkGrps<-data.frame(pkGrps,stringsAsFactors = F)
  pkGrpscolor=1:nrow(pkGrps)

  
  
  #axis(side = 2, at = c(1:ncol(pkGrps)), labels = colnames(pkGrps),cex.lab=0.5)
  
  ## Next we perform the alignment.
  res <- adjustRtime(res, param = fgp)
  
  ## Any grouping information was dropped
  #hasFeatures(res)
  
  ## Adjusterd retention times can be accessed using
  ## rtime(object, adjusted = TRUE) and adjustedRtime
  #all.equal(rtime(res), adjustedRtime(res))
  
  ## To get the raw, unadjusted retention times:
  #all.equal(rtime(res, adjusted = FALSE), rtime(raw_data))
  
  ## To extract the retention times grouped by sample/file:
  rts <- rtime(res, bySample = T)
  rtsraw <- rtime(raw_data)
  

   sampleftlength=sapply(rts,length) 
   rtsrawlist<-list()
   for (subset in 1: length(sampleftlength)){
     if (subset!=1){
      start=sum(sampleftlength[1:subset-1]) 
     }else{
      start=0
     }
     
     rtsrawlist[[as.character(subset)]]=rtsraw[(1+start):(sampleftlength[subset]+start)]
     
   }
  finalrts_sum<-function(rts){  
   finalrts<-data.frame()
  for (list in names(rts)){
    
    a<-rts[[list]]
    fid=names(a)
    fid=t(data.frame(str_split(fid,"\\.")))
    filenum=as.numeric(gsub("F","",fid[1,1]))
  a<-as.data.frame(rts[[list]]) 
  a$fID= fid[,2]
  colnames(a)=c(as.character(res@phenoData@data[["sampleNames"]][as.numeric(filenum)]),"fID")
  if (length(finalrts)==0){
  finalrts=a
  }else{
  finalrts=merge(finalrts,a,by="fID",all=TRUE) 
  }
  } 
   finalrts
     }
  finalrts=finalrts_sum(rts)
  finalrtsraw=finalrts_sum(rtsrawlist)
  pkGrpsadj=finalrts[(gsub("FT","",rownames(pkGrps))),]
  pkGrpsadj=pkGrpsadj[,-1]
  pkGrpsraw=finalrtsraw[(gsub("FT","",rownames(pkGrps))),]
  pkGrpsraw=pkGrpsraw[,-1]
  
  png("RTadjustment_QC.png",width = 15,height = 15,units = "in",res = 150)
  par(mfrow = c(2,2))
  
  plot(x = pkGrps[, 1], y = rep(1, nrow(pkGrps)), xlim =range(pkGrps,na.rm=T), col= pkGrpscolor,ylim=c(1,ncol(pkGrps)) ,xlab = "rt", ylab = "Data Files", yaxt = "n", main= "Representative feature RT")
  for (col in 2:ncol(pkGrps)){
    points(x = pkGrps[, col], y = rep(col, nrow(pkGrps)),col= pkGrpscolor)  
    segments(x0 = pkGrps[, col-1], x1 = pkGrps[, col],
             y0 = rep(col-1, nrow(pkGrps)), y1 = rep(col, nrow(pkGrps)),col= pkGrpscolor)  
  }
  grid()
  
  
  
  
  ## Plot the raw against the adjusted retention times.
  plot((rtime(raw_data)/60), (rtime(res)/60), pch = 16, cex = 0.25, col = fromFile(res),xlab = "Raw", ylab = "Adjusted",  main= "RT comparision before and after")
  
   plot(x = pkGrpsraw[, 1], y = rep(1, nrow(pkGrpsraw)), xlim =range(pkGrpsraw,na.rm=T), col= pkGrpscolor,ylim=c(1,ncol(pkGrpsraw)) ,xlab = "rt", ylab = "Data Files", yaxt = "n", main= "Before RT adjustment")
  for (col in 2:ncol(pkGrpsraw)){
    points(x = pkGrpsraw[, col], y = rep(col, nrow(pkGrpsraw)),col= pkGrpscolor)  
    segments(x0 = pkGrpsraw[, col-1], x1 = pkGrpsraw[, col],
             y0 = rep(col-1, nrow(pkGrpsraw)), y1 = rep(col, nrow(pkGrpsraw)),col= pkGrpscolor)  
  }
  grid()
  plot(x = pkGrpsadj[, 1], y = rep(1, nrow(pkGrpsadj)), xlim =range(pkGrpsadj,na.rm=T), col= pkGrpscolor,ylim=c(1,ncol(pkGrpsadj)) ,xlab = "rt", ylab = "Data Files", yaxt = "n", main= "After RT adjustment")
  for (col in 2:ncol(pkGrpsadj)){
    points(x = pkGrpsadj[, col], y = rep(col, nrow(pkGrpsadj)),col= pkGrpscolor)  
    segments(x0 = pkGrpsadj[, col-1], x1 = pkGrpsadj[, col],
             y0 = rep(col-1, nrow(pkGrpsadj)), y1 = rep(col, nrow(pkGrpsadj)),col= pkGrpscolor)  
  }
  grid()
  
 
  
  mtext("RTadjustment Quality control", side = 3, line = -2, outer = TRUE)
  dev.off()
  write.csv(finalrts,"finalRTadjustment_matrix.csv")
  return(list(finalrtsraw,finalrts))
}

Par_peakintegrate<-function(h,library_file,findScanTime,raw_data,intensity_type,Ion.bin){
  
  scanrangeL <- min(findScanTime@scantime)
  scanrangeU <- max(findScanTime@scantime)
  mzacqrangeL <- min(findScanTime@mzrange)
  mzacqrangeU <- max(findScanTime@mzrange)
  
  metabolite <-  library_file[h, ]
  R.Ion<-metabolite$Ref.ion
  RT.lowerset<-metabolite$RT.shfit.Lower
  RT.upperset<-metabolite$RT.shfit.upper
  R.Time<-metabolite$RT.median*60
  
  EIClower=R.Time-RT.lowerset
  EICupper=R.Time+RT.upperset
  EICmzlower=R.Ion-Ion.bin
  EICmzupper=R.Ion+Ion.bin
  if((R.Time-RT.lowerset)<=scanrangeL) EIClower=scanrangeL
  if((R.Time+RT.upperset)>=scanrangeU) EICupper=scanrangeU
  if((R.Ion-Ion.bin)<=mzacqrangeL) EICmzlower=mzacqrangeL
  if((R.Ion+Ion.bin)>=mzacqrangeU) EICmzupper=mzacqrangeU
  
  if (sum(EICupper<scanrangeL,EIClower>scanrangeU,(EICmzlower)>mzacqrangeU,(EICmzupper)<mzacqrangeL,na.rm = T)>=1) {
    Peakvalue=0
    return(Peakvalue)
  }else{
    
    IonExtract <-xcms::getEIC(raw_data,mzrange=cbind(EICmzlower,EICmzupper),rtrange=cbind(EIClower,EICupper))
    abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
    ifelse(intensity_type=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-flux::auc(abundance$rt, abundance$intensity*10))
    return(Peakvalue)
  }
  
}

Par_peakintegrate_slow<-function(h,library_file,findScanTime,raw_data,intensity_type,Ion.bin,name.file){
  
  
  scanrangeL <- min(findScanTime@scantime)
  scanrangeU <- max(findScanTime@scantime)
  mzacqrangeL <- min(findScanTime@mzrange)
  mzacqrangeU <- max(findScanTime@mzrange)
  
  metabolite <-  library_file[h, ]
  R.Ion<-metabolite$Ref.ion
  RT.lowerset<-metabolite$RT.shfit.Lower
  RT.upperset<-metabolite$RT.shfit.upper
  R.Time<-metabolite$RT.median*60
  
  EIClower=R.Time-RT.lowerset
  EICupper=R.Time+RT.upperset
  EICmzlower=R.Ion-Ion.bin
  EICmzupper=R.Ion+Ion.bin
  if((R.Time-RT.lowerset)<=scanrangeL) EIClower=scanrangeL
  if((R.Time+RT.upperset)>=scanrangeU) EICupper=scanrangeU
  if((R.Ion-Ion.bin)<=mzacqrangeL) EICmzlower=mzacqrangeL
  if((R.Ion+Ion.bin)>=mzacqrangeU) EICmzupper=mzacqrangeU
  
  if (sum(EICupper<scanrangeL,EIClower>scanrangeU,(EICmzlower)>mzacqrangeU,(EICmzupper)<mzacqrangeL,na.rm = T)>=1) {
    Peakvalue=0
    abundance=data.frame(rt=NA,intensity=NA)
  }else{
    
    IonExtract <-xcms::getEIC(raw_data,mzrange=cbind(EICmzlower,EICmzupper),rtrange=cbind(EIClower,EICupper))
    abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
    ifelse(intensity_type=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-flux::auc(abundance$rt, abundance$intensity*10))
    
  }
  
  result = list(Peakvalue=Peakvalue,Graphic.df=data.frame(Metabolite.Names=paste(metabolite$Name," (m/z:",R.Ion,")",sep="") ,Retention.Time=(abundance$rt)/60,
                                                          Intensity=abundance$intensity ,Sample.Names=name.file ))
  return(result)
}

testpar<-function(){
  for( h in 1: nrow(library_file)){
  
  scanrangeL <- min(findScanTime@scantime)
  scanrangeU <- max(findScanTime@scantime)
  mzacqrangeL <- min(findScanTime@mzrange)
  mzacqrangeU <- max(findScanTime@mzrange)
  
  metabolite <-  library_file[h, ]
  R.Ion<-metabolite$Ref.ion
  RT.lowerset<-metabolite$RT.shfit.Lower
  RT.upperset<-metabolite$RT.shfit.upper
  R.Time<-metabolite$RT.median*60
  
  EIClower=R.Time-RT.lowerset
  EICupper=R.Time+RT.upperset
  EICmzlower=R.Ion-Ion.bin
  EICmzupper=R.Ion+Ion.bin
  if((R.Time-RT.lowerset)<=scanrangeL) EIClower=scanrangeL
  if((R.Time+RT.upperset)>=scanrangeU) EICupper=scanrangeU
  if((R.Ion-Ion.bin)<=mzacqrangeL) EICmzlower=mzacqrangeL
  if((R.Ion+Ion.bin)>=mzacqrangeU) EICmzupper=mzacqrangeU
  
  if (sum(EICupper<scanrangeL,EIClower>scanrangeU,(EICmzlower)>mzacqrangeU,(EICmzupper)<mzacqrangeL,na.rm = T)>=1) {
    Peakvalue=0
    
  }else{
    
    IonExtract <-xcms::getEIC(raw_data,mzrange=cbind(EICmzlower,EICmzupper),rtrange=cbind(EIClower,EICupper))
    abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
    ifelse(intensity_type=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-flux::auc(abundance$rt, abundance$intensity*10))
    
  }
  
  result = list(Peakvalue=Peakvalue,Graphic.df=data.frame(Metabolite.Names=paste(metabolite$Name," (m/z:",R.Ion,")",sep="") ,Retention.Time=(abundance$rt)/60,
                                                          Intensity=abundance$intensity ,Sample.Names=name.file ))
  
}
  
}

adjustRtimePeakGroups <- function(object, param = PeakGroupsParam(),
                                  msLevel = 1L) {
  if (!is(object, "XCMSnExp"))
    stop("'object' has to be an 'XCMSnExp' object.")
  if (!hasFeatures(object))
    stop("No features present. Please run 'groupChromPeaks' first.")
  if (hasAdjustedRtime(object))
    warning("Alignment/retention time correction was already performed, ",
            "returning a matrix with adjusted retention times.")
  subs <- subset(param)
  if (!length(subs))
    subs <- seq_along(fileNames(object))
  nSamples <- length(subs)
  pkGrp <- .getPeakGroupsRtMatrix(
    peaks = chromPeaks(object, msLevel = msLevel),
    peakIndex = .peakIndex(
      .update_feature_definitions(featureDefinitions(object),
                                  rownames(chromPeaks(object)),
                                  rownames(chromPeaks(object,
                                                      msLevel = msLevel)))),
    sampleIndex = subs,
    missingSample = nSamples - (nSamples * minFraction(param)),
    extraPeaks = extraPeaks(param)
  )
  colnames(pkGrp) <- basename(fileNames(object))[subs]
  pkGrp
}

.getPeakGroupsRtMatrix <- function(peaks, peakIndex, sampleIndex,
                                   missingSample, extraPeaks) {
  ## For each feature:
  ## o extract the retention time of the peak with the highest intensity.
  ## o skip peak groups if they are not assigned a peak in at least a
  ##   minimum number of samples OR if have too many peaks from the same
  ##   sample assigned to it.
  nSamples <- length(sampleIndex)
  rt <- lapply(peakIndex, function(z) {
    cur_fts <- peaks[z, c("rt", "into", "sample"), drop = FALSE]
    ## Return NULL if we've got less samples that required or is the total
    ## number of peaks is larger than a certain threshold.
    ## Note that the original implementation is not completely correct!
    ## nsamp > nsamp + extraPeaks might be correct.
    nsamp <- length(unique(cur_fts[, "sample"]))
    if (nsamp < (nSamples - missingSample) |
        nrow(cur_fts) > (nsamp + extraPeaks))
      return(NULL)
    cur_fts[] <- cur_fts[order(cur_fts[, 2], decreasing = TRUE), ]
    cur_fts[match(sampleIndex, cur_fts[, 3]), 1]
  })
  rt <- do.call(rbind, (rt))
  ## Order them by median retention time. NOTE: this is different from the
  ## original code, in which the peak groups are ordered by the median
  ## retention time that is calculated over ALL peaks within the peak
  ## group, not only to one peak selected for each sample (for multi
  ## peak per sample assignments).
  ## Fix for issue #175
  if (is(rt, "matrix")) {
    rt <- rt[order(rowMedians(rt, na.rm = TRUE)), , drop = FALSE]
  }
  rt
}
.peakIndex <- function(object) {
  if (inherits(object, "DataFrame")) {
    idxs <- object$peakidx
    names(idxs) <- rownames(object)
  } else {
    if (!hasFeatures(object))
      stop("No feature definitions present. Please run groupChromPeaks first.")
    idxs <- featureDefinitions(object)$peakidx
    names(idxs) <- rownames(featureDefinitions(object))
  }
  idxs
}
.update_feature_definitions <- function(x, original_names, subset_names) {
  x$peakidx <- lapply(x$peakidx, function(z) {
    idx <- base::match(original_names[z], subset_names)
    idx[!is.na(idx)]
  })
  x[lengths(x$peakidx) > 0, ]
}


erah_pipeline<-function(workdir=getwd(),infopath=NULL){
  
  library(erah)
  if (is.null(infopath)){tk_choose.files(caption = "Select info file to normalized by SERRF")}
    #check if it is csv of xlsx
  if(grepl(".xls.?$", infopath)){
    pData <- openxlsx::read.xlsx(infopath, sheet = 1,colNames = T)
  }else if(grepl(".csv", infopath)){
    # file = "C:\\Users\\Sili Fan\\Downloads\\val (18).csv"
    if(nrow(d)<2){pData <- read.csv(infopath)}
    pData <- data.table::fread(infopath)
  }
    
  
  
}

xcmsRaw_dev<-function (filename, profstep = 1, profmethod = "bin", profparam = list(), 
          includeMSn = FALSE, mslevel = NULL, scanrange = NULL) 
{
  object <- new("xcmsRaw")
  object@env <- new.env(parent = .GlobalEnv)
  object@filepath <- xcmsSource(filename)
  rawdata <- loadRaw(object@filepath, includeMSn = includeMSn)
  rtdiff <- diff(rawdata$rt)
  if (any(rtdiff == 0)) 
    warning("There are identical scantimes.")
  if (any(rtdiff < 0)) {
    badtimes <- which(rtdiff < 0)
    stop(paste("Time for scan ", badtimes[1], " (", 
               rawdata$rt[[badtimes[1]]], ") greater than scan ", 
               badtimes[1] + 1, " (", rawdata$rt[[badtimes[1] + 
                                                    1]], ")", sep = ""))
  }
  object@scantime <- rawdata$rt
  object@tic <- rawdata$tic
  object@scanindex <- rawdata$scanindex
  object@env$mz <- rawdata$mz
  object@env$intensity <- rawdata$intensity
  if (length(scanrange) < 2) {
    scanrange <- c(1, length(object@scantime))
  }else {
    scanrange <- range(scanrange)
  }
  if (min(scanrange) < 1 | max(scanrange) > length(object@scantime)) {
    scanrange[1] <- max(1, scanrange[1])
    scanrange[2] <- min(length(object@scantime), scanrange[2])
    message("Provided scanrange was adjusted to ", 
            scanrange[1], " - ", scanrange[2])
  }
  if (!is.null(rawdata$acquisitionNum)) {
    object@acquisitionNum <- rawdata$acquisitionNum
  }
  if (!is.null(rawdata$polarity)) {
    object@polarity <- factor(rawdata$polarity, levels = c(0, 
                                                           1, -1), labels = c("negative", "positive", 
                                                                              "unknown"))
  }
  if (!is.null(rawdata$MSn)) {
    object@env$msnMz <- rawdata$MSn$mz
    object@env$msnIntensity <- rawdata$MSn$intensity
    object@msnScanindex <- rawdata$MSn$scanindex
    object@msnAcquisitionNum <- rawdata$MSn$acquisitionNum
    object@msnLevel <- rawdata$MSn$msLevel
    object@msnRt <- rawdata$MSn$rt
    object@msnPrecursorScan <- match(rawdata$MSn$precursorNum, 
                                     object@acquisitionNum)
    object@msnPrecursorMz <- rawdata$MSn$precursorMZ
    object@msnPrecursorIntensity <- rawdata$MSn$precursorIntensity
    object@msnPrecursorCharge <- rawdata$MSn$precursorCharge
    object@msnCollisionEnergy <- rawdata$MSn$collisionEnergy
  }
  xcms:::scanrange(object) <- as.numeric(scanrange)
  object <- object[scanrange[1]:scanrange[2]]
  xcms:::mslevel(object) <- as.numeric(mslevel)
  object@mzrange <- range(object@env$mz, na.rm = TRUE)
  object@profmethod <- profmethod
  object@profparam <- profparam
  if (profstep) 
    xcms:::profStep(object) <- profstep
  if (!missing(mslevel) & !is.null(mslevel)) {
    if (max(mslevel) > 1) {
      object <- msn2ms(object)
      object <- split(object, f = object@msnLevel == mslevel)$"TRUE"
    }
  }
  return(object)
}
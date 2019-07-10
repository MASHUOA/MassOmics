#' \tabular{ll}{
#' Package: \tab MassOmics \cr
#' License: \tab GPL (>= 2)\cr
#' URL: \tab 
#' }
#'
#' @docType package
#' @name MassOmics
#' @author George GUO \email{George.GUO@@auckland.ac.nz}
#' @references 
#' @keywords package
#' @import tcltk
#' @import tcltk2
#' @import Rcpp
#' @import openxlsx
#' @import tidyr
#' @import randomForest
#' @import rfviz
#' @import ggpubr
#' @import extrafont
#' @import mzR
#' @import xcms
#' @import svDialogs
#' @import svGUI
#' @import KEGGREST
#' @import RColorBrewer
#' @import rgl
#' @import mixOmics
#' @import plyr
#' @import flux
#' @import tkrplot
#' @import multtest
#' @import XML
#' @import CAMERA
#' @import qvalue
#' @import doParallel
#' @import BiocManager
#' @import pacman
#' @import backports
#' @import BiocParallel
#' @import pbapply
#' @importFrom stats na.omit runif
#' @importFrom utils download.file modifyList packageVersion read.table tail



######################################### clean.fix #########################################
clean.fix <- function (main.folder = tk_choose.dir(caption = "Select working directory"),
          amdis.report = read.csv(tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                  multi = FALSE, caption = "Select the AMDIS report in .TXT"),
                                  sep = "\t", colClass = "character"), ion.lib = read.csv(tk_choose.files(default = "Select the Idenditification report in .CSV",
                                                                                                          multi = FALSE, caption = "Select the Idenditification report in .CSV"),
                                                                                          sep = ","), save = TRUE, output = "Fianl Metab result")
{


  require(xcms)
  old.wd <- getwd()
  main.folder <- main.folder
  setwd(main.folder)
  names(ion.lib)[1]<-"Name"


  files <- c(dir())
  info <- file.info(files)
  isDir <- info$isdir
  conditions <- c(files[isDir == TRUE])
  largest <- 1

  for (i in 1:length(conditions)) {
    largest <- c(largest, length(grep(".CDF", c(list.files(conditions[i],
                                                           full.names = TRUE)), ignore.case = TRUE)))
  }
  largest <- max(largest)
  num.rep <- length(grep(".CDF", c(list.files(conditions[1],
                                              full.names = TRUE)), ignore.case = TRUE))
  samples <- grep(".CDF", c(list.files(conditions[1], full.names = TRUE)),
                  ignore.case = TRUE, value = TRUE)
  replicates <- matrix(c(1:length(samples)), nrow = 1, ncol = largest)
  if (length(conditions) > 1) {
    for (i in 2:length(conditions)) {
      num.rep <- c(num.rep, length(grep(".CDF", c(list.files(conditions[i],
                                                             full.names = TRUE)), ignore.case = TRUE)))
      samples <- c(samples, grep(".CDF", c(list.files(conditions[i],
                                                      full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      up.num <- length(grep(".CDF", c(list.files(conditions[i],
                                                 full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      if (i < 3) {
        low.num <- (replicates[length(replicates)] +
                      1)
      }
      else {
        low.num <- max(replicates[i - 1, !duplicated(replicates[i -
                                                                  1, ])]) + 1
      }
      replicates <- rbind(replicates, c((low.num):(low.num +
                                                     (up.num - 1))))
    }
  }
  if (is.data.frame(ion.lib) == FALSE) {
    ion.lib = read.csv(ion.lib, sep = ",")
  }
  library_file <- ion.lib
  library_file$Name <- gsub(" ", "", library_file$Name)
  if (is.data.frame(amdis.report) == TRUE) {
    amdis.report <- amdis.report
  }
  else {
    amdis.report = read.csv(amdis.report, sep = "\t")
  }
  for (i in 1:nrow(amdis.report)) {
    original <- as.character(amdis.report[i, 1])
    new <- unlist(strsplit(original, "\\", fixed = TRUE))
    new <- new[length(new)]
    amdis.report[i, 1] <- new
  }
  total_report <- amdis.report
  final.df <- 1
  for (c in 1:length(conditions)) {
    filenames <- dir(conditions[c])
    for (q in 1:length(filenames)) {
      file <- filenames[q]
      name.file <- file
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      file <- gsub("CDF", "FIN", file, ignore.case=TRUE)
      file <- which(total_report$FileName == file)
      if (length(file) > 0) {
        data.df <- total_report[file, ]
        Keep <- c("Name", "RT", "Scan", "Expec..RT",
                  "Base.Peak")
        data.df <- data.df[, Keep]
        data.df[2] <- as.numeric(data.df[, 2])
        data.df[3] <- as.numeric(data.df[, 3])
        data.df[4] <- as.numeric(data.df[, 4])
        data.df$deff <- data.df$RT - data.df$Expec..RT
        data.df <- subset(data.df, data.df$deff < 2.5 &
                            data.df$deff > -2.5)
        sure <- data.df[1, ]
        while (nrow(data.df) >= 1) {
          line <- data.df[1, ]
          sameRT <- subset(data.df, data.df$RT == line$RT)
          if (nrow(sameRT) > 1) {
            isoleucine <- grep("isoleucine", sameRT$Name,
                               ignore.case = TRUE)
            if (length(isoleucine) > 0) {
              isoleucine2 <- grep("isoleucine", sure$Name,
                                  ignore.case = TRUE)
              if (length(isoleucine2) > 0) {
                sameRT2 <- sameRT[-isoleucine, ]
                sure <- rbind(sure, sameRT2[which.min(which.min(abs(sameRT2$deff))),
                                            ])
                sure <- sure[!duplicated(sure$RT), ]
                data.df <- subset(data.df, row.names(data.df) !=
                                    row.names(sameRT))
              }
              else {
                sure <- rbind(sure, sameRT[isoleucine,
                                           ])
                sure <- sure[!duplicated(sure$RT), ]
                data.df <- subset(data.df, row.names(data.df) !=
                                    row.names(sameRT))
              }
            }
            else {
              sure <- rbind(sure, sameRT[which.min(which.min(abs(sameRT$deff))),
                                         ])
              sure <- sure[!duplicated(sure$RT), ]
              data.df <- subset(data.df, row.names(data.df) !=
                                  row.names(sameRT))
            }
          }
          else {
            sure <- rbind(sure, sameRT)
            sure <- sure[!duplicated(sure$RT), ]
            data.df <- subset(data.df, row.names(data.df) !=
                                row.names(sameRT))
          }
        }
        surefinal <- sure
        surefinal$Name <- gsub("? ", "", surefinal$Name,
                               fixed = TRUE)
        surefinal$Name <- gsub("?? ", "", surefinal$Name,
                               fixed = TRUE)
        surefinal$Name <- gsub("??? ", "", surefinal$Name,
                               fixed = TRUE)
        surefinal$Name <- gsub("?", "", surefinal$Name,
                               fixed = TRUE)
        surefinal$Name <- gsub("??", "", surefinal$Name,
                               fixed = TRUE)
        surefinal$Name <- gsub("???", "", surefinal$Name,
                               fixed = TRUE)
        raw_data <- xcmsRaw(filename = paste(conditions[c],
                                             filenames[q], sep = "\\"))
        dev.new.OS()
        plotChrom(raw_data)
        for (h in 1:nrow(surefinal)) {
          metabolite <- surefinal[h, ]
          amdis_scan <- (surefinal$Scan[h])
          ion.ref <- grep(gsub(" ", "", metabolite$Name),
                          library_file$Name, fixed = TRUE)
          if (length(ion.ref) == 0) {
            error.lib <- paste("There is no ion defined for metabolite",
                               metabolite$Name, "in the reference ion library",
                               sep = " ")
            print(error.lib)
          }
          metabolite <- library_file[ion.ref, ]
          if (is.na(metabolite$ref_ion)) {
            error.lib2 <- paste("The metabolite", metabolite$Name,
                                "was detected in the reference ion library, however, there is no ion fragment defined for that.",
                                sep = " ")
            print(error.lib2)
          }
          scan <- getScan(raw_data, scan = surefinal$Scan[h])
          scan <- data.frame(scan)
          scan$mz <- trunc(scan$mz)
          scan_ion <- (metabolite$ref_ion) - 2
          for (k in 1) {
            abundance <- subset(scan, scan$mz == (scan_ion +
                                                    k))
            maxabundance <- max(abundance$intensity)
            ion_set <- c(maxabundance)
          }
          for (k in 2:3) {
            abundance <- subset(scan, scan$mz == (scan_ion +
                                                    k))
            maxabundance <- max(abundance$intensity)
            ion_set <- c(ion_set, maxabundance)
          }
          correct_abundance <- max(ion_set)
          surefinal$Base.Peak[h] <- max(correct_abundance)
        }
        surefinal$RT <- NULL
        surefinal$Expec..RT <- NULL
        surefinal$deff <- NULL
        surefinal$Scan <- NULL
        surefinal <- surefinal[!duplicated(surefinal$Name),
                               ]
        names(surefinal)[2] <- name.file
        if (final.df == 1) {
          final.df <- surefinal
          confirmation <- paste("File", q, "(", name.file,
                                ")", "done!", sep = " ")
          print(confirmation)
        }
        else {
          final.df <- merge(final.df, surefinal, by.x = "Name",
                            by.y = "Name", all = TRUE)
          confirmation <- paste("File", q, "(", name.file,
                                ")", "done!", sep = " ")
          print(confirmation)
        }
      }
      else {
        confirmation <- paste("File", q, "(", name.file,
                              ")", "not analyzed! Probably there is no metabolites detected for this sample in the AMDIS report",
                              sep = " ")
        print(confirmation)
        if (final.df != 1) {
          final.df$error <- "Error"
          names(final.df)[ncol(final.df)] <- name.file
        }
      }
    }
  }
  for (r in 1:nrow(replicates)) {
    rep.row <- replicates[r, ]
    rep.row <- rep.row[!duplicated(rep.row)]
    rep.name <- rep(conditions[r], length(rep.row))
    if (r == 1) {
      rep.name.final <- c(rep.name)
      rep.name.final <- c("Replicates", rep.name.final)
    }
    else {
      rep.name.final <- c(rep.name.final, rep.name)
    }
  }
  final.df <- rbind(rep.name.final, final.df)
  final.df <- final.df[order(final.df[, 2], decreasing = T),
                       ]
  row.names(final.df) <- 1:nrow(final.df)
  if (save == TRUE) {
    sheet <- output
    store <- paste(main.folder, "\\", sheet, ".csv", sep = "")
    write.csv(final.df, file = store, row.names = FALSE)
  }
  setwd(old.wd)
  return(final.df)
}
######################################### clean.fix3 #########################################
clean.fix3 <- function (main.folder = tk_choose.dir(caption = "Select working directory"),
          amdis.report = read.csv(tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                  multi = FALSE, caption = "Select the AMDIS report in .TXT"),
                                  sep = "\t", stringsAsFactors=FALSE),
          save = TRUE,
          output = "Final Metab3 Results",
          Ion.bin= 0.5,
          RT.bin = 5,
          File.name = "Identification Report.csv",
          MS.L= tk_choose.files(caption="Select MS library (e.g. SVB_Totoal) in .msl")) {
  require(xcms)
  require(lattice)
  require(plyr)

  main.folder <- main.folder
  setwd(main.folder)

  ################# Step 1. Function transforms AMDIS library to spreadsheet of data for each metabolite

  getIonLib<-function(lib.fn = MS.L){
    lib.txt<-readLines(lib.fn)
    lib.txt<-lib.txt[-grep("^[ ]*$",lib.txt)]
    att.nms<-unique(sapply(strsplit(lib.txt[grep(":",lib.txt)],":"),function(x) x[1]))
    entry.nms<-sapply(strsplit(lib.txt[grep("NAME:",lib.txt)],":"), function(x) x[2])

    rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))

    starts=grep("NAME:", lib.txt)
    stops=c(starts[1:(length(starts)-1)]+diff(starts)-1,length(lib.txt))
    for (i in 1:length(starts)){
      #i=1
      tmp<-lib.txt[starts[i]:stops[i]]
      sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
    }
    rez<-data.frame(rez,stringsAsFactors=FALSE)
    rez=within(rez, {
      RT=as.numeric(RT)
      RSN=as.numeric(RSN)
      NUM.PEAKS=as.numeric(NUM.PEAKS)
    })
    return(rez)

  }
  libr<-getIonLib()## returns list of peaks


  ############ Step 2. Generate the list of metabolites and clean up
  AmRep<-amdis.report
  AmRep<-AmRep[!AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)+1.5 & AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)-1.5,]
  AmRep$Name <- gsub("?", "", AmRep$Name, fixed = TRUE)
  AmRep$Name <- gsub("^ ", "", AmRep$Name, perl=T)

  ## RT statistics from the AReport
  RT.stats<-t(sapply(split(AmRep$RT, AmRep$Name),function(x) c(RT.mean=round(mean(x,na.rm=TRUE),3),
                                                               RT.median=round(median(x,na.rm=T),3),
                                                               RT.sd=sd(x,na.rm=T),
                                                               Number.of.ID=length(x[!is.na(x)]))))

  AmRep.RT.stats <- as.data.frame(cbind(RT.stats,Exp.RT=0, ref_ion=0,Mol.mass=0, Diff.RT=0, CAS=0))##data frame for list of metabolites
  AmRep.RT.stats$ref_ion <- libr$RSN[match(rownames(RT.stats),libr$NAME)]
  AmRep.RT.stats$Mol.mass <- libr$FORM[match(rownames(RT.stats),libr$NAME)]
  AmRep.RT.stats$Exp.RT <- libr$RT[match(rownames(RT.stats),libr$NAME)]
  AmRep.RT.stats$Diff.RT <- AmRep.RT.stats$Exp.RT- AmRep.RT.stats$RT.median
  AmRep.RT.stats$CAS <- libr$CAS[match(rownames(RT.stats),libr$NAME)]

  AmRep.RT.stats <- AmRep.RT.stats[order(AmRep.RT.stats$RT.median, decreasing = F),]
  write.csv(AmRep.RT.stats, file = File.name)

  check <- tk_select.list(c("Good Identification Report (continue)",  "Bad Identificaiton Reprot (stop)","Load a corrected report"),  title = "Check Identification Report")


  if (check == "Bad Identificaiton Reprot (stop)"){
    stop("Data process terminated")

  }else if(check== "Load a corrected report" ){

    final.check.data = read.csv(tk_choose.files(caption = "Load a corrected report",  multi = FALSE),stringsAsFactors=FALSE)
  }else{

    final.check.data<-read.csv(File.name,stringsAsFactors=FALSE)
  }


  GGReport <- tk_select.list(c("Yes (Slower analysis speed)",  "No (Faster - but cann't perform the Diagnosis of Auto GC-Peak Integration )"),  title = "Turn on the Diagnosis function")



  names(final.check.data)[1]<-"Name"
  ion.lib<-final.check.data
  files <- c(dir())
  info <- file.info(files)
  isDir <- info$isdir
  conditions <- c(files[isDir == TRUE])
  largest <- 1

  for (i in 1:length(conditions)) {
    largest <- c(largest, length(grep(".CDF", c(list.files(conditions[i],
                                                           full.names = TRUE)), ignore.case = TRUE)))
  }
  largest <- max(largest)
  num.rep <- length(grep(".CDF", c(list.files(conditions[1],
                                              full.names = TRUE)), ignore.case = TRUE))

  samples <- grep(".CDF", c(list.files(conditions[1], full.names = TRUE)), ignore.case = TRUE, value = TRUE)

  replicates <- matrix(c(1:length(samples)), nrow = 1, ncol = largest)
  if (length(conditions) > 1) {
    for (i in 2:length(conditions)) {
      num.rep <- c(num.rep, length(grep(".CDF", c(list.files(conditions[i],
                                                             full.names = TRUE)), ignore.case = TRUE)))
      samples <- c(samples, grep(".CDF", c(list.files(conditions[i],
                                                      full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      up.num <- length(grep(".CDF", c(list.files(conditions[i],
                                                 full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      if (i < 3) {
        low.num <- (replicates[length(replicates)] +
                      1)
      }
      else {
        low.num <- max(replicates[i - 1, !duplicated(replicates[i -
                                                                  1, ])]) + 1
      }
      replicates <- rbind(replicates, c((low.num):(low.num +
                                                     (up.num - 1))))
    }
  }

  library_file <- ion.lib[,c("Name","ref_ion","RT.median")]
  final.df<-NULL
  remove(final.df)
  Graphic.df<-data.frame()

  if(length(conditions)==1) {
    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]
    for (q in 1:length(filenames)) {
      file1 <- filenames[q]
      name.file <- file1
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      surefinal <- data.frame(Name=ion.lib[,"Name"],Base.Peak=0)
      raw_data <- xcmsRaw(filename = paste(conditions,
                                           filenames[q], sep = "\\"))
      if (GGReport=="Yes (Slower analysis speed)"){
        dev.new.OS()
        plotChrom(raw_data) # plot TIC
      }
      for (h in 1:nrow(library_file)) {
        metabolite <-  library_file[h, ]
        R.Ion<-metabolite$ref_ion
        R.Time<-metabolite$RT.median*60
        IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(R.Time-RT.bin,R.Time+RT.bin))
        abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
        maxabundance <- max(abundance$intensity)
        surefinal$Base.Peak[h] <- maxabundance
        if (GGReport=="Yes (Slower analysis speed)"){
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
    final.df<-cbind(CAS=ion.lib$CAS,  Num.ID=ion.lib$Number.of.ID, Ref.Ion=ion.lib$ref_ion, Ret.Time=ion.lib$RT.median, final.df)
    final.df <- final.df[order(final.df$Ret.Time , decreasing = F),]
    row.names(final.df) <- 1:nrow(final.df)


  } else{


    for (c in 1:length(conditions)) {
      filenames <- dir(conditions[c])
      filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)] # remove the effect of fin and leu
      for (q in 1:length(filenames)) {
        file1 <- filenames[q]
        name.file <- file1
        name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
        surefinal <- data.frame(Name=ion.lib[,"Name"],Base.Peak=0)
        raw_data <- xcmsRaw(filename = paste(conditions[c],
                                             filenames[q], sep = "\\"))
        if (GGReport=="Yes (Slower analysis speed)"){
          dev.new.OS()
          plotChrom(raw_data) # plot TIC
        }
        for (h in 1:nrow(library_file)) {
          metabolite <-  library_file[h, ]
          R.Ion<-metabolite$ref_ion
          R.Time<-metabolite$RT.median*60
          IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(R.Time-RT.bin,R.Time+RT.bin))
          abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
          maxabundance <- max(abundance$intensity)
          surefinal$Base.Peak[h] <- maxabundance
          if (GGReport=="Yes (Slower analysis speed)"){
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
    }

    for (r in 1:nrow(replicates)) {
      rep.row <- replicates[r, ]
      rep.row <- rep.row[!duplicated(rep.row)]
      rep.name <- rep(conditions[r], length(rep.row))
      if (r == 1) {
        rep.name.final <- c(rep.name)
        rep.name.final <- c(rep("Replicates",4), rep.name.final)
      }
      else {
        rep.name.final <- c(rep.name.final, rep.name)
      }
    }

    final.df$Name<-as.character(final.df$Name)
    final.df<-cbind(CAS=ion.lib$CAS,  Num.ID=ion.lib$Number.of.ID, Ref.Ion=ion.lib$ref_ion, Ret.Time=ion.lib$RT.median, final.df)
    final.df$CAS<-as.character(final.df$CAS)
    final.df <- final.df[order(final.df$Ret.Time, decreasing = F),]
    final.df <- rbind(rep.name.final, final.df)
    row.names(final.df) <- 1:nrow(final.df)
  }


  if (save == TRUE) {
    sheet <- output
    store <- paste(main.folder, "\\", sheet, ".csv", sep = "")
    write.csv(final.df, file = store, row.names = FALSE)

    if (GGReport=="Yes (Slower analysis speed)"){
      write.csv(Graphic.df,file="GC-Peak Diagnosis Report.csv")
      cat(paste("Data process is completed and called",sheet, sep=" "))
    }
  }

  return(final.df)

}


######################################### IonExtractor_Window #########################################
IonExtractor_Window <- function(){
  require(RColorBrewer)
  require(flux)
  IonExtractor<- function(folder = tk_choose.dir(caption = "Select working directory"),
                          mz = 106,
                          RT.Time1 = 10.5,
                          RT.Time2 = 11,
                          Ion.bin = 0.5,
                          heightMAX= 50000,
                          Metabolite.Names = "D4-Methanol",
                          save = TRUE,
                          peak="Peak Height"
  ){

    require(xcms)
    old.wd <- getwd()
    main.folder <- folder
    folderRE <<- folder
    setwd(main.folder)
    files <- c(dir())
    info <- file.info(files)
    isDir <- info$isdir
    conditions_pre <- c(files[isDir == TRUE])
    conditions <- tk_select.list(conditions_pre, multiple = FALSE,
                              title = paste("Folder conatins all cdf files"))
    condtions1<<-conditions
    final.df<-NULL
    plot_1<-NULL
    rm(final.df)

    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]
    colors<-rainbow(length(filenames))
    # dev.new.OS()
    for (q in 1:length(filenames)) {
      raw_data <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "\\"))
      name.file <- filenames[q]
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      IonExtract <-getEIC(raw_data,mzrange=cbind(mz-Ion.bin,mz+Ion.bin),rtrange=cbind(RT.Time1*60,RT.Time2*60))
      abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
      ifelse(peak=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-auc(abundance$rt, abundance$intensity*10))

      if (!exists("final.df")){
        final.df <- data.frame(t(c(name.file,Peakvalue)))
      } else {
        final.df <- rbind(final.df, data.frame(t(c(name.file,Peakvalue))))
      }
      print(paste("File", q, name.file,"done!"))

      ## Plotting

      Intensity=abundance$intensity
      Retention.Time<-abundance$rt/60
      #CombinRT <- as.data.frame(cbind(Retention.Time,Intensity))
      if(q==1){
        dev.new.OS()
        ifelse(heightMAX<max(Intensity),heightMAX<-max(Intensity),heightMAX)
        plot(Intensity~Retention.Time,col= colors[q],bty="n", type="l",main=paste(Metabolite.Names," (",mz," m/z",") ","RT:",RT.Time1,"-",RT.Time2," min", sep=""),
             xlab="Retention Time (min)",ylab="Intensity", ylim=c(0,heightMAX))
        ## test new scrip
        #plot_1 <- ggplot(CombinRT , aes(x=Retention.Time, y=Intensity)) + geom_line(col= colors[q]) + labs(title=paste(Metabolite.Names," (",mz," m/z",") ","RT:",RT.Time1,"-",RT.Time2," min", sep=""), x = "Retention Time (min)", y ="Intensity")
      }else{
        lines(Intensity~Retention.Time,col= colors[q])
        #plot_1 <- plot_1 +  geom_line( aes(x=Retention.Time, y=Intensity*q), col=colors[q])
        #print(plot_1)
      }
    }

    final.df2<-NULL
    final.df2<- data.frame(t(final.df[2]))
    names(final.df2) <- t(final.df[1])
    row.names(final.df2) <- Metabolite.Names


    if (save == TRUE) {
      ifelse(peak=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
      store <- paste(main.folder, "\\", Metabolite.Names,"(",mz,")RT",RT.Time1,"-",RT.Time2,"_",ValueType,".csv", sep = "")
      write.csv(final.df2, file = store, row.names = TRUE)
      print("Complete!")
      print(paste("The file ",Metabolite.Names, ".csv", " was saved in the folder ", folder, sep=""))
    } else {
      print("No file was saved because the argument save was set as FALSE")
    }
    return(final.df2)
  }



  ReIonExtractor<- function( folder = tk_choose.dir(caption = "Select working directory"),
                             mz = 106,
                             RT.Time1 = 10.5,
                             RT.Time2 = 10.7,
                             Ion.bin = 0.5,
                             heightMAX= 50000,
                             Metabolite.Names = "D4-alanine",
                             conditions=1,
                             peak="Peak Height",
                             save = TRUE){
    require(xcms)


    old.wd <- getwd()
    main.folder <- folder
    setwd(main.folder)


    if(conditions==1){
      files <- c(dir())
      info <- file.info(files)
      isDir <- info$isdir
      conditions_pre <- c(files[isDir == TRUE])
      conditions <- tk_select.list(conditions_pre, multiple = FALSE, title = paste("Folder conatins all cdf files"))
    }


    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]

    final.df<-NULL
    rm(final.df)

    
    colors<-rainbow(length(filenames))
    
    for (q in 1:length(filenames)) {
      raw_data <- xcmsRaw(filename = paste(conditions, filenames[q], sep = "\\"))
      name.file <- filenames[q]
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      IonExtract <-getEIC(raw_data,mzrange=cbind(mz-Ion.bin,mz+Ion.bin),rtrange=cbind(RT.Time1*60,RT.Time2*60))
      abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
      ifelse(peak=="Peak Height", Peakvalue <- max(abundance$intensity) ,Peakvalue<-auc(abundance$rt, abundance$intensity*10))
      if (!exists("final.df")){
        final.df <- data.frame(t(c(name.file, Peakvalue)))
      } else {
        final.df <- rbind(final.df, data.frame(t(c(name.file, Peakvalue))))
      }
      print(paste("File",q, name.file,"done!"))

      ## Plotting

      Intensity=abundance$intensity
      Retention.Time<-abundance$rt/60

      if(q==1){
        dev.new.OS()
        ifelse(heightMAX<max(Intensity),heightMAX<-max(Intensity),heightMAX)
        plot(Intensity~Retention.Time,col= colors[q],bty="n", type="l",main=paste(Metabolite.Names," (",mz," m/z",") ","RT:",RT.Time1,"-",RT.Time2," min", sep=""),
             xlab="Retention Time (min)",ylab="Intensity", ylim=c(0,heightMAX))
      }else{
        lines(Intensity~Retention.Time,col= colors[q])
      }
    }

    final.df2<-NULL
    final.df2<- data.frame(t(final.df[2]))
    names(final.df2) <- t(final.df[1])
    row.names(final.df2) <- Metabolite.Names


    if (save == TRUE) {

      ifelse(peak=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
      store <- paste(main.folder, "\\", Metabolite.Names,"(",mz,")RT",RT.Time1,"-",RT.Time2,"_",ValueType,".csv", sep = "")
      write.csv(final.df2, file = store, row.names = TRUE)
      print("Complete!")
      print(paste("The file ",Metabolite.Names, ".csv", " was saved in the folder ", folder, sep=""))
    } else {
      print("No file was saved because the argument save was set as FALSE")
    }
    return(final.df2)
  }


  ### interphase

  require(tcltk)
  tclRequire("BWidget")
  kk <- tktoplevel()
  frame_All <- tkframe(kk)
  tkwm.title(kk,"GC-MS IonExtractor")

  xvar <- tclVar("106")
  yvar <- tclVar("10.5")
  zvar <- tclVar("10.9")
  wvar <- tclVar("D4-alanine")
  hvar <- tclVar("300000")
  Peakvar <- tclVar("Peak Height")

  x.entry <- tkentry(background="white",kk, textvariable=xvar, width=17)
  y.entry <- tkentry(background="white",kk, textvariable=yvar, width=17)
  z.entry <- tkentry(background="white",kk, textvariable=zvar, width=17)
  w.entry <- tkentry(background="white",kk, textvariable=wvar, width=17)
  h.entry <- tkentry(background="white",kk, textvariable=hvar, width=17)


  ## choose between peakheight or area
  Peak <- c("Peak Height","Peak Area")
  comboBox <- tkwidget(kk,"ComboBox",editable=FALSE,values=Peak,textvariable=tclVar("Peak Height"),width=15)


  ####################


  Re_analysis <-function () {
    w <- as.character(tclvalue(wvar))
    x <- as.numeric(tclvalue(xvar))
    y <- as.numeric(tclvalue(yvar))
    z <- as.numeric(tclvalue(zvar))
    h <- as.numeric(tclvalue(hvar))
    Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    if (!exists("folderRE")){c(folder<- tk_choose.dir(caption = "Select working directory"),conditions<-1)
    }else{c(folder<-folderRE, conditions <- condtions1)}
    ReIonExtractor(Metabolite.Names=w,mz=x,RT.Time1=y,RT.Time2=z, heightMAX=h, folder = folder, conditions =conditions ,peak= Peakvar)
  }

  re_analysis.but <- tkbutton(kk, text="Re-analysis", command=Re_analysis)

  reset <- function() {
    tclvalue(xvar)<-""
    tclvalue(yvar)<-""
    tclvalue(zvar)<-""
    tclvalue(wvar)<-""
    tclvalue(hvar)<-"300000"
  }

  reset.but <- tkbutton(kk, text="Reset", command=reset)
  submit <- function() {
    w <- as.character(tclvalue(wvar))
    x <- as.numeric(tclvalue(xvar))
    y <- as.numeric(tclvalue(yvar))
    z <- as.numeric(tclvalue(zvar))
    h <- as.numeric(tclvalue(hvar))
    Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    IonExtractor(Metabolite.Names=w,mz=x,RT.Time1=y,RT.Time2=z, heightMAX=h ,peak= Peakvar)
  }
  submit.but <- tkbutton(kk, text="Submit", command=submit)

  quit.but <- tkbutton(kk, text = "Close Session",
                       command = function() {
                         tkdestroy(kk)
                       }
  )
  # outline
  frameL<- tkframe(frame_All,borderwidth=2, relief="groove",padx=5,pady=5)
  tkgrid(tklabel(frameL,text="Metabolite Name"), w.entry, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameL,text="Reference Ion"), x.entry, tklabel(frameL, text="m/z"), pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameL,text="Initial Retention Time"), y.entry, tklabel(frameL, text="minutes"), pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameL,text="End Retention Time"), z.entry, tklabel(frameL, text="minutes"), pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameL,text="Min Intensitiy for y-axis"), h.entry, pady= 10, padx= 10, sticky="w")
  tkgrid(tklabel(frameL,text="Integration Types"),comboBox, sticky="w",pady= 10, padx= 10)

  frameButton<- tkframe(frame_All)
  tkgrid(tklabel(frameButton,text=""), submit.but, re_analysis.but, reset.but, quit.but, pady= 10, padx= 10)

  ## load images
  image.path_GC_IonExtractor<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCMS_IonExtractor_logo.gif", sep=""))
  GC_IonExtractor_logo <-tcl("image",  "create", "photo", "IonExtract_GC.image", file=image.path_GC_IonExtractor)

  # Interphase
  tkgrid(tklabel(frame_All,image=GC_IonExtractor_logo))
  tkgrid(tklabel(frame_All,text="ENTER PARAMETERS"),columnspan=4, pady = 10, sticky="w")
  tkgrid(frame_All)
  tkgrid(frameL,pady= 0, padx= 10)
  tkgrid(frameButton)
}
######################################### LcOrbitalTRap #########################################
LcOrbitalTRap <- function(){
  # Load all require library
  require(multtest)
  require(xcms)
  require(faahKO)
  require(CAMERA)
  require(flux)
  require(tcltk)

  require (mzmatch.R)
  mzmatch.init(version.1=FALSE)
  #require(RColorBrewer)

  ################# Setup working directory ######################

  WorkingDir <-function(){
    Main.Folder <- choose.dir(caption = "Select working directory")
    all.Main.Folder<<-Main.Folder
    setwd(all.Main.Folder)
  }

  ################# RUN LC Orbital Trap peak area integration ####

  Run.LCMS.Analysis<-function(ppm_b=2.5, peakwidth_b=c(5,20), AllScangRange=TRUE, RtRangeM=c(0,13), SignalNoiseR_b=50, Prefilter_b=c(3,5000), mzwid_b=0.015, bw_b=2, modeP="positive"){

    if(exists("all.Main.Folder")){
      setwd(all.Main.Folder)
    } else {
      all.Main.Folder<<-choose.dir(caption = "Select working directory")
      setwd(all.Main.Folder)
    }
    cdfpath<-all.Main.Folder

    ## Load all raw files
    print("STEP 1a: Upload all LC-orbitral raw data in .mzXML formate")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    cdffiles <- cdffiles[grep(c(".cdf$|.mzXML$"),cdffiles, ignore.case=TRUE)]
    print(cdffiles)

    ## Print parameters
    print("-----XCMS Parameter Settings----")
    print(paste("Mass Resolution:", ppm_b, "ppm"))
    print(paste("Peak width:", toString(peakwidth_b), "s"))
    print(paste("RT bin (bw):", bw_b, "s"))
    print(paste("mz varation (mzwid):", mzwid_b))
    print(paste("Signal-to-noise ratio:", SignalNoiseR_b))
    print(paste("Pre-filter:", toString(Prefilter_b)))
    print(paste("Polarity:", modeP))

    ## Peak identitication
    print("STEP 1b: Peak Identications (long waiting time: ~2 min per sample)")
    RtRangeS<-round(RtRangeM*60,0) # convert minute to second and no decemal place

    if(AllScangRange==TRUE) {xset <- xcmsSet(cdffiles, ppm = ppm_b, peakwidth = peakwidth_b, method = "centWave", nSlaves=4, snthresh=SignalNoiseR_b, prefilter= Prefilter_b, integrate = 1, fitgauss= TRUE, verbose.columns=TRUE) } else {
      xset <-xcmsSet(cdffiles, ppm = ppm_b, peakwidth = peakwidth_b, method = "centWave" , scanrange=RtRangeS, nSlaves=4, snthresh=SignalNoiseR_b, prefilter= Prefilter_b, integrate = 1, fitgauss= TRUE, verbose.columns=TRUE) }


    ## save the file
    # xset_RT<-xset
    #save(xset_RT,file="BackUp_ForRtFail.rda")

    ## Peak Match 1
    # print("STEP 3: Matching Peaks Across Samples")
    # xset <- group(xset, bw = 2, mzwid = 0.015)

    ## RT correction
    print("STEP 1c: Retention Time Correction")
    xset2 <- retcor(xset, method ="obiwarp", plottype = c("deviation"), profStep=0.5, center=3)

    ## Peak Match 2
    print("STEP 1d: Matching Peaks Across Samples")
    xset2 <- group(xset2, bw = bw_b , mzwid = mzwid_b)
    print(xset2)

    ## Fill up the missing values
    print("STEP 1e: Filling in Missing Peak Data")
    xset3 <- fillPeaks(xset2)

    ## Save the file in the working folder
    # peakTable(xset3, filebase="Final LC-MS result")
    # print("Data process is completed and save Final LC-MS result.tsv")
    ## save the process LC-MS data in the working folder and can be load for futher downstream analysis

    if(AllScangRange==TRUE) {RtRangFinal=" AllScan"} else { RtRangFinal=toString(RtRangeM) }
    LongFileName <-paste("Processed raw data(ppm_",ppm_b,")(peakwidth_",toString(peakwidth_b),")(RtRange_",RtRangFinal,")(bw_",bw_b,")(mzwid_",
                         mzwid_b,")(SignalNoise_",SignalNoiseR_b,")(Pre-filter_",toString(Prefilter_b),")(Polarity_",modeP,").rda",sep="")
    #rawdata_FileName <- tclvalue(tkgetSaveFile(initialfile= LongFileName))
    save(xset3,file= LongFileName)

    ## Annotation by CAMERA packages
    xset3_annotation <- xset3
    print("STEP 1f: Annotation of adduct fragments and isotopes in progress....")
    if(length(levels(xset3_annotation@phenoData$class))==2){
      peaklist<- annotateDiffreport(xset3_annotation,polarity=modeP)
      write.csv(peaklist, file=paste("Final LC-MS result(XCMS)_TTtest","_",modeP,".csv",sep=""))
    } else {
      an <- xsAnnotate(xset3_annotation)
      an <- groupFWHM(an)
      an <- findIsotopes(an)  # optional but recommended.
      an <- groupCorr(an) # optional but very recommended step
      an <- findAdducts(an,polarity=modeP)
      peaklist <- getPeaklist(an)
      write.csv(peaklist, file=paste("LC-MS result(XCMS)","_",modeP,".csv",sep=""))
    }
    print("LC-MS data extraction is completed and saved as LC-MS result(XCMS).csv")
    tkmessageBox(title = "MassOmics-LCMS processing", message = paste("LC-MS data extraction is completed and saved as LC-MS result(XCMS).csv"))
  }


  ########################## If the RT doesn't work

  #  RT_correct_error <- function(){
  #    Save_xsetRT<-choose.files(caption = "Select a processed raw data (e.g.BackUp_ForRtFail.rda")
  #    load(Save_xsetRT)
  #    print("STEP 3: Matching Peaks Across Samples")
  #    xset_RT <- group(xset_RT)
  #    xset_RT <- fillPeaks(xset_RT)
  #    xset3 <- xset_RT
  #    peakTable(xset_RT, filebase="Final LC-MS result(No RT correction)")
  ## need to save the files
  #    rawdata_FileName <- tclvalue(tkgetSaveFile(initialfile="Processed raw data(ppm)(peakwidth)(RrRange)_ForRTFail.rda"))
  #    save(xset3,file=rawdata_FileName)
  #  }

  ########################## mzMatch ####################################

  mzMatch<-function(mode_b = "positive", ppm_b = 5, rtwindow_b = 20, BlankName= "Blank", BlankFilter= TRUE, javaMemorySizeG=4){




    print("-----mzMatch Parameter Settings----")
    print(paste("Mass Resolution:", ppm_b, "ppm"))
    print(paste("Retention time window:", rtwindow_b, "s"))
    print(paste("Blank filter:", BlankFilter))
    print(paste("Name for blank class:", BlankName))
    print(paste("Polarity:", mode_b))
    print(paste("Memory usuage:", javaMemorySizeG, "GB"))
    javaMemorySize<-javaMemorySizeG*1000
    mzmatch.init(memorysize= javaMemorySize, version.1=FALSE)

    if(exists("all.Main.Folder")){
      setwd(all.Main.Folder)
    } else {
      all.Main.Folder<<-choose.dir(caption = "Select working directory")
      setwd(all.Main.Folder)
    }

    Save_xset3<-choose.files(caption = "Select a XCMS processed raw data (e.g.Processed raw data(ppm)(peakwidth)....rda)")
    load(Save_xset3)


    ## Convert xcms process file to .PeakML
    print("STEP 2a: Convert xcms process data to .peakML")
    PeakML.xcms.write(xset3, filename= "1_LCMS.peakML", ionisation= mode_b)

    ## Remove noise using blank samples
    print("STEP 2b:Remove noise using blank")
    if(BlankFilter==TRUE){
      PeakML.BlankFilter(filename="1_LCMS.peakML",
                         ionisation= mode_b,
                         outputfile="1b_LCMS_BlankFilter.peakml",
                         BlankSample=BlankName)
    }


    ## fill the gap
    print("STEP 2c: Gap Filling the missing values")
    if(BlankFilter==TRUE){
      PeakML.GapFiller(filename = "1b_LCMS_BlankFilter.peakml", ionisation = mode_b, outputfile = "2_LCMS_GapFilling.peakml", ppm = 0, rtwin = 0, fillAll=FALSE, Rawpath=NULL, nSlaves = 4)
    } else {
      PeakML.GapFiller(filename = "1_LCMS.peakML", ionisation = mode_b, outputfile = "2_LCMS_GapFilling.peakml", ppm = 0, rtwin = 0, fillAll=FALSE, Rawpath=NULL, nSlaves = 4)
    }

    ## Group peaks
    print('STEP 3c: Group realtive peak together')
    mzmatch.ipeak.sort.RelatedPeaks(i="2_LCMS_GapFilling.peakml", v=T,
                                    o="3_LCMS_GapFilling_GroupPeak.peakml",  ppm=ppm_b, rtwindow=rtwindow_b, JHeapSize= javaMemorySize)

    ## Identify peaks from databases
    print("STEP 4c: Peak annoatation by databases")
    DBS <- dir(paste(find.package("mzmatch.R"), "/dbs", sep=""), full.names=TRUE)
    DBS # print a list of library
    DBS <- paste(DBS[c(1:15)],collapse=",")

    if(mode_b=="positive"){
      mzmatch.ipeak.util.Identify(i="3_LCMS_GapFilling_GroupPeak.peakml", v=T,
                                  o="4_LCMS_GapFilling_GroupPeak_Annotation(Positive).peakml", ppm=ppm_b, databases=DBS, polarity= mode_b,
                                  adducts="M+H,M+ACN+Na,M+Na,M+K,M+ACN+H", JHeapSize= javaMemorySize) ## need to fix
      mzmatch.ipeak.convert.ConvertToText (
        i="4_LCMS_GapFilling_GroupPeak_Annotation(Positive).peakml",
        o="Final_annoated_LCMS_data(positive).txt", databases=DBS,
        annotations="identification,moleculeName,npeaks,ppm,adduct,relation.ship,charge", JHeapSize= javaMemorySize)
    } else {
      mzmatch.ipeak.util.Identify(i="3_LCMS_GapFilling_GroupPeak.peakml", v=T,
                                  o="4_LCMS_GapFilling_GroupPeak_Annotation(Negative).peakml", ppm=ppm_b, databases=DBS, polarity= mode_b,
                                  adducts="M-H", JHeapSize= javaMemorySize) ## need to fix
      mzmatch.ipeak.convert.ConvertToText (
        i="4_LCMS_GapFilling_GroupPeak_Annotation(Negative).peakml",
        o= "Final_annoated_LCMS_data(Negative).txt", databases=DBS,
        annotations="identification,moleculeName,npeaks,ppm,adduct,relation.ship,charge", JHeapSize= javaMemorySize)
    }
    tkmessageBox(title = "MassOmics-LCMS processing", message = paste("LC-MS data filter and annotation is completed and saved as Final_annoated_LCMS_data.txt"))
    print("Analysis is completed and saved as - Final_annoated_LCMS_data.txt")
    ## end if mzmatch
  }

  ## manually check data
  ManualCheck <-function(){
    require (tcltk)
    require (mzmatch.R)
    mzmatch.init(version.1=FALSE)
    PeakML.Viewer (JHeapSize=5000)
    tkmessageBox(title = "PeakML.Viewer", message = paste("May take long time to open PeakML.Viewer and open a file in .peakml format (e.g.4_LCMS_GapFilling_GroupPeak_Annotation.peakml"))
  }
  ExportToTXT <- function(javaMemorySizeG=4){
    require (mzmatch.R)
    require (tcltk)
    if(exists("all.Main.Folder")){
      setwd(all.Main.Folder)
    } else {
      all.Main.Folder<<-choose.dir(caption = "Select working directory")
      setwd(all.Main.Folder)
    }
    javaMemorySize<-javaMemorySizeG*1000
    mzmatch.init(memorysize= javaMemorySize, version.1=FALSE)
    DBS <- dir(paste(find.package("mzmatch.R"), "/dbs", sep=""), full.names=TRUE)
    DBS <- paste(DBS[c(1:15)],collapse=",")

    InputName <- basename(choose.files(caption = "Select a .peakmL file (e.g.manual check.peakml)"))
    OutputName <-  basename(tclvalue(tkgetSaveFile(initialfile= "LC-MS_manual_check.txt")))

    mzmatch.ipeak.convert.ConvertToText (
      i=InputName,
      o=OutputName, databases=DBS,
      annotations="identification,moleculeName,npeaks,ppm,adduct,relation.ship,charge",
      JHeapSize= javaMemorySize
    )
    print ("Done!")
  }



  ################# stataiscally analysis for 2 conditions ####
  LCMS_TTest<-function(nplot=100){
    if(exists("all.Main.Folder")){
      setwd(all.Main.Folder)
    } else {
      all.Main.Folder<<-choose.dir(caption = "Select working directory")
      setwd(all.Main.Folder)
    }

    #Auotmatically load the proces data if xset3 not present
    #if(!exists("xset3")){
    #  Save_xset3<-choose.files(caption = "Select Processed raw data")
    #  load(Save_xset3)
    #}

    Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
    load(Save_xset3)
    condition_1<- as.character(unique(xset3@ phenoData $class)[[1]])
    condition_2<- as.character(unique(xset3@ phenoData $class)[[2]])
    print("Statistical analysis in progress....")
    reporttab <- diffreport(xset3, condition_1, condition_2, "Final LC-MS result(T-Test)", nplot, metlin = 0.15, h=480, w=640)
    print("The TTest report is completed and save as Final LC-MS result(T-Test.tsv")
  }

  ################# Plot ion extracted chromatogram base on ID

  choose.xset3<-function(){
    Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
    load(Save_xset3)
    xset3_ID<<-xset3
  }

  plot.ID<-function(ID.NO=712){

    if(!exists("xset3_ID")){
      Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
      load(Save_xset3)
      xset3_ID<<-xset3
    }
    eiccor <- getEIC(xset3_ID, groupidx = ID.NO)
    # Mass range
    mzmin<-round(as.numeric(eiccor@ mzrange[1,1]),5)
    mzmax<-round(as.numeric(eiccor@ mzrange[1,2]),5)

    # RT range
    RTmin<-round(as.numeric(eiccor@ rtrange[1,1]),1)
    RTmax<-round(as.numeric(eiccor@ rtrange[1,2]),1)
   
    # plot and save image
    fig.names<-paste("ID(",ID.NO,")_Ion(",mzmin,"-",mzmax, "mz)_RT(", RTmin,"-",RTmax,"s).png", sep="")
    
    png(fig.names)
    plot(eiccor, xset3_ID)
    dev.off() ## Need to make it plot the figure.....>_<
  }

  ################# LC-MS Single ion extraction ##############

  RUN.LCMS.Ion.Extractor<-function(){
    require(xcms)
    require(flux)


    # Setup working dirctory

    WorkingDir <-function(){
      Main.Folder <- choose.dir(caption = "Select working directory")
      all.Main.Folder<<-Main.Folder
      setwd(all.Main.Folder)
    }

    # Located the process file
    choose.xset3_Ion<-function(){
      Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
      load(Save_xset3)
      xset3_Ion<<-xset3
    }

    # LC-MS IonExtractor
    LCMS.Ion.Extractor<-function(mz.1 = 210.1281,
                                 mz.2 = 210.1289,
                                 RT.Time1 = 190,
                                 RT.Time2 = 210,
                                 heightMAX= 1500000000,
                                 Metabolite.Names = "D5-tryptophan",
                                 peak="Peak Area"){

      if(exists("all.Main.Folder")){
        setwd(all.Main.Folder)
      } else {
        all.Main.Folder<<-choose.dir(caption = "Select working directory")
        setwd(all.Main.Folder)
      }

      if(!exists("xset3_Ion")){
        Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
        load(Save_xset3)
        xset3_Ion<<-xset3
      }

      # Extract specific ion and rentension time
      IonExtract <-getEIC(xset3_Ion,mzrange=cbind(mz.1,mz.2),rtrange=cbind(RT.Time1,RT.Time2))

      # Determine colour of figure
      #mycolors<-rainbow(length(t(unique(xset3@ phenoData))))
      #mycolors<-brewer.pal(length(t(unique(xset3@ phenoData))),"Set1")
      #names(mycolors) <- levels(xset3@ phenoData$class)
      #list.df<-t(xset3@ phenoData)
      #list.color<-mycolors[match(list.df,names(mycolors))]
      list.color<-rainbow(length(IonExtract@eic))

      IonExctor.df<-NULL

      for (i in 1:length(IonExtract@eic)){
        #Extraction
        Extra.Data<-data.frame(IonExtract@eic[i])
        ifelse(peak=="Peak Height", Peakvalue <- max(Extra.Data[,2]) ,Peakvalue<-auc(Extra.Data[,1], Extra.Data[,2]))
        name.file<-names(IonExtract@eic)[i]
        if (!exists("IonExctor.df")){
          IonExctor.df <- data.frame(t(c(name.file,Peakvalue)))
        } else {
          IonExctor.df <- rbind(IonExctor.df, data.frame(t(c(name.file,Peakvalue))))
        }
        print(paste("File", i, name.file,"done!"))

        # Plot
        Intensity<-Extra.Data[,2]
        Retention.Time<-Extra.Data[,1]
        if(i==1){
          dev.new.OS()
          ifelse(heightMAX<max(Intensity),heightMAX<-max(Intensity),heightMAX)
          plot(Intensity~Retention.Time,col= list.color[i],bty="n", type="l",main=paste(Metabolite.Names," (",mz.1,"-",mz.2, " m/z",") ","RT:",RT.Time1,"-",RT.Time2," min", sep=""),
               xlab="Retention Time (Sec)",ylab="Intensity", ylim=c(0,heightMAX))
        }else{
          lines(Intensity~Retention.Time,col= list.color[i])
        }
      }

      final.df2<- data.frame(t(IonExctor.df[2]))
      names(final.df2) <- t(IonExctor.df[1])
      row.names(final.df2) <- Metabolite.Names

      # Save files
      ifelse(peak=="Peak Height",ValueType<-"PeakHeight",ValueType<-"PeakArea")
      store <- paste(all.Main.Folder, "\\", Metabolite.Names,"(",mz.1,"-",mz.2," m.z",")RT",RT.Time1,"-",RT.Time2,"_",ValueType,".csv", sep = "")
      write.csv(final.df2, file =store, row.names = TRUE)
      print("Complete!")
      return(final.df2)
      print(paste("The file ",Metabolite.Names, ".csv", " was saved in ", all.Main.Folder, sep=""))
    }

    ## User interphse specific for LC-MS IonExtractor


    require(tcltk)
    tclRequire("BWidget")

    # Setup frame
    kk <- tktoplevel()
    frameAll <- tkframe(kk)
    tkwm.title(kk,"LCMS-IonExtractor")

    # Setup parameters
    xvar <- tclVar("210.1281")
    vvar <- tclVar("210.1289")
    yvar <- tclVar("190")
    zvar <- tclVar("210")
    wvar <- tclVar("D5-tryptophan")
    hvar <- tclVar("5000000")
    Peakvar <- tclVar("Peak Area")

    x.entry <- tkentry(background="white",kk, textvariable=xvar, width=17)
    v.entry <- tkentry(background="white",kk, textvariable=vvar, width=17)
    y.entry <- tkentry(background="white",kk, textvariable=yvar, width=17)
    z.entry <- tkentry(background="white",kk, textvariable=zvar, width=17)
    w.entry <- tkentry(background="white",kk, textvariable=wvar, width=17)
    h.entry <- tkentry(background="white",kk, textvariable=hvar, width=17)

    # Choose between peakheight or area
    Peak <- c("Peak Height","Peak Area")
    comboBox <- tkwidget(kk,"ComboBox",editable=FALSE,values=Peak,textvariable=tclVar("Peak Area"),width=15)

    # Reset function
    reset <- function() {
      tclvalue(xvar)<-""
      tclvalue(vvar)<-""
      tclvalue(yvar)<-""
      tclvalue(zvar)<-""
      tclvalue(wvar)<-""
      tclvalue(hvar)<-"5000000"
      Peakvar <- tclVar("Peak Area")
    }

    reset.but <- tkbutton(kk, text="Reset", command=reset, width=9)

    # Run function
    submit <- function() {
      w <- as.character(tclvalue(wvar))
      x <- as.numeric(tclvalue(xvar))
      v <- as.numeric(tclvalue(vvar))
      y <- as.numeric(tclvalue(yvar))
      z <- as.numeric(tclvalue(zvar))
      h <- as.numeric(tclvalue(hvar))
      Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]

      LCMS.Ion.Extractor(Metabolite.Names=w, mz.1=x, mz.2=v, RT.Time1=y, RT.Time2=z, heightMAX=h ,peak= Peakvar)
    }

    submit.but <- tkbutton(kk, text="Run", command=submit, width= 9)

    # Quit function
    quit.but <- tkbutton(kk, text = "Close Session",
                         command = function() {
                           tkdestroy(kk)
                         }
    )


    Wdir.but <- tkbutton(kk, text="Choose", command= WorkingDir ,width=9)
    Process.but <- tkbutton(kk, text="Choose", command= choose.xset3_Ion ,width=9)
    open.but <- tkbutton(kk,text="Open folder",width= 9,
                         command= function(){
                           dir = getwd()
                           shell.exec(dir)
                         }
    )

    ## load images
    image.path2<-c(paste(file.path(path.package(package="MassOmics")),"/R/LCMS_IonExtractor.gif", sep=""))
    logo_1<-tcl("image",  "create", "photo", "MaxC4.image", file=image.path2)

    # Upper frame
    frameUpper<- tkframe(frameAll,borderwidth=2, relief="groove",padx=5,pady=5)
    tkgrid(tklabel(frameUpper,text="Setup Working Directory        "),  Wdir.but, pady = 2, padx= 10, sticky="w")  #columnspan=4
    tkgrid(tklabel(frameUpper,text="Select processed LC-MS file     "),  Process.but, pady = 2, padx= 10, sticky="w")

    #middle frame
    frameMiddle<- tkframe(frameAll,borderwidth=2, relief="groove",padx=5,pady=5)
    tkgrid(tklabel(frameMiddle,text="Metabolite Name"), w.entry, pady= 10, padx= 10, sticky="w")
    tkgrid(tklabel(frameMiddle,text="Reference Ion (min,max)"), x.entry, v.entry, tklabel(frameMiddle,text="m/z"), pady= 8, padx= 10, sticky="w")
    tkgrid(tklabel(frameMiddle,text="Retention Time (min,max)"), y.entry, z.entry, tklabel(frameMiddle,text="Sec"), pady= 8, padx= 10, sticky="w")
    tkgrid(tklabel(frameMiddle,text="Min Intensitiy for y-axis"), h.entry, pady= 8, padx= 10, sticky="w")
    tkgrid(tklabel(frameMiddle,text="Integration type"),comboBox, pady= 8, padx= 10, sticky="w")

    # Button
    frameButton<- tkframe(frameAll)
    tkgrid(tklabel(frameButton,text=""), reset.but, submit.but, open.but, quit.but, pady= 10, padx= 10, sticky="w")

    # Interphase structure
    tkgrid(tklabel(kk,image=logo_1))
    tkgrid(tklabel(frameAll,text="1. SETUP"), pady = 10, padx= 20, sticky="w")
    tkgrid(frameAll)
    tkgrid(frameUpper, pady= 0, padx= 30, sticky="w")
    tkgrid(tklabel(frameAll,text="2. ENTER PARAMETERS"), pady = 10, padx= 20, sticky="w")  #columnspan=4
    tkgrid(frameMiddle, pady= 0, padx= 30)
    tkgrid(frameButton, pady= 10, padx= 30)
  }



  ## Annotation of Adduct and isotopes
  # LCMS_Adduct_isotpes<-function(mode="positive"){
  #    if(exists("all.Main.Folder")){
  #      setwd(all.Main.Folder)
  #    } else {
  #      all.Main.Folder<<-choose.dir(caption = "Select working directory")
  #      setwd(all.Main.Folder)
  #    }
  #
  #    Save_xset3<-choose.files(caption = "Select a processed raw data (e.g.Processed raw data(ppm)(peakwidth).rda)")
  #    load(Save_xset3)
  #    xset3_annotation <- xset3
  #
  #    print("Annotation of adduct fragments and isotopes in progress....")
  #
  #    if(length(levels(xset3_annotation@phenoData$class))==2){
  #      peaklist<- annotateDiffreport(xset3_annotation,polarity=mode)
  #      write.csv(peaklist, file=paste("Final LC-MS result(Adduct_Isotope)_TTtest","_",mode,".csv",sep=""))
  #    } else {
  #
  #      an <- xsAnnotate(xset3_annotation)
  #      an <- groupFWHM(an)
  #      an <- findIsotopes(an)  # optional but recommended.
  #      an <- groupCorr(an) # optional but very recommended step
  #      an <- findAdducts(an,polarity=mode)
  #      peaklist <- getPeaklist(an)
  #      write.csv(peaklist, file=paste("Final LC-MS result(Adduct_Isotope)","_",mode,".csv",sep=""))
  #    }



  ### Visilization ##
  #plotEICs(an,1)
  #plotPsSpectrum(an,1,maxlabel=10)


  ## find nature lost
  #xs.pseudo <- findNeutralLoss(an,mzdiff=18.01,mzabs=0.01)
  #xs.pseudo@peaks #WATER LOST


  ## Calcualte moleuclar composition
  #library(Rdisop)
  #isolist <- getIsotopeCluster(an)


  # test<-isolist[5][[1]]$peaks

  # molecules<- decomposeMass(488.158)
  # molecules<- decomposeIsotopes(test[,1],test[,2],ppm=10)

  # errors = MFerrors(molecules,test[,1],test[,2])

  # cbind(getFormula(molecules), getScore(molecules), getValid(molecules))

  #  print("Annotation of adduct fragments and isotopes are completed and save as Final LC-MS result(Adduct_Isotope).csv")
  ## save the process LC-MS data in the working folder and can be load for futher downstream analysis
  # rawdata_FileName <- tclvalue(tkgetSaveFile(initialfile=paste("Annotation raw data","(",mode, ").rda",sep="")))
  # save(an,file=rawdata_FileName)

  # }



  ################# User interphase ####################

  ## Load library
  require(tcltk)
  require(rgl)
  tclRequire("BWidget")

  ## Setup tk enviroment ##
  LC_Orbital <- tktoplevel()
  frameOverall0 <- tkframe(LC_Orbital)
  frameOverall <- tkframe(LC_Orbital)
  frameOverall2 <- tkframe(LC_Orbital)
  frameOverall3 <- tkframe(LC_Orbital)
  tkwm.title(LC_Orbital ,"LC-MS data analysis")

  ## Functions in the buttons
  workingdirctory<-function(){WorkingDir()
  }
  LC_MS_analysis<-function(){
    ppm_a <- as.numeric(tclvalue(ppm.var))
    peakwidth_a <- as.numeric(strsplit(tclvalue(peakwidth.var), ",",)[[1]])
    RtRangeMa <- as.numeric(strsplit(tclvalue(RtRangeM.var), ",",)[[1]])
    SignalNoiseR_a <-as.numeric(tclvalue(SignalNoiseR.var))
    mzwid_a <- as.numeric(tclvalue(mzwid.var))
    bw_a <- as.numeric(tclvalue(bw.var))
    Prefilter_a <- as.numeric(strsplit(tclvalue(Prefilter.var), ",",)[[1]])
    modeB  <- modeA[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    rtVal <- as.character(tclvalue(rtValue))
    if (rtVal=="1")
      AllScangRangeA<-TRUE
    if (rtVal=="0")
      AllScangRangeA<-FALSE
    Run.LCMS.Analysis(ppm_b=ppm_a,peakwidth_b=peakwidth_a, AllScangRange=AllScangRangeA, RtRangeM = RtRangeMa, SignalNoiseR_b=SignalNoiseR_a, mzwid_b=mzwid_a,  Prefilter_b= Prefilter_a,  bw_b= bw_a, modeP=modeB)
  }

  LCMS.hlep.page<-function(){
    TMS.help  <- tktoplevel()
    tkwm.title(TMS.help,"Help")
    scr <- tkscrollbar(TMS.help, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
    txt <- tktext(TMS.help,bg="white",font="courier",yscrollcommand=function(...)tkset(scr,...))
    tkgrid(txt,scr)
    tkgrid.configure(scr,sticky="ns")

    ## Insert text here
    tkinsert(txt,"end","Under construction\n")
    tkconfigure(txt, state="disabled")
    tkfocus(txt)
  }

  TT.Test<-function(){
    nplot_a <- as.numeric(tclvalue(nPlot.var))
    LCMS_TTest(nplot=nplot_a)
  }

  workingdirctoryID<-function(){choose.xset3()
  }

  ID.Ion.Extraction<-function(){
    ID.NO_a <- as.numeric(tclvalue(ID.NO.var))
    plot.ID(ID.NO=ID.NO_a)
  }

  #  adduct_isotope<-function(){
  #    modeB  <- modeA[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
  #    LCMS_Adduct_isotpes(mode=modeB)
  #  }

  # RT_Error<-function(){
  #    RT_correct_error()
  #  }

  xcms.Tips<-function () {
    require(tcltk)
    kk <- tktoplevel()
    tktitle(kk) <- "Suggested XCMS parameter settings"
    image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/XCMS_parameter_tips.gif", sep=""))
    FrontImage<-tcl("image",  "create", "photo", "setting.image", file=image.path)
    tkgrid(tklabel(kk,image=FrontImage))
  }

  MzMatch.ini <- function() {
    ppm_a <- as.numeric(tclvalue(ppm_mz.var))
    memory_a <-as.numeric(tclvalue( memory.var))
    rtwindow_a <- as.numeric(tclvalue(rtwindow_mz.var))
    BlankName_a <- as.character(tclvalue(BlankName_mz))
    BnaVal <- as.character(tclvalue(BlankName_apply_Value))
    if (BnaVal=="1")
      BlankFilter_a<-TRUE
    if (BnaVal=="0")
      BlankFilter_a<-FALSE
    mode_a  <-  mode_mz.var[[as.numeric(tclvalue(tcl(comboBox_MzMatch,"getvalue")))+1]]
    mzMatch(mode_b = mode_a, ppm_b = ppm_a, rtwindow_b = rtwindow_a, BlankName=  BlankName_a, BlankFilter= BlankFilter_a, javaMemorySizeG= memory_a)
  }

  manualCheck.ini <- function() {
    ManualCheck()
  }
  exportText.ini <- function() {
    memory_a <-as.numeric(tclvalue(memory_b.var))
    ExportToTXT(javaMemorySizeG= memory_a)
  }

  ## Creat buttons
  Wdir.but <- tkbutton(LC_Orbital, text="Choose", command=workingdirctory,width=6)
  integrate.but <- tkbutton(LC_Orbital, text="Run", command=LC_MS_analysis,width=6)
  help.but <-tkbutton(LC_Orbital, text="Help", command= LCMS.hlep.page,width=6)
  TT.Test.but <- tkbutton(LC_Orbital, text="Run", command=TT.Test, width=6)
  WdirID.but <- tkbutton(LC_Orbital, text="Choose", command=workingdirctoryID,width=6)
  ID.Ion.Extraction.but <- tkbutton(LC_Orbital, text="Run", command=ID.Ion.Extraction, width=6)
  Ion.Extraction.but <- tkbutton(LC_Orbital, text="Run", command=RUN.LCMS.Ion.Extractor, width=6)
  # adduct_isotope.but <- tkbutton(LC_Orbital, text="Run", command=adduct_isotope, width=6)
  # RT.Error.but <- tkbutton(LC_Orbital, text="Run", command=RT_Error, width=6)
  xcms.setting.but <- tkbutton(LC_Orbital, text="Tip", command=  xcms.Tips, width=6)
  MzMatch.but <- tkbutton(LC_Orbital, text="Run", command= MzMatch.ini, width=6)
  ManualCheck.but <- tkbutton(LC_Orbital, text="Run", command= manualCheck.ini, width=6)
  exportText.but <- tkbutton(LC_Orbital, text="Run", command= exportText.ini, width=6)



  ## create entry box and define their default values
  # PPM
  ppm.var <- tclVar("2.5")
  ppm.entry <- tkentry(background="white",LC_Orbital, textvariable=ppm.var, width=8)
  # +6Width
  peakwidth.var <- tclVar("5,20")
  peakwidth.entry <- tkentry(background="white",LC_Orbital, textvariable=peakwidth.var, width=8)
  # ScangeRange
  RtRangeM.var <- tclVar("0,13")
  RtRangeM.entry <- tkentry(background="white",LC_Orbital, textvariable=RtRangeM.var, width=8)
  rtx <- tkcheckbutton(LC_Orbital)
  rtValue <- tclVar("1")
  tkconfigure(rtx,variable=rtValue)
  # Signal noise ratio
  SignalNoiseR.var <- tclVar("50")
  SignalNoiseR.entry <- tkentry(background="white",LC_Orbital, textvariable=SignalNoiseR.var, width=8)
  # Prefilter
  Prefilter.var <- tclVar("3,5000")
  Prefilter.entry <- tkentry(background="white",LC_Orbital, textvariable= Prefilter.var, width=8)
  # mzwid
  mzwid.var <- tclVar("0.015")
  mzwid.entry <- tkentry(background="white",LC_Orbital, textvariable= mzwid.var, width=8)
  #bw
  bw.var <- tclVar("2")
  bw.entry <- tkentry(background="white",LC_Orbital, textvariable= bw.var, width=8)
  # Number of plot
  nPlot.var<- tclVar("10")
  nPlot.entry <- tkentry(background="white",LC_Orbital, textvariable=nPlot.var, width=8)
  # Chramtogram ID
  ID.NO.var<- tclVar("1")
  ID.NO.entry <- tkentry(background="white",LC_Orbital, textvariable=ID.NO.var, width=8)
  # Positive or negative for addcut determination
  modeA <- c("positive","negative")
  comboBox <- tkwidget(LC_Orbital,"ComboBox",editable=FALSE,values=modeA,textvariable=tclVar("positive"),width=7)



  ## variable for mzMatch

  ppm_mz.var <- tclVar("5")
  ppm_mz.entry <- tkentry(background="white",LC_Orbital, textvariable= ppm_mz.var, width=8)
  rtwindow_mz.var <- tclVar("10")
  rtwindow_mz.entry <- tkentry(background="white",LC_Orbital, textvariable= rtwindow_mz.var, width=8)
  BlankName_mz <- tclVar("Blank")
  BlankName_mz.entry <- tkentry(background="white",LC_Orbital, textvariable= BlankName_mz, width=8)
  BlankName_apply <- tkcheckbutton(LC_Orbital)
  BlankName_apply_Value <- tclVar("1")
  tkconfigure(BlankName_apply,variable=BlankName_apply_Value)
  mode_mz.var <- c("positive","negative")
  comboBox_MzMatch <- tkwidget(LC_Orbital,"ComboBox", editable=FALSE, values= mode_mz.var, textvariable=tclVar("positive"), width=7)
  ## memory
  memory.var <-tclVar("4")
  memory.entry <- tkentry(background="white",LC_Orbital, textvariable=  memory.var, width=8)

  memory_b.var <-tclVar("4")
  memory_b.entry <- tkentry(background="white",LC_Orbital, textvariable=  memory_b.var, width=8)


  ### Frame structure of user interphase
  ## UpperFrame
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper,text="Setup Working Directory            "),Wdir.but, tklabel(frameUpper,text="               "), sticky="w")
  ## MiddleFrame
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameMid,text="Mass Resolution"), ppm.entry, tklabel(frameMid,text="ppm"),sticky="w")
  tkgrid(tklabel(frameMid,text="Peak width (min,max)"), peakwidth.entry,tklabel(frameMid,text="second"),sticky="w")
  tkgrid(tklabel(frameMid,text="RT bin (bw)"), bw.entry, tklabel(frameMid,text="second"), sticky="w")
  tkgrid(tklabel(frameMid,text="mz varation (mzwid)"), mzwid.entry, sticky="w")
  tkgrid(tklabel(frameMid,text="Signal-to-noise ratio"), SignalNoiseR.entry, sticky="w")
  tkgrid(tklabel(frameMid,text="Pre-filter "), Prefilter.entry, sticky="w")
  tkgrid(tklabel(frameMid,text="All RT range (Ignore below)"), rtx, sticky="w")
  tkgrid(tklabel(frameMid,text="RT range (min,max)"), RtRangeM.entry, tklabel(frameMid,text="minute"),sticky="w")
  tkgrid(tklabel(frameMid,text="MS modes"), comboBox, tklabel(frameMid,text="      "),sticky="w")
  tkgrid(tklabel(frameMid,text="Suggested XCMS settings"),xcms.setting.but, sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic LC-peak integration"),integrate.but,sticky="w")

  ## MiddleFrame_b
  # frameRtError <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  # tkgrid(tklabel(frameRtError,text="No retention correction            "),RT.Error.but, tklabel(frameRtError,text="            "),sticky="w")
  ## Lower Frame_adduct&Isotope
  # frameLowerAdductIsotope <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  # tkgrid(tklabel(frameLowerAdductIsotope,text="MS modes"), comboBox, tklabel(frameLowerAdductIsotope,text="      "),sticky="w")
  # tkgrid(tklabel(frameLowerAdductIsotope,text="Annotate Adducts and Isotopes"), adduct_isotope.but ,sticky="w")


  # Mzmatch
  framLowerMzMatch <- tkframe(frameOverall2,relief="groove",borderwidth=2,padx=5,pady=5)
  #tkgrid(tklabel(framLowerMzMatch,text="Select processed XCMS raw data "), MzMatch.but, tklabel(framLowerMzMatch,text="            "), sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="Mass Resolution"), ppm_mz.entry, tklabel(framLowerMzMatch,text="ppm"),sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="RT bin to group related peak  "), rtwindow_mz.entry, tklabel(framLowerMzMatch,text="second"),sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="RAM memory reserve for java"),  memory.entry, tklabel(framLowerMzMatch,text="GB"),sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="Use blank filter"), BlankName_apply, sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="Name of blank folder"),  BlankName_mz.entry, sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="MS modes"), comboBox_MzMatch, tklabel(framLowerMzMatch,text="      "),sticky="w")
  tkgrid(tklabel(framLowerMzMatch,text="Peak annotations"),  MzMatch.but,sticky="w")

  ## Manually check data
  framLowerCheck <- tkframe(frameOverall2,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(framLowerCheck, text="Manually check peaks             "), ManualCheck.but, tklabel(framLowerCheck, text="               "), sticky="w")
  tkgrid(tklabel(framLowerCheck,text="RAM memory reserve for java"),  memory_b.entry, tklabel(framLowerCheck,text="GB"),sticky="w")
  tkgrid(tklabel(framLowerCheck, text="Export .peakml to .txt           "), exportText.but, tklabel(framLowerCheck, text="               "), sticky="w")



  ## lowerFrame_ID
  # frameLowerID <- tkframe(frameOverall2,relief="groove",borderwidth=2,padx=5,pady=5)
  # tkgrid(tklabel(frameLowerID,text="Select processed raw data"),WdirID.but, tklabel(frameLowerID,text="            "), sticky="w")
  # tkgrid(tklabel(frameLowerID ,text="Chromatrogram ID #"), ID.NO.entry, tklabel(frameLowerID ,text="      "),sticky="w")
  # tkgrid(tklabel(frameLowerID ,text="Ion Extraction based on ID        "), ID.Ion.Extraction.but ,sticky="w")

  ## lowerFrame_ionExtractor
  frameLowerIon <- tkframe(frameOverall2,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameLowerIon, text="Re-integrate a single peak       "), Ion.Extraction.but, tklabel(frameLowerIon, text="               "),sticky="w")

  ## lowerFrame_1
  frameLower1 <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameLower1,text="Dispaly most significant ion       "), nPlot.entry, tklabel(frameLower1,text="   plot     "),sticky="w")
  tkgrid(tklabel(frameLower1,text="TT-Test"), TT.Test.but ,sticky="w")

  # Quit button
  quit.but <- tkbutton(LC_Orbital, text = "Close Session",
                       command = function() {
                         tkdestroy(LC_Orbital)
                         rgl.close()
                       }
  )

  # Open working directory
  open.but <- tkbutton(LC_Orbital,text="Open folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )

  ## load images
  image.path1<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_SIlas_Gravida_logo.gif", sep=""))
  logo<-tcl("image",  "create", "photo", "LCMS.image", file=image.path1)
  ## 3 bottom buttons
  #button <- tkframe(frameOverall3,relief="groove", borderwidth=0,padx=5,pady=5)
  #tkgrid(tklabel(button ,text="                                  "),quit.but,tklabel(button ,text="  "),open.but, tklabel(button ,text="                    "),help.but)

  ## change to right
  button <- tkframe(frameOverall2,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text=""),quit.but,tklabel(button ,text=""),open.but, tklabel(button ,text=""),help.but)



  ## Add window scroll
  ## Userinterphase
  tkgrid(frameOverall0,columnspan=2)
  tkgrid(frameOverall,frameOverall2)
  # tkgrid(frameOverall3,columnspan=2)
  tkgrid(tklabel(frameOverall0,image=logo),pady= 10, padx= 10)
  # tkgrid(tklabel(frameOverall,text=""),sticky="w")
  tkgrid(tklabel(frameOverall,text="   STEP 1.  Setup"),sticky="w")
  tkgrid(frameUpper,pady= 10, padx= 10)

  #tkgrid(tklabel(frameOverall,text=""),sticky="w")
  tkgrid(tklabel(frameOverall,text="   STEP 2.  Peak extraction & RT correction (XCMS)"),sticky="w")
  tkgrid(frameMid,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall ,text="  STEP 2b.  Statsitcal anaylsis for 2 conditions (OPTION for XCMS)          "),sticky="w")
  tkgrid(frameLower1 ,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall,text=""),sticky="w")
  #tkgrid(tklabel(frameOverall,text="   STEP 2b. Use When RT Correction Fail"),sticky="w")
  #tkgrid(tklabel(frameOverall,text=""),sticky="w")
  #tkgrid(frameRtError,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall,text="   STEP 3.  Annotation of Adducts and Isotopes"),sticky="w")
  #tkgrid(frameLowerAdductIsotope ,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall2,text=""),sticky="w")
  #tkgrid(tklabel(frameOverall2,text="  STEP 4.  Singel ion chromatogram view based on ID"),sticky="w")
  #tkgrid(frameLowerID ,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall2,text="  STEP 3.  Peak filter, filling, annotation, and identifications"),sticky="w")
  tkgrid(framLowerMzMatch ,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall2,text="  STEP 4.  Manully check peaks using PeakMLViewer"),sticky="w")
  tkgrid(framLowerCheck ,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall2,text="  STEP #.  IonExtractor For LC-MS (OPTION)"),sticky="w")
  tkgrid(frameLowerIon ,pady= 10, padx= 30)
  tkgrid(tklabel(frameOverall2,text=""),sticky="w")



  #tkgrid(tklabel(frameOverall,text=""),sticky="w")
  #tkgrid(tklabel(frameOverall,text="                 STEP #.  Check LC-MS results"),sticky="w")
  #tkgrid(frameLower,pady= 10, padx= 10)
  tkgrid(button, pady= 10, padx= 10)
}

######################################### GCMS_analysis #########################################
GCMS_analysis <- function(workdir=tk_choose.dir(caption = "Select working directory")){

  require(xcms)
  require(lattice)
  require(plyr)
  require(flux)
  library(tcltk)
  library(stringr)

  WorkingDir <- function(workdir = tk_choose.dir(caption = "Select working directory")){
    workdir <<- workdir
    setwd(workdir)
  }



  # Peak integration





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
    tkinsert(txt,"end","      A MS library in .msl formate\n")
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
  tkwm.title(TMS,"GC-MS data processing")

  #All the functions
  workingdirctory<-function(){WorkingDir()}
  
  Id_report1<-function(){
    Ret.Time.Filter_s <- as.numeric(tclvalue(RT.filter.var))
    NbVal <- as.character(tclvalue(NbValue))
    generate_graph<-F
    generate_rt_shift_graph_s=as.numeric(tclvalue(generate_rt_shift_graph))
    rtchecks=tclvalue(RTcorrection)
    if (rtchecks=="0"){rtchecks<-FALSE}
    if (rtchecks=="1"){rtchecks<-TRUE}
    if (NbVal=="0"){MsLibrary_s<-"InHouse"}
    if (NbVal=="1"){MsLibrary_s<-"NIST"}
    if (generate_rt_shift_graph_s=="0"){generate_graph<-F}
    if (generate_rt_shift_graph_s=="1"){generate_graph<-T}
    RT.shift.limt_s <- as.numeric(tclvalue(multimode.var))
    mz_L_s=as.numeric(tclGetValue("mz_L"))
    mz_U_s=as.numeric(tclGetValue("mz_U"))
    amdis_id_Summary(workdir=workdir,
              MsLibrary=MsLibrary_s,
              Ret.Time.Filter=Ret.Time.Filter_s,
              RT.shift.limt=RT.shift.limt_s,
              mz_L=mz_L_s,
              mz_U=mz_U_s,
              generate_rt_shift_graph=generate_graph,
              RTcorrection=rtchecks
              )
    message(paste("Summary report generated:", MsLibrary_s,"library","\nRetention time:",Ret.Time.Filter_s,"Retention time shift:",RT.shift.limt_s,"\nmz range:",mz_L_s,"-",mz_U_s))
  }
  
  Tips1<-function(){Tips()}

  GCMS_integration1<-function(){
    # RT.bin_s <-as.numeric(tclvalue(SliderValue))
    Peakvar <- Peak[[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]]
    cbVal <- as.character(tclvalue(cbValue))
    if (cbVal=="0"){GGReport_s<-"Slow"}
    if (cbVal=="1"){GGReport_s<-"Fast"}

    
    Mzbin<-as.numeric(tclvalue(MzBinValue))
    GCMS_integration(workdir=workdir, intensity_type=Peakvar,GGReport=GGReport_s,Ion.bin=Mzbin)
  }
  
  merage_RT_standard1<-function(){merage_RT_standard(workdir=workdir)}
  
  TMS.hlep.page1<-function(){TMS.hlep.page()}

  MZrange_autodetect1<-function(){
    mzrange=MZrange_autodetect()
    try(tclSetValue("mz_L", mzrange[1]))
    try(tclSetValue("mz_U", mzrange[2]))
  }

  ### all the button function
  Wdir.but <- tkbutton(TMS, text="Choose", command=workingdirctory,width=6)
  ID_report.but <- tkbutton(TMS, text="Run", command=Id_report1,width=6)
  #Tips.but<-tkbutton(TMS, text="Tips", command=Tips1,width=6)
  integrate.but<-tkbutton(TMS, text="Run", command= GCMS_integration1,width=6)
  #mergeSM.but<-tkbutton(TMS, text="Run", command= merage_RT_standard1,width=6)
  help.but<-tkbutton(TMS, text="Help", command= TMS.hlep.page1,width=6)

  ## fast or slow mode
  cb <- tkcheckbutton(TMS)
  cbValue <- tclVar("1")
  tkconfigure(cb,variable=cbValue)

  ## NIST library mode
  Nb <- tkcheckbutton(TMS)
  NbValue <- tclVar("0")
  tkconfigure(Nb,variable=NbValue)
  
  ## generate_rt_shift_graph_check
  generate_rt_shift_graph_check <- tkcheckbutton(TMS)
  generate_rt_shift_graph<- tclVar("0")
  tkconfigure(generate_rt_shift_graph_check,variable=generate_rt_shift_graph)
  
  ## RTcorrection
  RTcorrection_check <- tkcheckbutton(TMS)
  RTcorrection<- tclVar("0")
  tkconfigure(RTcorrection_check,variable=RTcorrection)
  
  ## NIST library mode
  #MZrange_L <- tclVar(TMS)
 # MZrange_U <- tclVar(TMS)
  mz_L <- 0
  mz_U <- 1000
  
  #tkconfigure(MZrange_L,variable=mz_L)
  #tkconfigure(MZrange_U,variable=mz_U)
  #MZrange_L <- tkentry(background="white",TMS, textvariable=mz_L ,width=6)
  #MZrange_L <- tkentry(background="white",TMS,width="6",textvariable=mz_L,validate="key",validatecommand="string is double %P")
  tclSetValue("mz_L",mz_L)
  tclSetValue("mz_U",mz_U)
  MZrange_L <- tkentry(background="white",TMS,width="6",textvariable="mz_L",validate="key",validatecommand="string is double %P")
  MZrange_U <- tkentry(background="white",TMS,width="6",textvariable="mz_U",validate="key",validatecommand="string is double %P")
  MZrange_autodetect.but <- tkbutton(TMS, text="Detect", command=MZrange_autodetect1,width=6)
  
  ## RT filter
  RT.filter.var <- tclVar("2.5")
  RT.entry <- tkentry(background="white",TMS, textvariable=RT.filter.var,width=8)

  ## MZ bin
  MzBinValue <- tclVar("0.5")
  MzBin.entry <- tkentry(background="white",TMS, textvariable=MzBinValue ,width=8)

  ## find mutiple peak
  multimode.var <- tclVar("60")
  multimode.entry <- tkentry(background="white",TMS, textvariable=multimode.var ,width=8)

  ## Upper
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper,text="Setup Working Directory"),Wdir.but, sticky="w")
  tkgrid(tklabel(frameUpper,text="AutoSelect Reference Ions (NIST)"),Nb,sticky="w")
  
  tkgrid(tklabel(frameUpper,text="Specify the MZ range"),MZrange_L,MZrange_U,sticky="w")
  tkgrid(tklabel(frameUpper,text="Detect the MZ range"),MZrange_autodetect.but,sticky="w")
  #tkgrid(tklabel(frameUpper,text="Add StandardMix RT (option)"),mergeSM.but,sticky="w ")
  tkgrid(tklabel(frameUpper,text="Performe RT correction"), RTcorrection_check, sticky="w")
  tkgrid(tklabel(frameUpper,text="Detect Multiple Peaks - RT Range above"),  multimode.entry,tklabel(frameUpper,text="s"),sticky="w")
  tkgrid(tklabel(frameUpper,text="Retention Time filter by library"), RT.entry,tklabel(frameUpper,text="min  "),sticky="w")
  tkgrid(tklabel(frameUpper,text="Generate plots for RT shift"),generate_rt_shift_graph_check,sticky="w")
  tkgrid(tklabel(frameUpper,text="Create a Summary Report"),ID_report.but,sticky="w")

  ## Middle
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)


  ## Slider

  # SliderValue <- tclVar("2")
  # SliderValueLabel <- tklabel(TMS,text=as.character(tclvalue(SliderValue)),widt=2)
  # tkconfigure(SliderValueLabel,textvariable=SliderValue)
  # slider <- tkscale(frameMid, from=1, to=5,
  #                   showvalue=F, variable=SliderValue,
  #                   resolution=0.5, orient="horizontal")


  ## choose between peakheight or area

  Peak <- c("Peak Area","Peak Height")
  comboBox <- tkwidget(TMS,"ComboBox",editable=FALSE,values=Peak,textvariable=tclVar("Peak Height"),width=12)

  # tkgrid(tklabel(frameMid,text="Manually correct the Identification Report"),Tips.but,sticky="w")
  tkgrid(tklabel(frameMid,text="Fast Mode (No Chromatograms)"),cb,sticky="w")
  tkgrid(tklabel(frameMid,text="Peak Height or Peak Area        "),comboBox, tklabel(frameMid,text=" "), sticky="w")
  # tkgrid(tklabel(frameMid,text="Minimum RT Bin (Default = 2 s) "),slider,SliderValueLabel,tklabel(frameMid,text="s"),sticky="w")
  tkgrid(tklabel(frameMid,text="Ion Mass Bin (m/z)"),  MzBin.entry,sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic GC-Peak Integration              "),integrate.but,sticky="w")

  ## lower

  font.size.var <- tclVar("0.5")
  font.size.entry <- tkentry(background="white",TMS, textvariable=font.size.var,width=4)


  rb1 <- tkradiobutton(TMS)
  rb2 <- tkradiobutton(TMS)
  rbValue <- tclVar("PeakDiagnosis")
  tkconfigure(rb1,variable=rbValue,value="PeakDiagnosis")
  tkconfigure(rb2,variable=rbValue,value="IonExtractor_Window")

  OnOK <- function(){
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="PeakDiagnosis"){
      font.size<- as.numeric(tclvalue(font.size.var))
      ifelse(sum(dir()=="GC-Peak Diagnosis Report.csv")==1,PeakDiagnosis(font.size=font.size),stop("Please select correct working directory"))
    }
    if (rbVal=="IonExtractor_Window"){
      IonExtractor_Window()}
  }
  OK.but <- tkbutton(TMS, text="Run", command=OnOK, width=6)

  frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  #tkgrid(tklabel(frameLower,text="Check TMS results"))
  tkgrid(tklabel(frameLower,text="Overlay All Chromatograms                  "),rb1,tklabel(frameLower,text="Font Size "),font.size.entry, sticky="w")
  tkgrid(tklabel(frameLower,text="IonExtractor              "), rb2,OK.but,sticky="w")


  quit.but <- tkbutton(TMS, text = "Close Session",
                       command = function() {
                         tkdestroy(TMS)
                         rgl.close()
                       }
  )

  ## Open working directory

  open.but <- tkbutton(TMS,text="Open Folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )

  ## load images
  image.path1a<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_SIlas_logo.gif", sep=""))
  logo_2<-tcl("image",  "create", "photo", "MaxC3.image", file=image.path1a)




  button <- tkframe(frameOverall,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text="                              "),quit.but,open.but, tklabel(button ,text="                    "),help.but)


  ## Userinterphase

  tkgrid(tklabel(frameOverall,image=logo_2))
  tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="         STEP 1. Setup"),sticky="w")
  tkgrid(frameUpper,pady= 10, padx= 10)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="         STEP 2. Peak Integration"),sticky="w")
  tkgrid(frameMid,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="         STEP 3. Check GC-MS results"),sticky="w")
  tkgrid(frameLower,pady= 10, padx= 10)
  tkgrid(frameOverall)
  #tkgrid(quit.but,help.but, pady= 10, padx= 10)
  tkgrid(button,pady= 10, padx= 10)
  tkgrid(tklabel(frameOverall,text="    "))
 
  
}


######################################### raw.DIC #########################################
raw.DIC <- function (main.folder = choose.dir(caption = "Select Working Directory"),
                   save = TRUE, output = "DI findal data",ppm = 3)
{
  require(xcms)
  old.wd <- getwd()
  main.folder <- main.folder
  setwd(main.folder)
  files <- c(dir())
  info <- file.info(files)
  isDir <- info$isdir
  conditions <- c(files[isDir == TRUE])
  top<-100
  final.df<-NULL
  remove(final.df)

  if (length(conditions)==1){
    filenames <- dir(conditions)
    filenames <- filenames[grep(".cdf$",filenames, ignore.case=TRUE)]
    for (q in 1:length(filenames)) {
      file1 <- filenames[q]
      name.file <- file1
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      file.name <- paste(conditions,name.file, sep = "\\")
      spectrum<-xcmsRaw(file.name,includeMSn=F)
      scans<-lapply(1:length(spectrum@scantime), function(i) getScan(spectrum,i))
      lens<-sapply(scans, nrow)
      #scans<-scans[lens < median(lens)*1.2 ]
      tops<-lapply(scans, function(x) {
        tmp=x[order(x[,2],decreasing=T)[1:top],]
        tmp[order(tmp[,1]),]
      })
      topss<-matrix(unlist(tops), nrow=top)
      all.reads<-data.frame(
        mz=as.vector(topss[,ii<-seq(1,ncol(topss),by=2)]),
        I=as.vector(topss[,ii+1]),
        scN=rep(1:(ncol(topss)/2), each=top)
      )
      all.reads<-all.reads[order(all.reads$mz),]
      clusters<-vector()
      clusters<-unlist(c((all.reads[1,]), clust=1))
      curClust<-1
      for (i in 2:nrow(all.reads)){
        #   i<-2
        # print(i)
        temp<-all.reads[i,]
        if((temp$mz - all.reads[i-1, "mz"]) >=0.07){
          curClust<-curClust+1
          print(i)
        }
        clusters<-rbind(clusters, data.frame((temp), clust=curClust))
      }
      #######show cluster number "k". Show all.reads
      ############
      all.reads[120:150,]
      k<-5 ################
      clusters[clusters$clus==k,]
      #plot(clusters[clusters$clus==k,][,c(1,3)], pch="o", cex=1)
      ############################################################

      singleClust<-as.numeric(names(table(clusters$clust)[table(clusters$clust) < round(0.1*length(spectrum@scanindex),0)]))# removes clusters which have less than 10 per cent of elements
      clusters.clean<-clusters[!clusters$clust %in% singleClust,]
      length(unique(clusters$clust))
      length(unique(clusters.clean$clust))
      means<-vector()
      for(i in levels(factor(as.character((clusters.clean$clust))))){
        print(i)
        means<-rbind(means, c(i, colMeans(clusters.clean[clusters.clean$clust == i, c("mz", "I")])))

      }
      rez<-data.frame(means)[,-1]
      names(rez)[2] <- name.file
      rez[1]<-round(as.numeric(as.character(rez$mz)),ppm)

      if (!exists("final.df")){
        final.df <- rez
        confirmation <- paste("File", q, "(", name.file,")", "done!", sep = " ")
        print(confirmation)
      } else {
        final.df <- merge(final.df, rez, by.x = "mz",
                          by.y = "mz", all = TRUE)

        confirmation <- paste("File", q, "(", name.file,
                              ")", "done!", sep = " ")
        print(confirmation)
      }
    }
  }else {
    stop("place all the cdf files in a folder")
  }

  final.df <- final.df[order(final.df$mz, decreasing = F),]

  if (save == TRUE) {
    sheet <- output
    store <- paste(main.folder, "\\", sheet, ".csv", sep = "")
    write.csv(final.df, file = store, row.names = FALSE)
  }
  return(final.df)
}

######################################### raw.peaks #########################################
raw.peaks <- function (main.folder = choose.dir(caption = "Select Working Directory"),
                     correct.RT = TRUE, method = "loess", save = TRUE, output = "Total MassFragments Result")
{
  require(xcms)
  old.wd <- getwd()
  main.folder <- main.folder
  setwd(main.folder)
  files <- c(dir())
  info <- file.info(files)
  isDir <- info$isdir
  conditions <- c(files[isDir == TRUE])


  if (length(conditions)==1){
    num.rep <- length(grep(".CDF", c(list.files(conditions[1],
                                                full.names = TRUE)), ignore.case = TRUE))
    samples <- grep(".CDF", c(list.files(conditions, full.names = TRUE)),
                    ignore.case = TRUE, value = TRUE)
    replicates <- (c(1:length(samples)))

    xset <- xcmsSet(samples)
    xset <- group(xset)
    if (correct.RT == TRUE) {
      xset <- retcor(xset, method = method)
      xset <- group(xset)
    }
    mat.raw <- groupval(xset, "medret", "into")
    mat.raw <- mat.raw[!duplicated(row.names(mat.raw)), ]
    Names <- as.character(row.names(mat.raw))
    mat.raw <- cbind(Names, mat.raw)

  }else {
    num.rep <- length(grep(".CDF", c(list.files(conditions[1],
                                                full.names = TRUE)), ignore.case = TRUE))
    samples <- grep(".CDF", c(list.files(conditions[1], full.names = TRUE)),
                    ignore.case = TRUE, value = TRUE)
    replicates <- (c(1:length(samples)))


    for (i in 2:length(conditions)) {
      num.rep <- c(num.rep, length(grep(".CDF", c(list.files(conditions[i],
                                                             full.names = TRUE)), ignore.case = TRUE)))
      samples <- c(samples, grep(".CDF", c(list.files(conditions[i],
                                                      full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      up.num <- length(grep(".CDF", c(list.files(conditions[i],
                                                 full.names = TRUE)), ignore.case = TRUE, value = TRUE))
      if (i < 3) {
        low.num <- (replicates[length(replicates)] + 1)
      }
      else {
        low.num <- max(replicates[i - 1, !duplicated(replicates[i -
                                                                  1, ])]) + 1
      }
      replicates <- rbind(replicates, c((low.num):(low.num +
                                                     (up.num - 1))))
    }

    xset <- xcmsSet(samples)
    xset <- group(xset)
    if (correct.RT == TRUE) {
      xset <- retcor(xset, method = method)
      xset <- group(xset)
    }
    mat.raw <- groupval(xset, "medret", "into")
    mat.raw <- mat.raw[!duplicated(row.names(mat.raw)), ]
    Names <- as.character(row.names(mat.raw))
    mat.raw <- cbind(Names, mat.raw)
    for (r in 1:nrow(replicates)) {
      rep.row <- replicates[r, ]
      rep.row <- rep.row[!duplicated(rep.row)]
      rep.name <- rep(conditions[r], length(rep.row))
      if (r == 1) {
        rep.name.final <- c(rep.name)
        rep.name.final <- c("Replicates", rep.name.final)
      }
      else {
        rep.name.final <- c(rep.name.final, rep.name)
      }
    }
    mat.raw <- rbind(rep.name.final, mat.raw)
  }
  row.names(mat.raw) <- 1:nrow(mat.raw)
  mat.raw <- data.frame(mat.raw)
  if (save == TRUE) {
    store <- paste(main.folder, "\\", output,".csv", sep = "")
    write.csv(mat.raw, file = store, row.names = FALSE)
    cat(paste("The file", output, "was saved in the folder",
              main.folder, "\n"))
  }
  setwd(old.wd)
  return(mat.raw)
  print("Data processing is completed and saved as Meta Result!")
}

######################################### run #########################################
run <- function(){
  library("pacman")
  p_load(tcltk2)
  tt <- tktoplevel(width=600, height=200)
  tktitle(tt) <- "MassOmics"


  topMenu <- tk2menu(tt)
  tkconfigure(tt, menu = topMenu)
  fileMenu <- tk2menu(topMenu, tearoff = FALSE)
  software <- tk2menu(topMenu, tearoff = FALSE)  # Our cascaded menu
  installation  <- tk2menu(topMenu, tearoff = FALSE)
  Update<-tkmenu(topMenu, tearoff = FALSE)
  Statistics <- tk2menu(topMenu, tearoff = FALSE)
  PAPi<- tkmenu(topMenu, tearoff = FALSE)
  Graph<-tkmenu(topMenu, tearoff = FALSE)



  tkadd(fileMenu, "cascade", label = "Installation", menu = installation)
  tkadd(fileMenu, "command", label = "Quit", command = function() tkdestroy(tt))
  tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
  tkadd(topMenu, "cascade", label = "Softwares", menu = software)
  tkadd(topMenu, "cascade", label = "Statistics",menu=Statistics)
  tkadd(topMenu, "cascade", label = "Graph",menu=Graph)
  tkadd(topMenu, "cascade", label = "Update", menu=Update)


  #tkadd(software, "command", label = "Run IonExtractor (Batch)",command = function() IonE())

  #tkadd(software, "command", label = "Run Metab",command = function() clean.fix())
 # tkadd(software, "command", label = "Run Metab 3",command = function() clean.fix3())
 # tkadd(software, "command", label = "Run GC-MS data anlaysis (Liz version)",command = function() Liz_analysis())
  tkadd(software, "command", label = "GC-MS data processing",command = function() GCMS_analysis())

 # tkadd(software, "command", label = "Run TMS data anlaysis",command = function() TMS_analysis())
  #tkadd(software, "command", label = "IonExtractor Single Mode",command = function() IonExtractor_Window())
 # tkadd(software, "command", label = "IonExtractor Single Mode_5ppm",command = function() IonExtractor_Window_5ppm())


 # tkadd(software, "command", label = "Run Total MassFragment Analysis",command = function() raw.peaks())
 # tkadd(software, "command", label = "Run DI-TOF data analysis",command = function() raw.DIC())
 # tkadd(software, "command", label = "Remove false positive results",command = function() runFalsePostive())
#  tkadd(software, "command", label = "ChemStation data analysis I (OLD)",command = function() Chemstation_DataMining())
#  tkadd(software, "command", label = "ChemStation data analysis II (Advance)",command = function() Liz_analysis_chemstation())
  tkadd(software, "command", label = "LC-MS data processing",command = function() LcOrbitalTRap())
  tkadd(software, "command", label = "Create a sub-NIST library",command = function() create_sub_library())
  tkadd(software, "command", label = "Data clean-up and normalization",command = function() DataCorrection())
  tkadd(software, "command", label = "Gangliosides accurate mass extraction",command = function() RunGangliosides())


  tkadd(installation, "command", label = "Required R packages",command = function() install.MassOmics_dependency())
  tkadd(installation, "command", label = "Required R packages for MzMatch",command = function() install.mzmatch())
  tkadd(Update, "cascade", label = "Version 2.5 (11June2017)")

  ## Statistics
  tkadd(Statistics, "command", label = "Run ANOVA & T-test",command = function() Omics_htest())
  tkadd(Statistics, "command", label = "Run PCA analysis",command = function() Omics_PCA())
  tkadd(Statistics, "command", label = "Run sparse PCA analysis",command = function() Omics_sPCA())
  tkadd(Statistics, "command", label = "Run PLSDA analysis",command = function() Omics_PLSDA())
  
  # Graph
  tkadd(Graph, "command", label = "Diagnose GC-Peak Auto Integration",command = function() PeakDiagnosis())

  ##PAPi
  tkadd(software, "cascade", label = "Pathway analysis", menu = PAPi)
  tkadd(PAPi ,"command", label = "Building KEGG database in your PC", command = function() buildDatabase ())
  tkadd(PAPi ,"command", label = "Replace name by KeggCodes", command = function() addKeggCodes())
  tkadd(PAPi ,"command", label = "Run pathway analysis", command = function() papi())

  #image
  image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/LabLogo.gif", sep=""))
  FrontImage<-tcl("image",  "create", "photo", "MaxC1.image", file=image.path)
  tkgrid(tklabel(tt,image=FrontImage))
  tkgrid(tklabel(tt,text="Maintainer: George GUO (George.GUO@auckland.ac.nz)"))
  tkfocus(tt)
  ##Update 29April2016
}

######################################### runFalsePostive #########################################
runFalsePostive <- function(){

  FalsePostive<-function (true = 0.2,  conditions=1:3, save = TRUE,
                          output = "Remove false positive results")
  {
    report_folder <- choose.dir(caption = "Select the Working directory")
    setwd(report_folder)
    final.df<- read.csv(choose.files(default = "Select the CSV file containing the input data",
                                     multi = FALSE, caption = "Select the CSV file containing the input data"))

    final.stat <- final.df
    abc <- 1:100
    selected.samples <- numeric()

    for (q in 1:length(conditions)) {
      cat(paste("Select samples from condition:", q, "\n"))
      samples <- tk_select.list(names(final.stat), multiple = TRUE,
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
    data.df1<-rbind(group, data.df)
    group<-as.numeric(group)
    reps <- factor(group)
    cat(paste("Data clean up is in progress..", "\n"))


    for (z in 1:length(levels(reps))) {
      column <- which(as.character(levels(reps)[z]) ==
                        as.character( data.df1[1, ]))
      yMat <-  data.df1[column]
      pb <- txtProgressBar(min = 0, max = nrow(yMat), style = 3, width = 50)
      for (i in 1:nrow(yMat)) {
        setTxtProgressBar(pb, i)
        if (sum(!is.na(yMat[i, ]))/ncol(yMat) < true) {
          yMat[i, ] <- NA
        }
      }
      cat(paste(" Finish Condition",z, "\n"))
      data.df1[column] <- yMat
    }
    yMat <-  data.df1[-1,]
    missing <- t(apply(yMat, 1, function(y) tapply(as.numeric(y),
                                                   reps, function(x) sum(!is.na(x)))))
    missing2 <- apply(t(apply(missing, 1, function(x) x == 0)),
                      1, sum)
    missing2 <- as.data.frame(missing2)

    data.df2<-cbind(final.df[,-selected.samples],yMat)

    data.df2 <- merge( data.df2, missing2, by = 0)
    data.df2 <- subset( data.df2,  data.df2$missing2 < ((length(levels(reps)))-(length(levels(reps)))*true)) ## Change the
    data.df2$missing2 <- NULL
    data.df2$Row.names <- NULL



    if (save == TRUE) {
      sheet <- output
      store <- paste(report_folder, "\\", sheet,"(",true*100,"%",")", ".csv", sep = "")
      write.csv( data.df2, file = store, row.names = FALSE)
      cat(paste("Process Complete!",sheet ,"was saved in the folder",
                report_folder, "\n"))
    }
    return( data.df2)


  }

  ####### user interphase


  require(tcltk)
  ConditionN <- tclVar("")
  true <- tclVar("50")
  gg <- tktoplevel()
  tkwm.title(gg,"Remove false postive data")
  c.entry <- tkentry(background="white",gg, textvariable=ConditionN)
  p.entry <- tkentry(background="white",gg, textvariable=true)
  reset <- function() {
    tclvalue(ConditionN)<-""
    tclvalue(true)<-""
  }
  reset.but <- tkbutton(gg, text="Reset", command=reset)
  submit <- function() {
    c <- as.numeric(tclvalue(ConditionN))
    conditions<-1:c
    True<-as.numeric(tclvalue(true))
    FalsePostive(conditions=conditions,true=True/100)
  }
  submit.but <- tkbutton(gg, text="submit", command=submit)

  quit.but <- tkbutton(gg, text = "Close Session",
                       command = function() {
                         tkdestroy(gg)
                       }
  )
  tkgrid(tklabel(gg,text="Enter Parameters"),columnspan=3, pady = 10)
  tkgrid(tklabel(gg,text="Number of conditions"), c.entry, pady= 10, padx= 10)
  tkgrid(tklabel(gg,text="Exclude compounds less than"), p.entry, tklabel(gg,text="%"), pady= 10, padx= 10)
  tkgrid(submit.but, reset.but, quit.but, pady= 10, padx= 10)
}

######################################### RunGangliosides #########################################
RunGangliosides <- function (){

  runG <-function(          output = "Peak area for gangliosides",
                            Ion.bin= 0.01,
                            WorkingFolder= tk_choose.dir(caption = "Select working directory"),
                            heightMAX= 300000,
                            info = read.csv(tk_choose.files(caption = "Accurate mass list for gangliosides.csv",  multi = FALSE), stringsAsFactors=FALSE),
                            RTrange=c(8,14)

  ) {

    require(xcms)
    require(lattice)
    require(plyr)
    require(flux)
    library(tcltk)

    print(paste("Mass Resolution:", Ion.bin))
    print(paste("heightMAX:", heightMAX))
    print(paste("RT range:", RTrange))


    setwd(WorkingFolder)
    final.check.data <- info
    library_file <-final.check.data
    files <- c(dir())
    info <- file.info(files)
    isDir <- info$isdir
    conditions_pre1 <- c(files[isDir == TRUE])
    conditions_pre<-c(conditions_pre1 ,"Manually choose folder")
    conditions <- tk_select.list(conditions_pre, multiple = FALSE,
                              title = paste("Folder contains all cdf files"))

    if (conditions=="Manually choose folder"){ conditions<-tk_choose.dir(default = "", caption = "A folder contains all cdf files")}

    ###

    final.df<-NULL
    remove(final.df)
    Graphic.df<-data.frame()
    filenames <- dir(conditions)
    filenames <- filenames[grep(c(".cdf$|.mzXML$"),filenames, ignore.case=TRUE)]
    colors<-rainbow(length(library_file[,1]))

    for (q in 1:length(filenames)) {
      file1 <- filenames[q]
      name.file <- file1
      name.file <- gsub(".CDF", "", name.file, ignore.case=TRUE)
      surefinal <- data.frame(Name=paste(library_file[,"Gangliosides"],"_", library_file[,"Ion"],sep=""),Base.Peak=0)
      raw_data <- xcmsRaw(filename = paste(conditions,
                                           filenames[q], sep = "\\"))

      for (h in 1:nrow(library_file)) {
        metabolite <-  library_file[h, ]
        R.Ion<-metabolite$Ion
        R.Time<-metabolite$RT.min.*60
        RT.lowerset<- R.Time-metabolite$bin.min.*60
        RT.upperset<- R.Time+metabolite$bin.min.*60

        IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(RT.lowerset, RT.upperset))
        abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
        Peakvalue<-auc(abundance$rt, abundance$intensity*10)
        surefinal$Base.Peak[h] <-  Peakvalue


        ## Plotting

        Intensity <- abundance$intensity
        Retention.Time<-abundance$rt/60

        if(h==1){
          dev.new.OS()
          ifelse(heightMAX<max(Intensity),heightMAX<-max(Intensity),heightMAX)
          plot(Intensity~Retention.Time,col= colors[h],bty="n", type="l",main=filenames[q],
               xlab="Retention Time (min)",ylab="Intensity", ylim=c(0,heightMAX),xlim=RTrange)
        }else{
          lines(Intensity~Retention.Time,col= colors[h])

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


    finalF.df <- cbind( Ref.Ion=library_file$Ion, Ret.Time=library_file$RT.min., RT.bin=library_file$bin.min., final.df)
    finalF.df <- finalF.df[order(finalF.df$Ret.Time , decreasing = FALSE),]


    ## Save Files
    print("Save file")

    store <- paste(WorkingFolder, "\\", output, ".csv", sep = "")
    FileName <- tclvalue(tkgetSaveFile(initialfile=store))
    write.csv(finalF.df, file = FileName, row.names = FALSE)

    print(paste("Ganglioside data extraction is completed and saved as: ",output,".csv", sep=""))


  }


  ################################################################################################################




  ## Building an interfance
  library(tcltk)
  require(rgl)
  tclRequire("BWidget")


  #Setup
  TMS <- tktoplevel()
  frameOverall <- tkframe(TMS)
  tkwm.title(TMS,"Accurate mass extraction for gangliosides")


  # Parameters

  ## MZ bin
  MzBinValue <- tclVar("0.01")
  MzBin.entry <- tkentry(background="white",TMS, textvariable=MzBinValue ,width=8)
  heightValue <- tclVar("300000")
  height.entry <- tkentry(background="white",TMS, textvariable=heightValue ,width=8)
  RangeM.var <- tclVar("8,14")
  RangeM.entry <- tkentry(background="white",TMS, textvariable=RangeM.var, width=8)


  # All the functions
  GCMS_integration1<-function(){
    Mzbin<-as.numeric(tclvalue(MzBinValue))
    Range <- as.numeric(strsplit(tclvalue(RangeM.var),",")[[1]])
    height <- as.numeric(tclvalue(heightValue))
    runG(RTrange=Range, heightMAX=height ,Ion.bin=Mzbin)
  }


  ## the button function

  integrate.but<-tkbutton(TMS, text="Run", command= GCMS_integration1,width=6)
  quit.but <- tkbutton(TMS, text = "Close Session",
                       command = function() {
                         tkdestroy(TMS)
                         rgl.close()
                       }
  )



  ## Interphase
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)

  # tkgrid(tklabel(frameMid,text="Minimum RT Bin (Default = 2 s) "),slider,SliderValueLabel,tklabel(frameMid,text="s"),sticky="w")
  tkgrid(tklabel(frameMid,text="Mass resolution"),  MzBin.entry,tklabel(frameMid,text="(m/z)"), sticky="w")
  tkgrid(tklabel(frameMid,text="RT range (min,max)"), RangeM.entry, tklabel(frameMid,text="min"),sticky="w")
  tkgrid(tklabel(frameMid,text="Min Intensitiy for y-axis"), height.entry, sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic accurate mass extraction              "),integrate.but,sticky="w")


  ## Open working directory

  open.but <- tkbutton(TMS,text="Open Folder",
                       command= function(){
                         dir = getwd()
                         shell.exec(dir)
                       }
  )

  ## load images
  image.path1b<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_ganglioside_logogif.gif", sep=""))
  logo_3<-tcl("image",  "create", "photo", "MaxC3X.image", file=image.path1b)
  button <- tkframe(frameOverall,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text="      "),quit.but,open.but)

  ## Userinterphase
  tkgrid(tklabel(frameOverall,image=logo_3))
  tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="        Gangliosides accurate mass extraction"),sticky="w")
  tkgrid(frameMid,pady= 10, padx= 30)
  tkgrid(frameOverall)
  #tkgrid(quit.but,help.but, pady= 10, padx= 10)
  tkgrid(button,pady= 10, padx= 10)
  tkgrid(tklabel(frameOverall,text="    "))

}

######################################### TMS_analysis #########################################
TMS_analysis <- function(main.TMS.folder=tk_choose.dir(caption = "Select working directory")){

  require(xcms)
  require(lattice)
  require(plyr)
  require(MassOmics)

  WorkingDir <-function(main.TMS.folder = tk_choose.dir(caption = "Select working directory")){
    main.TMS.folder <<- main.TMS.folder
    setwd(main.TMS.folder)
  }


  Id_report<-function(      main.TMS.folder=main.TMS.folder,
                            MS.L= tk_choose.files(caption="Select MS library (e.g.AgilentGCMS_TMS_Library) in .msl"),
                            amdis.report = read.csv(tk_choose.files(default = "Select the AMDIS report in .TXT",
                                                                    multi = FALSE, caption = "Select the AMDIS report in .TXT"),
                                                    sep = "\t", stringsAsFactors=FALSE),
                            File.name = "Identification Report (TMS).csv"){

    main.TMS.folder <- main.TMS.folder
    setwd(main.TMS.folder)
    getIonLib<-function(lib.fn = MS.L){
      lib.txt<-readLines(lib.fn)
      lib.txt<-lib.txt[-grep("^[ ]*$",lib.txt)]
      att.nms<-unique(sapply(strsplit(lib.txt[grep(":",lib.txt)],":"),function(x) x[1]))
      entry.nms<-sapply(strsplit(lib.txt[grep("NAME:",lib.txt)],":"), function(x) x[2])

      rez<-matrix(NA, nrow=length(entry.nms), ncol=length(att.nms), dimnames=list(NULL,att.nms ))

      starts=grep("NAME:", lib.txt)
      stops=c(starts[1:(length(starts)-1)]+diff(starts)-1,length(lib.txt))
      for (i in 1:length(starts)){
        #i=1
        tmp<-lib.txt[starts[i]:stops[i]]
        sapply(strsplit(tmp[grep(":",tmp)],":"), function(x) rez[i,x[1]]<<-x[2])
      }
      rez<-data.frame(rez,stringsAsFactors=FALSE)
      rez=within(rez, {
        RT=as.numeric(RT)
        RSN=as.numeric(RSN)
        NUM.PEAKS=as.numeric(NUM.PEAKS)
      })
      return(rez)

    }
    libr<-getIonLib()## returns list of peaks


    ############ Step 2. Generate the list of metabolites and clean up
    AmRep<-amdis.report
    AmRep<-AmRep[!AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)+2 & AmRep$RT.RT.lib.>median(AmRep$RT.RT.lib.)-2,]## check this one ##
    AmRep$Name <- gsub("?", "", AmRep$Name, fixed = TRUE)
    AmRep$Name <- gsub("^ ", "", AmRep$Name, perl=T)

    ## RT statistics from the AReport
    RT.stats<-t(sapply(split(AmRep$RT, AmRep$Name),function(x) c(RT.median=round(median(x,na.rm=T),3),
                                                                 RT.mean=round(mean(x,na.rm=TRUE),3),
                                                                 RT.sd=sd(x,na.rm=T),
                                                                 Number.of.ID=length(x[!is.na(x)]))))

    AmRep.RT.stats <- as.data.frame(cbind(Name=0,RT.stats,Quality=0, Exp.RT=0,Diff.RT=0,Diff.RT.by.Ribitol=0,ref_ion=0))##data frame for list of metabolites
    AmRep.RT.stats$ref_ion <- libr$RSN[match(rownames(RT.stats),libr$NAME)]
    AmRep.RT.stats$Exp.RT <- libr$RT[match(rownames(RT.stats),libr$NAME)]
    AmRep.RT.stats$Diff.RT <- AmRep.RT.stats$Exp.RT- AmRep.RT.stats$RT.median
    AmRep.RT.stats$Diff.RT.by.Ribitol <-AmRep.RT.stats$Diff.RT- AmRep.RT.stats["16.1473 min RIBITOL","Diff.RT"]
    AmRep.RT.stats$Quality <- "Good/Bad"
    AmRep.RT.stats$Name <- rownames(RT.stats)
    AmRep.RT.stats <- AmRep.RT.stats[order(AmRep.RT.stats$RT.median, decreasing = F),]
    write.csv(AmRep.RT.stats, file = File.name,row.names = FALSE)
    print(paste("Identificaiton report was created and save in",main.TMS.folder, sep=" "))
  }



  ### MERAGE TWO FILE TOGETHER


  merage_RT_standard<-function(main.TMS.folder=main.TMS.folder,
                               File.name = "Identification Report (TMS) with StandardMix_RT.csv"
  ){
    main.TMS.folder <- main.TMS.folder
    setwd(main.TMS.folder)
    data.df=read.csv(tk_choose.files(caption = "Load Samples Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    standard.df=read.csv(tk_choose.files(caption = "Load Standard Mixture Identification Report.csv",  multi = FALSE),stringsAsFactors=FALSE)
    RT_standard_names <- tk_select.list(names(standard.df), multiple = FALSE, title = paste("Select Rentention time"))
    SelectStandard<-standard.df[,RT_standard_names]
    standard.df.final<-data.frame(cbind(Name=standard.df$Name,StandardMixture_RT=SelectStandard))
    Merge.data <- merge(data.df, standard.df.final, by.x = "Name",
                        by.y = "Name", all = TRUE)

    Merge.data  <- Merge.data [order(Merge.data $RT.median, decreasing = F),]

    write.csv(Merge.data, file = File.name,row.names = FALSE)
    print(paste("Identificaiton report with StandardMix RT was created and save in",main.TMS.folder))

  }



  Tips<-function () {
    require(tcltk)
    kk <- tktoplevel()
    tktitle(kk) <- "GC-Chromatogram of TMS standard mixtures"
    image.path<-c(paste(file.path(path.package(package="MassOmics")),"/R/TMSstandardMix.gif", sep=""))
    FrontImage<-tcl("image",  "create", "photo", "MaxC2.image", file=image.path)
    tkgrid(tklabel(kk,image=FrontImage))
  }


  RunTMS<-function(save = TRUE,
                   output = "Final TMS Results",
                   Ion.bin= 0.5,
                   main.TMS.folder=main.TMS.folder,
                   RT.bin = 5) {

    require(xcms)
    require(lattice)
    require(plyr)
    setwd(main.TMS.folder)
    final.check.data <- read.csv(tk_choose.files(caption = "Load a corrected Identification Report (TMS).csv",  multi = FALSE),stringsAsFactors=FALSE)


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
      dev.new.OS()
      plotChrom(raw_data) # plot TIC

      for (h in 1:nrow(library_file)) {
        metabolite <-  library_file[h, ]
        R.Ion<-metabolite$ref_ion
        R.Time<-metabolite$RT.median*60
        IonExtract <-getEIC(raw_data,mzrange=cbind(R.Ion-Ion.bin,R.Ion+Ion.bin),rtrange=cbind(R.Time-RT.bin,R.Time+RT.bin))
        abundance<-data.frame(IonExtract@eic$xcmsRaw[1])
        maxabundance <- max(abundance$intensity)
        surefinal$Base.Peak[h] <- maxabundance
        Graphic.df<-rbind.fill(Graphic.df,data.frame(Metabolite.Names=paste(metabolite$Name," (m/z:",R.Ion,")",sep="") ,Retention.Time=(abundance$rt)/60,
                                                     Intensity=abundance$intensity ,Sample.Names=name.file ))

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
    final.df<-cbind( Num.ID=ion.lib$Number.of.ID, Ref.Ion=ion.lib$ref_ion, Ret.Time=ion.lib$RT.median, Quality =ion.lib$Quality ,final.df)
    final.df <- final.df[order(final.df$Ret.Time , decreasing = F),]
    row.names(final.df) <- 1:nrow(final.df)



    if (save == TRUE) {
      sheet <- paste(output," (RT bin=", RT.bin,"s)", sep="")
      store <- paste(main.TMS.folder, "\\", sheet, ".csv", sep = "")
      write.csv(final.df, file = store, row.names = FALSE)


      write.csv(Graphic.df,file="GC-Peak Diagnosis Report.csv")
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



  ## Building an interfance for TMS anlaysis




  library(tcltk)
  require(rgl)

  TMS <- tktoplevel()
  frameOverall <- tkframe(TMS)
  tkwm.title(TMS,"TMS analysis")


  workingdirctory<-function(){WorkingDir()}
  Id_report1<-function(){Id_report(main.TMS.folder=main.TMS.folder)}
  Tips1<-function(){Tips()}
  RunTMS1<-function(){
    RT.bin_s <-as.numeric(tclvalue(SliderValue))
    RunTMS(RT.bin = RT.bin_s,main.TMS.folder=main.TMS.folder)
  }
  merage_RT_standard1<-function(){merage_RT_standard(main.TMS.folder=main.TMS.folder)}
  TMS.hlep.page1<-function(){TMS.hlep.page()}


  ### all the button function
  Wdir.but <- tkbutton(TMS, text="Choose", command=workingdirctory,width=6)
  ID_report.but <- tkbutton(TMS, text="Run", command=Id_report1,width=6)
  Tips.but<-tkbutton(TMS, text="Tips", command=Tips1,width=6)
  integrate.but<-tkbutton(TMS, text="Run", command= RunTMS1,width=6)
  mergeSM.but<-tkbutton(TMS, text="Run", command= merage_RT_standard1,width=6)
  help.but<-tkbutton(TMS, text="Help", command= TMS.hlep.page1,width=6)




  ## Upper
  frameUpper <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  tkgrid(tklabel(frameUpper,text="Setup Working Directory"),Wdir.but, sticky="w")
  tkgrid(tklabel(frameUpper,text="Create an Identification Report "),ID_report.but,sticky="w")
  tkgrid(tklabel(frameUpper,text="Add StandardMix RT (option)"),mergeSM.but,sticky="w")

  ## Middle
  frameMid <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)


  ## Slider

  SliderValue <- tclVar("5")
  SliderValueLabel <- tklabel(TMS,text=as.character(tclvalue(SliderValue)),width=2)
  tkconfigure(SliderValueLabel,textvariable=SliderValue)
  slider <- tkscale(frameMid, from=1, to=10,
                    showvalue=F, variable=SliderValue,
                    resolution=1, orient="horizontal")

  tkgrid(tklabel(frameMid,text="Manually correct the Identification Report"),Tips.but,sticky="w")
  tkgrid(tklabel(frameMid,text="Automatic GC-peak integration"),integrate.but,sticky="w")
  tkgrid(tklabel(frameMid,text="Rentention time bin (Default = 5 s) "),slider,SliderValueLabel,tklabel(frameMid,text="s"),sticky="w")

  ## lower

  rb1 <- tkradiobutton(TMS)
  rb2 <- tkradiobutton(TMS)
  rbValue <- tclVar("PeakDiagnosis")
  tkconfigure(rb1,variable=rbValue,value="PeakDiagnosis")
  tkconfigure(rb2,variable=rbValue,value="IonExtractor_Window")

  OnOK <- function(){
    rbVal <- as.character(tclvalue(rbValue))
    if (rbVal=="PeakDiagnosis")
      ifelse(sum(dir()=="GC-Peak Diagnosis Report.csv")==1,PeakDiagnosis(),stop("Please select correct working directory"))
    if (rbVal=="IonExtractor_Window")
      IonExtractor_Window()
  }
  OK.but <- tkbutton(TMS,text="RUN",command=OnOK,width=6)

  frameLower <- tkframe(frameOverall,relief="groove",borderwidth=2,padx=5,pady=5)
  #tkgrid(tklabel(frameLower,text="Check TMS results"))
  tkgrid(tklabel(frameLower,text="Diagnose GC-Peak Auto Integration "),rb1, OK.but,sticky="w")
  tkgrid(tklabel(frameLower,text="IonExtractor Single Mode"),rb2,sticky="w")


  quit.but <- tkbutton(TMS, text = "Close Session",
                       command = function() {
                         tkdestroy(TMS)
                         rgl.close()
                       }
  )


  image.path1<-c(paste(file.path(path.package(package="MassOmics")),"/R/GCPM_SIlas_logo.gif", sep=""))
  logo<-tcl("image",  "create", "photo", "MaxC3.image", file=image.path1)




  button <- tkframe(frameOverall,relief="groove", borderwidth=0,padx=5,pady=5)
  tkgrid(tklabel(button ,text="                                         "),quit.but,tklabel(button ,text="                              "),help.but)


  ## Userinterphase

  tkgrid(tklabel(frameOverall,image=logo))
  tkgrid(tklabel(frameOverall,text="Step 1: Setup"))
  tkgrid(frameUpper,pady= 10, padx= 10)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="Step 2: Integration"))
  tkgrid(frameMid,pady= 10, padx= 30)
  #tkgrid(tklabel(frameOverall,text="    "))
  tkgrid(tklabel(frameOverall,text="Step 3: Check TMS results"))
  tkgrid(frameLower,pady= 10, padx= 10)
  tkgrid(frameOverall)
  #tkgrid(quit.but,help.but, pady= 10, padx= 10)
  tkgrid(button,pady= 10, padx= 10)
}





###############################MassOmics.GUI##################
#' MassOmics.GUI
#'
#' This function initiate a graphical user interface for MassOmics. To initiate the GG-MS data processing software, open the R console and type the codes in the example to initate the user interface.
#' @return None
#'
#' @examples
#' MassOmics.GUI()
#'
#' @export
MassOmics.GUI<-function(){
 library(tcltk2) 
  run()
  
}








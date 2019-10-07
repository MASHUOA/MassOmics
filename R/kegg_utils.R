######################################### addKeggCodes #########################################
addKeggCodes <- function (inputData, keggCodes, folder, save = TRUE, output = "Results with kegg codes",
                          addCodes = TRUE, formating_before_import=T,smpdb=load_smpdb())
{
  require(KEGGREST)
  require(svDialogs)
  require(svGUI)
  require(tcltk)
  require(tcltk2)
  isCSVdlg <- function(titleMSG, errorMSG) {
    t = 0
    while (t == 0) {
      checkIfCsv <- dlgOpen(title = titleMSG, multiple = FALSE)$res
      checkIfCsv2 <- basename(checkIfCsv)
      checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
      checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
      if (checkIfCsv4 %in% c("CSV", "csv")) {
        t = 1
        return(checkIfCsv)
      }      else {
        dlgMessage(errorMSG)
      }
    }
  }
  isCSV <- function(pathFile, errorMSG) {
    t = 0
    checkIfCsv <- pathFile
    checkIfCsv2 <- basename(checkIfCsv)
    checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
    checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
    if (checkIfCsv4 %in% c("CSV", "csv")) {
      t = 1
      return(t)
    }    else {
      return(t)
    }
  }
  OSsystem <- Sys.info()["sysname"]
  if (OSsystem == "Windows") {
    FolderDivisor <- "\\"
  }  else {
    FolderDivisor <- "/"
  }
  if (missing(inputData)) {
    inputData <- isCSVdlg("Select the CSV file containing the input data",
                          "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
    
    dataKegg <- read.csv(inputData, colClasses = "character")
    print("Input file loaded...")
  }  else {
    if (is.data.frame(inputData) == TRUE) {
      dataKegg <- inputData
      print("Data frame loaded...")
    }    else {
      if (is.character(inputData)) {
        checkIfCsv <- isCSV(inputData)
        if (checkIfCsv == 1) {
          inputTest <- file.access(inputData, 0)
          if (inputTest == 0) {
            dataKegg = read.csv(inputData, colClasses = "character")
          }          else {
            dlgMessage("The input file specified is not accessible. Please, choose a valid CSV file to be used as input data.")
            inputData <- isCSVdlg("Select the CSV file containing the input data",
                                  "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
            dataKegg <- read.csv(inputData, colClasses = "character")
            print("Input file loaded...")
          }
        }        else {
          dlgMessage("The input file specified is not in CSV format. Please, choose a valid CSV file to be used as input data.")
          inputData <- isCSVdlg("Select the CSV file containing the input data",
                                "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
          dataKegg <- read.csv(inputData, colClasses = "character")
          print("Input file loaded...")
        }
      }      else {
        dlgMessage("The path to the input data must be specified as character string. Please, choose a valid CSV file to be used as input data.")
        inputData <- isCSVdlg("Select the CSV file containing the input data",
                              "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
        dataKegg <- read.csv(inputData, colClasses = "character")
        print("Input file loaded...")
      }
    }
  }
  if (missing(keggCodes)) {
    keggCodes <- isCSVdlg("Select the CSV file containing the KEGG codes",
                          "The KEGG code library MUST be in the format of comma-separated value (csv). Please, choose a KEGG code library showing the extension .csv.")
    codesKegg <- read.csv(keggCodes, colClasses = "character")
    print("keggCodes - csv file loaded...")
  }  else {
    if (is.data.frame(keggCodes) == TRUE) {
      codesKegg <- keggCodes
      print("keggCodes - Data frame loaded...")
    }    else {
      if (is.character(keggCodes)) {
        checkIfCsv <- isCSV(keggCodes)
        if (checkIfCsv == 1) {
          inputTest <- file.access(keggCodes, 0)
          if (inputTest == 0) {
            codesKegg = read.csv(keggCodes, colClasses = "character")
          }          else {
            dlgMessage("The KEGG code library specified is not accessible. Please, choose a valid CSV file to be used as KEGG code library.")
            keggCodes <- isCSVdlg("Select the CSV file containing KEGG codes",
                                  "The KEGG code library MUST be in the format of comma-separated value (csv). Please, choose a KEGG code library showing the extension .csv.")
            codesKegg <- read.csv(keggCodes, colClasses = "character")
            print("keggCodes - csv file loaded...")
          }
        }        else {
          dlgMessage("The KEGG code library specified is not in CSV format. Please, choose a valid CSV file to be used as KEGG code library.")
          keggCodes <- isCSVdlg("Select the CSV file containing KEGG code library",
                                "The KEGG code library MUST be in the format of comma-separated value (csv). Please, choose a KEGG code library showing the extension .csv.")
          codesKegg <- read.csv(keggCodes, colClasses = "character")
          print("keggCodes - csv file loaded...")
        }
      }      else {
        dlgMessage("The path to the KEGG code library must be specified as character string. Please, choose a valid CSV file to be used as KEGG code library.")
        keggCodes <- isCSVdlg("Select the CSV file containing the KEGG code library",
                              "The KEGG code library MUST be in the format of comma-separated value (csv). Please, choose a KEGG code library showing the extension .csv.")
        codesKegg <- read.csv(keggCodes, colClasses = "character")
        print("keggCodes - csv file loaded...")
      }
    }
  }
  if (save == TRUE) {
    if (missing(folder)) {
      print("No folder was defined to save the results.")
      print("Please, point to the folder where the results should be saved.")
      folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
    }
    else {
      if (is.character(inputData)) {
        isFolder <- file.access(as.character(folder),
                                0)
        if (isFolder == 0) {
          isFolder <- file.info(folder)
          if (isFolder$isdir != TRUE) {
            print("The folder defined to save the results is not a valid path.")
            print("Please, point to the folder where the results should be saved.")
            folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
          }
        }
        else {
          print("The folder defined to save the results is not a valid path.")
          print("Please, point to the folder where the results should be saved.")
          folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
        }
      }
      else {
        print("The path to the folder where the results will be saved must be specified as character.")
        print("Please, point to the folder where the results should be saved.")
        folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
      }
    }
  }
  
  dataKegg=kegg_format(dataKegg)
  
  if (dataKegg[1, 1] %in% c("Replicates", "Replicate", "replicates",
                            "replicate")) {
    replicates <- as.character(dataKegg[1, ])
    dataKegg <- dataKegg[-1, ]
    rep <- 1
  }  else{
    rep <- 0
  }
  if (names(codesKegg)[1]=="COMPbase"){
    codesKegg=codesKegg[,1:3]
    names(codesKegg)=c("Name", "kegg",
                       "Pathway")
    
  }
  if (!is.null(dataKegg$Name)){
    dataKegg= dataKegg[,c("Name",colnames(dataKegg)[colnames(dataKegg)!="Name"])]
    
  }
  
  
  codesKegg=codesKegg[,c("kegg","Name", 
                         "Pathway")]
  original.data <- dataKegg
  original.kegglib <- codesKegg
  original.comps <- dataKegg
  
  names(dataKegg)[1] <- "Name"
  dataKegg[1] <- gsub(" ", "", dataKegg[, 1])
  #names(codesKegg)[c(1, 2)] <- c("kegg", "Name")
  #names(codesKegg)[c(1, 2)] <- c("Name", "kegg")
  #codesKegg[2] <- gsub(" ", "", codesKegg[, 2])
  codesKegg[["Name"]] <- gsub(" ", "", codesKegg[["Name"]])
  if (addCodes == TRUE) {
    pre.data <- merge(dataKegg, codesKegg, by.x = "Name",
                      by.y = "Name", all.x = TRUE)
    missingCpd <- pre.data[is.na(pre.data$kegg), ]
    #uniquesmissingcpd=NULL
    missingCpdname=unique(missingCpd[,1])
    missingCpd<-NULL
    missingCpd$Name=missingCpdname
    missingCpd=data.frame(missingCpd,stringsAsFactors = F)
    #missingCpd <- missingCpd[missingCpd$Name %in% unique(missingCpd$Name),]
    #missingCpd <- unique(missingCpd[,1])
    if (nrow(missingCpd) == 0) {
      print("Every compound in the inputData was found in the keggCodes library.")
    } else {
      message(paste("Found new cpd, ",nrow(missingCpd),"in total\n",paste(head(original.comps[which(gsub(" ", "", original.comps[,
                                                              1]) %in% missingCpd[, 1]), ]),collapse = "\n")))
      addCPD <- dlgMessage("Some compounds of your input data are not listed in the kegg code library used. Some of these compounds are listed in the R console. Would you like to include them into your kegg library right now?",
                           type = "yesno")$res
      if (addCPD == "yes") {
        if (OSsystem == "Windows") {
        } else {
          dlgMessage("This is  is going to happen now: KEGG database will be searched for the missing compounds. For each missing compound, a list of potential matches will be presented to you. If the missing compound is part of the list, select it and its respective KEGG code will be added to your KEGG code library. If the desired compound is not listed, you can go to the end of the list and click on OTHER to manually add its KEGG code, or you can click on SKIP if you have no KEGG code for this compound. Compounds showing no KEGG codes will not be analyzed by PAPi.")
        }
        autocpdChosen=F
        
        for (i in 1:nrow(missingCpd)) {
          
          compToSearch = original.comps[which(gsub(" ","", original.comps[, 1]) == missingCpd[i,1]), "CAS"]
          compToSearch =unique(compToSearch)
          #compToSearch=gsub("-","",compToSearch)
          compToSearch=gsub("-","",compToSearch)
          compToSearch=gsub("\\(+-\\)","",compToSearch)
          compToSearch=gsub("\\(+\\)","",compToSearch)
          compToSearch=gsub("\\(-\\)","",compToSearch)
          compToSearch=gsub("^-","",compToSearch)
          #compToSearch=gsub(")","",compToSearch)
          potentialCPDS=smpdb$CAS[which(smpdb$CAS==compToSearch)]
          potentialCPDS=character(0)
          
          potentialCPDS <- try(keggFind("compound", as.character(compToSearch)))
          
          if (length(potentialCPDS)==0 | class(potentialCPDS)=="try-error"){
            possible_cas=retieve_cpd(compToSearch)
            if(!is.na(possible_cas)){
             potentialCPDS=try(keggFind("compound",possible_cas[[1]][1]))
             
            }
            
          }
          
          if (length(potentialCPDS)!=0 && class(potentialCPDS)!="try-error"){
          if (autocpdChosen){
            #Sys.sleep(runif(1,0.01, 0.02))
            cpdChosen <- potentialCPDS[1]
            cpdChosen <- names(cpdChosen)
            cpdChosen <- gsub("cpd:", "", cpdChosen)
            
            original.kegglib <- rbind(original.kegglib,
                                      c(kegg=as.character(compToSearch), Name=cpdChosen,Pathways=getPathways(cpdChosen)))
            
          }else{
          
          cpdChosen <- tk_select.list(c("AUOT-SELECT TOP HIT FOR THE REST NEW CPDS", "OTHER",
                                        "SKIP", potentialCPDS),preselect="AUOT-SELECT TOP HIT FOR THE REST NEW CPDS", title = as.character(compToSearch))
          if (cpdChosen == "OTHER") {
            cpdChosen <- dlgInput(paste("Enter the KEGG code for ",
                                        as.character(compToSearch), " or leave it as absent.",
                                        sep = ""), default = "absent")$res
            if (!length(cpdChosen)) {
              cpdChosen <- "absent"
            }
            original.kegglib <- rbind(original.kegglib,
                                      c(cpdChosen, as.character(compToSearch)))
          } else {
            if (cpdChosen == "SKIP") {
            }else if (cpdChosen == "AUOT-SELECT TOP HIT FOR THE REST NEW CPDS"){
              autocpdChosen=T
              
            } else {
              cpdChosen <- potentialCPDS[which(potentialCPDS ==
                                                 cpdChosen)]
              cpdChosen <- names(cpdChosen)
              cpdChosen <- gsub("cpd:", "", cpdChosen)
              
              original.kegglib <- rbind(original.kegglib,
                                        c(kegg=as.character(compToSearch), Name=cpdChosen,Pathways=getPathways(cpdChosen)))
            }
          }
        }
        }
      }
          

        saveLib <- dlgMessage("Would you like to save your new KEGG library to a CSV file?",
                              type = c("yesno"))$res
        if (saveLib == "yes") {
          placeToSave <- dlgDir(title = "Select the folder where the library will be saved.")$res
          dateForFile <- Sys.time()
          dateForFile <- gsub(" ", "", dateForFile)
          dateForFile <- gsub(":", "", dateForFile)
          fileName <- dlgInput(message = "What is the name of the new library?",
                               default = paste("KEGGLibrary", dateForFile,
                                               sep = ""))$res
          fileName <- paste(fileName, ".csv", sep = "")
          placeToSave <- paste(placeToSave, FolderDivisor,
                               fileName, sep = "")
          write.csv(original.kegglib, file = placeToSave,
                    row.names = FALSE)
          print("Your new KEGG code library has been saved.")
        }
        codesKegg <- original.kegglib
      }
    }
  }

  codesKegg<-codesKegg[,c("kegg", "Name","Pathway")]
  codesKegg[1] <- gsub("cpd:", "", codesKegg[, 1])
  codesKegg[2] <- gsub(" ", "", codesKegg[, 2])
  FinalData <- merge(dataKegg, codesKegg, by.x = "Name", by.y = "Name",
                     all.x = TRUE)
  FinalData <- subset(FinalData, !(FinalData$kegg == "absent"))
  FinalData$Name <- FinalData$kegg
  FinalData$kegg <- NULL
  if (rep == 1) {
    FinalData <- rbind(c(replicates), FinalData)
  }
  if (save == TRUE) {
    sheet <- output
    store <- paste(folder, FolderDivisor, sheet, ".csv",
                   sep = "")
    write.csv(FinalData, file = store, row.names = FALSE)
    print(paste("The file ", output, ".csv", " was saved in the folder ",
                folder, sep = ""))
  }
  else {
    print("The final data frame was not saved because the argument save was set as FALSE")
  }
  return(FinalData)
}

######################################### build KEGG Database #########################################
buildDatabase_old <- function (save = FALSE, folder, saveAs="default")
{ 
  require(KEGGREST)
  require(svDialogs)
  require(svGUI)
  library(tcltk2)
  checkAdminstration <- tk_select.list(c("Yes (continue)",  "No (stop)"),  title = "Run R as administrator")
  
  
  if (checkAdminstration  == "No (Stop)"){
    stop("KEGG library won't able to save in you computer")
  }
  
  
  testInternet <- try(keggList("pathway"), TRUE)
  if ("try-error" %in% class(testInternet)) {
    stop("We could not connect to KEGG database. You probably have no internet connection.")
  }
  OSsystem <- Sys.info()["sysname"]
  if (OSsystem == "Windows") {
    FolderDivisor <- "\\"
  }  else {
    FolderDivisor <- "/"
  }
  apply_pb <- function(X, MARGIN, FUN, ...) {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    wrapper <- function(...) {
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1, envir = env)
      setTxtProgressBar(get("pb", envir = env), curVal +
                          1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  NameToSave <- function() {
    dateForFile <- Sys.time()
    dateForFile <- gsub(" ", "", dateForFile)
    dateForFile <- gsub(":", "", dateForFile)
    NameFile <- dlgInput(message = "What is the name of the new database?",
                         default = saveAs)
    NameFile=(NameFile$res)
    NameFile <- windows_filename(NameFile)
    return(NameFile)
  }
  error1 <- "The folder defined to save the database is not a valid path."
  error2 <- "Please, point to the folder where the database should be saved."
  titleDLG <- "Select the folder where the new database will be saved."
  if (save == TRUE) {
    if (missing(folder)) {
      print(error1)
      print(error2)
      folder = dlgDir(title = titleDLG)$res
      if (missing(saveAs)) {
        fileName <- NameToSave()
      } else {
        fileName <- saveAs
      }
    }else {
      if (is.character(folder)) {
        isFolder <- file.access(as.character(folder), 0)
        if (isFolder == 0) {
          isFolder <- file.info(folder)
          if (isFolder$isdir != TRUE) {
            print(error1)
            print(error2)
            folder = dlgDir(title = titleDLG)$res
            if (missing(saveAs)) {
              fileName <- NameToSave()
            }
            else {
              fileName <- saveAs
            }
          }
        } else {
          print(error1)
          print(error2)
          folder = dlgDir(title = titleDLG)$res
          if (missing(saveAs)) {
            fileName <- NameToSave()
          } else {
            fileName <- saveAs
          }
        }
      }
      else {
        print("The path to the folder where the results will be saved must be specified as character.")
        print(error2)
        folder = dlgDir(title = titleDLG)$res
      }
    }
  }  else {
    if (missing(saveAs)) {
      fileName <- NameToSave()
    }
    else {
      fileName <- saveAs
    }
  }

  listOfComps <- keggList("compound")
  listOfComps <- data.frame(listOfComps)
  listOfComps[2] <- row.names(listOfComps)
  row.names(listOfComps) <- 1:nrow(listOfComps)
  #BPPARAM=bpparam()
  #bpprogressbar(BPPARAM)=T
  
  #listOfCompstest <- bplapply(as.list(t(listOfComps[2])),getPathways,BPPARAM = BPPARAM)
  #listOfComps[3] <-unlist(listOfCompstest)
  listOfComps[3] <- apply_pb(listOfComps[2],1,getPathways)
  FileName <- "COMPbase"
  dateForFile <- Sys.time()
  dateForFile <- gsub(" ", "", dateForFile)
  dateForFile <- gsub(":", "", dateForFile)
  userName <- Sys.info()
  userName <- userName[[7]]
  names(listOfComps)[c(1, 2, 3)] <- c("Name", "kegg",
                                      "Pathway")
  numbcomp <- function(x) {
    compounds <- keggGet(x)
    compounds <- data.frame(compounds[[1]]$COMPOUND)
    compounds <- nrow(compounds)
    if (compounds == 0) {
      x <- gsub("path:map", "ko", x, fixed = TRUE)
      compounds <- try(keggGet(x), silent = TRUE)
      if ("try-error" %in% class(compounds)) {
        compounds <- 0
      }
      else {
        compounds <- data.frame(compounds[[1]]$COMPOUND)
        compounds <- nrow(compounds)
      }
    }
    return(compounds)
  }
  pathname <- data.frame(keggList("pathway"))
  pathname[2] <- row.names(pathname)
  row.names(pathname) <- 1:nrow(pathname)
  pathname[3] <- apply_pb(pathname[2], 1, function(x) numbcomp(x))
  FileName <- "PATHbase"
  dateForFile <- Sys.time()
  dateForFile <- gsub(" ", "", dateForFile)
  dateForFile <- gsub(":", "", dateForFile)
  userName <- Sys.info()
  userName <- userName[[7]]
  names(pathname)[c(1, 2, 3)] <- c(FileName, dateForFile, userName)
  dir.create(paste(R.home("library/PAPi/databases"), fileName,
                   sep = FolderDivisor), recursive = TRUE)
  placeToSave <- paste(R.home("library/PAPi/databases"), fileName,
                       "/COMPbase.csv", sep = FolderDivisor)
  write.csv(listOfComps, file = placeToSave, row.names = FALSE)
  placeToSave <- paste(R.home("library/PAPi/databases"), fileName,
                       "/PATHbase.csv", sep = FolderDivisor)
  write.csv(pathname, file = placeToSave, row.names = FALSE)
  print("Your new KEGG database has been installed.")
  if (save == TRUE) {
    dir.create(paste(folder, fileName, sep = FolderDivisor))
    placeToSave <- paste(folder, fileName,
                         "COMPbase.csv", sep =FolderDivisor)
    write.csv(listOfComps, file = placeToSave, row.names = FALSE)
    placeToSave <- paste(folder, fileName, 
                         "PATHbase.csv", sep = FolderDivisor)
    write.csv(pathname, file = placeToSave, row.names = FALSE)
    print("Your new KEGG database has been saved.")
  }
}


getPathways <- function(x) {
  TotalReport <- KEGGREST::keggGet(x)
  pathways <- TotalReport[[1]]$PATHWAY
  pathways <- data.frame(pathways)
  pathways <- row.names(pathways)
  pathways <- paste(pathways, collapse = ";")
  return(pathways)
}





######################################### papi #########################################
papi <- function (inputData, save = TRUE, folder, output = "Papi results",
                  offline = TRUE, localDatabase = "default")
{
  require(KEGGREST)
  require(svDialogs)
  require(svGUI)
  require(tcltk)
  require(tcltk2)
  
  isCSVdlg <- function(titleMSG, errorMSG) {
    t = 0
    while (t == 0) {
      checkIfCsv <- dlgOpen(title = titleMSG, multiple = FALSE)$res
      checkIfCsv2 <- basename(checkIfCsv)
      checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
      checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
      if (checkIfCsv4 %in% c("CSV", "csv")) {
        t = 1
        return(checkIfCsv)
      }
      else {
        dlgMessage(errorMSG)
      }
    }
  }
  isCSV <- function(pathFile, errorMSG) {
    t = 0
    checkIfCsv <- pathFile
    checkIfCsv2 <- basename(checkIfCsv)
    checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
    checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
    if (checkIfCsv4 %in% c("CSV", "csv")) {
      t = 1
      return(t)
    }
    else {
      return(t)
    }
  }
  OSsystem <- Sys.info()["sysname"]
  if (OSsystem == "Windows") {
    FolderDivisor <- "\\"
  }
  else {
    FolderDivisor <- "/"
  }
  online <- 0
  if (offline == TRUE && localDatabase == "default") {
    listOfComps <- try(read.csv(R.home("library/PAPi/databases/default/COMPbase.csv"),
                                colClasses = "character"), TRUE)
    if ("try-error" %in% class(listOfComps)) {
      print(class(listOfComps))
      dlgMessage("The compound's database could not be loaded. We will look for additional databases to be used.")$res
      listOfDatabases <- try(list.files(R.home("library/PAPi/databases/")),
                             TRUE)
      if ("try-error" %in% class(listOfDatabases)) {
        print(class(listOfDatabases))
        decision1 <- dlgMessage("We could not search for additional databases. The related error is printed on your screen. Press OK to apply PAPi online or CANCEL to stop the function.",
                                type = "okcancel")$res
        if (decision1 == "cancel") {
          stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
        }
        else {
          online <- 1
        }
      }
      else {
        baseToUse <- tk_select.list(listOfDatabases)
        baseToUse <- paste(R.home("library/PAPi/databases/"),
                           baseToUse, sep = "")
        listOfComps <- try(read.csv(paste(baseToUse,
                                          "/", "COMPbase.csv", sep = ""), colClasses = "character"),
                           TRUE)
        if ("try-error" %in% class(listOfComps)) {
          print(class(listOfComps))
          decision2 <- dlgMessage("The compound's database could not be loaded. Would you like to do this analysis online ?.",
                                  type = "okcancel")$res
          if (decision2 == "cancel") {
            stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
          }
          else {
            online <- 1
          }
        }
        else {
          pathname <- try(read.csv(paste(baseToUse, "/",
                                         "PATHbase.csv", sep = ""), colClasses = "character"),
                          TRUE)
          if ("try-error" %in% class(pathname)) {
            print(class(pathname))
            decision3 <- dlgMessage("The pathways' database could not be loaded. Would you like to do this analysis online ?.",
                                    type = "okcancel")$res
            if (decision3 == "cancel") {
              stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
            }
            else {
              online <- 1
            }
          }
          pathname[2] <- gsub("path:map", "ko", pathname[,
                                                         2])
        }
      }
    }
    else {
      pathname <- try(read.csv(R.home("library/PAPi/databases/default/PATHbase.csv"),
                               colClasses = "character"), TRUE)
      if ("try-error" %in% class(pathname)) {
        print(class(listOfComps))
        dlgMessage("The pathways' database could not be loaded. We will look for additional databases to be used.")$res
        listOfDatabases <- try(list.files(R.home("library/PAPi/databases/")),
                               TRUE)
        if ("try-error" %in% class(listOfDatabases)) {
          print(class(listOfDatabases))
          decision1 <- dlgMessage("We could not search for additional databases. The related error is printed on your screen. Press OK to apply PAPi online or CANCEL to stop the function.",
                                  type = "okcancel")$res
          if (decision1 == "cancel") {
            stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
          }
          else {
            online <- 1
          }
        }
        else {
          baseToUse <- tk_select.list(listOfDatabases)
          baseToUse <- paste(R.home("library/PAPi/databases/"),
                             baseToUse, sep = "")
          listOfComps <- try(read.csv(paste(baseToUse,
                                            "/", "COMPbase.csv", sep = ""), colClasses = "character"),
                             TRUE)
          if ("try-error" %in% class(listOfComps)) {
            print(class(listOfComps))
            decision2 <- dlgMessage("The compound's database could not be loaded. Would you like to do this analysis online ?.",
                                    type = "okcancel")$res
            if (decision2 == "cancel") {
              stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
            }
            else {
              online <- 1
            }
          }
          else {
            pathname <- try(read.csv(paste(baseToUse,
                                           "/", "PATHbase.csv", sep = ""), colClasses = "character"),
                            TRUE)
            if ("try-error" %in% class(pathname)) {
              print(class(pathname))
              decision3 <- dlgMessage("The pathways' database could not be loaded. Would you like to do this analysis online ?.",
                                      type = "okcancel")$res
              if (decision3 == "cancel") {
                stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
              }
              else {
                online <- 1
              }
            }
            pathname[2] <- gsub("path:map", "ko", pathname[,
                                                           2])
          }
        }
      }
      pathname[2] <- gsub("path:map", "ko", pathname[,
                                                     2])
    }
  }
  if (offline == TRUE && localDatabase == "choose") {
    listOfDatabases <- try(list.files(R.home("library/PAPi/databases/")),
                           TRUE)
    if ("try-error" %in% class(listOfDatabases)) {
      print(class(listOfDatabases))
      decision1 <- dlgMessage("We could not search for databases. The related error is printed on your screen. Press OK to apply PAPi online or CANCEL to stop the function.",
                              type = "okcancel")$res
      if (decision1 == "cancel") {
        stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
      }
      else {
        online <- 1
      }
    }
    else {
      baseToUse <- tk_select.list(listOfDatabases, title = "Which database do you want to use?")
      baseToUse <- paste(R.home("library/PAPi/databases/"),
                         baseToUse, sep = "")
      listOfComps <- try(read.csv(paste(baseToUse, "/",
                                        "COMPbase.csv", sep = ""), colClasses = "character"),
                         TRUE)
      if ("try-error" %in% class(listOfComps)) {
        print(class(listOfComps))
        decision2 <- dlgMessage("The compound's database could not be loaded. Would you like to do this analysis online ?.",
                                type = "okcancel")$res
        if (decision2 == "cancel") {
          stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
        }
        else {
          online <- 1
        }
      }
      else {
        pathname <- try(read.csv(paste(baseToUse, "/",
                                       "PATHbase.csv", sep = ""), colClasses = "character"),
                        TRUE)
        if ("try-error" %in% class(pathname)) {
          print(class(pathname))
          decision3 <- dlgMessage("The pathways' database could not be loaded. Would you like to do this analysis online ?.",
                                  type = "okcancel")$res
          if (decision3 == "cancel") {
            stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
          }
          else {
            online <- 1
          }
        }
        pathname[2] <- gsub("path:map", "ko", pathname[,
                                                       2])
      }
    }
  }
  else {
    if (offline == TRUE) {
      listOfDatabases <- try(list.files(R.home("library/PAPi/databases/")),
                             TRUE)
      if (localDatabase %in% listOfDatabases) {
        baseToUse <- paste(R.home("library/PAPi/databases/"),
                           localDatabase, sep = "")
        listOfComps <- try(read.csv(paste(baseToUse,
                                          "/", "COMPbase.csv", sep = ""), colClasses = "character"),
                           TRUE)
        if ("try-error" %in% class(listOfComps)) {
          print(class(listOfComps))
          decision2 <- dlgMessage("The compound's database could not be loaded. Would you like to do this analysis online ?.",
                                  type = "okcancel")$res
          if (decision2 == "cancel") {
            stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
          }
          else {
            online <- 1
          }
        }
        else {
          pathname <- try(read.csv(paste(baseToUse, "/",
                                         "PATHbase.csv", sep = ""), colClasses = "character"),
                          TRUE)
          if ("try-error" %in% class(pathname)) {
            print(class(pathname))
            decision3 <- dlgMessage("The pathways' database could not be loaded. Would you like to do this analysis online ?.",
                                    type = "okcancel")$res
            if (decision3 == "cancel") {
              stop("The process was stopped by the user. See ?buildDatabase for creating a local database.")
            }
            else {
              online <- 1
            }
          }
          pathname[2] <- gsub("path:map", "ko", pathname[,
                                                         2])
        }
      }
      else {
        stop("The selected database is not installed. Try to reapply the function using localDatabase = \"choose\". See ?buildDatabase for creating a local database.")
      }
    }
  }
  if (offline == FALSE) {
    testInternet <- try(keggList("pathway"), TRUE)
    if ("try-error" %in% class(testInternet)) {
      stop("We could not connect to KEGG database. You probably have no internet connection. Try to use a local database by reapplying papi using localDatabase = \"choose\".")
    }
    online <- 1
  }
  if (online == 1) {
    offline <- FALSE
    print("Online analysis being performed...")
  }
  else {
    print("Offline analysis being performed...")
  }
  print("PAPi in progess...")
  if (missing(inputData)) {
    inputData <- isCSVdlg("Select the CSV file containing the input data",
                          "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
    omics.data.frame <- read.csv(inputData, colClasses = "character")
    print("Input file loaded...")
  }
  else {
    if (is.data.frame(inputData) == TRUE) {
      omics.data.frame <- inputData
      print("Data frame loaded...")
    }
    else {
      if (is.character(inputData)) {
        checkIfCsv <- isCSV(inputData)
        if (checkIfCsv == 1) {
          inputTest <- file.access(inputData, 0)
          if (inputTest == 0) {
            omics.data.frame = read.csv(inputData, colClasses = "character")
          }
          else {
            print("The input file specified is not accessible. Please, choose a valid CSV file to be used as input data.")
            inputData <- isCSVdlg("Select the CSV file containing the input data",
                                  "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
            omics.data.frame <- read.csv(inputData, colClasses = "character")
            print("Input file loaded...")
          }
        }
        else {
          print("The input file specified is not in CSV format. Please, choose a valid CSV file to be used as input data.")
          inputData <- isCSVdlg("Select the CSV file containing the input data",
                                "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
          omics.data.frame <- read.csv(inputData, colClasses = "character")
          print("Input file loaded...")
        }
      }
      else {
        print("The path to the input data must be specified as character string. Please, choose a valid CSV file to be used as input data.")
        inputData <- isCSVdlg("Select the CSV file containing the input data",
                              "The input file MUST be in the format of comma-separated value (csv). Please, choose an input file showing the extension .csv.")
        omics.data.frame <- read.csv(inputData, colClasses = "character")
        print("Input file loaded...")
      }
    }
  }
  if (save == TRUE) {
    if (missing(folder)) {
      print("No folder was defined to save the results.")
      print("Please, point to the folder where the results should be saved.")
      folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
    }
    else {
      if (is.character(inputData)) {
        isFolder <- file.access(as.character(folder),
                                0)
        if (isFolder == 0) {
          isFolder <- file.info(folder)
          if (isFolder$isdir != TRUE) {
            print("The folder defined to save the results is not a valid path.")
            print("Please, point to the folder where the results should be saved.")
            folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
          }
        }
        else {
          print("The folder defined to save the results is not a valid path.")
          print("Please, point to the folder where the results should be saved.")
          folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
        }
      }
      else {
        print("The path to the folder where the results will be saved must be specified as character.")
        print("Please, point to the folder where the results should be saved.")
        folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
      }
    }
  }
  if (omics.data.frame[1, 1] %in% c("Replicates", "Replicate",
                                    "replicates", "replicate")) {
    replicates <- as.character(omics.data.frame[1, ])
    omics.data.frame <- omics.data.frame[-1, ]
    rep <- 1
    reps <- factor(replicates[-1])
  }
  else {
    rep <- 0
  }
  if (online == 1) {
    getcomp <- function(x) {
      comp <- try(keggGet(x), TRUE)
      if ("try-error" %in% class(comp)) {
        comp <- 0
      }
      else {
        comp <- comp[[1]]$PATHWAY
        comp <- data.frame(comp)
        names(comp)[1] <- "value"
      }
      return(comp)
    }
    numbcomp <- function(x) {
      compounds <- try(keggGet(x), TRUE)
      if ("try-error" %in% class(compounds)) {
        compounds <- 0
      }
      else {
        compounds <- data.frame(compounds[[1]]$COMPOUND)
        names(compounds)[1] <- "value"
        compounds <- nrow(compounds)
      }
      return(compounds)
    }
    pathname <- data.frame(keggList("pathway"))
    pathname[2] <- row.names(pathname)
    row.names(pathname) <- 1:nrow(pathname)
    names(pathname)[c(1, 2)] <- c("pathwayname", "idpathway")
    pathname[2] <- gsub("path:map", "ko", pathname[, 2])
  }
  else {
    getcomp <- function(x) {
      comp <- which(listOfComps[2] == x)
      if (length(comp) > 0) {
        comp <- listOfComps[comp[1], 3]
        comp <- gsub("map", "ko", comp)
        comp <- data.frame(unlist(strsplit(comp, ";")))
        pathwaysFound <- which(pathname[, 2] %in% comp[,
                                                       1])
        if (length(pathwaysFound) > 0) {
          pathwaysFound <- pathname[pathwaysFound, ]
          row.names(pathwaysFound) <- pathwaysFound[,
                                                    2]
          pathwaysFound <- pathwaysFound[1]
          names(pathwaysFound)[1] <- "value"
        }
        else {
          pathwaysFound <- 0
        }
      }
      else {
        pathwaysFound <- 0
      }
      return(pathwaysFound)
    }
    numbcomp <- function(x) {
      compounds <- which(pathname[2] == x)
      if (length(compounds) > 0) {
        compounds <- as.numeric(pathname[compounds[1],
                                         3])
      }
      else {
        compounds <- 0
      }
      return(compounds)
    }
  }
  pb <- txtProgressBar(min = 0, max = ncol(omics.data.frame),
                       style = 3, width = 50)
  papi.frame <- numeric()
  for (j in 2:ncol(omics.data.frame)) {
    Sys.sleep(0.01)
    setTxtProgressBar(pb, j)
    data.df <- omics.data.frame[, c(1, j)]
    data.df <- subset(data.df, !is.na(data.df[2]))
    complist <- as.character(data.df[, 1])
    complist <- gsub("C", "cpd:C", complist, ignore.case = TRUE)
    getpath.final <- numeric()
    for (i in 1:length(complist)) {
      getpath <- getcomp(complist[i])
      if (getpath != 0) {
        if (length(getpath.final) == 0) {
          getpath$rate <- as.numeric(data.df[i, 2])
          getpath.final <- getpath
        }
        else {
          getpath$rate <- as.numeric(data.df[i, 2])
          getpath.final <- rbind(getpath, getpath.final)
        }
      }
      else {
      }
    }
    getpath.final$rate <- as.numeric(as.character(getpath.final$rate))
    res.arr <- with(getpath.final, tapply(rate, value, sum))
    res.df <- data.frame(pathwayname = names(res.arr), rate = res.arr)
    rownames(res.df) <- NULL
    pathwayfreq <- data.frame(table(getpath.final$value))
    names(pathwayfreq)[1] <- "pathwayname"
    names(pathname)[c(1, 2)] <- c("pathwayname", "idpathway")
    freqname <- merge(pathwayfreq, pathname)
    freqname$idpathway <- gsub("path:map", "ko", freqname$idpathway)
    total <- nrow(pathwayfreq)
    selecrow.final <- numeric()
    for (i in 1:nrow(res.df)) {
      selecrow <- freqname[i, ]
      numberCompounds <- try(numbcomp(selecrow$idpathway),
                             silent = TRUE)
      if ("try-error" %in% class(numberCompounds)) {
        selecrow$percentage <- 1
      }
      else {
        if (numberCompounds != 0) {
          selecrow$percentage <- (selecrow$Freq/numberCompounds)
        }
        else {
          selecrow$percentage <- 1
        }
        if (i == 1) {
          selecrow.final <- selecrow
        }
        else {
          selecrow.final <- rbind(selecrow.final, selecrow)
        }
      }
    }
    if (online == 1) {
      selecrow.final <- selecrow.final[, -c(2, 3)]
    }
    else {
      selecrow.final <- selecrow.final[, -c(2, 3, 4)]
    }
    final2.df <- merge(selecrow.final, res.df)
    final2.df[2] <- final2.df$rate/final2.df$percentage
    names(final2.df)[2] <- names(data.df)[2]
    final2.df[3] <- NULL
    if (j == 2) {
      papi.frame <- final2.df
    }
    else {
      papi.frame <- merge(papi.frame, final2.df, all = TRUE)
    }
  }
  close(pb)
  papi.frame <- papi.frame[order(papi.frame[, 2], decreasing = T),
                           ]
  papi.frame[1] <- as.character(papi.frame[, 1])
  if (rep == 1) {
    papi.frame <- rbind(c(replicates), papi.frame)
  }
  if (save == TRUE) {
    sheet <- output
    store <- paste(folder, FolderDivisor, sheet, ".csv",
                   sep = "")
    write.csv(papi.frame, file = store, row.names = TRUE)
    print(paste("The file ", output, ".csv", " was saved in the folder ",
                folder, sep = ""))
  }
  else {
    print("No file was saved because the argument save was set as FALSE")
  }
  return(papi.frame)
}


######################################### format the integration result table for pathway analysis#######
kegg_format<-function(inputData,save=T){
  
  library(data.table)
  library(stringr)
  
  if (class(inputData)=="character"){
   inputData.df=fread(inputData) 
  }else if("data.frame" %in% class(inputData)) {
    inputData.df=(inputData)
  }
  
  
  inputData.df[,1]=gsub("\\(split peak .\\)","",inputData.df[,1])
  
  inputData.df[,1]=gsub(" ","",inputData.df[,1])
  
  replicates=data.frame(t(data.frame(as.character(colnames(inputData.df)),stringsAsFactors = F)),stringsAsFactors = F)
  
  
  
  colnames(replicates)=replicates[1,]
  
  replicates[1,1]="Replicates"
  
  inputData.df=rbind(replicates,inputData.df)
  
  rownames(inputData.df)=1:nrow(inputData.df)
  
  if (save && class(inputData)=="character"){
    
    write.csv(inputData.df,paste0(gsub(".csv$","",inputData),"_pathway_formated.csv"))
    
  }
  
  return(inputData.df)
  
}

retieve_cpd<-function(x){
 library("webchem") 
  
  cir_query(x,"names")
}

######################################### build KEGG Database #########################################
buildDatabase_old_2 <- function (save = FALSE, folder, saveAs="default")
{ 
  require(KEGGREST)
  require(svDialogs)
  require(svGUI)
  library(tcltk2)
  #library(RbioRXN)
  checkAdminstration <- tk_select.list(c("Yes (continue)",  "No (stop)"),  title = "Run R as administrator")
  
  
  if (checkAdminstration  == "No (Stop)"){
    stop("KEGG library won't able to save in you computer")
  }
  
  
  testInternet <- try(keggList("pathway"), TRUE)
  if ("try-error" %in% class(testInternet)) {
    stop("We could not connect to KEGG database. You probably have no internet connection.")
  }
  OSsystem <- Sys.info()["sysname"]
  if (OSsystem == "Windows") {
    FolderDivisor <- "\\"
  }  else {
    FolderDivisor <- "/"
  }
  apply_pb <- function(X, MARGIN, FUN, ...) {
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
    wrapper <- function(...) {
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1, envir = env)
      setTxtProgressBar(get("pb", envir = env), curVal +
                          1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  NameToSave <- function() {
    dateForFile <- Sys.time()
    dateForFile <- gsub(" ", "", dateForFile)
    dateForFile <- gsub(":", "", dateForFile)
    NameFile <- dlgInput(message = "What is the name of the new database?",
                         default = saveAs)
    NameFile=(NameFile$res)
    NameFile <- windows_filename(NameFile)
    return(NameFile)
  }
  error1 <- "The folder defined to save the database is not a valid path."
  error2 <- "Please, point to the folder where the database should be saved."
  titleDLG <- "Select the folder where the new database will be saved."
  if (save == TRUE) {
    if (missing(folder)) {
      print(error1)
      print(error2)
      folder = dlgDir(title = titleDLG)$res
      if (missing(saveAs)) {
        fileName <- NameToSave()
      } else {
        fileName <- saveAs
      }
    }else {
      if (is.character(folder)) {
        isFolder <- file.access(as.character(folder), 0)
        if (isFolder == 0) {
          isFolder <- file.info(folder)
          if (isFolder$isdir != TRUE) {
            print(error1)
            print(error2)
            folder = dlgDir(title = titleDLG)$res
            if (missing(saveAs)) {
              fileName <- NameToSave()
            }
            else {
              fileName <- saveAs
            }
          }
        } else {
          print(error1)
          print(error2)
          folder = dlgDir(title = titleDLG)$res
          if (missing(saveAs)) {
            fileName <- NameToSave()
          } else {
            fileName <- saveAs
          }
        }
      }
      else {
        print("The path to the folder where the results will be saved must be specified as character.")
        print(error2)
        folder = dlgDir(title = titleDLG)$res
      }
    }
  }  else {
    if (missing(saveAs)) {
      fileName <- NameToSave()
    }
    else {
      fileName <- saveAs
    }
  }
  
  listOfComps <- keggList("compound")
  listOfComps <- data.frame(listOfComps)
  listOfComps[2] <- row.names(listOfComps)
  row.names(listOfComps) <- 1:nrow(listOfComps)
  #BPPARAM=bpparam()
  #bpprogressbar(BPPARAM)=T
  
  #listOfCompstest <- bplapply(as.list(t(listOfComps[2])),getPathways,BPPARAM = BPPARAM)
  #listOfComps[3] <-unlist(listOfCompstest)
  listOfComps[3] <- apply_pb(listOfComps[2],1,getPathways)
  FileName <- "COMPbase"
  dateForFile <- Sys.time()
  dateForFile <- gsub(" ", "", dateForFile)
  dateForFile <- gsub(":", "", dateForFile)
  userName <- Sys.info()
  userName <- userName[[7]]
  names(listOfComps)[c(1, 2, 3)] <- c("Name", "kegg",
                                      "Pathway")
  numbcomp <- function(x) {
    compounds <- keggGet(x)
    compounds <- data.frame(compounds[[1]]$COMPOUND)
    compounds <- nrow(compounds)
    if (compounds == 0) {
      x <- gsub("path:map", "ko", x, fixed = TRUE)
      compounds <- try(keggGet(x), silent = TRUE)
      if ("try-error" %in% class(compounds)) {
        compounds <- 0
      }
      else {
        compounds <- data.frame(compounds[[1]]$COMPOUND)
        compounds <- nrow(compounds)
      }
    }
    return(compounds)
  }
  pathname <- data.frame(keggList("pathway"))
  pathname[2] <- row.names(pathname)
  row.names(pathname) <- 1:nrow(pathname)
  pathname[3] <- apply_pb(pathname[2], 1, function(x) numbcomp(x))
  FileName <- "PATHbase"
  dateForFile <- Sys.time()
  dateForFile <- gsub(" ", "", dateForFile)
  dateForFile <- gsub(":", "", dateForFile)
  userName <- Sys.info()
  userName <- userName[[7]]
  names(pathname)[c(1, 2, 3)] <- c(FileName, dateForFile, userName)
  dir.create(paste(R.home("library/PAPi/databases"), fileName,
                   sep = FolderDivisor), recursive = TRUE)
  placeToSave <- paste(R.home("library/PAPi/databases"), fileName,
                       "/COMPbase.csv", sep = FolderDivisor)
  write.csv(listOfComps, file = placeToSave, row.names = FALSE)
  placeToSave <- paste(R.home("library/PAPi/databases"), fileName,
                       "/PATHbase.csv", sep = FolderDivisor)
  write.csv(pathname, file = placeToSave, row.names = FALSE)
  print("Your new KEGG database has been installed.")
  if (save == TRUE) {
    dir.create(paste(folder, fileName, sep = FolderDivisor))
    placeToSave <- paste(folder, fileName,
                         "COMPbase.csv", sep =FolderDivisor)
    write.csv(listOfComps, file = placeToSave, row.names = FALSE)
    placeToSave <- paste(folder, fileName, 
                         "PATHbase.csv", sep = FolderDivisor)
    write.csv(pathname, file = placeToSave, row.names = FALSE)
    print("Your new KEGG database has been saved.")
  }
}
buildDatabase<-function(){
  library(PAPi)
  buildDatabase()
}

getPathways <- function(x) {
  TotalReport <- KEGGREST::keggGet(x)
  pathways <- TotalReport[[1]]$PATHWAY
  pathways <- data.frame(pathways)
  pathways <- row.names(pathways)
  pathways <- paste(pathways, collapse = ";")
  return(pathways)
}


load_smpdb<-function(path="Z:\\George skyline results\\maldiimaging\\DB\\SMPdb.csv"){
  smpdb<-read.csv(path)
  
  return(smpdb)
  
  
  
}

convert_PAPi<-function(metabolomicsData,localDatabase = "default",matchingscore=0.6,top_N=10,derivatisation=c("TMS","derivative","methyl","ester")){
  library(stringr)
  library(PAPi)
  library(stringdist)
  library(dplyr)
  library(BiocParallel)
  
  library(tcltk2)
  if (missing(metabolomicsData)) metabolomicsData=read_table_generic(tk_choose.files(caption = "Select formatted data for KEGG code conversion"))
  

  keggLibrary=try(read.csv(paste0(R.home(),"/library/PAPi/databasesdefault/COMPbase.csv"),stringsAsFactors = F))
  if (grep("Error",keggLibrary[1],ignore.case = T)==T){
    keggLibrary=read_table_generic(tk_choose.files(caption = "Select COMPbase.csv for KEGG code conversion"))
  }
  
    metabolomicsData$trimedNames=metabolomicsData$Names
  
  for (derivate in derivatisation){
    metabolomicsData$trimedNames=gsub(derivate,"",metabolomicsData$trimedNames,ignore.case = T)
  }
  
  metabolomicsData$trimedNames=gsub("\\(split peak .\\)","",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub("\\(split peak ..\\)","",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub("\\(peak .\\)","",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub("\\(peak ..\\)","",metabolomicsData$trimedNames)
  #metabolomicsData$trimedNames=gsub("methyl ester"," ",metabolomicsData$trimedNames)
  #metabolomicsData$trimedNames=gsub(" TMS derivative"," ",metabolomicsData$trimedNames)
  #metabolomicsData$trimedNames=gsub("methyl"," ",metabolomicsData$trimedNames)
  #metabolomicsData$trimedNames=gsub("Methyl"," ",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub("-"," ",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub(","," ",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub(" {2,}"," ",metabolomicsData$trimedNames)
  metabolomicsData$trimedNames=gsub("[0-9]","",metabolomicsData$trimedNames)
  #metabolomicsData_Names=gsub("-"," ",metabolomicsData_Names)
  
  metabolomicsData_Names=metabolomicsData$trimedNames[2:length(metabolomicsData$trimedNames)]
  keggLibraryentry=as.vector(keggLibrary[,1])
  metabolomicsData_Names=unique(metabolomicsData_Names)
  #stringsim(metabolomicsData_Names,keggLibraryentry)
  #keggLibraryentry_split=stringr::str_split(keggLibraryentry,"; ")
  
    keggLibraryentry_split_index=lapply(1:length(keggLibraryentry),function(x,keggLibraryentry){
      keggLibraryentry_split_index_name=unlist(stringr::str_split(keggLibraryentry[x],"; ") )
      keggLibraryentry_split_index_index=length(keggLibraryentry_split_index_name)
      return(data.frame(Name=keggLibraryentry_split_index_name,index=base::rep(x,keggLibraryentry_split_index_index),stringsAsFactors = F))
    },keggLibraryentry) 
    keggLibraryentry_split_index=do.call(rbind,keggLibraryentry_split_index)
    #a=as.matrix(a,nrow=2)
    #a=as.data.frame(t(as.matrix(a,nrow=2)),stringsAsFactors = F)
    
    matchentry=(sapply(sapply(metabolomicsData_Names, tolower), stringsim,sapply(keggLibraryentry_split_index$Name, tolower)))
    
  stringsim_list<-function(x,list,t){
    max(stringsim(t,list[[x]]))
  }
  
  BPPARAM=bpparam()
  bpworkers(BPPARAM)=4
  bpprogressbar(BPPARAM)=T
  #matchentry=BiocParallel::bplapply(metabolomicsData_Names,function(x,keggLibraryentry_split){
  #  library(stringdist)
  #  stringsim_list<-function(x,list,t){
  #    max(stringsim(t,list[[x]]))
  #  }
    
   # unlist(lapply(1:length(keggLibraryentry_split),stringsim_list,list=keggLibraryentry_split,t=x))
    
   # },keggLibraryentry_split,BPPARAM = BPPARAM)
  
  #rownames(matchentry)=keggLibraryentry
  matchentry=as.data.frame(matchentry,stringsAsFactors = F)
  
  
  matchentry_topN<-lapply(metabolomicsData_Names, function(x,matchentry,keggLibraryentry,top_N){
    a=as.numeric(as.character(matchentry[,x]))
    as=tail(sort(a),top_N)
    ass=sapply(as, function(as,a){
     which(a==as)},a) 
    assunlist=unlist(ass)
    assunlist=tail(unique(assunlist),top_N)
    
    #ass=a %in% as
    #which(ass)
    asss=data.frame(original=x,retrievedID=keggLibraryentry[assunlist],SimScore=a[assunlist],stringsAsFactors = F)
 
  },matchentry,keggLibraryentry_split_index$Name,top_N)
  
  matchentry_topNbind=do.call(rbind,matchentry_topN)
  
  autoselect_top_1<-function(x,metabolomicsData_Names,matchentry_topN,matchingscore=0.6){
    
    name=metabolomicsData_Names[x]
    df=matchentry_topN[[x]]
    if (max(df$SimScore)>=matchingscore){
    return(data.frame(origin=name,covert=df[head(which(df$SimScore==max(df$SimScore)),1),2],stringsAsFactors = F))
    }else{
    return(data.frame(origin=name,covert="Not_found",stringsAsFactors = F))
    }
    
  }
  
  #metabolomicsData_Names_autoconvert=base::lapply(1:length(metabolomicsData_Names), autoselect_top_1,metabolomicsData_Names,matchentry_topN)
  metabolomicsData_Names_autoconvert=base::lapply(1:length(metabolomicsData_Names), autoselect_top_1,metabolomicsData_Names,matchentry_topN,matchingscore=matchingscore) %>% unlist
  metabolomicsData_Names_autoconvert=matrix(data=metabolomicsData_Names_autoconvert,nrow = 2)
  metabolomicsData_Names_autoconvert=as.data.frame(t(metabolomicsData_Names_autoconvert),stringsAsFactors = F)
  names(metabolomicsData_Names_autoconvert)=c("trimedNames","converted")
  metabolomicsData_Names_autoconvert=merge(metabolomicsData_Names_autoconvert,metabolomicsData,by="trimedNames",all.y = T)
  keggLibraryentry_split_index$converted=keggLibraryentry_split_index$Name
 #metabolomicsData_Names_autoconvert$origin=metabolomicsData$Names[metabolomicsData$trimedNames==metabolomicsData_Names_autoconvert$trimed]
  metabolomicsData_Names_autoconvert=merge(keggLibraryentry_split_index[,c("converted","index")],metabolomicsData_Names_autoconvert,by= "converted",all.y = T)
  metabolomicsData_Names_autoconvert$Keggcode<-keggLibrary[metabolomicsData_Names_autoconvert$index,2]
  papiData=metabolomicsData_Names_autoconvert[,c("Keggcode",colnames(metabolomicsData_Names_autoconvert)[colnames(metabolomicsData_Names_autoconvert) %in% colnames(metabolomicsData)])]
  papiData=rbind(papiData[grepl("Replicates",papiData$Names,ignore.case = T),],papiData[!grepl("Replicates",papiData$Names,ignore.case = T),])
  papiData$Keggcode[1]="Replicates"	
  colnames(papiData)[1]="Name"
  papiData$Name=gsub("^cpd:","",unlist(unname(papiData$Name)))
  #sapply(papiData$Name,gsub,)
  rownames(papiData)=1:nrow(papiData)
  papiData$trimedNames=NULL
  papiData$Names=NULL
  papiData$V1=NULL
  workdir=getwd()
  if (dir.exists(paste(workdir,"/Pathway",sep=""))==FALSE){dir.create(paste(workdir,"/Pathway",sep=""))}
  
  write.csv(matchentry_topNbind,"Pathway/match_entry_topN.csv",row.names = F)
  write.csv(papiData,"Pathway/papiData.csv",row.names = F)
  message(paste("Result data have been saved to",paste(workdir,"/Pathway",sep="")))
return(papiData)
  }

run_PAPi<-function(papiData,localDatabase = "default"){
library(PAPi)
  
  library(tcltk2)
  if (missing(papiData)) papiData=read_table_generic(tk_choose.files(caption = "Select papiData.csv for PAPi analysis"))
  papiData=papiData[!is.na(papiData$Name),]
  #data(papiData)
  papiData=as.data.frame(papiData)
  papiResults <- PAPi::papi(papiData, save = FALSE, offline = TRUE, localDatabase = localDatabase)
  
  #data=matchentry %>% top_n( 10,1)
standardcond=unique(as.character(papiData[1,]))

if (sum(grepl("qc",standardcond,ignore.case = T))>=1) {Ref.cond=standardcond[grepl("qc",standardcond,ignore.case = T)]
}else {Ref.cond=standardcond[23]}
  
  

workdir=getwd()
if (dir.exists(paste(workdir,"/Pathway",sep=""))==FALSE){dir.create(paste(workdir,"/Pathway",sep=""))}
plot_papiResults=cbind(papiResults[,1],papiResults[,(papiResults[1,] %in% c("C","IC","OC"))])
plot_papiResults=plot_papiResults[plot_papiResults$C1!=0,]
png(paste0("Pathway/PAPi_graph_",Ref.cond,".png"),width = 1600,height = 900)
papiLine(
  plot_papiResults,
   relative = TRUE,setRef.cond = T, Ref.cond=Ref.cond,
   save = FALSE,
  yscale=c(-5,3)
   )

dev.off()
write.csv(papiResults,"Pathway/papiResults.csv",row.names = F)
message(paste("Result data have been saved to",paste(workdir,"/Pathway",sep="")))

}

PAPi_formatting<-function(selectfile,
                          selectinfofile,
                          output=TRUE){
  library(tcltk)
  library(tcltk2)
  if (missing(selectfile)) selectfile=tk_choose.files(caption = "Select area.csv files to be converted for PAPi analysis")
  if (missing(selectinfofile)) selectinfofile=tk_choose.files(caption = "Select info file for PAPi analysis")
  
  rootdir=dirname(selectfile)
  fileName=gsub(paste0(rootdir,"/"),"",selectfile)
  file<-read_table_generic(selectfile)
  file<-as.data.frame(file)
  #file<-openxlsx::read.xlsx(selectfile)
  infofile<-read_table_generic(selectinfofile)
  file<-data_test_rename(c("CAS","Name"),file)
  infofile<-data_test_rename(c("Name","Condition"),infofile)
  condition<-colnames(file) %in% infofile$Name
   #file<-data.frame(file,stringsAsFactors = F)
  #infofile$Name[sort(order(infofile$Name)[which(condition)])]
   existfile<-infofile$Name %in% colnames(file)
    newfile<-file[,infofile$Name[existfile]]
   
   newfile<-cbind((file[,"Name"]),newfile)
   colnames(newfile)[1]="Names"
   conditionnames<-matrix(nrow=length(which(condition)),ncol = 2)
   i=1
   for (conditionname in which(condition)){
     conditionnames[i,1]=colnames(file)[conditionname]
     conditionnames[i,2]=infofile$Condition[infofile$Name==colnames(file)[conditionname]]
     i=i+1
   }
   newfile$Names=as.character(newfile$Names)
   #conditionnames<-infofile$Condition[infofile$Name==colnames(file)]
   newfile1<-rbind(as.character(c("Replicates",infofile$Condition[existfile])),newfile)
   if (output){
     write.csv(newfile1,paste0(rootdir,"\\KEGG_formated_",fileName),row.names = F)
   }
   message(paste("Result data have been saved to",paste(rootdir,sep="")))
   return(newfile1)
}


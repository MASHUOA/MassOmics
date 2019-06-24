buildLib <-
  function(AmdisLib, folder, save = TRUE, output = "ion_lib", verbose = TRUE, mz_L=0,mz_U=3000) {
    
    ######### Check if a .msl file ###########################
    isMSLdlg <- function(titleMSG, errorMSG) {
      t = 0
      while (t == 0) {
        checkIfCsv <- dlgOpen(title = titleMSG, multiple = FALSE)$res
        checkIfCsv2 <- basename(checkIfCsv)
        checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
        checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
        if (toupper(checkIfCsv4) == "MSL") {
          t = 1
          return(checkIfCsv)
        } else {
          dlgMessage(errorMSG)
        }
      }
    }
    
    ######### Check if a .msl file  NO DIALOG BOX ############
    isMSL <- function(pathfile, errorMSG) {
      t = 0
      checkIfCsv2 <- basename(pathfile)
      checkIfCsv3 <- unlist(strsplit(checkIfCsv2, "\\."))
      checkIfCsv4 <- checkIfCsv3[length(checkIfCsv3)]
      if (toupper(checkIfCsv4) == "MSL") {
        t = 1
        return(t)
      } else {
        t = 0
        return(t)
      }
    }
    
    ########### is.odd function #######
    is.odd <- function(x) x%%2 != 0
    ###################################
    
    ####### Change language so it accepts weird characters #######
    Sys.setlocale("LC_ALL", "C")
    ##############################################################
    
    ############# Load AMDIS' library #######################################################################################
    if (missing(AmdisLib)) {
      AmdisLib <- isMSLdlg("Select the .msl file of the AMDIS library in use.", 
                           "The selected file is not a msl file")
      AmdisLib <- data.frame(read.csv(AmdisLib, sep = "\t", header = FALSE))
    } else {
      if (is.data.frame(AmdisLib)) {
        AmdisLib <- AmdisLib
      } else {
        if (is.character(AmdisLib)) {
          checkIfCsv <- isMSL(AmdisLib)
          if (checkIfCsv == 1) {
            inputTest <- file.access(AmdisLib, 0)
            if (inputTest == 0) {
              AmdisLib = data.frame(read.csv(AmdisLib, sep = "\t", header = FALSE))
            } else {
              dlgMessage("The AmdisLib specified is not accessible. Please, choose a msl file of the AMDIS library in use.")
              AmdisLib <- isMSLdlg("Select the .msl file of the AMDIS library in use.", 
                                   "The selected file is not a msl file")
              AmdisLib = data.frame(read.csv(AmdisLib, sep = "\t", header = FALSE))
            }
          } else {
            dlgMessage("The AmdisLib specified is not a msl file. Please, choose a msl file of the AMDIS library in use.")
            AmdisLib <- isMSLdlg("Select the .msl file of the AMDIS library in use.", 
                                 "The selected file is not a msl file")
            AmdisLib = data.frame(read.csv(AmdisLib, sep = "\t", header = FALSE))
          }
        } else {
          dlgMessage("The AmdisLib must be specified as character. Please, choose a msl file of the AMDIS library in use.")
          AmdisLib <- isMSLdlg("Select the .msl file of the AMDIS library in use.", 
                               "The selected file is not a msl file")
          AmdisLib = data.frame(read.csv(AmdisLib, sep = "\t", header = FALSE))
        }
      }
    }
    
    ############################### Set folder if SAVE == TRUE ###################################
    if (save) {
      if (missing(folder)) {
        folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
      } else {
        if (is.character(folder)) {
          isFolder <- file.access(as.character(folder), 0)
          if (isFolder == 0) {
            isFolder <- file.info(folder)
            if (!isFolder$isdir) {
              message("The folder defined to save the results is not a valid path.")
              folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
            }
          } else {
            message("The folder defined to save the results is not a valid path.")
            folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
          }
        } else {
          message("The path to the folder where the results will be saved must be specified as character.")
          folder = dlgDir(title = "Select the folder where the output file will be saved.")$res
        }
      }
    }
    
    
    lib <- AmdisLib
    finalLib <- c()
    
    ########## Set progress bar #####################################
    if (verbose) {
      bar <- 0
      toBar <- nrow(lib)
      pb <- txtProgressBar(min = 0, max = toBar, style = 3, width = 50)
    }
    #################################################################
    
    ########## Start collecting info from the MSL file ###############
    while (nrow(lib) > 0) {
      i <- 1
      ## Get the row containing the name of the compound ###
      getname <- c()
      while (length(getname) == 0) {
        getname <- grep("NAME:", lib[i, ], ignore.case = FALSE)
        i <- i + 1
      }
      ## Keep just the name of the compound ###
      getname <- lib[i - 1, ]
      getname <- gsub("NAME:", "", getname, ignore.case = FALSE)
      nameLoc <- i - 1
      ## Get Retention time ###
      getRT <- c()
      i0 <- i
      while (length(getRT) == 0) {
        getRT <- grep("RT:", lib[i, ], ignore.case = FALSE)
        i <- i + 1
        
        
        if ((i - i0) == 20) {
          getRT<-"NA"
          i=i0
          #message("\nThere is no information about retention time (RT) in the selected MSL file or AMDIS library.")
        }
      }
      ## Keep just the name of the compound ###
      #getRT <- lib[i - 1, ]
      getRT <- gsub("RT:", "", getRT, ignore.case = FALSE)
      ## Get rows containing ions and intensities ###
      getion <- c()
      while (length(getion) == 0) {
        getion <- grep("NUM PEAKS:", lib[i, ], ignore.case = FALSE, )
        i <- i + 1
      }
      startpeak <- i
      getname2 <- c()
      checkNames <- grep("NAME:", lib[, 1], ignore.case = FALSE)
      if (length(checkNames) > 1) {
        while (length(getname2) == 0) {
          getname2 <- grep("NAME:", lib[i, ], ignore.case = FALSE)
          i <- i + 1
        }
        endpeak <- i - 2
      } else {
        endpeak <- nrow(lib)
      }
      peak_rows <- lib[startpeak:endpeak, ]
      peak_rows <- data.frame(as.character(peak_rows))
      for (k in 1:nrow(peak_rows)) {
        row1 <- data.frame(unlist(strsplit(as.character(peak_rows[k, ]), 
                                           ") \\(")))
        row1 <- data.frame(gsub(")", "", row1[, 1]))
        row1 <- data.frame(gsub("( ", "", row1[, 1], fixed = TRUE))
        row1 <- data.frame(gsub("(", "", row1[, 1], fixed = TRUE))
        row1 <- apply(row1[1], 1, function(x) data.frame(unlist(strsplit(x, 
                                                                         "[[:blank:]]"))))
        row1 <- lapply(row1, function(x) x[x != ""])
        row1 <- unlist(row1)
        frags <- row1[is.odd(1:length(row1))]
        if (length(frags) > 0) {
          int <- row1[!is.odd(1:length(row1))]
          oneGroup <- data.frame(cbind(frags, int))
          if (k == 1) {
            totalPeak <- oneGroup
          } else {
            totalPeak <- rbind(totalPeak, oneGroup)
          }
        }
      }
      totalPeak[1] <- as.numeric(as.character(totalPeak[, 1]))
      totalPeak[2] <- as.numeric(as.character(totalPeak[, 2]))
      totalPeak=totalPeak['&'(totalPeak$frags>=mz_L,totalPeak$frags<=mz_U),]
      totalPeak <- totalPeak[order(totalPeak[, 2], decreasing = T), ]
      if (nrow(totalPeak) >= 4) {
        ions <- totalPeak[1:4, ]
      } else {
        ions <- totalPeak
        while (nrow(ions) < 4) {
          ions <- rbind(ions, c(NA, NA), TRUE)
        }
      }
      
      ### Now, prepare the data frame with name of compounds, ref_ion1, ref_ion2, etc...
      if (length(finalLib) == 0) {
        finalLib <- data.frame(Name = getname, RT = getRT, ref_ion1 = ions[1, 
                                                                           1], ref_ion2 = ions[2, 1], ref_ion3 = ions[3, 1], ref_ion4 = ions[4, 
                                                                                                                                             1], ion2to1 = as.numeric(ions[2, 2])/ions[1, 2], ion3to1 = as.numeric(ions[3, 
                                                                                                                                                                                                                        2])/ions[1, 2], ion4to1 = as.numeric(ions[4, 2])/ions[1, 2])
        
      } else {
        finalLib2 <- data.frame(Name = getname, RT = getRT, ref_ion1 = ions[1, 
                                                                            1], ref_ion2 = ions[2, 1], ref_ion3 = ions[3, 1], ref_ion4 = ions[4, 
                                                                                                                                              1], ion2to1 = as.numeric(ions[2, 2])/ions[1, 2], ion3to1 = as.numeric(ions[3, 
                                                                                                                                                                                                                         2])/ions[1, 2], ion4to1 = as.numeric(ions[4, 2])/ions[1, 2])
        finalLib <- rbind(finalLib, finalLib2)
        
      }
      lib <- data.frame(lib[-c(nameLoc:endpeak), ])
      
      ############ Update progress bar #########
      if (verbose) {
        bar <- bar + length(c(nameLoc:endpeak))
        ## Sys.sleep(0.01)
        setTxtProgressBar(pb, bar)
      }
      ##########################################
    }
    if (verbose)
      close(pb)
    
    ##### Select compounds showing less than 1 min of difference in their RTs #####
    SameRTandIons <- data.frame()
    finalLibRec <- finalLib
    finalLibRec[2] <- as.numeric(as.character(finalLibRec[, 2]))
    finalLibRec <- finalLibRec[order(finalLibRec[3]), ]
    while (nrow(finalLibRec) > 1) {
      row1 <- finalLibRec[1, ]
      sameIon <- finalLibRec[finalLibRec[3] == row1[1, 3], ]
      finalLibRec <- finalLibRec[-c(which(finalLibRec[3] == row1[1, 3])), ]
      if (nrow(sameIon) > 1) {
        sameIon <- sameIon[order(sameIon[, 2]), ]
        origSameIon <- sameIon
        sameIon[2] <- sameIon[, 2] * 2
        sameIon[2] <- round(sameIon[, 2])
        sameIon[2] <- sameIon[, 2]/2
        sameIon[2] <- round(sameIon[, 2])
        sameIon <- origSameIon[duplicated(sameIon[, 2]) | duplicated(sameIon[nrow(sameIon):1, 
                                                                             2])[nrow(sameIon):1], ] ## Select compounds showing less than 1 minute difference between compounds.
        if (nrow(sameIon) > 1) {
          SameRTandIons <- rbind(SameRTandIons, sameIon)
        }
      }
    }
    
    ######## Present compounds to the user ######################################
    if (verbose){
    if (nrow(SameRTandIons) > 0) {
      message("The following compounds have similar retention times (less than one minute difference) and they also share the same ion mass fragment as reference. For these compounds, we suggest the use of different ion mass fragments in order to avoid a strong coelution effect. This is not an error, just a warning.")
      pandoc.table(data.frame(Names = SameRTandIons[, 1], RT = SameRTandIons[, 
                                                                             2], Ion = SameRTandIons[, 3]), justify = "left", emphasize.cols = 3)
    }
    }
    ########## Reorder results accroding to RT ###################################
    finalLib <- finalLib[order(as.numeric(as.character(finalLib[, 2]))), ]
    
    ########## Save results if save == TRUE ######################################    
    if (save) {
      sheet <- output
      store <- file.path(folder, paste(sheet, ".csv", sep = ""))
      inputTest <- file.access(store, 0)
      if (inputTest == 0) {
        addFile <- 1
        while (inputTest == 0) {
          store <- file.path(folder, paste(sheet, addFile, ".csv", sep = ""))
          inputTest <- file.access(store, 0)
          addFile <- addFile + 1
        }
      }
      write.csv(finalLib, file = store, row.names = FALSE)
      message("File saved: ", store, "\n")
    }
    ######## Set SysLocale to default ######
    Sys.setlocale("LC_ALL", "")
    ########################################
    return(finalLib)
  }

parse_msp<- function(peak_rows){
  

peak_rows <- data.frame(as.character(peak_rows))
for (k in 1:nrow(peak_rows)) {
  row1 <- data.frame(unlist(strsplit(as.character(peak_rows[k, ]), 
                                     ") \\(")))
  row1 <- data.frame(gsub(")", "", row1[, 1]))
  row1 <- data.frame(gsub("( ", "", row1[, 1], fixed = TRUE))
  row1 <- data.frame(gsub("(", "", row1[, 1], fixed = TRUE))
  row1 <- apply(row1[1], 1, function(x) data.frame(unlist(strsplit(x, 
                                                                   "[[:blank:]]"))))
  row1 <- lapply(row1, function(x) x[x != ""])
  row1 <- unlist(row1)
  frags <- row1[is.odd(1:length(row1))]
  if (length(frags) > 0) {
    int <- row1[!is.odd(1:length(row1))]
    oneGroup <- data.frame(cbind(frags, int))
    if (k == 1) {
      totalPeak <- oneGroup
    } else {
      totalPeak <- rbind(totalPeak, oneGroup)
    }
  }
}
totalPeak$frags=as.numeric(as.character(totalPeak$frags))
totalPeak$int=as.numeric(as.character(totalPeak$int))
return(totalPeak)
}


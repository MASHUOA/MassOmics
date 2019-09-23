read_table_generic<-function(path,header=T){
  
  if(grepl(".xlsx", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1,colNames = header)
  }else if(grepl(".xls", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1,colNames = header)
  }else if(grepl(".xlsm", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1,colNames = header)
  }else if(grepl(".csv", path)){
    d <- data.table::fread(path,header = header)
    if(nrow(d)<2){d <- read.csv(path,header = header)}
  }
  return(d)
}



windows_filename<- function(stringX){
  stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\|]")
  stringX<-gsub("\"", "", stringX)
  if (nchar(stringX)>=100){stringX=strtrim(stringX,100)}
  return(stringX)
  
}
windows_filename2<- function(stringX){
  library(stringr)
  stringX<-str_replace_all(stringX, "[^[:alnum:]]", " ")
  stringX<-str_replace_all(stringX, ",", " ")
  return(stringX)
  
}

Advanced_building<-function(){
  install.packages("ggplot2")
  install.packages("pryr")
  install.packages("devtools")
  devtools::install_github("hadley/lineprof")
  library(pryr)
  object_size(1:10)
  object_size(mean)
  sizes <- sapply(0:50, function(n) object_size(seq_len(n)))
  plot(0:50, sizes, xlab = "Length", ylab = "Size (bytes)", 
       type = "s")
  object_size(numeric())
  plot(0:50, sizes - 40, xlab = "Length", 
       ylab = "Bytes excluding overhead", type = "n")
  abline(h = 0, col = "grey80")
  abline(h = c(8, 16, 32, 48, 64, 128), col = "grey80")
  abline(a = 0, b = 4, col = "grey90", lwd = 4)
  lines(sizes - 40, type = "s")
  mem_used()
  
  read_delim <- function(file, header = TRUE, sep = ",") {
    # Determine number of fields by reading first line
    first <- scan(file, what = character(1), nlines = 1,
                  sep = sep, quiet = TRUE)
    p <- length(first)
    
    # Load all fields as character vectors
    all <- scan(file, what = as.list(rep("character", p)),
                sep = sep, skip = if (header) 1 else 0, quiet = TRUE)
    
    # Convert from strings to appropriate types (never to factors)
    all[] <- lapply(all, type.convert, as.is = TRUE)
    
    # Set column names
    if (header) {
      names(all) <- first
    } else {
      names(all) <- paste0("V", seq_along(all))
    }
    
    # Convert list into data frame
    as.data.frame(all)
  }
  
  library(ggplot2)
  write.csv(diamonds, "diamonds.csv", row.names = FALSE)
  
  library(lineprof)
  
  source("Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy/R/read-delim.R")
  prof <- lineprof(read_delim("diamonds.csv"))
  shine(prof)
}
string_slice <- function(string, size) {
  pat <- paste0('(?<=.{',size,'})')
  strsplit(string, pat, perl=TRUE)
}



getMonomass <- function(formular){
  Rdisop::getMolecule(formular)[["exactmass"]][1]}
getMonomass_para <- function(i,formularlist){
  return(Rdisop::getMolecule(formularlist[i])$exactmass)}
setworkdir<-function(workdir){
  if (dir.exists(workdir)==FALSE){dir.create(workdir)}
  setwd(workdir)
}




parse_msl_par<-function(i,lib.txt,starts,stops,mz_L,mz_U){
  is.odd <- function(x) x%%2 != 0
  tmp<-lib.txt[starts[i]:stops[i]]
  tmp<-tmp[grep("^\\(",tmp)]
  tmp<-paste(tmp,collapse = "",sep="")
  
  #tmp<-parse_msp(tmp)
  tmp <- data.frame((tmp),stringsAsFactors = F)
  
  row1 <- data.frame(unlist(strsplit(as.character(tmp[,1]),") \\(")))
  row1 <- data.frame(gsub(")", "", row1[, 1]))
  row1 <- data.frame(gsub("( ", "", row1[, 1], fixed = TRUE))
  row1 <- data.frame(gsub("(", "", row1[, 1], fixed = TRUE))
  row1 <- apply(row1[1], 1, function(x) data.frame(unlist(strsplit(x, 
                                                                   "[[:blank:]]"))))
  row1 <- unlist(lapply(row1, function(x) x[x != ""]))
  
  frags <- row1[is.odd(1:length(row1))]
  int <- row1[!is.odd(1:length(row1))]
  totalPeak <- data.frame(cbind(frags, int))
  #totalPeak <- oneGroup  
  
  totalPeak$frags=as.numeric(as.character(totalPeak$frags))
  totalPeak$int=as.numeric(as.character(totalPeak$int))
  totalPeak<-totalPeak['&'(totalPeak$frags>=mz_L,totalPeak$frags<=mz_U),]
  
  max(totalPeak[which(totalPeak$int==max(totalPeak$int)),"frags"])
  #totalPeak[which(totalPeak$int==max(totalPeak$int)),"frags"]
}


dev.new.OS<-function(){
switch(Sys.info()[['sysname']],
       Windows= {windows()},
       Linux  = {X11()},
       Darwin = {quartz()})  
    
  
}

std <- function(x) sd(x)/sqrt(length(x))

std.outliner.rm <- function(x) {
  x<-x[!('|'(x>(3*sd(x)+median(x)),x<(-3*sd(x)+median(x))))]
  sd(x)/sqrt(length(x))
}

Combine_result<-function(e,f,p){
  
  finaltable<-rbind(t(p),e)
  
  return(finaltable)
  
  
}
Combine_result_file<-function(path){
  
  e=read.csv(paste0(path,"/","e.csv"),stringsAsFactors = F)
  
  enorm=read.csv(paste0(path,"/","SERRF_normalized.csv"),stringsAsFactors = F)
  
  f=read.csv(paste0(path,"/","f.csv"),stringsAsFactors = F)
  
  p=read.csv(paste0(path,"/","p.csv"),header = F,stringsAsFactors = F)
  
  p=p[,2:5]
  
  
  
  rownames(p)=as.character(p[,1])
  rownames(p)[1]="label"
  tp=as.data.frame(t(p),stringsAsFactors=F)
  colnames(e)=as.character(p[,1])
  finaltable<-rbind(tp,e)
  colnames(enorm)=as.character(p[,1])
  finaltablenorm<-rbind(tp,enorm)
  
  return(list(raw=finaltable,norm=finaltablenorm))
  
  
}

data_test_rename<-function(required_col,df){
  
  testcolumeresult = testcolume(df,required_col)
  
  if (length(testcolumeresult$failcol)>0){
    
    stop(paste(testcolumeresult$failcol,"column is missing in the input file, please check your datafile"))  
    
  }
  
  if (length(testcolumeresult$renamecol)>0){
    
    lapply(testcolumeresult$renamecol,gsub,testcolumeresult$renamecol,ignore.case = T,x=names(df))
    
    for (i in length(testcolumeresult$renamecol)){
      
      names(df)=gsub(testcolumeresult$renamecol,testcolumeresult$renamecol,ignore.case = T,x=names(df))
      
    }
    
  }
  
  if (length(testcolumeresult$duplicatecol)>0){
    
    message(paste("Found ambigous or duplicate column: ",testcolumeresult$duplicatecol))
    
  }
  
  return(df)
  
}

testcolume<-function(df,testcolnames,match_exact=T){
  
  library(stringr)
  
  if (match_exact){testcolnames=paste0("^",testcolnames,"$")}
  
  testcolnames=unique(testcolnames)
  
  test=sapply(testcolnames,FUN = grepl,names(df))
  
  testresult=sapply(colnames(test),FUN = function(x,df){sum(df[,x])},test) == 0
  
  case_test=sapply(testcolnames,FUN = grepl,names(df),ignore.case = T)
  
  case_testresult=sapply(testcolnames,FUN = function(x,df){sum(df[,x])},case_test) > 1
  
  duplicatecol=names(case_testresult)[case_testresult==T]
  
  if (length(duplicatecol)==0){duplicatecol=NULL}
  
  testfail=names(testresult)[testresult==T]
  
  passcol=names(testresult)[testresult==F]
  
  renamecol=NULL
  
  if (length(passcol)<length(testcolnames)){
    
    case_test=sapply(testfail,FUN = grepl,names(df),ignore.case = T)
    
    case_testresult=sapply(colnames(case_test),FUN = function(x,df){sum(df[,x])},case_test) == 1
    
    renamecol=names(case_testresult)[case_testresult==T]
    
    if (length(renamecol)>0){
      
      failcol=testfail[!grep(renamecol,testfail)]  
      
    }else{
      
      failcol=NULL
      
    }
    
    
    
  } else{
    
    failcol=NULL
    
  }
  
  
  
  return(list(passcol=passcol,renamecol=renamecol,failcol=failcol,duplicatecol=duplicatecol))
}


Creat_short_cut<-function(){
  
  int_script=paste(file.path(path.package(package="MassOmics")),"/R/runMassOmics.R", sep="")
  int_script=shortPathName(int_script)
  r_sript_path=paste0(R.home("bin"),"/Rscript.exe")
  
  switch(Sys.info()[['sysname']],
         Windows= {
           batpath=paste(file.path(path.package(package="MassOmics")),"/R/runMassOmics.bat", sep="")
           batpath=shortPathName(batpath)
           short_cut_icon=shortPathName(paste(file.path(path.package(package="MassOmics")),"/R/ico_eyw_icon.ico", sep=""))
           short_cut_path=shortPathName(paste(file.path(path.package(package="MassOmics")),"/R/shortcut.bat", sep=""))
           #short_cut_path_bk=paste(file.path(path.package(package="MassOmics")),"/R/shortcut_bk.bat", sep="")
           #short_cut_path_line=readLines(con=short_cut_path)
           #short_cut_path_line_bk=readLines(con=short_cut_path_bk)
           #identical(short_cut_path_line,short_cut_path_line_bk)
           writeLines(text=c(paste(r_sript_path,"--ess",int_script),"pause"),con=batpath)
           writeLines(text=c(
                             "@echo off",
                             quote('echo Set oWS = WScript.CreateObject("WScript.Shell") > CreateShortcut.vbs'),
                             'echo sLinkFile = "%HOMEDRIVE%%HOMEPATH%\\Desktop\\MassOmics.lnk" >> CreateShortcut.vbs',
                             'echo Set oLink = oWS.CreateShortcut(sLinkFile) >> CreateShortcut.vbs',
                             paste0('echo oLink.TargetPath = "',batpath,'" >> CreateShortcut.vbs'),
                             paste0('echo oLink.IconLocation = "',short_cut_icon,'" >> CreateShortcut.vbs'),
                             'echo oLink.Save >> CreateShortcut.vbs',
                             'cscript CreateShortcut.vbs',
                             'del CreateShortcut.vbs'
                             ),con=short_cut_path)
           shell(short_cut_path)
         },
         Linux  = {message("Sorry the Massomics shortcut only available for windows now.")},
         Darwin = {message("Sorry the Massomics shortcut only available for windows now.")}) 
  
}
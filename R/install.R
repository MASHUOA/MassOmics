######################################### install #########################################
install.MassOmics_dependency <- function(){
  if (length(.libPaths())!=1){librarypath=.libPaths()[2]}else{librarypath=.libPaths()[1]}
  
  message(paste0("Installing packages to ",librarypath," that may take a while!"))
  install.packages(c("pacman","devtools","backports"))
  library("pacman")
  library("devtools")
  library("backports")
  #devtools::install_github("tidyverse")
  BiocManager::install(c("KEGGREST", "mixOmics", "CAMERA", "qvalue"))
  p_install(c("openxlsx", "randomForest", "rfviz", "extrafont", "svDialogs", "svGUI", "KEGGREST", "rgl", "mixOmics", "tkrplot", "CAMERA", "qvalue"),try.bioconductor = T)
  
  p_install(c("tcltk","tcltk2","stringdist","pbapply","PAPi","Rcpp","openxlsx","tidyr","faahKO","randomForest","rfviz",
              "ggpubr","extrafont","mzR","Rcpp","xcms","svDialogs","svGUI","KEGGREST",
              "RColorBrewer","rgl","mixOmics","plyr","flux","tkrplot","multtest","XML",
              "CAMERA","qvalue", "doParallel","ggplot2","ggplot","MALDIquant","CAMERA"),try.bioconductor = TRUE,lib=librarypath,force = T,character.only = F)
  #font_import(prompt=F)
  #loadfonts(device = "win")
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
  BiocManager::install(c("mzR","Rcpp","xcms","svDialogs","svGUI","KEGGREST","RColorBrewer","rgl","mixOmics","plyr","flux","tkrplot","multtest","CAMERA","qvalue", "doParallel"))
  message("Installation Complete!")
}
######################################### install.mzmatch #########################################
install.mzmatch <- function(){
  require (tcltk)
  library("pacman")
  if (length(.libPaths())!=1){librarypath=.libPaths()[2]}else{librarypath=.libPaths()[1]}
  tkmessageBox(title = "Installation of MzMatch(R)", message = paste("Both the 32-bit and the 64 bit Java must be installed"))
  source("http://bioconductor.org/biocLite.R")
  if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
  BiocManager::install(pkgs=c("Rcpp","multtest", "rJava", "XML", "snow", "caTools","bitops", "ptw", "gplots", "tcltk2", "R.utils"),
                       ask = F,lib=librarypath)
  library(XML)
  source ("http://puma.ibls.gla.ac.uk/mzmatch.R/install_mzmatch.R")
  require (mzmatch.R)
  PeakML.Viewer()
  mzmatch.init(version.1=FALSE)
  print("Installation is completed!")
}
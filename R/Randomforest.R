
SERRF_norm<-function (e, f, p, batch = define_batch(e, f, p), QC.index, time = "time", 
                      cl = automakeCluster()) 
{
  library("randomForest")
  library("pbapply")
  pboptions(type="timer")
  #library("caret")
  #library("e1071")
  qc = rep(F, nrow(p))
  qc[QC.index] = T
  e. = e
  if (sum(qc[QC.index] == T) < length(qc)) {
    for (i in 1:nrow(e.)) {
      e.[i, qc] = unlist(by(data.frame(e.[i, ], qc), factor(batch[1, ]), function(x) {
        diff = (median(x[x[, 2], 1]) - median(x[!x[, 2], 1]))
        if (is.na(diff)) diff=0
      x[x[, 2], 1] - diff}))
      }
    #message("QC level")
  }
  message(c("predict start"))
  #cl=autoStopCluster(parallel::makeCluster(try(detectCores()/2)))
  #doParallel::registerDoParallel(cl)
  
  message(paste(nrow(f),"molecules to be classified.", ncol(e),"observations to be considered"))
  #message(f[1,])
  subsetqc <- function(X,limit=500) 
  { X[ ifelse(X<=limit, TRUE,FALSE)]
  }
  randomForest_par<-function(j, e., batch, QC.index, time) {
    set.seed(1)
    #message(c("randomForest for ", rownames(f)[j]))
    data = data.frame(y = e.[j, ], t(e.[-j, ]), time = as.numeric(time), batch = as.numeric(batch[1, ]))
    colnames(data) = c("y", paste0("X", 1:nrow(e.))[-j], "time","batch")
    #data=data[QC.index,]
    #model = randomForest.default(y=data$y, x = data, importance = F)
    #model = randomForest.default(y ~ ., data = data, subset = QC.index, importance = F)
    #data=data[1:500,]
    model = randomForest::randomForest(y ~ ., data = data, subset = QC.index, importance = T,keep.inbag=T)
    newdata = data.frame(t(e.[-j, ]), time = as.numeric(time), batch = as.numeric(batch[1, ]))
    colnames(newdata) = c(paste0("X", 1:nrow(e.))[-j], "time","batch")
    #newdata=newdata[1:500,]
    new = (e.[j, ]/predict(model,newdata=newdata)) * median(e.[j, ])
    return(new)
    }
  #res1pbcl <- pblapply(1:B, function(i) fun(bid[,i]), cl = cl)
  pred = pblapply(  1:nrow(f), FUN=randomForest_par,cl = cl, e.=e., batch=batch, QC.index=QC.index, time=p[["time"]])
  
  #pred = parLapply(cl = cl,  1:nrow(f), fun=randomForest_par, e.=e., batch=batch, QC.index=QC.index, time=p[[time]])
  pred= matrix(unlist(t(pred)),ncol = length(pred))
  message("predict finished")
  e_SERRF_pred = t(pred)
  for (i in 1:nrow(e_SERRF_pred)) {
    e_SERRF_pred[i, p$sampleType == "QC"] = e_SERRF_pred[i,p$sampleType == "QC"] + (median(e[i, p$sampleType == "QC"], na.rm = T) - median(e[i, p$sampleType != "QC"], 
                                                                                                                        na.rm = T))
    e_SERRF_pred[i, e_SERRF_pred[i, ] < 0] = 0.5 * min(e_SERRF_pred[i,e_SERRF_pred[i, ] > 0])
  }
  return(list(e = e_SERRF_pred, p = p, f = f))
}

readData = function(path =  "G:\\data\\D\\data project D.xlsx",data=NULL){
  #check if it is csv of xlsx
  message("Data loading")
  if (is.null(data)){
  if(grepl(".xlsx", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1,colNames = FALSE)
  }else if(grepl(".csv", path)){
    # file = "C:\\Users\\Sili Fan\\Downloads\\val (18).csv"
    d <- data.table::fread(path)
  } 
  }else{d=data}

  
  # make "" as NA
  d[d==""] <- NA
  d[d=="Inf"] <- NA
  d[d=="NaN"] <- NA
  #### fData
  d.df<-as.data.frame(d)
  fData <- d.df[!is.na(d.df[,1]),c(which(is.na(d.df[1,])),sum(is.na(d.df[1,]))+1)] # The first row and column is critical of formating the data.
  colnames(fData) = as.character(fData[1,]); fData = data.frame(fData[-1,],stringsAsFactors = F,check.names = FALSE);rownames(fData) = 1:nrow(fData);
  # following steps keeps the column type.
  fData.=lapply(fData,function(x){
    if(sum(!is.na(suppressWarnings(as.numeric(x)))) == length(x)){
      as.numeric(x)
    }else{
      x
    }
  })
  fData. = do.call(cbind, lapply(fData., data.frame, stringsAsFactors=FALSE))
  colnames(fData.) = colnames(fData)
  fData = fData.
  
  fData = fData[,c(ncol(fData),2:ncol(fData)-1)]
  fData[[1]] = make.unique(fData[[1]], sep = '_')
  
  #### pData
  pData <- d.df[c(which(is.na(d[,1])),max(which(is.na(d[,1])))+1) ,!is.na(d[1,])]
  pData <- t(pData); colnames(pData) = pData[1,]; pData = data.frame(pData[-1,],stringsAsFactors = F,check.names = FALSE)
  # following steps keeps the column type.
  pData.=lapply(pData,function(x){
    if(sum(!is.na(suppressWarnings(as.numeric(x)))) == length(x)){
      as.numeric(x)
    }else{
      x
    }
  })
  pData. = do.call(cbind, lapply(pData., data.frame, stringsAsFactors=FALSE))
  colnames(pData.) = colnames(pData)
  pData = pData.
  
  pData = pData[,c(ncol(pData),2:ncol(pData)-1)]
  pData[[1]] = make.unique(make.names(pData[[1]]), sep = '_')
  
  #### eData
  eData <- d.df[!is.na(d[,1]),!is.na(d[1,])][-1,-1]
  eData <- sapply(eData, as.numeric)
  eData <- data.frame(eData,stringsAsFactors = F)
  colnames(eData) = pData[[1]]; rownames(eData) = fData[[1]]
  
  # # remove any unwanted character in columns of eData, fData and pData to _.
  # colnames(eData) = gsub("([_])|[[:punct:]]", "_", colnames(eData))
  # colnames(fData) = gsub("([_])|[[:punct:]]", "_", colnames(fData))
  # colnames(pData) = gsub("([_])|[[:punct:]]", "_", colnames(pData))
  
  # remove all the NA. And replace NA with "NA" Otherwise DataTables will give error.datatables warning requested unknown parameter
  # eData[is.na(eData)]="NA"
  # fData[is.na(fData)]="NA"
  # pData[is.na(pData)]="NA"
  
  # remove unwanted character in p.
  # for(i in 1:nrow(pData)){
  #   for(j in 1:ncol(pData)){
  #     pData[i,j] = gsub("\\+|~|-", " ", pData[i,j])
  #   }
  # }
  
  
  if(nrow(eData)>0){message("Data loaded")}
  return(list(e = eData, f = fData, p = pData))
}
#' @export
readData_serrf_native = function(path =  "G:\\data\\D\\data project D.xlsx",infopath="G:\\data\\D\\data project D.xlsx"){
  #check if it is csv of xlsx
  message("MassOmics native Data loading for SERRF")
  if(grepl(".xlsx", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1)
  }else if(grepl(".xls", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1)
  }else if(grepl(".xlsm", path)){
    d <- openxlsx::read.xlsx(path, sheet = 1)
  }else if(grepl(".csv", path)){
    # file = "C:\\Users\\Sili Fan\\Downloads\\val (18).csv"
    d <- data.table::fread(path)
    if(nrow(d)<2){d <- read.csv(path)}
  }
  
  # make "" as NA
  d[d==""] <- NA
  
  #### fData
  d.df<-as.data.frame(d)
  fData <- as.data.frame(d.df[,"Name"],stringasfactor=F)
  fData$No=rownames(d.df)
  colnames(fData) = c("label","No")
  
  #### pData
  if(grepl(".xls.?$", infopath)){
    pData <- openxlsx::read.xlsx(infopath, sheet = 1,colNames = T)
  }else if(grepl(".csv", infopath)){
    # file = "C:\\Users\\Sili Fan\\Downloads\\val (18).csv"
    if(nrow(d)<2){pData <- read.csv(infopath)}
    pData <- data.table::fread(infopath)
  }
  
  if (is.null(pData$time)){pData$time=rownames(pData)}
  pData<-pData[,c("Name","Batch", "Type","time")]
  colnames(pData)=c("label","batch","sampleType","time")
  sampleID=as.data.frame(colnames(d.df))
  #sampleID
  pData=pData[pData$label %in% colnames(d.df),]
  #### eData
  eData <- d.df
  rownames(eData)=d.df$Name
  
  eData <- eData[as.character(fData$label),as.character(pData$label)]
 #eData <- data.frame(eData,stringsAsFactors = F)
  
  
  # # remove any unwanted character in columns of eData, fData and pData to _.
  # colnames(eData) = gsub("([_])|[[:punct:]]", "_", colnames(eData))
  # colnames(fData) = gsub("([_])|[[:punct:]]", "_", colnames(fData))
  # colnames(pData) = gsub("([_])|[[:punct:]]", "_", colnames(pData))
  
  # remove all the NA. And replace NA with "NA" Otherwise DataTables will give error.datatables warning requested unknown parameter
  # eData[is.na(eData)]="NA"
  # fData[is.na(fData)]="NA"
  # pData[is.na(pData)]="NA"
  
  # remove unwanted character in p.
  # for(i in 1:nrow(pData)){
  #   for(j in 1:ncol(pData)){
  #     pData[i,j] = gsub("\\+|~|-", " ", pData[i,j])
  #   }
  # }
  
  
  if(nrow(eData)>0){message("Data loaded")}
  return(list(e = eData, f = fData, p = pData))
}
#' @export
Data_prep_area_SERRF<- function(selectfile=tk_choose.files(caption = "Select area.csv files to normalized by SERRF"),
                                selectinfofile=tk_choose.files(caption = "Select info file to normalized by SERRF"),
                                ourput=TRUE){
  file<-read.table(selectfile, header = F, fill = TRUE,sep = "\t",as.is=T)
  infofile<-openxlsx::read.xlsx(selectinfofile, sheet = 1,colNames = T, fill = TRUE)
  #file[file[,1]=="Alignment ID",file[file[,1]=="Alignment ID",]=="Metabolite name"]
  newfile<-NULL
  newfile$No<-as.data.frame(file[,1])
  newfile$"Metabolite name"<-as.data.frame(file[,file[file[,1]=="Alignment ID",]=="Metabolite name"])
  newfile<-as.data.frame(newfile)
  newfile[[1]]<-as.character(newfile[[1]])
  newfile[[2]]<-as.character(newfile[[2]])
  newfile[,3:(length(file)-20+3)]<-1
  newfile[,3:(length(file)-20+3)]<-file[,20:length(file)]
  #newfile[,3:(length(file)-20+3)]<-as.character(as.data.frame(newfile[,3:(length(file)-20+3)]))
  newfile[4,1]<-"No"
  
  #newfile<-rbind((lapply(newfile[4,],function(x){infofile[infofile[,"ReName"]==x,"Batch"]})),(newfile))
  newfile[1,3:ncol(newfile)]=infofile[infofile[,"ReName"] %in% newfile[4,3:ncol(newfile)],"Batch"]
  newfile[1,1]<-""
  newfile[1,2]<-"batch"
  newfile[2,2]<-"sampleType"
  newfile[3,2]="time"
  newfile[4,2]="label"
  newfile <- data.frame(lapply(newfile, as.character), stringsAsFactors=FALSE)
  a<-lapply((newfile[2,]),as.character)
  newfile[2,]<-a
  newfile<-newfile[,newfile[2,]!="Blank"]
  if (ourput==TRUE){write.table(newfile, file = paste(selectfile,"_output.csv",sep=""),row.names=FALSE, na="",col.names=FALSE, sep=",")
  }
  message(paste(selectfile,"_output.csv",sep=""))
  return(paste(selectfile,"_output.csv",sep=""))
}





#' Systematic error removal using random forest (SERRF) normalization
#'
#' This systematic error removal method using random forest algorithm to detect molecule cluster based on sequential intensities profile and use median intensity of this cluster to normalize the data and eliminat the unwanted systematic variations in large sample sets. 
#' 
#' @param input locate the molecule intensities data
#' @param Predict_level specify the sample class that will used for cluster prediction
#' @return None
#'
#' @examples
#' SERRF(input = "Area.csv",Predict_level=c("QC","Sample"))
#'
#' @export
#' 
#' 
#' 
#' 
#' 
SERRF <- function(input = "Area.csv",Predict_level="QC",data=NULL,datatype=c("MSDIAL","MASSOMICS"),infopath=NULL,Log_trans=F,zero_imputaion=T){
  library(tcltk)
  if ('&'(!file.exists(input) , is.null(data))){input =tk_choose.files(caption = "Select area.csv files to normalized by SERRF")}
  if ('&'(datatype=="MASSOMICS" , is.null(data))){infopath =tk_choose.files(caption = "Select caseinfo files")}
  
  SampleType=Predict_level
  setwd(base::dirname(input))
  # library(rgeolocate)
  #message(input)
  pacman::p_load(ggplot2, randomForest, parallel, officer, dplyr, rvg,data.table,ggpubr,extrafont)
  
  cl <- autoStopCluster(makeCluster(try(detectCores())))
  
  start = Sys.time()
  message(start)
  message(paste("Log transformation beforehand:",Log_trans))
  message(paste("Prediction level:",paste(SampleType,sep = "+",collapse = "+")))
 
  
  # info = read.csv(paste0("http://localhost:5984/serrf/info/info.csv"), stringsAsFactors = F,na.strings = "")
  #
  #
  # file <- system.file("extdata","ip2_sample.bin", package = "rgeolocate")
  #
  # code = ip2location(ip, file, c("country_code"))[[1]]
  # info$num[info$code == ip2location(ip, file, c("country_code"))[[1]]] =   info$num[info$code == code]+1
  #
  # put_att_csv = function(projectID = 'tryThu.Aug.17.14.53.35.2017', attname = 'test.csv', att = data.table::fread("G:\\initialize MetDA\\user_active.csv")){
  #   projectUrl <- paste0("http://localhost:5984/serrf/",projectID)
  #   projectList <- jsonlite::fromJSON(projectUrl)
  #
  #   new_att = projectList[["_attachments"]]
  #   new_att = new_att[!names(new_att)%in%attname]
  #   new_att[[attname]] = list(content_type="text/csv", data = RCurl::base64(
  #     paste0(R.utils::captureOutput(write.csv(att,stdout(), row.names=F)),collapse = "\n")
  #   ))
  #   projectList[["_attachments"]] = new_att
  #   result = RCurl::getURL(paste0("http://localhost:5984/serrf/",projectID),customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(projectList,auto_unbox = T, force = T))
  #   while(grepl("error",result)) {
  #     projectList <- jsonlite::fromJSON(projectUrl)
  #     new_att = projectList[["_attachments"]]
  #     new_att = new_att[!names(new_att)%in%attname]
  #     new_att[[attname]] = list(content_type="text/csv", data = RCurl::base64(
  #       paste0(R.utils::captureOutput(write.csv(att,stdout(), row.names=F)),collapse = "\n")
  #     ))
  #     projectList[["_attachments"]] = new_att
  #     result = RCurl::getURL(paste0("http://localhost:5984/serrf/",projectID),customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(projectList,auto_unbox = T, force = T))
  #     if(grepl("ok",result)){
  #       break;
  #     }
  #   }
  #   return(result)
  # }
  #
  # put_att_csv("info", "info.csv", info)
  
  
  if (and(file.exists(input) ,datatype=="MSDIAL")) {
    data = readData(path=input)
  }else  if (and(file.exists(input) ,datatype=="MASSOMICS")) {
    data = readData_serrf_native(path=input,infopath = infopath)
  }else if(!is.null(data) ){
    data = readData(data=data)
  
  }else{stop("Critical data input is missing!")}
  
  e = data$e
  f = data$f
  p = data$p
  
  p$sampleType=gsub("qc","QC",p$sampleType)
  
  p$sampleType=gsub("sample","Sample",p$sampleType)
  
  for (col in colnames(e)){
  e[,col]=as.numeric(as.character(e[,col]))  
  }
  
  # impute missing value.
    e[e==0]<-NA
    e[e=="Inf"]<-NA
    e[e=="-Inf"]<-NA
    e[e==""]<-NA
    e[e=="NaN"]<-NA
    
    #e[e<0]<-NA
  if (Log_trans) {
      e<-log(e)
      e[(e==-Inf)] <- NA
    }
  
    
    
  if(sum(is.na(e))>0){
    missing_compounds = unique(which(is.na(e), arr.ind = T)[,1])
    
    for(i in missing_compounds){
      
      if (zero_imputaion){
      
         e[i, is.na(e[i,])] = 1/2 * min(e[i,!is.na(e[i,])]) 
      } else  {
        e[i, is.na(e[i,])] = 0 
      }
      
    }
    
  }

  e = data.matrix(e)
  
  
  theme.scatter = theme(
    plot.title = element_text(size = rel(2), hjust = 0.5,face = 'bold',family = "Arial"),#title size.
    # axis.title = element_text(size = rel(2)),
    axis.text	 = element_text(colour = 'black'),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.key = element_rect(fill = "white",colour = "white"),
    legend.title = element_text(face = 'bold'),
    text=element_text(family="Arial")
  )
  
  # define batch.
  batch = matrix(rep(p$batch, nrow(f)), nrow = nrow(f), byrow = T)
  #QC.index = which(p$sampleType == Norm_datatype)SampleType
  
  QC.index = which(p$sampleType %in% SampleType)
  
  # parallel computing.
  # check if it is example.xlsx, if so, give result directly.
  # SERRF normalization
  
  norm = SERRF_norm(e, f, p, batch, QC.index, time = "time",cl=cl)
  
  
  
  # evaluation methods.
  # RSD
  RSD = function(e,f,p,robust = F,cl){
    if(robust){
      result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,remove_outlier,e){
        x = remove_outlier(e[i,])[[1]]
        sd(x,na.rm=T)/mean(x,na.rm=T)
      },remove_outlier,e)
    }else{
      result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,e){
        x = e[i,]
        sd(x,na.rm=T)/mean(x,na.rm=T)
      },e)
    }
    
    
    return(result)
  }
  SERRF.validate = RSD(norm$e[,p$sampleType=="Sample"],f,p[p$sampleType=="Sample",],cl=cl)
  raw.validate = RSD(e[,p$sampleType=="Sample"],f,p[p$sampleType=="Sample",],cl=cl)
  
  if (sum(is.na(SERRF.validate))==length(SERRF.validate)){
    SERRF.validate = RSD(norm$e[,p$sampleType=="QC"],f,p[p$sampleType=="QC",],cl=cl)
    raw.validate = RSD(e[,p$sampleType=="QC"],f,p[p$sampleType=="QC",],cl=cl)
  }
  # stopCluster(cl)
  # PCA
  # generate PCA plot.
  generate_PCA = function(e, f, p, QC.index, batch, method){
    rownames(e)=f$label
    colnames(e)=p$label
    
    for(i in 1:nrow(e)){
      if (base::max(e[i,!is.na(e[i,])])!=0){
        e[i,is.na(e[i,])] = min(e[i,!is.na(e[i,])])}}
    
    sds = apply(e,1,sd,na.rm = T)
    sd_pos = ifelse(and(sds>0 ,!is.na(sds)),T,F) 
    
    if (!is.null(sd_pos)){
      pca = prcomp(t(e[sd_pos,]), center = T, scale. = T)
      #ggbiplot(pca)
      variance = pca$sdev^2/sum(pca$sdev^2)
      pca.data = data.frame(pca$x,batch = batch[1,],order = 1:nrow(pca$x))
      batch.QC = batch[1,];
      batch.QC[p$sampleType=="QC"] = "QC"
      qc = rep(F, nrow(p))
      qc[p$sampleType=="QC"] = TRUE
      PC1_PC2=ggplot(pca.data, aes(PC1, PC2, color = batch.QC,group=batch.QC, size = ifelse(qc==T,2,1), order = order)) +
        geom_point(alpha = 0.8) + scale_size(guide="none",breaks =c(1, 2),range =c(1, 2)) +
        stat_ellipse( linetype = 2, size = 0.5) +
        labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"),
             title = method) +
        xlim(-25, 25) +
        ylim(-25, 25) +
        theme.scatter
    plot_pcs<-function(e){
      pc.df=data.frame(e,stringsAsFactors = F)
      pc.df[,"median"]=0
      for (row in rownames(pc.df)){
        pc.df[row,"median"]=median(as.numeric(pc.df[row,])) 
        
      }
      pc.df=head(pc.df[rev(order(pc.df$median)),],3)
      
      pc.df=pc.df[,-ncol(pc.df)]
      
      pc.df=as.data.frame(t(pc.df))
      
      pcplot.df=data.frame()
      
      for (col in 1:(ncol(pc.df))){
        
        pcplot.df=rbind(pcplot.df,data.frame(data=pc.df[[col]],time=1:nrow(pc.df),cpd=rep(colnames(pc.df)[col],nrow(pc.df)),type=factor(p$sampleType,levels=(c("qc","QC","Sample","sample"))),size=ifelse(qc ==FALSE,1,1.5),batch.QC=batch.QC))
        
      }
      #pcplot.df$size=as.factor(as.character(pcplot.df$size))
      #pc.df$timepoint=1:nrow(pc.df)
      outliers <- boxplot(pcplot.df$data, plot=FALSE)$out
      if (length(outliers)!=0){
        range=pcplot.df$data[-which(pcplot.df$data %in% outliers)]
        range=c(min(range),max(range))
      }else{
        range=c(min(pcplot.df$data),max(pcplot.df$data))
      }
      
      pcplot=ggplot(pcplot.df %>% arrange(type),aes(x=time,y=data,color=cpd,group=batch.QC, size = ifelse(pcplot.df$type %in% c("qc","QC"),2,1), fill =type)) + geom_point() +scale_size(guide="none",breaks =c(1, 2),range =c(1, 2)) +
            ylim(range)
      
      
      
      (pcplot)
    }
      
    #plot_pcs(e)
  
      
      
      
      return(list(PC1_PC2,plot_pcs(e)))
      
    }
  }
  SERRFpc = generate_PCA(norm$e[!is.na(SERRF.validate),],f[!is.na(SERRF.validate),],p,QC.index, batch , "SERRF")
  rawpc = generate_PCA(e,f,p,QC.index, batch , "raw")
  
  SERRFpca=SERRFpc[[1]]
  rawpca=rawpc[[1]]
  SERRFpcs=SERRFpc[[2]]
  rawpcs=rawpc[[2]]
  # par(mfrow = c(1,2))
  # max = max(c(e[i,], norm$e[i,]), na.rm = T)
  # min = min(c(e[i,], norm$e[i,]), na.rm = T)
  # plot(e[i,], col = factor(qc), ylim= c(min, max))
  # plot(norm$e[i,], col = factor(qc), ylim= c(min, max))
  # i = i + 1
  #
  #
  newwd=paste0(input,"SERRF_",paste(SampleType,collapse="+"))
  if (dir.exists(newwd)==FALSE){dir.create(newwd)}
  setwd(newwd)
  
  # save results.
  rownames(norm$e) = f$label
  
  write.csv(norm$e, "SERRF_normalized.csv")
  write.csv(data.frame(label=f$label,rawValidateRSD = raw.validate, SERRFValidateRSD = SERRF.validate), 'SERRF_validateRSD.csv')
  
  
  
  
  
  combine_figure <- ggarrange(rawpca, SERRFpca,
                              ncol = 2, nrow = 1,legend="right",common.legend = TRUE)
  ggsave("combine_figure.png",combine_figure, width = 20, height = 10,limitsize = FALSE)
  combine_figures <- ggarrange(rawpcs, SERRFpcs,
                              ncol = 2, nrow = 1,legend="bottom",common.legend = TRUE)
  ggsave("combine_figure_pcs.png",combine_figures, width = 20, height = 10,limitsize = FALSE)
  write.csv(e,"e.csv")
  write.csv(f,"f.csv")
  write.csv(p,"p.csv")
  Norm_finaltable<-Combine_result_file(newwd)
  write.table(Norm_finaltable[["raw"]], file="raw finaltb.csv",row.names=FALSE, col.names=FALSE, sep=",")
  write.table(Norm_finaltable[["norm"]], file="norm finaltb.csv",row.names=FALSE, col.names=FALSE, sep=",")
  #read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
  #  ph_with_vg(code = print(rawpca), type = "body", width = 10,
  #             height = 8, offx = 0, offy = 0) %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
  #  ph_with_vg(code = print(SERRFpca), type = "body", width = 10,
  #             height = 8, offx = 0, offy = 0) %>% print(target = "PCAs.pptx") %>%
  #  invisible()
  if (datatype=="MASSOMICS"){
    inputmassomics<-read.csv(input)
    inputmassomics<-inputmassomics[,grep("CAS",colnames(inputmassomics)):grep("Name",colnames(inputmassomics))]
    massomicsoupt<-merge(inputmassomics,Norm_finaltable[["norm"]],by.x="Name",by.y="label")[, union(names(inputmassomics), names(Norm_finaltable[["norm"]])[2:ncol(Norm_finaltable[["norm"]])])]
    dirname(input)
    write.table(massomicsoupt, file=paste0(dirname(input),"/Serrf_normed_as_",Predict_level,".csv"),row.names=FALSE, col.names=T, sep=",")
  }
  
  
  
  return(Norm_finaltable)
  
  
  # doc = pptx( )
  # doc = addSlide(doc, slide.layout = "Title and Content")
  # doc = addPlot(doc, fun = function() print(rawpca),
  #               vector.graphic = TRUE, width = 6, height = 6)
  # doc = addSlide(doc, slide.layout = "Title and Content")
  # doc = addPlot(doc, fun = function() print(SERRFpca),
  #               vector.graphic = TRUE, width = 6, height = 6)
  #
  # # write the document to a file
  # writeDoc(doc, file = "PCAs.pptx")
  
  #  zip(files = c(
  #    "SERRF_normalized.csv",
  #    'performance -  validateRSD.csv',
  #    "SERRFpca.png",
  #    "rawpca.png",
  #    "PCAs.pptx"
  #  ), zipfile = "SERRF - results.zip")
  #end = Sys.time()
  
  
  
  
  
  
  
  # get ip summary
  # ip_summ = t(info$num)
  # colnames(ip_summ) = info$code
  #
  # ip_summ = data.frame(ip_summ)
  #
  # ip_summ = jsonlite::toJSON(ip_summ)
  #ip_summ = NA
  
  
  
  
  
  #return(list(validateSERRF = signif(median(SERRF.validate),3)*100, validateraw = signif(median(raw.validate),3)*100, count_reduced = sum(SERRF.validate<raw.validate, na.rm = T), perc_reduced = sum(SERRF.validate<raw.validate, na.rm = T)/nrow(f), count_less_20_raw = sum(raw.validate<.2, na.rm = T), count_less_20_SERRF = sum(SERRF.validate<.2, na.rm = T), perc_less_20_raw = signif(sum(raw.validate<.2, na.rm = T)/nrow(f),3) * 100, perc_less_20_SERRF = signif(sum(SERRF.validate<.2, na.rm = T)/nrow(f),3)*100, runtime = signif(as.numeric(end - start)/60,3),ip_summ=ip_summ))
  
}

randomForest_norm_par<-function(j, data, batch, QC.index, time) {
  set.seed(1)
  #message(c("randomForest for ", rownames(f)[j]))
  data = data.frame(y = data[j, ], t(data[-j, ]), batch = batch, time = time)
  colnames(data) = c("y", paste0("X", 1:nrow(data))[-j],"batch", "time")
  model = randomForest::randomForest(y ~ ., data = data, subset = QC.index, importance = F)
  newdata = data.frame(t(data[-j, ]), batch = batch,time = time)
  colnames(newdata) = c(paste0("X", 1:nrow(data))[-j], "batch", "time")
  new = (data[j, ]/predict(model, newdata = newdata)) * median(data[j, ])
  return(new)
}

SERRF_native <- function(input = "Area.csv",Predict_level="QC",data=NULL){
  library(tcltk)
  if (and(!file.exists(input) , is.null(data))){input =tk_choose.files(caption = "Select area.csv files to normalized by SERRF")}
  SampleType=Predict_level
  setwd(base::dirname(input))
  # library(rgeolocate)
  #message(input)
  pacman::p_load(ggplot2, randomForest, parallel, officer, dplyr, rvg,data.table,ggpubr,extrafont)
  
  cl <- autoStopCluster(makeCluster(try(detectCores())))
  
  start = Sys.time()
  message(start)
  message(paste("Prediction level:",paste(SampleType,sep = "+",collapse = "+")))
  
  
  # info = read.csv(paste0("http://localhost:5984/serrf/info/info.csv"), stringsAsFactors = F,na.strings = "")
  #
  #
  # file <- system.file("extdata","ip2_sample.bin", package = "rgeolocate")
  #
  # code = ip2location(ip, file, c("country_code"))[[1]]
  # info$num[info$code == ip2location(ip, file, c("country_code"))[[1]]] =   info$num[info$code == code]+1
  #
  # put_att_csv = function(projectID = 'tryThu.Aug.17.14.53.35.2017', attname = 'test.csv', att = data.table::fread("G:\\initialize MetDA\\user_active.csv")){
  #   projectUrl <- paste0("http://localhost:5984/serrf/",projectID)
  #   projectList <- jsonlite::fromJSON(projectUrl)
  #
  #   new_att = projectList[["_attachments"]]
  #   new_att = new_att[!names(new_att)%in%attname]
  #   new_att[[attname]] = list(content_type="text/csv", data = RCurl::base64(
  #     paste0(R.utils::captureOutput(write.csv(att,stdout(), row.names=F)),collapse = "\n")
  #   ))
  #   projectList[["_attachments"]] = new_att
  #   result = RCurl::getURL(paste0("http://localhost:5984/serrf/",projectID),customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(projectList,auto_unbox = T, force = T))
  #   while(grepl("error",result)) {
  #     projectList <- jsonlite::fromJSON(projectUrl)
  #     new_att = projectList[["_attachments"]]
  #     new_att = new_att[!names(new_att)%in%attname]
  #     new_att[[attname]] = list(content_type="text/csv", data = RCurl::base64(
  #       paste0(R.utils::captureOutput(write.csv(att,stdout(), row.names=F)),collapse = "\n")
  #     ))
  #     projectList[["_attachments"]] = new_att
  #     result = RCurl::getURL(paste0("http://localhost:5984/serrf/",projectID),customrequest='PUT',httpheader=c('Content-Type'='application/json'),postfields= jsonlite::toJSON(projectList,auto_unbox = T, force = T))
  #     if(grepl("ok",result)){
  #       break;
  #     }
  #   }
  #   return(result)
  # }
  #
  # put_att_csv("info", "info.csv", info)
  
  if (and(file.exists(input) , is.null(data))) {data = readData(path=input)}
  
  e = data$e
  f = data$f
  p = data$p
  
  # impute missing value.
  if(sum(is.na(e))>0){
    missing_compounds = which(is.na(e), arr.ind = T)[,1]
    for(i in missing_compounds){
      e[i, is.na(e[i,])] = 1/2 * min(e[i,!is.na(e[i,])])
    }
  }
  e = data.matrix(e)
  
  
  theme.scatter = theme(
    plot.title = element_text(size = rel(2), hjust = 0.5,face = 'bold',family = "Arial"),#title size.
    # axis.title = element_text(size = rel(2)),
    axis.text	 = element_text(colour = 'black'),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.key = element_rect(fill = "white",colour = "white"),
    legend.title = element_text(face = 'bold'),
    text=element_text(family="Arial")
  )
  
  # define batch.
  batch = matrix(rep(p$batch, nrow(f)), nrow = nrow(f), byrow = T)
  #QC.index = which(p$sampleType == Norm_datatype)SampleType
  
  QC.index = which(p$sampleType %in% SampleType)
  
  # parallel computing.
  # check if it is example.xlsx, if so, give result directly.
  # SERRF normalization
  
  norm = SERRF_norm(e, f, p, batch, QC.index, time = "time",cl=cl)
  
  
  
  # evaluation methods.
  # RSD
  RSD = function(e,f,p,robust = F,cl){
    if(robust){
      result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,remove_outlier,e){
        x = remove_outlier(e[i,])[[1]]
        sd(x,na.rm=T)/mean(x,na.rm=T)
      },remove_outlier,e)
    }else{
      result=parSapply(cl=cl,X=1:nrow(e),FUN = function(i,e){
        x = e[i,]
        sd(x,na.rm=T)/mean(x,na.rm=T)
      },e)
    }
    
    
    return(result)
  }
  SERRF.validate = RSD(norm$e[,p$sampleType=="Sample"],f,p[p$sampleType=="Sample",],cl=cl,robust = T)
  raw.validate = RSD(e[,p$sampleType=="Sample"],f,p[p$sampleType=="Sample",],cl=cl)
  # stopCluster(cl)
  # PCA
  # generate PCA plot.
  generate_PCA = function(e, f, p, QC.index, batch, method){
    
    
    for(i in 1:nrow(e)){
      if (base::max(e[i,!is.na(e[i,])])!=0){
        e[i,is.na(e[i,])] = min(e[i,!is.na(e[i,])])}}
    
    sds = apply(e,1,sd,na.rm = T)
    sd_pos = sds>0 
    if (!is.null(sd_pos)){
      pca = prcomp(t(e[sd_pos,]), center = T, scale. = T)
      variance = pca$sdev^2/sum(pca$sdev^2)
      pca.data = data.frame(pca$x,batch = batch[1,],order = 1:nrow(pca$x))
      batch.QC = batch[1,];
      batch.QC[p$sampleType=="QC"] = "QC"
      qc = rep(F, nrow(p))
      qc[p$sampleType=="QC"] = TRUE
      ggplot(pca.data, aes(PC1, PC2, color = batch.QC,group=batch.QC, size = qc, order = order)) +
        geom_point(alpha = 0.8) +
        stat_ellipse( linetype = 2, size = 0.5) +
        labs(x = paste0("PC1: ",signif(variance[1]*100,3),"%"), y = paste0("PC2: ",signif(variance[2]*100,3),"%"),
             title = method) +
        xlim(-25, 25) +
        ylim(-25, 25) +
        theme.scatter
    }
  }
  SERRFpca = generate_PCA(norm$e[!is.na(SERRF.validate),],f,p,QC.index, batch , "SERRF")
  rawpca = generate_PCA(e,f,p,QC.index, batch , "raw")
  
  
  
  
  # par(mfrow = c(1,2))
  # max = max(c(e[i,], norm$e[i,]), na.rm = T)
  # min = min(c(e[i,], norm$e[i,]), na.rm = T)
  # plot(e[i,], col = factor(qc), ylim= c(min, max))
  # plot(norm$e[i,], col = factor(qc), ylim= c(min, max))
  # i = i + 1
  #
  #
  newwd=paste0(input,"SERRF_",paste(SampleType,collapse="+"))
  if (dir.exists(newwd)==FALSE){dir.create(newwd)}
  setwd(newwd)
  
  # save results.
  rownames(norm$e) = f$label
  write.csv(norm$e, "SERRF_normalized.csv")
  write.csv(data.frame(label=f$label,rawValidateRSD = raw.validate, SERRFValidateRSD = SERRF.validate), 'SERRF_validateRSD.csv')
  ggsave("SERRFpca.png",SERRFpca, width = 8, height = 8)
  ggsave("rawpca.png",rawpca, width = 8, height = 8)
  combine_figure <- ggarrange(rawpca, SERRFpca, 
                              labels = c("RAW", "SERRF"),
                              ncol = 2, nrow = 1)
  ggsave("combine_figure.png",combine_figure, width = 16, height = 8)
  #read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
  #  ph_with_vg(code = print(rawpca), type = "body", width = 10,
  #             height = 8, offx = 0, offy = 0) %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
  #  ph_with_vg(code = print(SERRFpca), type = "body", width = 10,
  #             height = 8, offx = 0, offy = 0) %>% print(target = "PCAs.pptx") %>%
  #  invisible()
  
  
  
  
  # doc = pptx( )
  # doc = addSlide(doc, slide.layout = "Title and Content")
  # doc = addPlot(doc, fun = function() print(rawpca),
  #               vector.graphic = TRUE, width = 6, height = 6)
  # doc = addSlide(doc, slide.layout = "Title and Content")
  # doc = addPlot(doc, fun = function() print(SERRFpca),
  #               vector.graphic = TRUE, width = 6, height = 6)
  #
  # # write the document to a file
  # writeDoc(doc, file = "PCAs.pptx")
  
  #  zip(files = c(
  #    "SERRF_normalized.csv",
  #    'performance -  validateRSD.csv',
  #    "SERRFpca.png",
  #    "rawpca.png",
  #    "PCAs.pptx"
  #  ), zipfile = "SERRF - results.zip")
  end = Sys.time()
  
  
  
  
  
  
  
  # get ip summary
  # ip_summ = t(info$num)
  # colnames(ip_summ) = info$code
  #
  # ip_summ = data.frame(ip_summ)
  #
  # ip_summ = jsonlite::toJSON(ip_summ)
  ip_summ = NA
  
  
  
  
  
  #return(list(validateSERRF = signif(median(SERRF.validate),3)*100, validateraw = signif(median(raw.validate),3)*100, count_reduced = sum(SERRF.validate<raw.validate, na.rm = T), perc_reduced = sum(SERRF.validate<raw.validate, na.rm = T)/nrow(f), count_less_20_raw = sum(raw.validate<.2, na.rm = T), count_less_20_SERRF = sum(SERRF.validate<.2, na.rm = T), perc_less_20_raw = signif(sum(raw.validate<.2, na.rm = T)/nrow(f),3) * 100, perc_less_20_SERRF = signif(sum(SERRF.validate<.2, na.rm = T)/nrow(f),3)*100, runtime = signif(as.numeric(end - start)/60,3),ip_summ=ip_summ))
  
}

predict.randomForest <-  function (object, newdata, type = "response", norm.votes = TRUE,
            predict.all=FALSE, proximity = FALSE, nodes=FALSE, cutoff, ...)
  {
    if (!inherits(object, "randomForest"))      stop("object not of class randomForest")
    if (is.null(object$forest)) stop("No forest component in the object")
    out.type <- charmatch(tolower(type),
                          c("response", "prob", "vote", "class"))
    if (is.na(out.type))
      stop("type must be one of 'response', 'prob', 'vote'")
    if (out.type == 4) out.type <- 1
    if (out.type != 1 && object$type == "regression")
      stop("'prob' or 'vote' not meaningful for regression")
    if (out.type == 2) norm.votes <- TRUE
    if (missing(newdata)) {
      p <- if (! is.null(object$na.action)) {
        napredict(object$na.action, object$predicted)
      } else {
        object$predicted
      }
      if (object$type == "regression") return(p)
      if (proximity & is.null(object$proximity))
        warning("cannot return proximity without new data if random forest object does not already have proximity")
      if (out.type == 1) {
        if (proximity) {
          return(list(pred = p,
                      proximity = object$proximity))
        } else return(p)
      }
      v <- object$votes
      if (!is.null(object$na.action)) v <- napredict(object$na.action, v)
      if (norm.votes) {
        t1 <- t(apply(v, 1, function(x) { x/sum(x) }))
        class(t1) <- c(class(t1), "votes")
        if (proximity) return(list(pred = t1, proximity = object$proximity))
        else return(t1)
      } else {
        if (proximity) return(list(pred = v, proximity = object$proximity))
        else return(v)
      }
    }
    if (missing(cutoff)) {
      cutoff <- object$forest$cutoff
    } else {
      if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
          length(cutoff) != length(object$classes)) {
        stop("Incorrect cutoff specified.")
      }
      if (!is.null(names(cutoff))) {
        if (!all(names(cutoff) %in% object$classes)) {
          stop("Wrong name(s) for cutoff")
        }
        cutoff <- cutoff[object$classes]
      }
    }
    
    if (object$type == "unsupervised")
      stop("Can't predict unsupervised forest.")
    
    if (inherits(object, "randomForest.formula")) {
      newdata <- as.data.frame(newdata)
      rn <- row.names(newdata)
      Terms <- delete.response(object$terms)
      x <- model.frame(Terms, newdata, na.action = na.omit)
      keep <- match(row.names(x), rn)
    } else {
      if (is.null(dim(newdata)))
        dim(newdata) <- c(1, length(newdata))
      x <- newdata
      if (nrow(x) == 0)
        stop("newdata has 0 rows")
      if (any(is.na(x)))
        stop("missing values in newdata")
      keep <- 1:nrow(x)
      rn <- rownames(x)
      if (is.null(rn)) rn <- keep
    }
    vname <- if (is.null(dim(object$importance))) {
      names(object$importance)
    } else {
      rownames(object$importance)
    }
    if (is.null(colnames(x))) {
      if (ncol(x) != length(vname)) {
        stop("number of variables in newdata does not match that in the training data")
      }
    } else {
      if (any(! vname %in% colnames(x)))
        stop("variables in the training data missing in newdata")
      x <- x[, vname, drop=FALSE]
    }
    if (is.data.frame(x)) {
      isFactor <- function(x) is.factor(x) & ! is.ordered(x)
      xfactor <- which(sapply(x, isFactor))
      if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
        for (i in xfactor) {
          if (any(! levels(x[[i]]) %in% object$forest$xlevels[[i]]))
            stop("New factor levels not present in the training data")
          x[[i]] <-
            factor(x[[i]],
                   levels=levels(x[[i]])[match(levels(x[[i]]), object$forest$xlevels[[i]])])
        }
      }
      cat.new <- sapply(x, function(x) if (is.factor(x) && !is.ordered(x))
        length(levels(x)) else 1)
      if (!all(object$forest$ncat == cat.new))
        stop("Type of predictors in new data do not match that of the training data.")
    }
    mdim <- ncol(x)
    ntest <- nrow(x)
    ntree <- object$forest$ntree
    maxcat <- max(object$forest$ncat)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))
    x <- t(data.matrix(x))
    
    if (predict.all) {
      treepred <- if (object$type == "regression") {
        matrix(double(ntest * ntree), ncol=ntree)
      } else {
        matrix(integer(ntest * ntree), ncol=ntree)
      }
    } else {
      treepred <- numeric(ntest)
    }
    proxmatrix <- if (proximity) matrix(0, ntest, ntest) else numeric(1)
    nodexts <- if (nodes) integer(ntest * ntree) else integer(ntest)
    
    if (object$type == "regression") {
      if (!is.null(object$forest$treemap)) {
        object$forest$leftDaughter <-
          object$forest$treemap[,1,, drop=FALSE]
        object$forest$rightDaughter <-
          object$forest$treemap[,2,, drop=FALSE]
        object$forest$treemap <- NULL
      }
      
      keepIndex <- "ypred"
      if (predict.all) keepIndex <- c(keepIndex, "treepred")
      if (proximity) keepIndex <- c(keepIndex, "proximity")
      if (nodes) keepIndex <- c(keepIndex, "nodexts")
      ## Ensure storage mode is what is expected in C.
      if (! is.integer(object$forest$leftDaughter))
        storage.mode(object$forest$leftDaughter) <- "integer"
      if (! is.integer(object$forest$rightDaughter))
        storage.mode(object$forest$rightDaughter) <- "integer"
      if (! is.integer(object$forest$nodestatus))
        storage.mode(object$forest$nodestatus) <- "integer"
      if (! is.double(object$forest$xbestsplit))
        storage.mode(object$forest$xbestsplit) <- "double"
      if (! is.double(object$forest$nodepred))
        storage.mode(object$forest$nodepred) <- "double"
      if (! is.integer(object$forest$bestvar))
        storage.mode(object$forest$bestvar) <- "integer"
      if (! is.integer(object$forest$ndbigtree))
        storage.mode(object$forest$ndbigtree) <- "integer"
      if (! is.integer(object$forest$ncat))
        storage.mode(object$forest$ncat) <- "integer"
      
      ans <- .C("regForest",
                as.double(x),
                ypred = double(ntest),
                as.integer(mdim),
                as.integer(ntest),
                as.integer(ntree),
                object$forest$leftDaughter,
                object$forest$rightDaughter,
                object$forest$nodestatus,
                nrnodes,
                object$forest$xbestsplit,
                object$forest$nodepred,
                object$forest$bestvar,
                object$forest$ndbigtree,
                object$forest$ncat,
                as.integer(maxcat),
                as.integer(predict.all),
                treepred = as.double(treepred),
                as.integer(proximity),
                proximity = as.double(proxmatrix),
                nodes = as.integer(nodes),
                nodexts = as.integer(nodexts),
                #DUP=FALSE,
                PACKAGE = "randomForest")[keepIndex]
      ## Apply bias correction if needed.
      yhat <- rep(NA, length(rn))
      names(yhat) <- rn
      if (!is.null(object$coefs)) {
        yhat[keep] <- object$coefs[1] + object$coefs[2] * ans$ypred
      } else {
        yhat[keep] <- ans$ypred
      }
      if (predict.all) {
        treepred <- matrix(NA, length(rn), ntree,
                           dimnames=list(rn, NULL))
        treepred[keep,] <- ans$treepred
      }
      if (!proximity) {
        res <- if (predict.all)
          list(aggregate=yhat, individual=treepred) else yhat
      } else {
        res <- list(predicted = yhat,
                    proximity = structure(ans$proximity,
                                          dim=c(ntest, ntest), dimnames=list(rn, rn)))
      }
      if (nodes) {
        attr(res, "nodes") <- matrix(ans$nodexts, ntest, ntree,
                                     dimnames=list(rn[keep], 1:ntree))
      }
    } else {
      countts <- matrix(0, ntest, nclass)
      t1 <- .C("classForest",
               mdim = as.integer(mdim),
               ntest = as.integer(ntest),
               nclass = as.integer(object$forest$nclass),
               maxcat = as.integer(maxcat),
               nrnodes = as.integer(nrnodes),
               jbt = as.integer(ntree),
               xts = as.double(x),
               xbestsplit = as.double(object$forest$xbestsplit),
               pid = object$forest$pid,
               cutoff = as.double(cutoff),
               countts = as.double(countts),
               treemap = as.integer(aperm(object$forest$treemap,
                                          c(2, 1, 3))),
               nodestatus = as.integer(object$forest$nodestatus),
               cat = as.integer(object$forest$ncat),
               nodepred = as.integer(object$forest$nodepred),
               treepred = as.integer(treepred),
               jet = as.integer(numeric(ntest)),
               bestvar = as.integer(object$forest$bestvar),
               nodexts = as.integer(nodexts),
               ndbigtree = as.integer(object$forest$ndbigtree),
               predict.all = as.integer(predict.all),
               prox = as.integer(proximity),
               proxmatrix = as.double(proxmatrix),
               nodes = as.integer(nodes),
               #DUP=FALSE,
               PACKAGE = "randomForest")
      if (out.type > 1) {
        out.class.votes <- t(matrix(t1$countts, nrow = nclass, ncol = ntest))
        if (norm.votes)
          out.class.votes <-
            sweep(out.class.votes, 1, rowSums(out.class.votes), "/")
        z <- matrix(NA, length(rn), nclass,
                    dimnames=list(rn, object$classes))
        z[keep, ] <- out.class.votes
        class(z) <- c(class(z), "votes")
        res <- z
      } else {
        out.class <- factor(rep(NA, length(rn)),
                            levels=1:length(object$classes),
                            labels=object$classes)
        out.class[keep] <- object$classes[t1$jet]
        names(out.class)[keep] <- rn[keep]
        res <- out.class
      }
      if (predict.all) {
        treepred <- matrix(object$classes[t1$treepred],
                           nrow=length(keep), dimnames=list(rn[keep], NULL))
        res <- list(aggregate=res, individual=treepred)
      }
      if (proximity)
        res <- list(predicted = res, proximity = structure(t1$proxmatrix,
                                                           dim = c(ntest, ntest),
                                                           dimnames = list(rn[keep], rn[keep])))
      if (nodes) attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
                                              dimnames=list(rn[keep], 1:ntree))
    }
    res
}
mylevels <- function(x) if (is.factor(x)) levels(x) else 0
randomForest.default <- function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree=500,
           mtry=if (!is.null(y) && !is.factor(y))
             max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
           replace=TRUE, classwt=NULL, cutoff, strata,
           sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, ...) {
  RFpars<-function(){
    xtest=NULL
    ytest=NULL
    ntree=500
    mtry=if (!is.null(y) && !is.factor(y)) max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))
    replace=TRUE
    classwt=NULL
    cutoff
    strata
    sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x))
    nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1
    maxnodes=NULL
    importance=FALSE
    localImp=FALSE
    nPerm=1
    proximity
    oob.prox=proximity
    norm.votes=TRUE
    do.trace=FALSE
    keep.forest=!is.null(y) && is.null(xtest)
    corr.bias=FALSE
    keep.inbag=FALSE
    
    
  }
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest
    
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(x) != ncol(xtest))
        stop("x and xtest must have same number of columns")
      ntest <- nrow(xtest)
      xts.row.names <- rownames(xtest)
    }
    
    ## Make sure mtry is in reasonable range.
    if (mtry < 1 || mtry > p)
      warning("invalid mtry: reset to within valid range")
    mtry <- max(1, min(p, round(mtry)))
    if (!is.null(y)) {
      if (length(y) != n) stop("length of response must be the same as predictors")
      addclass <- FALSE
    } else {
      if (!addclass) addclass <- TRUE
      y <- factor(c(rep(1, n), rep(2, n)))
      x <- rbind(x, x)
    }
    
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
    
    if (is.data.frame(x)) {
      xlevels <- lapply(x, mylevels)
      ncat <- sapply(xlevels, length)
      ## Treat ordered factors as numerics.
      ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
      x <- data.matrix(x)
      if(testdat) {
        if(!is.data.frame(xtest))
          stop("xtest must be data frame if x is")
        xfactor <- which(sapply(xtest, is.factor))
        if (length(xfactor) > 0) {
          for (i in xfactor) {
            if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
              stop("New factor levels in xtest not present in x")
            xtest[[i]] <-
              factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])],
                     levels=xlevels[[i]])
          }
        }
        xtest <- data.matrix(xtest)
      }
    } else {
      ncat <- rep(1, p)
      names(ncat) <- colnames(x)
      xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 53)
      stop("Can not handle categorical predictors with more than 53 categories.")
    
    if (classRF) {
      nclass <- length(levels(y))
      ## Check for empty classes:
      if (any(table(y) == 0)) stop("Can't have empty classes in y.")
      if (!is.null(ytest)) {
        if (!is.factor(ytest)) stop("ytest must be a factor")
        if (!all(levels(y) == levels(ytest)))
          stop("y and ytest must have the same levels")
      }
      if (missing(cutoff)) {
        cutoff <- rep(1 / nclass, nclass)
      } else {
        if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
            length(cutoff) != nclass) {
          stop("Incorrect cutoff specified.")
        }
        if (!is.null(names(cutoff))) {
          if (!all(names(cutoff) %in% levels(y))) {
            stop("Wrong name(s) for cutoff")
          }
          cutoff <- cutoff[levels(y)]
        }
      }
      if (!is.null(classwt)) {
        if (length(classwt) != nclass)
          stop("length of classwt not equal to number of classes")
        ## If classwt has names, match to class labels.
        if (!is.null(names(classwt))) {
          if (!all(names(classwt) %in% levels(y))) {
            stop("Wrong name(s) for classwt")
          }
          classwt <- classwt[levels(y)]
        }
        if (any(classwt <= 0)) stop("classwt must be positive")
        ipi <- 1
      } else {
        classwt <- rep(1, nclass)
        ipi <- 0
      }
    } else addclass <- FALSE
    
    if (missing(proximity)) proximity <- addclass
    if (proximity) {
      prox <- matrix(0.0, n, n)
      proxts <- if (testdat) matrix(0, ntest, ntest + n) else double(1)
    } else {
      prox <- proxts <- double(1)
    }
    
    if (localImp) {
      importance <- TRUE
      impmat <- matrix(0, p, n)
    } else impmat <- double(1)
    
    if (importance) {
      if (nPerm < 1) nPerm <- as.integer(1) else nPerm <- as.integer(nPerm)
      if (classRF) {
        impout <- matrix(0.0, p, nclass + 2)
        impSD <- matrix(0.0, p, nclass + 1)
      } else {
        impout <- matrix(0.0, p, 2)
        impSD <- double(p)
        names(impSD) <- x.col.names
      }
    } else {
      impout <- double(p)
      impSD <- double(1)
    }
    
    nsample <- if (addclass) 2 * n else n
    Stratify <- length(sampsize) > 1
    if ((!Stratify) && sampsize > nrow(x)) stop("sampsize too large")
    if (Stratify && (!classRF)) stop("sampsize should be of length one")
    if (classRF) {
      if (Stratify) {
        if (missing(strata)) strata <- y
        if (!is.factor(strata)) strata <- as.factor(strata)
        nsum <- sum(sampsize)
        if (length(sampsize) > nlevels(strata))
          stop("sampsize has too many elements.")
        if (any(sampsize <= 0) || nsum == 0)
          stop("Bad sampsize specification")
        ## If sampsize has names, match to class labels.
        if (!is.null(names(sampsize))) {
          sampsize <- sampsize[levels(strata)]
        }
        if (any(sampsize > table(strata)))
          stop("sampsize can not be larger than class frequency")
      } else {
        nsum <- sampsize
      }
      nrnodes <- 2 * nsum + 1
    } else {
      ## For regression trees, need to do this to get maximal trees.
      nrnodes <- 2 * sampsize + 1
    }
    if (!is.null(maxnodes)) {
      ## convert # of terminal nodes to total # of nodes
      maxnodes <- 2 * maxnodes - 1
      if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
      nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    ## Compiled code expects variables in rows and observations in columns.
    x <- t(x)
    storage.mode(x) <- "double"
    if (testdat) {
      xtest <- t(xtest)
      storage.mode(xtest) <- "double"
      if (is.null(ytest)) {
        ytest <- labelts <- 0
      } else {
        labelts <- TRUE
      }
    } else {
      xtest <- double(1)
      ytest <- double(1)
      ntest <- 1
      labelts <- FALSE
    }
    nt <- if (keep.forest) ntree else 1
    
    if (classRF) {
      cwt <- classwt
      threshold <- cutoff
      error.test <- if (labelts) double((nclass+1) * ntree) else double(1)
      rfout <- .C("classRF",
                  x = x,
                  xdim = as.integer(c(p, n)),
                  y = as.integer(y),
                  nclass = as.integer(nclass),
                  ncat = as.integer(ncat),
                  maxcat = as.integer(maxcat),
                  sampsize = as.integer(sampsize),
                  strata = if (Stratify) as.integer(strata) else integer(1),
                  Options = as.integer(c(addclass,
                                         importance,
                                         localImp,
                                         proximity,
                                         oob.prox,
                                         do.trace,
                                         keep.forest,
                                         replace,
                                         Stratify,
                                         keep.inbag)),
                  ntree = as.integer(ntree),
                  mtry = as.integer(mtry),
                  ipi = as.integer(ipi),
                  classwt = as.double(cwt),
                  cutoff = as.double(threshold),
                  nodesize = as.integer(nodesize),
                  outcl = integer(nsample),
                  counttr = integer(nclass * nsample),
                  prox = prox,
                  impout = impout,
                  impSD = impSD,
                  impmat = impmat,
                  nrnodes = as.integer(nrnodes),
                  ndbigtree = integer(ntree),
                  nodestatus = integer(nt * nrnodes),
                  bestvar = integer(nt * nrnodes),
                  treemap = integer(nt * 2 * nrnodes),
                  nodepred = integer(nt * nrnodes),
                  xbestsplit = double(nt * nrnodes),
                  errtr = double((nclass+1) * ntree),
                  testdat = as.integer(testdat),
                  xts = as.double(xtest),
                  clts = as.integer(ytest),
                  nts = as.integer(ntest),
                  countts = double(nclass * ntest),
                  outclts = as.integer(numeric(ntest)),
                  labelts = as.integer(labelts),
                  proxts = proxts,
                  errts = error.test,
                  inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(n),
                  #DUP=FALSE,
                  PACKAGE="randomForest")[-1]
      if (keep.forest) {
        ## deal with the random forest outputs
        max.nodes <- max(rfout$ndbigtree)
        treemap <- aperm(array(rfout$treemap, dim = c(2, nrnodes, ntree)),
                         c(2, 1, 3))[1:max.nodes, , , drop=FALSE]
      }
      if (!addclass) {
        ## Turn the predicted class into a factor like y.
        out.class <- factor(rfout$outcl, levels=1:nclass,
                            labels=levels(y))
        names(out.class) <- x.row.names
        con <- table(observed = y,
                     predicted = out.class)[levels(y), levels(y)]
        con <- cbind(con, class.error = 1 - diag(con)/rowSums(con))
      }
      out.votes <- t(matrix(rfout$counttr, nclass, nsample))[1:n, ]
      oob.times <- rowSums(out.votes)
      if (norm.votes)
        out.votes <- t(apply(out.votes, 1, function(x) x/sum(x)))
      dimnames(out.votes) <- list(x.row.names, levels(y))
      class(out.votes) <- c(class(out.votes), "votes")
      if (testdat) {
        out.class.ts <- factor(rfout$outclts, levels=1:nclass,
                               labels=levels(y))
        names(out.class.ts) <- xts.row.names
        out.votes.ts <- t(matrix(rfout$countts, nclass, ntest))
        dimnames(out.votes.ts) <- list(xts.row.names, levels(y))
        if (norm.votes)
          out.votes.ts <- t(apply(out.votes.ts, 1,
                                  function(x) x/sum(x)))
        class(out.votes.ts) <- c(class(out.votes.ts), "votes")
        if (labelts) {
          testcon <- table(observed = ytest,
                           predicted = out.class.ts)[levels(y), levels(y)]
          testcon <- cbind(testcon,
                           class.error = 1 - diag(testcon)/rowSums(testcon))
        }
      }
      cl <- match.call()
      cl[[1]] <- as.name("randomForest")
      out <- list(call = cl,
                  type = if (addclass) "unsupervised" else "classification",
                  predicted = if (addclass) NULL else out.class,
                  err.rate = if (addclass) NULL else t(matrix(rfout$errtr,
                                                              nclass+1, ntree,
                                                              dimnames=list(c("OOB", levels(y)), NULL))),
                  confusion = if (addclass) NULL else con,
                  votes = out.votes,
                  oob.times = oob.times,
                  classes = levels(y),
                  importance = if (importance)
                    matrix(rfout$impout, p, nclass+2,
                           dimnames = list(x.col.names,
                                           c(levels(y), "MeanDecreaseAccuracy",
                                             "MeanDecreaseGini")))
                  else matrix(rfout$impout, ncol=1,
                              dimnames=list(x.col.names, "MeanDecreaseGini")),
                  importanceSD = if (importance)
                    matrix(rfout$impSD, p, nclass + 1,
                           dimnames = list(x.col.names,
                                           c(levels(y), "MeanDecreaseAccuracy")))
                  else NULL,
                  localImportance = if (localImp)
                    matrix(rfout$impmat, p, n,
                           dimnames = list(x.col.names,x.row.names)) else NULL,
                  proximity = if (proximity) matrix(rfout$prox, n, n,
                                                    dimnames = list(x.row.names, x.row.names)) else NULL,
                  ntree = ntree,
                  mtry = mtry,
                  forest = if (!keep.forest) NULL else {
                    list(ndbigtree = rfout$ndbigtree,
                         nodestatus = matrix(rfout$nodestatus,
                                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                         bestvar = matrix(rfout$bestvar, ncol = ntree)[1:max.nodes,, drop=FALSE],
                         treemap = treemap,
                         nodepred = matrix(rfout$nodepred,
                                           ncol = ntree)[1:max.nodes,, drop=FALSE],
                         xbestsplit = matrix(rfout$xbestsplit,
                                             ncol = ntree)[1:max.nodes,, drop=FALSE],
                         pid = rfout$classwt, cutoff=cutoff, ncat=ncat,
                         maxcat = maxcat,
                         nrnodes = max.nodes, ntree = ntree,
                         nclass = nclass, xlevels=xlevels)
                  },
                  y = if (addclass) NULL else y,
                  test = if(!testdat) NULL else list(
                    predicted = out.class.ts,
                    err.rate = if (labelts) t(matrix(rfout$errts, nclass+1,
                                                     ntree,
                                                     dimnames=list(c("Test", levels(y)), NULL))) else NULL,
                    confusion = if (labelts) testcon else NULL,
                    votes = out.votes.ts,
                    proximity = if(proximity) matrix(rfout$proxts, nrow=ntest,
                                                     dimnames = list(xts.row.names, c(xts.row.names,
                                                                                      x.row.names))) else NULL),
                  inbag = if (keep.inbag) matrix(rfout$inbag, nrow=nrow(rfout$inbag),
                                                 dimnames=list(x.row.names, NULL)) else NULL)
    } else {
      ymean <- mean(y)
      y <- y - ymean
      ytest <- ytest - ymean
      rfout <- .C("regRF",
                  x,
                  as.double(y),
                  as.integer(c(n, p)),
                  as.integer(sampsize),
                  as.integer(nodesize),
                  as.integer(nrnodes),
                  as.integer(ntree),
                  as.integer(mtry),
                  as.integer(c(importance, localImp, nPerm)),
                  as.integer(ncat),
                  as.integer(maxcat),
                  as.integer(do.trace),
                  as.integer(proximity),
                  as.integer(oob.prox),
                  as.integer(corr.bias),
                  ypred = double(n),
                  impout = impout,
                  impmat = impmat,
                  impSD = impSD,
                  prox = prox,
                  ndbigtree = integer(ntree),
                  nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                  leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                  rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                  nodepred = matrix(double(nrnodes * nt), ncol=nt),
                  bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                  xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                  mse = double(ntree),
                  keep = as.integer(c(keep.forest, keep.inbag)),
                  replace = as.integer(replace),
                  testdat = as.integer(testdat),
                  xts = xtest,
                  ntest = as.integer(ntest),
                  yts = as.double(ytest),
                  labelts = as.integer(labelts),
                  ytestpred = double(ntest),
                  proxts = proxts,
                  msets = double(if (labelts) ntree else 1),
                  coef = double(2),
                  oob.times = integer(n),
                  inbag = if (keep.inbag)
                    matrix(integer(n * ntree), n) else integer(1),
                  #DUP=FALSE,
                  PACKAGE="randomForest")[c(16:28, 36:41)]
      ## Format the forest component, if present.
      if (keep.forest) {
        max.nodes <- max(rfout$ndbigtree)
        rfout$nodestatus <-
          rfout$nodestatus[1:max.nodes, , drop=FALSE]
        rfout$bestvar <-
          rfout$bestvar[1:max.nodes, , drop=FALSE]
        rfout$nodepred <-
          rfout$nodepred[1:max.nodes, , drop=FALSE] + ymean
        rfout$xbestsplit <-
          rfout$xbestsplit[1:max.nodes, , drop=FALSE]
        rfout$leftDaughter <-
          rfout$leftDaughter[1:max.nodes, , drop=FALSE]
        rfout$rightDaughter <-
          rfout$rightDaughter[1:max.nodes, , drop=FALSE]
      }
      cl <- match.call()
      cl[[1]] <- as.name("randomForest")
      ## Make sure those obs. that have not been OOB get NA as prediction.
      ypred <- rfout$ypred
      if (any(rfout$oob.times < 1)) {
        ypred[rfout$oob.times == 0] <- NA
      }
      out <- list(call = cl,
                  type = "regression",
                  predicted = structure(ypred + ymean, names=x.row.names),
                  mse = rfout$mse,
                  rsq = 1 - rfout$mse / (var(y) * (n-1) / n),
                  oob.times = rfout$oob.times,
                  importance = if (importance) matrix(rfout$impout, p, 2,
                                                      dimnames=list(x.col.names,
                                                                    c("%IncMSE","IncNodePurity"))) else
                                                                      matrix(rfout$impout, ncol=1,
                                                                             dimnames=list(x.col.names, "IncNodePurity")),
                  importanceSD=if (importance) rfout$impSD else NULL,
                  localImportance = if (localImp)
                    matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                             x.row.names)) else NULL,
                  proximity = if (proximity) matrix(rfout$prox, n, n,
                                                    dimnames = list(x.row.names, x.row.names)) else NULL,
                  ntree = ntree,
                  mtry = mtry,
                  forest = if (keep.forest)
                    c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                              "rightDaughter", "nodepred", "bestvar",
                              "xbestsplit")],
                      list(ncat = ncat), list(nrnodes=max.nodes),
                      list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                  coefs = if (corr.bias) rfout$coef else NULL,
                  y = y + ymean,
                  test = if(testdat) {
                    list(predicted = structure(rfout$ytestpred + ymean,
                                               names=xts.row.names),
                         mse = if(labelts) rfout$msets else NULL,
                         rsq = if(labelts) 1 - rfout$msets /
                           (var(ytest) * (n-1) / n) else NULL,
                         proximity = if (proximity)
                           matrix(rfout$proxts / ntree, nrow = ntest,
                                  dimnames = list(xts.row.names,
                                                  c(xts.row.names,
                                                    x.row.names))) else NULL)
                  } else NULL,
                  inbag = if (keep.inbag)
                    matrix(rfout$inbag, nrow(rfout$inbag),
                           dimnames=list(x.row.names, NULL)) else NULL)
    }
    class(out) <- "randomForest"
    return(out)
}

RFtutorial<-function(){
  
  #### Margin plot
  set.seed(1)
  data(iris)
  iris.rf <- randomForest(Species ~ ., iris, keep.forest=FALSE)
  plot(margin(iris.rf))
  
  #### MDS plot
  iris.rf <- randomForest(Species ~ ., iris, proximity=TRUE,
                          keep.forest=FALSE)
  MDSplot(iris.rf, iris$Species)
  ## Using different symbols for the classes:
  MDSplot(iris.rf, iris$Species, palette=rep(1, 3), pch=as.numeric(iris$Species))
  
  #### partialPlot
  data(iris)
  set.seed(543)
  iris.rf <- randomForest(Species~., iris)
  partialPlot(iris.rf, iris, Petal.Width, "versicolor")
  ## Looping over variables ranked by importance:
  data(airquality)
  airquality <- na.omit(airquality)
  set.seed(131)
  ozone.rf <- randomForest(Ozone ~ ., airquality, importance=TRUE)
  imp <- importance(ozone.rf)
  impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
  op <- par(mfrow=c(2, 3))
  for (i in seq_along(impvar)) {
    partialPlot(ozone.rf, airquality, impvar[i], xlab=impvar[i],
                main=paste("Partial Dependence on", impvar[i]),
                ylim=c(30, 70))
  }
  par(op)
  
  ### error plot
  data(mtcars)
  plot(randomForest(mpg ~ ., mtcars, keep.forest=FALSE, ntree=100), log="y")
  
  ### classification
  set.seed(71)
  iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE,
                          proximity=TRUE)
  print(iris.rf)
  ## Look at variable importance:
  round(importance(iris.rf), 2)
  ## Do MDS on 1 - proximity:
  iris.mds <- cmdscale(1 - iris.rf$proximity, eig=TRUE)
  op <- par(pty="s")
  pairs(cbind(iris[,1:4], iris.mds$points), cex=0.6, gap=0,
        col=c("red", "green", "blue")[as.numeric(iris$Species)],
        main="Iris Data: Predictors and MDS of Proximity Based on RandomForest")
  par(op)
  print(iris.mds$GOF)
  
  ### test
  ### classification
  model = randomForest::randomForest(batch ~ ., data = data, importance = T,keep.inbag=T,proximity=TRUE,keep.forest=F)
  
  model = randomForest::randomForest(y ~ ., data = data, subset = QC.index, importance = T,keep.inbag=T,proximity=TRUE,keep.forest=F)
  imp <- importance(model)
  roundimp=round(importance(model), 2)
  model.mds <- cmdscale(1 - model$proximity, eig=TRUE)
  
  op <- par(pty="s")
  pairs(cbind(data[,1:4], model.mds$points), cex=0.6, gap=0,
        col=c("red", "green", "blue")[as.numeric(data$y)],
        main="Model Data: Predictors and MDS of Proximity Based on RandomForest")
  par(op)
  print(model.mds$GOF)
  
  
  
  impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
  op <- par(mfrow=c(2, 3))
  for (i in 1:6) {
    partialPlot(model, data, impvar[i], xlab=impvar[i],
                main=paste("Partial Dependence on", impvar[i]))
  }
  
  MDSplot(model, as.factor(data$y))
  plot(model,log="y")
  
  for (i in 1:6) {
    ggplot(data=data, aes(x=))
  }
  
}

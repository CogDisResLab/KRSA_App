


list_of_rows <- function(fullDataSet) {
  sampleNameCol=grep("Sample name", fullDataSet)
  sampleNameRow=which(fullDataSet[,sampleNameCol]=="Sample name")
  cycleCol=grep("Cycle", fullDataSet)
  cycleRow=which(fullDataSet[,cycleCol]=="Cycle")
  exposureTimeCol=grep("Exposure time", fullDataSet)
  exposureTimeRow=which(fullDataSet[,exposureTimeCol]=="Exposure time")
  sequenceCol=grep("Sequence", fullDataSet)
  sequenceRow=which(fullDataSet[,sequenceCol]=="Sequence")
  substrateCol=unlist(as.numeric(which(apply(fullDataSet, 2, grep, pattern="^ID$", value=FALSE, perl=TRUE)!=0)))
  
  
  #Clean the data
  goodRows=which(fullDataSet[,substrateCol]!="#REF" & fullDataSet[,substrateCol]!="ART_025_CXGLRRWSLGGLRRWSL" & 
                   fullDataSet[,substrateCol]!="pVASP_150_164" & fullDataSet[,substrateCol]!="pTY3H_64_78" &
                   fullDataSet[,substrateCol]!="" & fullDataSet[,substrateCol]!="ID" & fullDataSet[,substrateCol]!="UniprotAccession")
  goodRowsTable=fullDataSet[goodRows,(substrateCol):ncol(fullDataSet)]
  goodRowsNames<-goodRowsTable[,1]
  
  return(list(cycleRow = cycleRow,
              exposureTimeRow = exposureTimeRow,
              sequenceRow = sequenceRow,
              goodRows = goodRows,
              goodRowsTable = goodRowsTable,
              goodRowsNames = goodRowsNames
              
              ))
  
}



make_graphs <- function(controlSamplesTest, experimentalSamplesTest, fullDataSet, opt1, opt2, lower, upper){
 
  #Get list of unique control samples (if the user selected more than one sample type)
  controlSamples_u=unique(as.character(controlSamplesTest[sampleNameRow, ]))
  experimentalSamples_u=unique(as.character(experimentalSamplesTest[sampleNameRow, ]))
  allSamples_u=union(controlSamples_u, experimentalSamples_u)
  
  par(mfrow=c(1,1))
  
  #for ggplot
  graphs=list()
  for(sub in goodRows){
    graphs[[sub-sequenceRow]]=ggplot(as.data.frame(matrix(NA, ncol=2, nrow=5)))
  }
  
  barcodeColumns5=as.data.frame(matrix(nrow=nrow(controlSamplesTest), ncol=0))
  
  #For each unique sample (control)
  results_table_c=NULL
  for(samp in 1:length(controlSamples_u)){
    
    #Pull out the correct samples and reduce this to the post wash phase
    barcodeColumns=sapply(controlSamples_u[samp], grep, controlSamplesTest)
    barcodeColumns2=controlSamplesTest[,barcodeColumns]
    postWashCol=min(which(sapply(1:(length(as.numeric(barcodeColumns2[cycleRow,]))-1), temp_fun, as.numeric(barcodeColumns2[cycleRow,]))==0))
    barcodeColumns3=barcodeColumns2[,postWashCol:ncol(barcodeColumns2)]
    barcodeColumns5=cbind(barcodeColumns5,barcodeColumns3)
    
    #Some empty vectors to store the results of the linear regression
    r2s=NULL
    slopes=NULL
    #For each substrate
    for(sub in goodRows){
      
      #Do and save linear regression results (maybe plot)

      temp_graph=graphs[[sub-sequenceRow]] + geom_line(aes(x=V1, y=V2), as.data.frame(cbind(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,])))) + labs( x= "Time", y= "Intensity")
      graphs[[sub-sequenceRow]]=temp_graph
      df=as.data.frame(cbind(as.numeric(barcodeColumns3[exposureTimeRow,]), as.numeric(barcodeColumns3[sub,])))
      linear=lm(V2~V1, df)
      #linear=lm(y~0+.,df)
      r2s=c(r2s,summary(linear)$r.squared)
      slopes=c(slopes, coef(linear)[["V1"]])
    }
    
    #If there is more than one sample then average the r2s and slopes
    results_table_c=cbind(results_table_c,goodRowsTable[,1], r2s, slopes)
    
  }
  kitty<<-cbind(results_table_c[,1], apply(results_table_c[,seq(3,ncol(results_table_c), 3)],2,as.numeric))
  for(x in 2:ncol(kitty)){
    colnames(kitty)[x]=paste0("control_", (x-1))
  }
  results_table_c=cbind(results_table_c[,1], apply(apply(results_table_c[,seq(2,ncol(results_table_c), 3)],2,as.numeric), 1, mean),apply(apply(results_table_c[,seq(3,ncol(results_table_c), 3)],2,as.numeric), 1, mean))
  
  #For each unique sample (experimental)
  results_table_e=NULL
  for(samp in 1:length(experimentalSamples_u)){
    
    #Pull out the correct samples and reduce this to the post wash phase
    barcodeColumns=sapply(experimentalSamples_u[samp], grep, experimentalSamplesTest)
    barcodeColumns2=experimentalSamplesTest[,barcodeColumns]
    postWashCol=min(which(sapply(1:(length(as.numeric(barcodeColumns2[cycleRow,]))-1), temp_fun, as.numeric(barcodeColumns2[cycleRow,]))==0))
    barcodeColumns4=barcodeColumns2[,postWashCol:ncol(barcodeColumns2)]
    barcodeColumns5=cbind(barcodeColumns5,barcodeColumns4)
    
    #Some empty vectors to store the results of the linear regression
    r2s=NULL
    slopes=NULL
    
    #For each substrate
    for(sub in goodRows){
      
      #Do and save linear regression results (maybe plot)
      
      temp_graph=graphs[[sub-sequenceRow]] + geom_line(color="red", aes(x=V1, y=V2), as.data.frame(cbind(as.numeric(barcodeColumns4[exposureTimeRow,]), as.numeric(barcodeColumns4[sub,])))) + labs( x= "Time", y= "Intensity")
      graphs[[sub-sequenceRow]]=temp_graph
      df=as.data.frame(cbind(as.numeric(barcodeColumns4[exposureTimeRow,]), as.numeric(barcodeColumns4[sub,])))
      linear=lm(V2~V1, df)
      r2s=c(r2s,summary(linear)$r.squared)
      slopes=c(slopes, coef(linear)[["V1"]])
    }
    
    #If there is more than one sample then average the r2s and slopes
    # if(samp == 1){
    results_table_e=cbind(results_table_e,goodRowsTable[,1], r2s, slopes)

  }
  kitty2<<-cbind(results_table_e[,1], apply(results_table_e[,seq(3,ncol(results_table_e), 3)],2,as.numeric))
  for(x in 2:ncol(kitty2)){
    colnames(kitty2)[x]=paste0("experimental_", (x-1))
  }
  
  results_table_e=cbind(results_table_e[,1], apply(apply(results_table_e[,seq(2,ncol(results_table_e), 3)],2,as.numeric), 1, mean),apply(apply(results_table_e[,seq(3,ncol(results_table_e), 3)],2,as.numeric), 1, mean))
  
  results_table=inner_join(as.data.frame(results_table_c), as.data.frame(results_table_e), by="V1")
  new_kitty=inner_join(as.data.frame(kitty), as.data.frame(kitty2), by="V1")
  new_kitty[,ncol(new_kitty)+1]=as.numeric(as.character(results_table[,3]))/as.numeric(as.character(results_table[,5]))
  results_table_download2=cbind(new_kitty, results_table[,2:5])
  colnames(results_table_download2)=c("Substrate", colnames(new_kitty)[2:(ncol(new_kitty)-1)], "Fold_Change", "R2_Control", "Average_Control", "R2_Experimental", "Average_Experimental")
  results_table_download2<<-results_table_download2
  
  print(barcodeColumns5)
  if(opt1==T){ #Remove substrates without a max value of above 2 for each sample
    print("in opt 1")
    temp=matrix(ncol=length(allSamples_u), nrow=length(goodRows))
    for(i in 1:(ncol(barcodeColumns5)/5)){
      #for(i in 1:length(allSamples_u)){ #this is where i am testing a fix
      print(i)
      temp[,i]=as.numeric(as.character(apply(barcodeColumns5[goodRows,((1:5)+(i-1)*5)],1,max)))
    }
    results_table=results_table[which(apply(temp, 1, min)>=2),]
    x<<-length(which(apply(temp, 1, min)<2))

  }
  if(opt2==T){ #Remove substrates without an average R2 of greater than 0.9
    print("in opt 2")

    y<<-length(which(as.numeric(as.character(results_table[,2]))<0.9 | as.numeric(as.character(results_table[,4]))<0.9))
    results_table=results_table[which(as.numeric(as.character(results_table[,2]))>=0.9 & as.numeric(as.character(results_table[,4]))>=0.9),]
    
  }

  results_table[,6]=as.numeric(as.character(results_table[,3]))/as.numeric(as.character(results_table[,5]))
  results_table_prime=new_kitty
  num_good<<-nrow(results_table)
  results_table_download1=results_table
  colnames(results_table_download1)=c("Substrate", "R2_Control", "Average_Control", "R2_Experimental", "Average_Experimental", "Fold_Change")
  results_table_download1<<-results_table_download1
  
  #TODO: save previous results table for heatmap
  results_table=results_table[which(results_table[,6]<=lower | results_table[,6]>=upper),]
  
  return(list(results_table, graphs, results_table_prime))
}



#This function is used in the post wash phase identification to find when the cycle is no longer increasing
temp_fun <- function(ind,list){
  dist(c(list[ind], list[ind+1]))
}

#To test if a value is non-zero
is.non.zero<-function(x){
  if(is.numeric(x)){
    return(x>0)
  }
}
# heatmap code

plotInputH <- function(mx, srt, end, low, mid, high, clust){
  #For the breaks
  breaks=seq(srt, #start point of color key
             end,  #end point of color key
             by=(end-srt)/50) #length of sub-division
  
  #heatmap colors
  mycol=colorpanel(n=length(breaks)-1,low=low,mid=mid,high=high)
  
  #Make heatmap
    heatmap.2(mx, #the matrix
            scale="row", 
            Colv=F, # No clustering of columns
            Rowv = F, #no clustering of rows
            hclustfun = function(x) hclust(x,method = clust), #clustering function
            col=mycol, #colors used in heatmap
            breaks=breaks, #color key details
            trace="none", #no trace on map
            na.rm=TRUE, #ignore missing values
            margins = c(15,10), # size and layout of heatmap window
            xlab = "Conditions", #x axis title
            ylab =  "Substrates",
            srtCol=90) # y axis title
}
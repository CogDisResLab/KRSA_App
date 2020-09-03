
library(tidyverse)

reading_data <- function(df) {
  
  meta_rows <- which(df[,1] == "")
  epmty_col <- which(df[meta_rows[1],] == "")
  
  meta_info <- df[meta_rows,(epmty_col[length(epmty_col)]+1):ncol(df)]
  
  rows_end <- nrow(df)
  rows_start <- meta_rows[(length(meta_rows))]+2
  
  data <- df[rows_start:rows_end,]
  data <- data[,c(1,(epmty_col[length(epmty_col)]+2):ncol(data))]
  data <- data %>% filter(! V1 %in% c("#REF", "ART_025_CXGLRRWSLGGLRRWSL", "pVASP_150_164", "pTY3H_64_78"))
  
  df2 <- data.frame(t(data[-1]), stringsAsFactors = F)
  colnames(df2) <- data[, 1]
  colnames(df2) <- make.unique(colnames(df2)) 
  
  df3 <- data.frame(t(meta_info[-1]), stringsAsFactors = F)
  colnames(df3) <- meta_info[, 1]
  cbind(df3,df2) -> df4
  
  df4 %>% gather((length(meta_rows)+1):ncol(.), key = "Peptide", value  = "Signal") -> tidydata
  
  tidydata$`Exposure time` <- as.numeric(tidydata$`Exposure time`)
  tidydata$Cycle <- as.numeric(tidydata$Cycle)
  tidydata$Signal <- as.numeric(tidydata$Signal)
  
  return(tidydata)
}







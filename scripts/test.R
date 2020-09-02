
fulldataset <- read.csv("data/datasets/DPLFC_MvsF_STK.txt", header=T, 
                        sep="\t", quote='"', stringsAsFactors = F)

source("KRSA2/R/main.R")

rows_list <- list_of_rows(fulldataset)

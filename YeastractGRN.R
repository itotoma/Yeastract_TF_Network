library(ggraph)
library(igraph)
library(tidygraph)
if(!exists("GRNDATADir.PATH")){cat("PLEASE SET THE PATH to GRN DATA in GRNDATADir.PATH\n")}

YRGRNBoth <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_both.tsv"), sep=";", header=F)
YRGRNAct <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_activator.tsv"), sep=";", header=F)
YRGRNInh <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_inhibitor.tsv"), sep=";", header=F) 

NameSerialDF <- c(YRGRNBoth$V1, YRGRNBoth$V2) %>% unique() %>% as.data.frame() %>% set_names("Name")
NameSerialDF <- NameSerialDF %>% mutate(SerialNum = c(1:nrow(NameSerialDF)))

NameTable_To_SerialTable <- function(TwoColumnNameTable){
  TwoColumnNameTable %>% left_join(NameSerialDF, by = c("V1" = "Name")) %>% 
    left_join(NameSerialDF, by = c("V2" = "Name")) %>% select(SerialNum.x, SerialNum.y) %>% 
    return()
}

tbl_YRGRNBoth = YRGRNBoth %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)
tbl_YRGRNAct = YRGRNAct %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)
tbl_YRGRNInh = YRGRNInh %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)



YRGRNDegDF = degree(tbl_YRGRNBoth, mode= "all") %>% as.data.frame() %>% 
  rownames_to_column() %>% setNames(c("SerialNum", "degree"))
YRGRNDegDF$SerialNum = as.numeric(YRGRNDegDF$SerialNum)
YRGRNDegDF = left_join(NameSerialDF, YRGRNDegDF, by="SerialNum")

Deg1Genes = YRGRNDegDF %>% filter(degree == 1)




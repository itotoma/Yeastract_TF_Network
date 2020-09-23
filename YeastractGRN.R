library(ggraph)
library(igraph)
library(tidygraph)
if(!exists("GRNDATADir.PATH")){cat("PLEASE SET THE PATH to GRN DATA in GRNDATADir.PATH\n")}

yeast_gene_name = read.csv("/Users/itoutouma/Lab_Analysis/GRN/orftogene.txt", sep="\t", header=T)
yeast_gene_name %>% dim()

Convert_GeneTwoColum_To_ORFTwoColumn <- function(DF){
  DF %>% left_join(yeast_gene_name, by=c("V1" = "Name")) %>% 
    left_join(yeast_gene_name, by=c("V2" = "Name")) %>% 
    mutate(ORF.y = ifelse(is.na(ORF.y), V2, ORF.y)) %>% 
    mutate(ORF.x = ifelse(is.na(ORF.x), V1, ORF.x)) %>% select(ORF.x, ORF.y) %>% 
    setNames(c("From", "To")) %>% return()
}

YRGRNBoth <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_both.tsv"), sep=";", header=F) %>% 
  Convert_GeneTwoColum_To_ORFTwoColumn()

YRGRNAct <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_activator.tsv"), sep=";", header=F) %>% 
  Convert_GeneTwoColum_To_ORFTwoColumn()

YRGRNInh <- read.csv(file.path(GRNDATADir.PATH, "RegulationTwoColumnTable_inhibitor.tsv"), sep=";", header=F) %>% 
  Convert_GeneTwoColum_To_ORFTwoColumn()

NameSerialDF <- data.frame(ORF = unique(c(YRGRNBoth$From, YRGRNBoth$To)), 
                           SerialNum = 1:length(unique(c(YRGRNBoth$From, YRGRNBoth$To))))

NameTable_To_SerialTable <- function(TwoColumnNameTable){
  TwoColumnNameTable %>% left_join(NameSerialDF, by = c("From" = "ORF")) %>% 
    left_join(NameSerialDF, by = c("To" = "ORF")) %>% select(SerialNum.x, SerialNum.y) %>% 
    return()
}

tbl_YRGRNBoth = YRGRNBoth %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)
tbl_YRGRNAct = YRGRNAct %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)
tbl_YRGRNInh = YRGRNInh %>% NameTable_To_SerialTable() %>% as_tbl_graph(directed = TRUE)


#Find_Feedback

#YDR463W
FindFeedBackOf <- function(target){
  
  #Narrow down candidate genes
  RegulatedByTarget = YRGRNBothORF %>% filter(From == target)
  RegulateTarget = YRGRNBothORF %>% filter(is.element(From, RegulatedByTarget$To)) %>% filter(To == target)
  FromGenes = YRGRNBothORF %>% filter(From == target) %>% filter(is.element(To, RegulateTarget$From ))
  
  #Search the interaction type
  ActivatedByTarget = YRGRNActORF %>% filter(From == target) %>% filter(is.element(To, FromGenes$To))
  InhibitedByTarget = YRGRNInhORF %>% filter(From == target) %>% filter(is.element(To, FromGenes$To))
  ActivateTarget = YRGRNActORF %>% filter(To == target) %>% filter(is.element(From, RegulateTarget$From))
  InhibitTarget =YRGRNInhORF %>% filter(To == target) %>% filter(is.element(From, RegulateTarget$From))
  
  AA = intersect(ActivatedByTarget$To, ActivateTarget$From)
  AI = intersect(ActivatedByTarget$To, InhibitTarget$From)
  IA = intersect(InhibitedByTarget$To, ActivateTarget$From)
  II = intersect(InhibitedByTarget$To, InhibitTarget$From)
  
  print(list(FromGenes, RegulateTarget))
  
  list(AA, AI, IA, II) %>% setNames(c("A_Back:A", "A_Back:I", "I_Back:A", "I_Back:I")) %>% 
    return()
}


# %>% select(ORF.x, ORF.y)


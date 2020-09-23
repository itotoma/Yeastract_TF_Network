setwd("/Users/itoutouma/Lab_Analysis/")

library(readr)
library(tm)
library(tidyverse)
library(ggraph)
library(igraph)
library(tidygraph)
library(gridExtra)

validate <- function(DF){
  if(dim(DF)[1]!=6564){cat("There is a lack of information!\n")}else{
    cat("Successfully Processed\n")
  }
}

### Build TF network object

GRNDATADir.PATH <- "/Users/itoutouma/Lab_Analysis/GRN"
source("./GRN/YeastractGRN.R")

YRGRNDegDF = degree(tbl_YRGRNBoth, mode= "all") %>% cbind(degree(tbl_YRGRNBoth, mode= "in")) %>% 
  cbind(degree(tbl_YRGRNBoth, mode= "out")) %>% cbind(NameSerialDF) %>% 
  set_names(c("alldeg", "indeg", "outdeg", "ORF", "SerialNum"))
validate(YRGRNDegDF)

### Build scDD object

library(scDD)
library(SingleCellExperiment)

raw_count_data <- read.table("./elife-51254-code1/data/103118_SS_Data.tsv", header=T, sep="\t")
YPD_cond = raw_count_data %>% filter(Condition == "YPD") %>% group_by(Genotype_Group)
conds = YPD_cond %>% group_keys(Genotype_Group) %>% pull(1)
count_by_cond <- YPD_cond %>% group_split() %>% set_names(conds) %>% 
  lapply(function(x){return(select(x, -KANMX, -NATMX, -Genotype, -Genotype_Group, -Replicate, -Condition, -tenXBarcode))})

build_sceObj <- function(geno1, geno2){
  n1 = nrow(count_by_cond[[geno1]])
  n2 = nrow(count_by_cond[[geno2]])
  condition <- c(rep(1,n1), rep(2, n2)) %>% 
    as.matrix() %>% t() %>% as.data.frame() %>%  
    data.table::setnames(c( paste(geno1, c(1:n1), sep=""), paste(geno2, c(1:n2), sep=""))) %>% 
    t()
  wt_targ_counts = cbind(t(count_by_cond[[geno1]]), t(count_by_cond[[geno2]]))
  wt_targ_sce <- SingleCellExperiment(assays=list(counts=wt_targ_counts), colData=data.frame(condition))
  return(wt_targ_sce)
}
detect_DD_from_normsceObj <- function(norm_sce){
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  scDatExSim <- scDD(norm_sce, prior_param=prior_param, testZeroes=TRUE)
  return(scDatExSim)
}
scDDpipeline <- function(genotype, zero_th){
  scDatEx.scran = build_sceObj("WT(ho)", genotype) %>% 
    preprocess(zero.thresh=zero_th, scran_norm=TRUE) %>% 
    detect_DD_from_normsceObj() %>% list() %>% return()
}

#This normalization makes the value scale of condition1 and condition2 same
scDatEx.scran = build_sceObj("WT(ho)","stp1") %>% 
  preprocess(zero.thresh=1, scran_norm=TRUE)

stp1DD <- detect_DD_from_normsceObj(scDatEx.scran) %>% results("Genes") %>% 
  select(gene, DDcategory)

## Start analysis from here!


#Most of genes that are regulated by STP1 don't shows difference expression 
#in the stp1 knock out group.
nth(filter(YRGRNBoth, From == "YDR463W"),2) %>% as.data.frame() %>% 
  setNames("Gene") %>% left_join(stp1DD, by=c("Gene"="gene")) %>% 
  filter(DDcategory == "NS")

#Genes that regulate STP1 don't show DE except YMR043W that feed back to STP1.
# DE: YMR043W,  DZ: STP1,  NS: Other(10)
nth(filter(YRGRNBoth, To == "YDR463W"),1) %>% as.data.frame() %>% 
  setNames("Gene") %>% left_join(stp1DD, by=c("Gene"="gene")) %>% 
  #filter(DDcategory == "DZ")
  select(DDcategory) %>% table()

YRGRNDegDD = YRGRNDegDF %>% left_join(stp1DD, by=c("ORF"="gene"))
validate(YRGRNDegDD)


#I wonder why the genes that are regulated by STP1 are classified as NS.
#Maybe the upstream genes compensate the lack of STP1 like below diagram.
#This hypothesis is supported if upstream gene ( Gene X) is classified as DE.

# /---| Gene X
# |      |
# |     \|/
#STP1 -> NS (RegBySTP1NS) ?

#STP1から抑制を受けるもの
InhBySTP1 <- YRGRNInh %>% filter(From == "YDR463W") %>% select(To)
RegbyXwithDD <- function(geneX, DD, Data){
  Data %>% filter(From == geneX) %>% 
    select(To) %>% left_join(stp1DD, by=c("To"="gene")) %>% 
    filter(DDcategory == DD) %>% return()
}

#STP1の下流遺伝子にも関わらずSTP1 KO で発現量が変わらなかった遺伝子
RegBySTP1NS <- RegbyXwithDD("YDR463W", "NS", YRGRNBoth)

#RegBySTP1NSの上流遺伝子
GenesX <- YRGRNBoth[is.element(YRGRNBoth$To, RegBySTP1NS$To),]$From %>% 
  intersect(InhBySTP1$To)

#Are they DE genes?  
stp1DD[is.element(stp1DD$gene, GenesX),]
#No, they are all NS genes. Above hypothesis is not suppoerted


#Genes which is regulated from both STP1 and STP2

#The half of genes that are regulated by STP1 are also reguated by STP2, which is STP1's paralog
RegBySTP1 = YRGRNBoth %>% filter(From == "YDR463W")
RegBySTP1And2 = YRGRNBoth[is.element(YRGRNBoth$To, RegBySTP1$To),] %>% 
  filter(From == "YHR006W")
nrow(RegBySTP1And2)/nrow(RegBySTP1) #0.4993271

RegBySTP1And2NS = YRGRNDegDD[is.element(YRGRNDegDD$ORF, RegBySTP1And2$To),] %>%
  filter(DDcategory == "NS")
RegBySTP1And2NS %>% nrow()
nrow(RegBySTP1And2NS)/nrow(RegBySTP1And2)
nrow(RegBySTP1NS)/nrow(RegBySTP1)

#Proportion of genes which show DE and have a interaction with STP1 and STP2
DEgenes = stp1DD %>% filter(DDcategory == "DE") %>% nth(1)
RegBySTP1 = RegBySTP1[is.element(RegBySTP1$To, DEgenes),]
RegBySTP1And2DE = RegBySTP1And2[is.element(RegBySTP1And2$To,DEgenes),]
nrow(RegBySTP1And2DE)/length(DEgenes)
# 0.03804348  

#But, I realize that most of DE genes aren't regulated by STP1
nrow(RegBySTP1)/length(DEgenes)
# 0.1413043

#This result is a natural because there must be a lot of unknown interaction.
#However, why the gene that have interaction with STP1 is classfied as NS?


#NS genes that regulated by STP1, named, RegBySTP1NS.
#RegBySTP1NS might has a DE upstream genes.
RBSNS_upstream = YRGRNBoth[is.element(YRGRNBoth$To, RegBySTP1NS$To),] %>% 
  left_join(stp1DD, by=c("From"="gene")) %>% group_by(To) %>% group_split()
names(RBSNS_upstream) = lapply(RBSNS_upstream, function(x){select(x, To) %>% unlist() %>% return()}) %>% 
  unlist() %>% unique()

Upstream_DE = sapply(RBSNS_upstream, function(data){ 
         data %>% filter(DDcategory == "DE") %>% nth(1) %>% return() })
#The hypothesis is not supported. Few of the upstream gene show DE.
#Interastingly, all DE upstream gene is YMR043W

#RegbySTP1NS may be actually suspicious interaction.
#Compare "YRGRNBoth" with more reliable data "YRGRNBothStrict"
reliableRegBySTP1NS = YRGRNBothStrict %>% filter(From == "YDR463W") %>% dplyr::nth(2)
stp1DD[is.element(stp1DD$gene, reliableRegBySTP1NS),]
#Even the genes that have reliable interaction with STP1, shows NS. 


#Maybe STP1 is not neccessary for YPD condition.
#The circuit may be off.
#if you set Zero.th 1, most of genes are dropped.
#So, in conclusion, most of interaction is not activated.



FindFeedBackOf("YDR463W") #STP1
#YMR043W #YEL009C #YML007W
#YHR006W

FindFeedBackOf("YHR006W") #STP2











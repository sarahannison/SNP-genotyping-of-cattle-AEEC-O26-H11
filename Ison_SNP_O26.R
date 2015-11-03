
#############Sarah's Data##################
getwd()
setwd("/Users/sarahison/Desktop/SNP")
search()
library(gdata)
library(plyr)
library(reshape2)
library(openxlsx)

SNP <- read.xlsx("O26 SNP Results.xlsx",4)
#Subset by Zygosity
SNP1 <- SNP[which(SNP$Zygosity == 'Homozygous'),] 

SNP1$SN <- "P"
SNP1$Region <- paste(SNP1$SN, SNP1$Region, sep="_")

#Reshape data to wide format
SNP2 <- dcast(SNP1, Id ~ Region, value.var = "Allele")

#Subset table Reference Genome
SNP3 <- SNP2[c("Id", "P_80", "P_264", "P_450", "P_620",
               "P_820",
               "P_940",
               "P_1144",
               "P_1319",
               "P_1507",
               "P_1720",
               "P_1891",
               "P_2061",
               "P_2224",
               "P_2391",
               "P_2549",
               "P_2781",
               "P_2900",
               "P_3147",
               "P_3270",
               "P_3445",
               "P_3645",
               "P_3832",
               "P_4032",
               "P_4251",
               "P_4436",
               "P_4653",
               #"P_4819", #Not Present
               "P_4970",
               "P_5183",
               "P_5332",
               "P_5514",
               "P_5648",
               "P_5764",
               "P_5967",
               "P_6145",
               "P_6333",
               "P_6515",
               "P_6725",
               "P_6852",
               "P_7060",
               "P_7221",
               "P_7379",
               "P_7600",
               "P_7785",
               "P_7952",
               "P_8092",
               "P_8281",
               "P_8435",
               "P_8591",
               "P_8755",
               "P_8897")]

SNP.5 <- read.xlsx("O26 SNP Results.xlsx",8)
SNP5 <-SNP.5

SNP5 <- transform(SNP5,concat1=paste0(P_80, P_264, P_450, P_620,P_820,P_940,sep = '')) 
SNP5 <- transform(SNP5,concat2=paste0(P_1144,P_1319,P_1507,P_1720,P_1891,P_2061,P_2224,sep = ''))
SNP5 <- transform(SNP5,concat3=paste0(P_2391,P_2549,P_2781,P_2900,P_3147,P_3270,P_3445,sep=''))
SNP5 <- transform(SNP5,concat4=paste0(P_3645,P_3832,P_4032,P_4251,P_4436,P_4653,P_4819,sep='')) 
SNP5 <- transform(SNP5,concat5=paste0(P_4970,P_5183,P_5332,P_5514,P_5648,P_5764,P_5967,sep = ''))
SNP5 <- transform(SNP5,concat6=paste0(P_6145,P_6333,P_6515,P_6725,P_6852,P_7060,P_7221,sep = ''))
SNP5 <- transform(SNP5,concat7=paste0(P_7379,P_7600,P_7785,P_7952,P_8092,P_8281,P_8435, sep = ''))
SNP5 <- transform(SNP5,concat=paste0(concat1,concat2,concat3,concat4,concat5,concat6,concat7,sep =''))

#Output Bletz Clonal Complex references
SNPCC <- read.xlsx("O26 SNP Results.xlsx",3)
SNPCC_B <-SNPCC

SNPCC_B <- transform(SNPCC_B,concat1=paste0(P_80, P_264, P_450, P_620,P_820,P_940,sep = '')) 
SNPCC_B <- transform(SNPCC_B,concat2=paste0(P_1144,P_1319,P_1507,P_1720,P_1891,P_2061,P_2224,sep = ''))
SNPCC_B <- transform(SNPCC_B,concat3=paste0(P_2391,P_2549,P_2781,P_2900,P_3147,P_3270,P_3445,sep=''))
SNPCC_B <- transform(SNPCC_B,concat4=paste0(P_3645,P_3832,P_4032,P_4251,P_4436,P_4653,P_4819,sep='')) 
SNPCC_B <- transform(SNPCC_B,concat5=paste0(P_4970,P_5183,P_5332,P_5514,P_5648,P_5764,P_5967,sep = ''))
SNPCC_B <- transform(SNPCC_B,concat6=paste0(P_6145,P_6333,P_6515,P_6725,P_6852,P_7060,P_7221,sep = ''))
SNPCC_B <- transform(SNPCC_B,concat7=paste0(P_7379,P_7600,P_7785,P_7952,P_8092,P_8281,P_8435,sep = ''))
SNPCC_B <- transform(SNPCC_B,concat=paste0(concat1,concat2,concat3,sep =''))

#Compare O26:H11 concatenated SNP output to BLETZ reference strains 
SNP5$match <- revalue(SNP5$concat, c("GGGCGGGTTGCGCTCCGATGCAAGGTGTCTGTGGCCGTCCCACAGCCG"="HUSE020", 
                                     "GTGCAGCTAGCGCTCCGATGCCAGGTGCATGCGGCCGTACCACAGCCG"="HUSE018",
                                     "GTGCAGGTTGCGCTCCGATGCAAGGTGTCTGTGGCCGTCCCACAGCCG"="HUSE019",
                                     "GTTTAGCTACCGCTCCGAAGGCAGGTGCATGCGGATATACCACAGCCG"="HUSE015", 
                                     "TTTTAGCTACCGCTCCGAAGGCAGGTGCATGCGGATATACCACAGCCG"="HUSE013",
                                     "TTTTATCGACAACCCTGGAAGCAAAAGCAAACGGATAGATACCGATGT"="CC4",
                                     "TTTTATCTACAAACTTAGAAGCCAAAGCAAACAAATAGATACTGATGT"="11368",
                                     "TTTTATCTACAACCCTGGAAGCAAAAGCAAACGGATAGACACCAATGT"="1226/65"))
#Output CSV files
write.csv(SNP1, file="SNP1.csv")
write.csv(SNP2, file="SNP2.csv")
write.csv(SNP3, file="SNP3.csv")

#Norman SNPs
SNP5 <- transform(SNP5,concatNorman=paste0(P_8591, P_8755, P_8897,sep = ''))

write.csv(SNP5, file="SNP5.csv")

SNP.5 <- read.xlsx("O26 SNP Results.xlsx",3)
SNPBletz <-SNP.5

SNPBletz <- transform(SNPBletz,concat1=paste0(P_80, P_264, P_450, P_620,P_820,P_940,sep = '')) 
SNPBletz <- transform(SNPBletz,concat2=paste0(P_1144,P_1319,P_1507,P_1720,P_1891,P_2061,P_2224,sep = ''))
SNPBletz <- transform(SNPBletz,concat3=paste0(P_2391,P_2549,P_2781,P_2900,P_3147,P_3270,P_3445,sep=''))
SNPBletz <- transform(SNPBletz,concat4=paste0(P_3645,P_3832,P_4032,P_4251,P_4436,P_4653,P_4819,sep='')) 
SNPBletz <- transform(SNPBletz,concat5=paste0(P_4970,P_5183,P_5332,P_5514,P_5648,P_5764,P_5967,sep = ''))
SNPBletz <- transform(SNPBletz,concat6=paste0(P_6145,P_6333,P_6515,P_6725,P_6852,P_7060,P_7221,sep = ''))
SNPBletz <- transform(SNPBletz,concat7=paste0(P_7379,P_7600,P_7785,P_7952,P_8092,P_8281,P_8435, sep = ''))
SNPBletz <- transform(SNPBletz,concat=paste0(concat1,concat2,concat3,concat4,concat5,concat6,concat7,sep =''))

write.csv(SNPBletz, file="SNPBletz.csv")

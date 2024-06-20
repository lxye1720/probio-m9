(A) the Simpson diversity index of the fecal microbiome 
#利用菌种数据计算其α多样性指数
library(vegan)
data_lxy <- read.table("clipboard", row.names=1, header=T,sep = "\t")
Shannon.Wiener <- diversity(data_lxy, indx = "shannon")
Simpson <- diversity(data_lxy, index = "simpson")
Inverse.Simpson <- diversity(data_lxy, index = "inv")
S <- specnumber(data_lxy)
plot(S)
J <- Shannon.Wiener/log(S)
div<-data.frame(Shannon.Wiener,Simpson,Inverse.Simpson,S,J)
write.csv(div,"Data.csv")



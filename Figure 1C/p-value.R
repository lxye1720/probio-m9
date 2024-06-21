#（C）Boxplots showing the intramammary levels of six cytokines, interferon (IFN)-γ, interleukin (IL)-4, IL-10, IL-1β, IL-2, and IL-6
#IFN_γ
library(ggpubr)
library(ggplot2)
setwd("path/to/your/directory")  #设置工作路径
df <-read.table("IFN_γ.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)

#IL_1β
df <-read.table("IL_1β.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)

#IL_2
df <-read.table("IL_2.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)

#IL_4
df <-read.table("IL_4.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)

#IL_6
df <-read.table("IL_6.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)

#IL_10
df <-read.table("IL_10.txt",header=T,check.names=FALSE,sep="\t")
p <- ggboxplot(df, x="group", y="len", color = "group", palette = "jco", add = "jitter",facet.by = "dose")     
p
compare_means(len~group, data=df, method = "wilcox.test", paired = F,  group.by = NULL)
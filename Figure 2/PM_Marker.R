#################################################################
# Function: Biomarker Test
# Call: Rscript PM_Marker -m map_file -i abund_file -o out_path -p prefix
# Rscript PM_Marker.R -i phylum.txt -m mapping_file.txt -P F -o phylum_marker -t 40
# R packages used: optparse ,utils
# Last update: 2016-03-29
#################################################################

## install necessary libraries
p <- c("optparse","utils")
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
    suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)

# make option list and parse command line
option_list <- list(  
    make_option(c("-i", "--abund_file"), type="character", help="Input feature table with Relative Abundance (*.Abd) [Required]."),
    make_option(c("-m", "--meta_data"), type="character", help="Input meta-data file [Required]."),
    make_option(c("-o", "--outfile"), type="character", default='Marker', help="Output path [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default='OUT',help="Output file prefix [Optional, default %default]"), #小写
    make_option(c("-P", "--Paired"), type="logical",default=FALSE,help="If paired samples [Optional, default %default]"), #大写
    make_option(c("-t", "--threshold"), type="double", default=0.01, help="Threshold of significance [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$meta_data)) stop('Please supply a mapping file.')
if(is.null(opts$abund_file)) stop('Please supply an abundance table.')

filename<-opts$abund_file
metadata.filename<-opts$meta_data

outpath<-opts$outfile
outpath<-paste(outpath,"/",sep="")

dir.create(outpath)

outprefix<-opts$prefix
if (outprefix != ""){
outperfix<-paste(outprefix,".",sep="")
}

MeanAB.Cutoff=0 # To remove variables whose mean relative abundance less than "MeanAB.Cutoff"
Zero.p=0.01          # The minimum percentage of zero value in each variable
alpha=opts$threshold         # The confidence level for each comparison

#--------------------------------------------------Fucntion
BetweenGroup.test <-function(data, group, p.adj.method="bonferroni",paired=opts$Paired){
    # p.adjust.methods
    # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    
    n_group<-nlevels(group)
    if(!is.numeric(n_group) | n_group==1)
        stop("group must be a numeric and up to two levels\n")
    if(n_group==2){
        output1<-matrix(NA,ncol=9,nrow=ncol(data))
        rownames(output1)<-colnames(data)
        colnames(output1)<-c(paste("mean_",levels(group)[1],sep=""),paste("mean_",levels(group)[2],sep=""),
                             paste("sd_",levels(group)[1],sep=""),paste("sd_",levels(group)[2],sep=""),"Var.test","T-test","Wilcoxon.test",
                             paste("T-test_",p.adj.method,sep=""),paste("Wilcoxon.test_",p.adj.method,sep=""))
        for(i in 1:ncol(data))
        {
            output1[i,1]<-mean(data[which(group==levels(group)[1]),i])
            output1[i,2]<-mean(data[which(group==levels(group)[2]),i])
            output1[i,3]<-sd(data[which(group==levels(group)[1]),i])
            output1[i,4]<-sd(data[which(group==levels(group)[2]),i])
            output1[i,5]<-var.test(data[,i]~group)$p.value
            if(output1[i,5]<0.01 | output1[i,5]=="NaN")
                output1[i,6]<-t.test(data[,i]~group,var.equal = FALSE,paired=paired)$p.value
            else
                output1[i,6]<-t.test(data[,i]~group, var.equal=TRUE,paired=paired)$p.value
            output1[i,7]<-wilcox.test(data[,i]~group, paired=paired, conf.int=TRUE, exact=FALSE, correct=FALSE)$p.value
            output1[i,8]<-NA
            output1[i,9]<-NA
        }
        output1[,8]<-p.adjust(output1[,6], method = p.adj.method, n = ncol(data))
        output1[,9]<-p.adjust(output1[,7], method = p.adj.method, n = ncol(data))
        
        return(data.frame(output1))
    }else{
        output2<-matrix(NA,ncol=n_group+5,nrow=ncol(data))
        rownames(output2)<-colnames(data)
        colnames.output2<-array(NA)
        for(j in 1:ncol(output2)){
            if(j<=n_group){
                colnames.output2[j]<-c(paste("mean_",levels(group)[j],sep=""))
            }else{
                colnames.output2[(n_group+1):(n_group+5)]<-c("Var.test","Oneway-test","Kruskal.test",
                                                             paste("Oneway-test_",p.adj.method,sep=""),paste("Kruskal.test_",p.adj.method,sep=""))
            }
        }
        colnames(output2)<-colnames.output2
        for(i in 1:ncol(data))
        {
            for(j in 1:n_group)
            {
                output2[i,j]<-mean(data[which(group==levels(group)[j]),i])
            }
            output2[i,(n_group+1)]<-bartlett.test(data[,i]~group)$p.value
            if(output2[i,(n_group+1)]<0.01)
            {output2[i,(n_group+2)]<-oneway.test(data[,i]~group)$p.value
            }else{
                output2[i,(n_group+2)]<-oneway.test(data[,i]~group, var.equal=T)$p.value}
            output2[i,(n_group+3)]<-kruskal.test(data[,i]~group)$p.value
            output2[i,(n_group+4)]<-NA
            output2[i,(n_group+5)]<-NA
        }
        output2[ ,(n_group+4)]<-p.adjust(output2[,(n_group+2)], method = p.adj.method, n = ncol(data))
        output2[ ,(n_group+5)]<-p.adjust(output2[,(n_group+3)], method = p.adj.method, n = ncol(data))
        return(data.frame(output2))
    }
    
    
}

#-------------------------------
con <- file(paste(outpath,outprefix,"GroupMeanComps.log",sep=""))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")
cat("
    -------------------------------
    Input  1. Data matrix: 
    row.names	Sample_id
    col.names	Varibles
    For example, data should be organized like this:
    Sample_id	group	V1	V2	etc...
    sample_0001	A	6	25
    sample_0002	B	9	32
    etc...        
    2. Metadata:
    row.names	Sample_id
    col.names	Sample_categories
    -------------------------------   
    Output 1. Table 1: Statistical results of all taxa
    Table 2: Statistical summary of significant taxa
    2. Logfile
    
    Last update: 2015-06-16
    
    
    ")
cat("
    Input Parameters: \n",
    "                  filename=",filename,"\n",
    "                  metadata.filename=",metadata.filename,"\n",
    "                  The cutoff of Mean AB of each taxon=",MeanAB.Cutoff,"\n",
    "                  The cutoff of percentage of Zero value for each taxon=",Zero.p, "\n",
    "                  The confidence level of statistical tests for taxa=", alpha, "\n",
    "\n"
)         

#g<-read.table(paste(filename,".txt",sep=""),header=T,row.names=1)
g<-read.table(filename,header=T,row.names=1,sep="\t")
g<-t(g)
#allmetadata<-read.table(paste(metadata.filename,".txt",sep=""),header=T,sep="\t",row.names=1)
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
#print(allmetadata)
cat(paste("    The number of all input variables : ", ncol(g) ,sep=""), "\n")
#-------------------------------filtering taxa with X% zero
Zero.p<-1
g<-g[,which(colSums(g==0)<Zero.p*nrow(g))][order(rownames(g)),] #删除80%样品中都是0的菌

cat(paste("    The number of variables (removed variables containing over ",Zero.p*100,"% zero): ", ncol(g) ,sep=""), "\n")
#-------------------------------
g<-g[,colMeans(data.matrix(g))>MeanAB.Cutoff]
cat(paste("    The number of variables (removed variables whose mean relative abundance less than ", MeanAB.Cutoff," :", ncol(g) ,sep=""), "\n")
gmat<-data.matrix(g)
#print(gmat)
#-------------------------------
#metadata<-allmetadata[sapply(allmetadata,class)=="factor"][rownames(gmat),]
metadata = allmetadata[rownames(gmat), ]
#print(metadata)
if(length(metadata)==nrow(gmat)){ 
    group<-metadata
    metadata.num<-1}else{
        metadata.num<-ncol(metadata)}

for(i in 1:metadata.num){
    if(metadata.num>1){
        group<-factor(metadata[,i],levels=rev(as.vector(unique(metadata[,i]))))
        group.name<-colnames(metadata)[i]}else{
            group<-factor(metadata,levels=rev(as.vector(unique(allmetadata[,1]))))
            #group<-factor(metadata,levels=rev(as.vector(unique(metadata[,1]))))
            group.name<-names(allmetadata)
        }
    #-------------------------------
    # BetweenGroup.test
    #-------------------------------
    if(nlevels(group)<3){
        #print(head(gamt), head(group))
        result<- BetweenGroup.test(gmat,group,p.adj.method="BH",paired=opts$Paired)
        result.sig<-result[which(result$Wilcoxon.test<alpha),]
        cat(paste("    The number of variables significantly different between ",levels(group)[1]," and ",levels(group)[2]," =",nrow(result.sig)," (alpha=" ,alpha,")",sep=""), "\n")
    }else{
        #-------------------------------
        #print(head(gmat))
        #print(metadata)
        result<- BetweenGroup.test(gmat,group,p.adj.method="BH",paired=opts$Paired)
        result.sig<-result[which(result$Kruskal.test<alpha),]
        cat(paste("    The number of variables significantly different among ",paste(levels(group),collapse=", ")," =",nrow(result.sig)," (alpha=" ,alpha,")",sep=""), "\n")
    }
    #-------------------------------
    sink(paste(outpath,outprefix,".",group.name,".all.meanTests.xls",sep="")); cat("\t"); write.table(result,quote=FALSE,sep="\t");sink()
    sink(paste(outpath,outprefix,".",group.name,".sig.meanTests.xls",sep="")); cat("\t"); write.table(result.sig,quote=FALSE,sep="\t");sink()
    #-------------------------------
    cat("    The number of differentially abundant taxa between ",group.name," is: ", nrow(result.sig)," (alpha=" ,alpha,")","\n",sep="")
    #-------------------------------
    if(nrow(result.sig)>1){
        gmat.sig<-gmat[,rownames(result.sig)]
        pdf(paste(outpath,outprefix,".", group.name,".sig.boxplot.pdf",sep=""),40,20)
        coul.group<-rainbow(length(levels(group)))
        par(mfrow = c(6,8))
        for(m in 1:ncol(gmat.sig)){
            NumMaxLabel<-max(nchar(as.character(group)))
            LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
            par(mar = c(4,LeftSpace,4,3))
            boxplot(gmat.sig[,m] ~ group, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                    main=colnames(gmat.sig)[m],cex.main=2,cex.axis=1.5,font.main=3,
                    xlab="Relative abundance (%)",cex.lab=1.5,
                    boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                    whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
            #title(ylab=group.name, line = LeftSpace-2)
            stripchart(gmat.sig[,m] ~ group,method="jitter",jitter=.2,add=T,pch=20) 
            
        }
        invisible(dev.off())
    }
    
    if(nrow(result.sig)==1){
        gmat.sig<-cbind(gmat[,rownames(result.sig)],1)
        colnames(gmat.sig)<-c(rownames(result.sig),"Temp")
        
        pdf(paste(outpath,outprefix,".", group.name,".sig.boxplot.pdf",sep=""),40,20)
        coul.group<-rainbow(length(levels(group)))
        par(mfrow = c(6,8))
        NumMaxLabel<-max(nchar(as.character(group)))
        LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
        par(mar = c(4,LeftSpace,4,3))
        boxplot(gmat.sig[,1] ~ group, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                main=colnames(gmat.sig)[1],cex.main=2,cex.axis=1.5,font.main=3,
                xlab="Relative abundance (%)",cex.lab=1.5,
                boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
        #title(ylab=group.name, line = LeftSpace-2)
        stripchart(gmat.sig[,1] ~ group,method="jitter",jitter=.2,add=T,pch=20) 
        
        invisible(dev.off())
    }
    
    if(nlevels(group)>2){
        if (metadata.num == 1){
            metadata<-cbind(allmetadata,"additional")
        }
        combination <- combn(na.omit(unique(group)),2)
        Num_combination <- ncol(combination)
        for(j in 1: Num_combination) {
            vec <- as.vector(combination[,j])
            gmat_combn <- rbind(gmat[row.names(metadata[which(metadata[,i]==vec[1]),]),],gmat[row.names(metadata[which(metadata[,i]==vec[2]),]),])
            gmat_combn<-gmat_combn[,which(colSums(gmat_combn==0)<Zero.p*nrow(gmat_combn))]  #删除80%样品中都是0的菌
            metadata_combn <- as.factor(as.vector(metadata[row.names(gmat_combn),i]))            
            try_test<- try(BetweenGroup.test(gmat_combn,metadata_combn,p.adj.method="BH",paired=opts$Paired),silent=TRUE) #判断程序是否有报错
            if('try-error' %in% class(try_test)){
            cat(paste(vec[1],vec[2],sep="_VS_"),":","unpaired","\n",sep=" ")
            result<-BetweenGroup.test(gmat_combn,metadata_combn,p.adj.method="BH",paired=F)
            }else{
            cat(paste(vec[1],vec[2],sep="_VS_"),":",opts$Paired,"\n",sep=" ")
            result<-BetweenGroup.test(gmat_combn,metadata_combn,p.adj.method="BH",paired=opts$Paired)
            }
            result.sig<-result[which(result$Wilcoxon.test<alpha),]
            sink(paste(outpath,outprefix,".",group.name,".",vec[1],".VS.",vec[2],".all.meanTests.xls",sep="")); cat("\t"); write.table(result,quote=FALSE,sep="\t");sink()
            sink(paste(outpath,outprefix,".",group.name,".",vec[1],".VS.",vec[2],".sig.meanTests.xls",sep="")); cat("\t"); write.table(result.sig,quote=FALSE,sep="\t");sink()
            #-------------------------------
            cat("    The number of differentially abundant taxa between ",paste(group.name,".",vec[1],".VS.",vec[2],sep="")," is: ", nrow(result.sig)," (alpha=" ,alpha,")", "\n",sep="")
            #-------------------------------
            if(nrow(result.sig)>1){
                gmat.sig<-gmat_combn[,rownames(result.sig)]
                pdf(paste(outpath,outprefix,".", group.name,".",vec[1],".VS.",vec[2],".sig.boxplot.pdf",sep=""),40,20)
                coul.group<-rainbow(length(levels(metadata_combn)))
                par(mfrow = c(6,8))
                for(m in 1:ncol(gmat.sig)){
                    NumMaxLabel<-max(nchar(as.character(metadata_combn)))
                    LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
                    par(mar = c(4,LeftSpace,4,3))
                    boxplot(gmat.sig[,m] ~ metadata_combn, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                            main=colnames(gmat.sig)[m],cex.main=2,cex.axis=1.5,font.main=3,
                            xlab="Relative abundance (%)",cex.lab=1.5,
                            boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                            whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
                    #title(ylab=group.name, line = LeftSpace-2)
                    stripchart(gmat.sig[,m] ~ metadata_combn,method="jitter",jitter=.2,add=T,pch=20) 
                }
                invisible(dev.off())
                
                gmat.sig<-gmat[,rownames(result.sig)]
                pdf(paste(outpath,outprefix,".", group.name,".",vec[1],".VS.",vec[2],".all_samples.sig.boxplot.pdf",sep=""),40,20)
                coul.group<-rainbow(length(levels(group)))
                par(mfrow = c(6,8))
                for(m in 1:ncol(gmat.sig)){
                    NumMaxLabel<-max(nchar(as.character(group)))
                    LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
                    par(mar = c(4,LeftSpace,4,3))
                    boxplot(gmat.sig[,m] ~ group, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                            main=colnames(gmat.sig)[m],cex.main=2,cex.axis=1.5,font.main=3,
                            xlab="Relative abundance (%)",cex.lab=1.5,
                            boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                            whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
                    #title(ylab=group.name, line = LeftSpace-2)
                    stripchart(gmat.sig[,m] ~ group,method="jitter",jitter=.2,add=T,pch=20) 
                    
                }
                invisible(dev.off())
                
            }
            
            if(nrow(result.sig)==1){
                gmat.sig<-cbind(gmat_combn[,rownames(result.sig)],1)
                colnames(gmat.sig)<-c(rownames(result.sig),"Temp")
                pdf(paste(outpath,outprefix,".", group.name,".",vec[1],".VS.",vec[2],".sig.boxplot.pdf",sep=""),40,20)
                coul.group<-rainbow(length(levels(metadata_combn)))
                par(mfrow = c(6,8))
                NumMaxLabel<-max(nchar(as.character(metadata_combn)))
                LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
                par(mar = c(4,LeftSpace,4,3))
                boxplot(gmat.sig[,1] ~ metadata_combn, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                        main=colnames(gmat.sig)[1],cex.main=2,cex.axis=1.5,font.main=3,
                        xlab="Relative abundance (%)",cex.lab=1.5,
                        boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                        whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
                #title(ylab=group.name, line = LeftSpace-2)
                stripchart(gmat.sig[,1] ~ metadata_combn,method="jitter",jitter=.2,add=T,pch=20) 
                invisible(dev.off())
                
                gmat.sig<-cbind(gmat[,rownames(result.sig)],1)
                colnames(gmat.sig)<-c(rownames(result.sig),"Temp")
                
                pdf(paste(outpath,outprefix,".", group.name,".",vec[1],".VS.",vec[2],".all_samples.sig.boxplot.pdf",sep=""),40,20)
                coul.group<-rainbow(length(levels(group)))
                par(mfrow = c(6,8))
                NumMaxLabel<-max(nchar(as.character(group)))
                LeftSpace<-ifelse(NumMaxLabel>4,NumMaxLabel,4)
                par(mar = c(4,LeftSpace,4,3))
                boxplot(gmat.sig[,1] ~ group, data=gmat.sig,col=coul.group,horizontal = TRUE, las=1,
                        main=colnames(gmat.sig)[1],cex.main=2,cex.axis=1.5,font.main=3,
                        xlab="Relative abundance (%)",cex.lab=1.5,
                        boxcol=coul.group,medcol="white",medlwd=2,outcex=4,outpch=20,outcol=coul.group,outline=FALSE,
                        whisklwd=2,whisklty=5,staplelwd=2,staplecol=coul.group,whiskcol=coul.group) # ordered by median
                #title(ylab=group.name, line = LeftSpace-2)
                stripchart(gmat.sig[,1] ~ group,method="jitter",jitter=.2,add=T,pch=20) 
                invisible(dev.off())
            }
        }
    }
}

#-------------------------------
sink() 
sink(type="message")
#-------------------------------

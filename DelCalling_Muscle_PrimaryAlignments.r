#load packages 
library(GenomicAlignments)
library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(data.table)

#load mitoanno, make gRange for plotting
locs=read.csv("/Users/amyvandiver/Box/Nanopore/mitogenes.csv")
gr= GRanges(seqnames=rep("M",nrow(locs)), IRanges(start=locs$Start,end=locs$End),strand=rep('*',nrow(locs)),name=locs$Name,type=locs$Type,complex=locs$Complex)
gr$names1=gr$name
gr$names1[which(gr$type=="tRNA")]=NA
gr$names1[which(gr$type=="D-Loop")]=NA
gr$names1[which(gr$name=="ATP8")]=NA
gr$names1[which(gr$name=="ATP6")]="ATP8-ATP6"

seqlengths(gr)<-("M"=max(end(gr)))

setwd("/Users/amyvandiver/Box/Nanopore/Wanagat_Muscle/")

#function for plotting 
plot_data= function (data, sample) {
    age=anno$Age[which(anno$Sample==sample)]
    if(length(gr_Dels[which(gr_Dels$sample==sample)])>0){
        ggbio() + circle(gr,geom='rect',aes(fill=type),space.skip=0.001)+ 
        circle(gr_Dels[which(gr_Dels$sample==sample)],geom = "link",aes(color=col),linked.to ="to.gr",radius=35) + 
        scale_color_manual(values=cols) +
        labs(title=paste0(age," Years"))
        }
    else{
         ggbio() + circle(gr,geom='rect',aes(fill=type),space.skip=0.001)+ 
        labs(title=paste0(age," Years"))
    }
}

#read in anotation file
anno=read.csv("anno.csv")


#define mode function to get summary stats
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#loop through each bam and get summary stats
anno$chrMcount=rep(NA,nrow(anno))
anno$chrMmean=rep(NA,nrow(anno))
anno$chrMmedian=rep(NA,nrow(anno))
anno$chrMmode=rep(NA,nrow(anno))
anno$large=rep(NA,nrow(anno))
anno$strand=rep(NA,nrow(anno))

for (i in 1:nrow(anno)){
    
bam=paste0("bams/",anno$Run[i],".rot.bam")
reads1=scanBam(BamFile(bam))
meta=reads1[[1]]
widths=meta$qwidth[which(meta$flag%in%c(0,16))]
    anno$chrMcount[i]=length(widths)
    anno$chrMmean[i]=mean(widths)
    anno$chrMmedian[i]=median(widths)
    anno$chrMmode[i]=getmode(widths)
    anno$large[i]=length(which(widths>15000))/length(widths)
    allstrands[i]=summary(factor(strand(datM),levels=c("+","-")))
    anno$strand[i]=as.numeric(allstrands[1]/(allstrands[1]+allstrands[2])*100)
}


write.csv(anno,"Anno_sum.csv")

anno=read.csv("Anno_sum.csv")

#stats for ChrM counts, enrichment
anno$enrich=anno$chrMcount/anno$reads
summary(anno$reads)
summary(anno$enrich)


#write.csv(anno,"Anno_sum.csv")

anno=read.csv("Anno_sum.csv")

##Call deletions in primary reads by parsing Bam 

#make empty df to append dels to
dels=data.frame(sample=NA,read=NA,start=NA,len=NA,end=NA,strand=NA)

#loop through annotation, read in bams, get chrM reads
for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    dat=readGAlignments(bam,use.names=TRUE)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    
    st=rep(NA,length(datM))
    len=rep(NA,length(datM))
 
#for each chrM read, make cigar into table to be queried for dels >100 bp    
    for (j in 1:length(datM)){
        tab=data.frame(letter=unlist(explodeCigarOps(cigar(datM)[j])),
                   number=as.numeric(unlist(explodeCigarOpLengths(cigar(datM)[j]))))
#first look at reads with 1 del >100 bp, pull out location of del
        if (length(which(tab$letter=="D"&tab$number>100))==1){
        idx=which(tab$letter=="D"&tab$number>100)
        p=tab[1:idx-1,]
        start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
        len=tab$number[idx[1]]
        end=start+len
        del1=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end,strand=as.character(strand(datM)[j]))
        dels=rbind(dels,del1)    
    }
#now for reads with >1 del over 100 bp
        else{
        if (length(which(tab$letter=="D"&tab$number>100))>1){
            idx=which(tab$letter=="D"&tab$number>100)
        #first get start, length, end for first deletion IDed in this read    
                p=tab[1:idx[1]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[1]]
                end=start+len
                del1=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end,strand=as.character(strand(datM)[j]))

        #now look at subsequent deletions, see if they overlap    
            for (k in 2:length(idx)){
                p=tab[1:idx[k]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[k]]
                end1=start+len
                if(start-del1$end<301){
                    del1$len=end1-del1$start
                    del1$end=end1
                    dels=rbind(dels,del1)
                }
                else{
                del2=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end1,strand=as.character(strand(datM)[j]))
                dels=rbind(dels,del1,del2)   
               }     
        }
        }
 
}
}
}            
#remove blank row at the beginning
dels=dels[-1,]
write.csv(dels,paste0("deletions_100bp.csv"))

#look at strand distribution of deletions 
strands=summary(factor(dels$strand,levels=c("+","-")),)
strands[1]/(strands[1]+strands[2])*100


allstrands=summary(factor(strand(datM),levels=c("+","-")))
allstrands[1]/(allstrands[1]+allstrands[2])*100


Dels=read.csv("deletions_100bp.csv")

summary(Dels)

summary(factor(Dels$sample))
summary(factor(Dels$sample[which(Dels$len>2000)]))

#Repeat above deletion finding but on random subset of 30k reads per samples 

set.seed(1234)

dels=data.frame(sample=NA,read=NA,start=NA,len=NA,end=NA)

#difference is subsetting ChrM reads to 300000 randomly chosen 
for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    dat=readGAlignments(bam,use.names=TRUE)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    datM=datM[sample(1:length(datM),30000,replace=FALSE),]
    
    st=rep(NA,length(datM))
    len=rep(NA,length(datM))
    for (j in 1:length(datM)){
        tab=data.frame(letter=unlist(explodeCigarOps(cigar(datM)[j])),
                   number=as.numeric(unlist(explodeCigarOpLengths(cigar(datM)[j]))))
#need to separately process reads with just 1 deletion vs more
        if (length(which(tab$letter=="D"&tab$number>100))==1){
        idx=which(tab$letter=="D"&tab$number>100)
        p=tab[1:idx-1,]
        start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
        len=tab$number[idx[1]]
        end=start+len
        del1=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end)
        dels=rbind(dels,del1)    
    }
        else{
        if (length(which(tab$letter=="D"&tab$number>100))>1){
            idx=which(tab$letter=="D"&tab$number>100)
        #first get start, length, end for first deletion IDed in this read    
                p=tab[1:idx[1]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[1]]
                end=start+len
                del1=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end)

        #now look at subsequent deletions, see if they overlap    
            for (k in 2:length(idx)){
                p=tab[1:idx[k]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[k]]
                end1=start+len
                if(start-del1$end<301){
                    del1$len=end1-del1$start
                    del1$end=end1
                    dels=rbind(dels,del1)
                }
                else{
                del2=data.frame(sample=anno$Sample[i],read=names(datM)[j],start=start,len=len,end=end1)
                dels=rbind(dels,del1,del2)   
               }     
        }
        }
 
}
}
}            

dels=dels[-1,]
write.csv(dels,paste0("deletions_100bp_30k.csv"))

Dels=read.csv("deletions_100bp_30k.csv")
summary(factor(Dels$sample))
summary(factor(Dels$sample[which(Dels$len>2000)]))

# Prep to plot deletions from 30k subset
#first unrotate coordinates (cause aligned to rotated genome)
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#make gRanges 
gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

#subset to deletions >2kb
gr_Dels=gr_Dels[which(gr_Dels$size>2000)]

#identify deletions that are the common deletion based on start and length, flag these as different color for plotting

Dels$sample=factor(Dels$sample,levels=anno$Sample)
gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8460 & start(gr_Dels) <8480 & gr_Dels$size<5000)]="red"
gr_Dels[which(start(gr_Dels)>8460 & start(gr_Dels) <8480 & gr_Dels$size<5000)]$size
cols <- c("grey33","red")

#plot using plot function defined at top
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)

pdf("primaryread_dels_30k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)
dev.off()

##Now back to big set from all reads for further analysis...
Dels=read.csv("deletions_100bp.csv")

#unrotate coordinates for this one 
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#remove deletions within 300 bp of cut site, but there aren't any, so ignore
#Dels=Dels[-which(Dels$start>1247 & Dels$start<1847),]
write.csv(Dels,"prim_unrot.csv")

Dels=read.csv("prim_unrot.csv")

#count deletions that are likely common deletion
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"

Dels$len[which(Dels$cat=="common")]

summary(factor(Dels$cat))
summary(factor(Dels$cat[which(Dels$len>2000)]))

#make empty df to calculate correlation with deletions of each size
tab=data.frame(start=as.numeric(c("0","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")),end=as.numeric(c("1000","2000","3000",
                                        "4000","5000","6000","7000","8000","9000","10000","16000")))
tab$name=paste0(tab$start,"-",tab$end)
tab$R_age=NA
tab$R_ddpcr=NA


#loop through size cutoffs, to calculate correlation to age and ddPCR for each size
for (i in 1:nrow(tab)){
    count=rep(NA,nrow(anno))
    for (j in 1:nrow(anno)){
    count[j]=(length(which(Dels$sample==anno$Sample[j] & Dels$len>tab$end[i]))+1)/anno$chrMcount[j]
    }
    tab$R_age[i]=summary(lm((log10(count)~anno$Age)))$r.squared
    tab$R_ddPCR[i]=summary(lm((log10(count)~log10(anno$ddPCR))))$r.squared

}

#plot correlations

pdf("Dels_prim_sizeopt.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(tab$R_age~tab$end,pch=19,cex=1.5,col="darkmagenta",main="R2 to Age vs Del Size \n Primary Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
plot(tab$R_ddPCR~tab$end,pch=19,cex=1.5,col="darkgreen",main="R2 to ddPCR vs Del Size \n Primary Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
dev.off()

#looking at dels greater than 3kbp, calculate frequency for each sample
for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$sample==anno$Sample[i] & Dels$len>3000))+1
}

#plot correlations for 4kbp cutoff
pdf("Dels_prim_corr.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(log10(anno$Dels/anno$chrMcount)~anno$Age,pch=19,col="darkmagenta",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",
     main=paste0("Primary Dels >3kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~anno$Age))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~anno$Age),lty=2)

plot(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR),pch=19,col="darkgreen",cex=1.5,xlab="ddPCR Dels",ylab="Dels/Read",
     xlim=c(-5.5,-2.5),ylim=c(-5.5,-2.5),main=paste0("Primary Dels >3kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)),lty=2)
abline(0,1,col="red",lty=2)

dev.off()

#plot ddPCR to age corr for supplemental 
pdf("ddPCR_age.pdf",width=5,height=5)
plot(log10(anno$ddPCR)~anno$Age,pch=19,col="green",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/mtDNA",
     main=paste0("ddPCR Deletion Frequency vs Age \n R2=",round(summary(lm(log10(anno$ddPCR)~anno$Age))$r.squared,2)))
abline(lm(log10(anno$ddPCR)~anno$Age),lty=2)
dev.off()

sessionInfo()



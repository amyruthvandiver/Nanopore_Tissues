#load all required packages 

library(GenomicAlignments)
#library(ggbio)
#library(GenomicRanges)
#library(ggplot2)
library(gridExtra)
library(data.table)
library(rBLAST)
library("seqinr")

#load annotation for plotting 

locs=read.csv("/Users/amyvandiver/Box/Nanopore/mitogenes.csv")
gr= GRanges(seqnames=rep("M",nrow(locs)), IRanges(start=locs$Start,end=locs$End),strand=rep('*',nrow(locs)),name=locs$Name,type=locs$Type,complex=locs$Complex)
gr$names1=gr$name
gr$names1[which(gr$type=="tRNA")]=NA
gr$names1[which(gr$type=="D-Loop")]=NA
gr$names1[which(gr$name=="ATP8")]=NA
gr$names1[which(gr$name=="ATP6")]="ATP8-ATP6"

seqlengths(gr)<-("M"=max(end(gr)))

setwd("/Users/amyvandiver/Box/Nanopore/Wanagat_Muscle/")

#load previously generated annotation 
anno=read.csv("Anno_sum.csv")

#plot function for circle plots
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

#load blast database
bl <- blast(db="/Users/amyvandiver/Library/CloudStorage/Box-Box/Nanopore/Analysis/simulated/mydb")

bl

#create empty data frame to add deletion info to
all=data.frame(a1=NA,a2=NA,b1=NA,b2=NA,start=NA,end=NA,len=NA, qa1=NA,qa2=NA,qb1=NA,qb2=NA,
               sample=NA,width=NA,qstart=NA,qend=NA,sep=NA)

#loop through each sample
for (i in 1:nrow(anno)){
    reads=scanBam(BamFile(paste0("bams/",anno$Run[i],".rot.bam")),param=ScanBamParam(what=scanBamWhat(),flag = scanBamFlag(isSupplementaryAlignment = FALSE)))
    meta=reads[[1]]

    dat=meta$seq
    qwidth=meta$qwidth
    cig=meta$cigar
    sample=anno$Sample[i]
    
#pull out reads for which read length is very different from cigar length 
    seqs=dat[which(abs(qwidth-cigarWidthAlongReferenceSpace(cig))>500)]
    widths=qwidth[which(abs(qwidth-cigarWidthAlongReferenceSpace(cig))>500)]
#for these reads, re-align using rBLAST
    for (j in 1:length(seqs)){
        cl=predict(bl,seqs[j])
#BLAST alignmentsi n right orientation only        
        cl=cl[which(cl$S.start<cl$S.end),]
        qwidth=widths[j]
        name=names[j]
#for those with two blast alignments        
        if (nrow(cl)==2){
             hits=IRanges(start=cl$S.start,end=cl$S.end)
#determine if they overlap, if not move ahead            
            if (length(findOverlaps(hits[1],hits[2]))==0){
                breaks=data.frame(a1=NA,a2=NA,b1=NA,b2=NA,start=NA, end=NA,len=NA)
                breaks$a1=cl$S.start[1]
                breaks$a2=cl$S.end[1]
                breaks$qa1=cl$Q.start[1]
                breaks$qa2=cl$Q.end[1]
                breaks$b1=cl$S.start[2]
                breaks$b2=cl$S.end[2]
                breaks$qb1=cl$Q.start[2]
                breaks$qb2=cl$Q.end[2]
                breaks$sample=sample
                #identify deletion as location between alignments 
                breaks$start=breaks$a2                
                breaks$start[which(breaks$a2>breaks$b2)]=breaks$b2[which(breaks$a2>breaks$b2)]
                breaks$end=breaks$b1
                breaks$end[which(breaks$a2>breaks$b2)]=breaks$a1[which(breaks$a2>breaks$b2)]
                breaks$len=breaks$end-breaks$start
                breaks$width=qwidth
                #identify locations on query, distance between alignments on query
                breaks$qstart=breaks$qa2                
                breaks$qstart[which(breaks$a2>breaks$b2)]=breaks$qb2[which(breaks$a2>breaks$b2)]
                breaks$qend=breaks$qb1
                breaks$qend[which(breaks$a2>breaks$b2)]=breaks$qa1[which(breaks$a2>breaks$b2)]
                breaks$sep=breaks$qend-breaks$qstart
            all=rbind(all,breaks)
    }
}
}
}
all=all[-1,]

write.csv(all, all_sep_muscle.csv)

#write.csv(all, "all_sep_muscle.csv")
all=read.csv("all_sep_muscle.csv")
all=all[-which(all$width>17000),]

length(which(all$len>2000))
length(which(all$len>2000 & abs(all$sep)<300))

#remove putative deletions that are too far apart on query to be real
all=all[which(abs(all$sep)<300),]

all$sample=factor(all$sample,levels=anno$Sample)
Dels=all
Dels$len=Dels$end-Dels$start

#unrotate coordinates for downstream analysis 

Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#remove those with start very close to cut site
Dels=Dels[-which(Dels$start>1247 & Dels$start<1847),]

Dels=Dels[which(Dels$len>100),]



#save unrotated list for downstream 
write.csv(Dels,"all_blast_unrot.csv")

all=read.csv("all_blast_unrot.csv")

#next will work with subset of 30k reads, set. seed so this is consistent between runs 
set.seed(1234)

##Same deletion finding as above, but on subset of 30k reads to normalize for plotting purpose

all=data.frame(a1=NA,a2=NA,b1=NA,b2=NA,start=NA,end=NA,len=NA, qa1=NA,qa2=NA,qb1=NA,qb2=NA,
               sample=NA,width=NA,qstart=NA,qend=NA,sep=NA)

for (i in 1:nrow(anno)){
    reads=scanBam(BamFile(paste0("bams/",anno$Run[i],".rot.bam")),param=ScanBamParam(what=scanBamWhat(),flag = scanBamFlag(isSupplementaryAlignment = FALSE)))
    meta=reads[[1]]

    ind=sample(1:length(meta$seq),30000,replace=F)
    dat=meta$seq[ind]
    qwidth=meta$qwidth[ind]
    cig=meta$cigar[ind]
    sample=anno$Sample[i]

    seqs=dat[which(abs(qwidth-cigarWidthAlongReferenceSpace(cig))>500)]
    widths=qwidth[which(abs(qwidth-cigarWidthAlongReferenceSpace(cig))>500)]

    for (j in 1:length(seqs)){
        cl=predict(bl,seqs[j])
        cl=cl[which(cl$S.start<cl$S.end),]
        qwidth=widths[j]
        if (nrow(cl)==2){
             hits=IRanges(start=cl$S.start,end=cl$S.end)
            if (length(findOverlaps(hits[1],hits[2]))==0){
                breaks=data.frame(a1=NA,a2=NA,b1=NA,b2=NA,start=NA, end=NA,len=NA)
                breaks$a1=cl$S.start[1]
                breaks$a2=cl$S.end[1]
                breaks$qa1=cl$Q.start[1]
                breaks$qa2=cl$Q.end[1]
                breaks$b1=cl$S.start[2]
                breaks$b2=cl$S.end[2]
                breaks$qb1=cl$Q.start[2]
                breaks$qb2=cl$Q.end[2]
                breaks$sample=sample
                #identify deletion as location between alignments 
                breaks$start=breaks$a2                
                breaks$start[which(breaks$a2>breaks$b2)]=breaks$b2[which(breaks$a2>breaks$b2)]
                breaks$end=breaks$b1
                breaks$end[which(breaks$a2>breaks$b2)]=breaks$a1[which(breaks$a2>breaks$b2)]
                breaks$len=breaks$end-breaks$start
                breaks$width=qwidth
                #identify locations on query, distance between alignments on query
                breaks$qstart=breaks$qa2                
                breaks$qstart[which(breaks$a2>breaks$b2)]=breaks$qb2[which(breaks$a2>breaks$b2)]
                breaks$qend=breaks$qb1
                breaks$qend[which(breaks$a2>breaks$b2)]=breaks$qa1[which(breaks$a2>breaks$b2)]
                breaks$sep=breaks$qend-breaks$qstart
            all=rbind(all,breaks)
    }
}
}
}

all=all[-1,]     


write.csv(all,"all_blast_30k_sep.csv")


#all=all[-which(all$width>17000)]
#write.csv(all,"all_blast_30k_sep.csv")
all=read.csv("all_blast_30k_sep.csv")

#load packages for plotting 
library(ggbio)
library(GenomicRanges)
library(ggplot2)

Dels=read.csv("all_blast_30k_sep.csv")

#remove deletions which are too far apart on query to make sense 
Dels=Dels[(which(abs(Dels$sep)<300)),]

dim(Dels)



#Dels=read.csv("all_blast_30k.csv")

Dels$sample=factor(Dels$sample,levels=anno$Sample)

#unrotate coordinates of subset dels for plotting 

Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#remove those that are close to cut site 
Dels=Dels[-which(Dels$start>1247 & Dels$start<1847),]



library("ggbio")

#make GRanages of dels for plotting 
gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

#subset to size of interest for plots 
gr_Dels=gr_Dels[which(gr_Dels$size>2000)]

#flag common deletion for plotting 
gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size<5100 & gr_Dels$size>4800 )]="red"
gr_Dels[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size<5000 & gr_Dels$size>4000)]$size
cols <- c("grey33","red")

#plot using function above
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)



pdf("blast_dels_30k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)
dev.off()

#back to big deletions list for rest of analysis 
Dels=read.csv("all_blast_unrot.csv")


length(which(Dels$len>2000))


#how many deletions are common deletion?  based on start and length 
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4000 & Dels$len<5000)]="common"
Dels$cat[which(Dels$end>8400 & Dels$end <8500 & Dels$len>4000 & Dels$len<5000)]="common"

Dels$len[which(Dels$cat=="common")]

summary(factor(Dels$cat))
summary(factor(Dels$cat[which(Dels$len>2000)]))

#make empty df to calculate correlations with age and ddpcr for different sized deletions 
tab=data.frame(start=as.numeric(c("0","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")),end=as.numeric(c("1000","2000","3000",
                                        "4000","5000","6000","7000","8000","9000","10000","16000")))
tab$name=paste0(tab$start,"-",tab$end)
tab$R=NA
tab$R_ddPCR=NA


#calculate correlation to age and ddPCR for each size del
for (i in 1:nrow(tab)){
    count=rep(NA,nrow(anno))
    for (j in 1:nrow(anno)){
    count[j]=(length(which(Dels$sample==anno$Sample[j] & Dels$len>tab$end[i]))+1)/anno$chrMcount[j]
    }
    tab$R[i]=summary(lm((log10(count)~anno$Age)))$r.squared
    tab$R_ddPCR[i]=summary(lm((log10(count)~log10(anno$ddPCR))))$r.squared

}

#plot correlations 
pdf("Rsq_vs_Size_Blast.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(tab$R~tab$end,pch=19,cex=1.5,col="darkmagenta",main="R2 to Age vs Del Size \n Blast Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
plot(tab$R_ddPCR~tab$end,pch=19,cex=1.5,col="darkgreen",main="R2 to ddPCR vs Del Size \n Blast Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
dev.off()

#calculate frequency of reads with del of 5kbp for each sample 

for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$sample==anno$Sample[i] & Dels$len>2000))+1
}



pdf("Dels_blast_corr.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(log10(anno$Dels/anno$chrMcount)~anno$Age,pch=19,col="darkmagenta",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",
     main=paste0("Blast Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~anno$Age))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~anno$Age),lty=2)

plot(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR),pch=19,col="darkgreen",cex=1.5,xlab="Log 10 ddPCR Dels",ylab="Log 10 Dels/Read",
     main=paste0("Blast Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)),lty=2)
dev.off()

sessionInfo()



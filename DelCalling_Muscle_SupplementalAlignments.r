#load packages 
library(GenomicAlignments)
library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
library(data.table)

#load annotation file for plotting 
locs=read.csv("/Users/amyvandiver/Box/Nanopore/mitogenes.csv")
gr= GRanges(seqnames=rep("M",nrow(locs)), IRanges(start=locs$Start,end=locs$End),strand=rep('*',nrow(locs)),name=locs$Name,type=locs$Type,complex=locs$Complex)
gr$names1=gr$name
gr$names1[which(gr$type=="tRNA")]=NA
gr$names1[which(gr$type=="D-Loop")]=NA
gr$names1[which(gr$name=="ATP8")]=NA
gr$names1[which(gr$name=="ATP6")]="ATP8-ATP6"

seqlengths(gr)<-("M"=max(end(gr)))

#plot function
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

setwd("/Users/amyvandiver/Box/Nanopore/Wanagat_Muscle/")

#load annotation generated in Prim Del callign
anno=read.csv("Anno_sum.csv")

##Call deletions from reads with both primary and supplemental alignments 
#make empty df to append dels to
Dels1=data.frame('SuppName'=NA,'SuppPos'=NA,'SuppWidth'=NA,'SuppCig'=NA,'flag'=NA,'strand'=NA,'qwidth'=NA,'SuppEnd'=NA,'PrimPos'=NA,'PrimWidth'=NA,'PrimName'=NA,'PrimCig'=NA,'PrimEnd'=NA,'Sample'=NA)

#loop through reads
for (j in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[j],".rot.bam")
    reads1=scanBam(BamFile(bam))
    meta=reads1[[1]]
#index reads which have supplemental alingment on chrM, pull out position and cigar from supp alignment , use cigar to calculate end
    idx=which(meta$flag%in%c(2048,2064))
    Dels=data.frame(SuppName=meta$qname[idx],SuppPos=meta$pos[idx],
                SuppWidth=cigarWidthAlongReferenceSpace(meta$cigar[idx]),SuppCig=meta$cigar[idx],flag=meta$flag[idx],strand=meta$strand[idx],qwidth=meta$qwidth[idx])
    Dels$SuppEnd=Dels$SuppPos+Dels$SuppWidth

#empty columns for stuff we will calculate   
    Dels$PrimPos=NA
    Dels$PrimWidth=NA
    Dels$PrimName=NA
    Dels$PrimCig=NA
    Dels$qwidth=NA

    for (i in 1:nrow(Dels)){
#index to match primary alignment with the supplemental
        ind=which(meta$qname==Dels$SuppName[i]& meta$flag%in%c(0,16))
#pull out primary alignment name, position and cigar. use cigar to calculate end
        Dels$PrimName[i]=meta$qname[ind]
        Dels$PrimPos[i]=meta$pos[ind]
        Dels$PrimWidth[i]=cigarWidthAlongReferenceSpace(meta$cigar[ind])
        Dels$PrimCig[i]=meta$cigar[ind]
        Dels$PrimEnd[i]=Dels$PrimPos[i]+Dels$PrimWidth[i]
#pull out length of query 
        Dels$qwidth[i]=meta$qwidth[ind]
}
#double check that primary and supplemental were from same read
identical(Dels$SuppName,Dels$PrimName)

#get rid of weird ones with NA locations
#Dels=Dels[-which(is.na(Dels$DelStart)),]
    
#get rid of Dels which start/end at cut site, these are likely fake    
#Dels=Dels[-c(which(Dels$SuppEnd>16568 & Dels$PrimPos<50), which(Dels$PrimEnd>16568 & Dels$SuppPos<50)),]

#add which sample the deletion came from, bind to larger df 
Dels$Sample=anno$Sample[j]
Dels1=rbind(Dels1,Dels)
    
#get rid of empty row at beginning    
Dels1=Dels1[-1,]
}
   
Dels=Dels1

#get rid of extra long reads, these are likely concatemers 
Dels=Dels[-which(Dels$qwidth>17000),]

#want to pull out locations on query sequence to look at arrangement of alignments, make blank lists for this 
sstart=rep(NA,nrow(Dels))
ssend=rep(NA,nrow(Dels))
#get query locations supp alignments 
for (i in 1:nrow(Dels)){
if(unlist(explodeCigarOps(Dels$SuppCig[i]))[1]%in%c("H","S")){
    sstart[i]=unlist(explodeCigarOpLengths(Dels$SuppCig[i]))[1]
}
else{
    sstart[i]=1
}
    ssend[i]=sstart[i]+cigarWidthAlongQuerySpace(Dels$SuppCig[i],after.soft.clipping=TRUE)
    }

#get query locations for primary alignments
pstart=rep(NA,nrow(Dels))
pend=rep(NA,nrow(Dels))
for (i in 1:nrow(Dels)){
if(unlist(explodeCigarOps(Dels$PrimCig[i]))[1]%in%c("H","S")){
    pstart[i]=unlist(explodeCigarOpLengths(Dels$PrimCig[i]))[1]
}
else{
    pstart[i]=1
}
    pend[i]=pstart[i]+cigarWidthAlongQuerySpace(Dels$PrimCig[i],after.soft.clipping=TRUE)
    }

#now want to look at order of alignments on query and whether alignments overlap on query or on reference to filter out overlapping ones
overlap=rep(NA,nrow(Dels))
overlap1=rep(NA,nrow(Dels))
order=rep("Prim",nrow(Dels))

for (i in 1:nrow(Dels)){
hits1=IRanges(start=sstart[i],end=ssend[i])
hits2=IRanges(start=pstart[i],end=pend[i])
if (length(findOverlaps(hits1,hits2))==0){
    overlap[i]=0    
}
    else{
    overlap[i]=width(intersect(hits1,hits2))}
    
if (sstart[i]<pstart[i]){
    order[i]="Supp"
}    
}

#look at whether alignments overlap on reference 
for (i in 1:nrow(Dels)){
hits1=IRanges(start=Dels$SuppPos[i],end=Dels$SuppEnd[i])
hits2=IRanges(start=Dels$PrimPos[i],end=Dels$PrimEnd[i])
if (length(findOverlaps(hits1,hits2))==0){
    overlap1[i]=0    
}
    else{
    overlap1[i]=width(intersect(hits1,hits2))}
}    
    
#calculate distance between alignments, start and end is based on order determined above
dist=NA
dist[which(order=="Supp")]=pstart[which(order=="Supp")]-ssend[which(order=="Supp")]

dist[which(order=="Prim")]=sstart[which(order=="Prim")]-pend[which(order=="Prim")]

#make new data frame with just pertinent info for deletions 
new=data.frame(name=Dels$SuppName,strand=Dels$strand,qwidth=Dels$qwidth,SuppPos=Dels$SuppPos,SuppWidth=Dels$SuppWidth,SuppEnd=Dels$SuppEnd,PrimPos=Dels$PrimPos,
               PrimWidth=Dels$PrimWidth,PrimEnd=Dels$PrimEnd,sstart,ssend,pstart,pend,overlap,overlap1,order,dist,Sample=Dels$Sample)    

#filter out potential dels which: overlap too much, are too far apart on query
new=new[which(new$overlap<50 & new$overlap1<50 & abs(new$dist)<300),]
#filter out potential dels for which one alignment is shorter than 200 bp
new=new[which(new$SuppWidth>200 & new$PrimWidth>200),]
    
#now that we have filtered list, calculate deletion location as space between alignments based on order determined above
new$DelStart=new$PrimEnd
new$DelEnd=new$SuppPos
new$DelStart[which(new$order=="Supp")]=new$SuppEnd[which(new$order=="Supp")]
new$DelEnd[which(new$order=="Supp")]=new$PrimPos[which(new$order=="Supp")]
new$newwidth=new$DelEnd-new$DelStart

#for those which alignments split near circle end (so start is much greater than end), correct for circularity
new$newwidth[which((new$DelStart-new$DelEnd)>50)]=new$DelEnd[which((new$DelStart-new$DelEnd)>50)]-(new$DelStart[which((new$DelStart-new$DelEnd)>50)]-16569)
    
#save both original list and filtered list    
write.csv(Dels1,"Deletions_supplemental_base.csv")
write.csv(new,"Deletions_supplemental_filt.csv")

Dels=read.csv("Deletions_supplemental_filt.csv")
#unrotate coordinates for deletions for plotting and downstream
Dels$DelStart[which(Dels$DelStart > 15022)]=Dels$DelStart[which(Dels$DelStart > 15022)]-15022
Dels$DelStart[which(Dels$DelStart < 15022)]=Dels$DelStart[which(Dels$DelStart < 15022)]+1547

Dels$DelEnd[which(Dels$DelEnd > 15022)]=Dels$DelEnd[which(Dels$DelEnd > 15022)]-15022
Dels$DelEnd[which(Dels$DelEnd < 15022)]=Dels$DelEnd[which(Dels$DelEnd < 15022)]+1547

Dels$end=Dels$DelEnd

Dels$start=Dels$DelStart



#save unrotated dels for downstream 
write.csv(Dels,"supplemental_new_unrot.csv")

Dels=read.csv("supplemental_new_unrot.csv")
Dels$len=Dels$newwidth
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]

summary(factor(Dels$cat))
summary(factor(Dels$cat[which(Dels$len>2000)]))

#make empty df for calculating R2s with idfferent size cutoffs 
tab=data.frame(start=as.numeric(c("0","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")),end=as.numeric(c("1000","2000","3000",
                                        "4000","5000","6000","7000","8000","9000","10000","16000")))
tab$name=paste0(tab$start,"-",tab$end)
tab$R_age=NA
tab$R_ddpcr=NA

#loop through each size cutoff and calculate R2 for age and ddPCR 
for (i in 1:nrow(tab)){
    count=rep(NA,nrow(anno))
    for (j in 1:nrow(anno)){
    count[j]=(length(which(Dels$Sample==anno$Sample[j] & Dels$len>tab$end[i]))+1)/anno$chrMcount[j]
    }
    tab$R_age[i]=summary(lm((log10(count)~anno$Age)))$r.squared
    tab$R_ddPCR[i]=summary(lm((log10(count)~log10(anno$ddPCR))))$r.squared

}

#look at results
tab

#plot correlations 

pdf("Dels_supp_sizeopt.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(tab$R_age~tab$end,pch=19,cex=1.5,col="darkmagenta",main="R2 to Age vs Del Size \n Supp Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
plot(tab$R_ddPCR~tab$end,pch=19,cex=1.5,col="darkgreen",main="R2 to ddPCR vs Del Size \n Supp Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
dev.off()

#using best cutoff, calculate Del frequency for each sample, add to annotation
for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$Sample==anno$Sample[i] & Dels$len>2000))+1
}

#plot correlations for frequency of >2kbp dels 
#pdf("Dels_supp_corr.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(log10(anno$Dels/anno$chrMcount)~anno$Age,pch=19,col="darkmagenta",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",
     main=paste0("Supplemental Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~anno$Age))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~anno$Age),lty=2)

plot(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR),pch=19,col="darkgreen",cex=1.5,xlab="Log 10 ddPCR Dels",ylab="Log 10 Dels/Read",
     main=paste0("Supplemental Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)),lty=2)
abline(0,1,col="red",lty=2)
#dev.off()

#now combine deletions called in supplemental alignments with those called in primary alignments 
Dels=read.csv("supplemental_new_unrot.csv")
prim=read.csv("prim_unrot.csv")
new=data.frame(Sample=c(prim$sample,Dels$Sample),start=c(prim$start,Dels$start),len=c(prim$len,Dels$newwidth))
Dels=new

Dels=Dels[-which(Dels$len<100),]


dim(Dels)
dim(Dels[which(Dels$len>2000),])
summary(Dels)
summary(prim)


#flag the common deletion in all deletions file to count 
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]

summary(factor(Dels$cat))
summary(factor(Dels$cat[which(Dels$len>2000)]))

#what is the % common deletion
100/2390*100

#make empty df to calculate correlations for big deletions list for each deletion size cutoff
tab=data.frame(start=as.numeric(c("0","1000","2000","3000","4000","5000","6000","7000","8000","9000","10000")),end=as.numeric(c("1000","2000","3000",
                                        "4000","5000","6000","7000","8000","9000","10000","16000")))
tab$name=paste0(tab$start,"-",tab$end)
tab$R_age=NA
tab$R_ddpcr=NA

#loop through and calculate R2 for each cutoff
for (i in 1:nrow(tab)){
    count=rep(NA,nrow(anno))
    for (j in 1:nrow(anno)){
    count[j]=(length(which(Dels$Sample==anno$Sample[j] & Dels$len>tab$end[i]))+1)/anno$chrMcount[j]
    }
    tab$R_age[i]=summary(lm((log10(count)~anno$Age)))$r.squared
    tab$R_ddPCR[i]=summary(lm((log10(count)~log10(anno$ddPCR))))$r.squared

}

#Plot correlations 
pdf("Dels_combined_sizeopt.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(tab$R_age~tab$end,pch=19,cex=1.5,col="darkmagenta",main="R2 to Age vs Del Size \n Supp Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
plot(tab$R_ddPCR~tab$end,pch=19,cex=1.5,col="darkgreen",main="R2 to ddPCR vs Del Size \n Supp Dels",
     ylab="R Squared",xlab="Del Min Size (bp)")
dev.off()

#calculation deletion frequency for each sample from combined list, using 2kbp cutoff
for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$Sample==anno$Sample[i] & Dels$len>2000))
}

#pdf("Dels_combined_corr.pdf",width=5,height=10)
par(mfrow=c(2,1))
plot(log10(anno$Dels/anno$chrMcount)~anno$Age,pch=19,col="darkmagenta",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",
     main=paste0("Combined Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~anno$Age))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~anno$Age),lty=2)

plot(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR),pch=19,col="darkgreen",cex=1.5,xlab="Log 10 ddPCR Dels",ylab="Log 10 Dels/Read",
      xlim=c(-5.5,-2.0),ylim=c(-5.5,-2.0),main=paste0("Combined Dels >2kb \n R2=",round(summary(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)))$r.squared,2)))
abline(lm(log10(anno$Dels/anno$chrMcount)~log10(anno$ddPCR)),lty=2)
abline(0,1,col="red",lty=2)

#dev.off()



#what is the average deletion frequency for specific ages? how much does this change?  
fit=lm(log10(anno$Dels/anno$chrMcount)~anno$Age)
10^(fit$coef[1]+fit$coef[2]*25)
10^(fit$coef[1]+fit$coef[2]*75)
10^(fit$coef[1]+fit$coef[2]*75)/10^(fit$coef[1]+fit$coef[2]*25)

#save file with frequencies for further plotting 
write.csv(anno,"muscle_freqs.csv")

#flag which deletions involve major arc vs those totally within minor
Dels$loc="major"
Dels$end=Dels$start+Dels$len
Dels$loc[which(Dels$start>407 &Dels$start <5746 &Dels$end>407 &Dels$end<5746)]="minor"

#what are these numbers?
summary(factor(Dels$loc))

anno$DelsAll=NA
anno$DelsMinor=NA

#count number of deletions in minor vs major for all deletions vs just ones >2kbp
for (i in 1:nrow(anno)){
anno$DelsAll[i]=length(which(Dels$Sample==anno$Sample[i]))+1    
anno$DelsMinor[i]=length(which(Dels$Sample==anno$Sample[i]  & Dels$loc=="minor"))+1
anno$DelsMajor[i]=length(which(Dels$Sample==anno$Sample[i]  & Dels$loc=="major"))+1
anno$DelsMinorBig[i]=length(which(Dels$Sample==anno$Sample[i]  & Dels$len>2000 & Dels$loc=="minor"))+1
anno$DelsMajorBig[i]=length(which(Dels$Sample==anno$Sample[i]  & Dels$len>2000 & Dels$loc=="major"))+1
}

#how does deletion frequency from each location correlate with age?
anno$DelsMinorBig
anno$DelsMajorBig
summary(lm(log10(anno$DelsMajor/anno$chrMcount)~anno$Age))
summary(lm(log10(anno$DelsMinor/anno$chrMcount)~anno$Age))
summary(lm(log10(anno$DelsMajorBig/anno$chrMcount)~anno$Age))
summary(lm(log10(anno$DelsMinorBig/anno$chrMcount)~anno$Age))

#plot major vs minor arc deletion frequency for each size cutoff 
pdf("Dels_minor.pdf",width=10,height=10)
par(mfrow=c(2,2))
plot(log10(anno$DelsMinor/anno$chrMcount)~anno$Age,pch=19,col="darkblue",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",ylim=c(-3,-1.8))
abline(lm(log10(anno$DelsMinor/anno$chrMcount)~anno$Age),lty=2,col="darkblue")
points(log10(anno$DelsMajor/anno$chrMcount)~anno$Age,pch=19,cex=1.5,col="darkorange")
abline(lm(log10(anno$DelsMajor/anno$chrMcount)~anno$Age),lty=2,col="darkorange")
Rmaj=round(summary(lm(log10(anno$DelsMajor/anno$chrMcount)~anno$Age))$r.squared,2)
Rmin=round(summary(lm(log10(anno$DelsMinor/anno$chrMcount)~anno$Age))$r.squared,2)
legend("topright",legend=c(paste0('Minor Arc R2=',Rmin),paste0('Major Arc R2=',Rmaj)),
                          pch=c(19, 19), col=c('darkblue', 'darkorange'))
plot(log10(anno$DelsMinorBig/anno$chrMcount)~anno$Age,pch=19,col="darkblue",cex=1.5,xlab="Age in years",ylab="Log 10 Dels/Read",ylim=c(-6,-2))
abline(lm(log10(anno$DelsMinorBig/anno$chrMcount)~anno$Age),lty=2,col="darkblue")
points(log10(anno$DelsMajorBig/anno$chrMcount)~anno$Age,pch=19,cex=1.5,col="darkorange")
abline(lm(log10(anno$DelsMajorBig/anno$chrMcount)~anno$Age),lty=2,col="darkorange")
Rmaj=round(summary(lm(log10(anno$DelsMajorBig/anno$chrMcount)~anno$Age))$r.squared,2)
Rmin=round(summary(lm(log10(anno$DelsMinorBig/anno$chrMcount)~anno$Age))$r.squared,2)
legend("topleft",legend=c(paste0('Minor Arc R2=',Rmin),paste0('Major Arc R2=',Rmaj)),
                          pch=c(19, 19), col=c('darkblue', 'darkorange'))

dev.off()

#calculate mean deletion size for each sample
anno$mean=NA

for (i in 1:nrow(anno)){
dat=Dels[which(Dels$Sample==anno$Sample[i]),] 
anno$mean[i]=mean(dat$len)  
}

#add age annotation to deletion df for plotting 
Dels$Age=anno$Age[match(Dels$Sample,anno$Sample)]


#whats the relationship between mean size and age? 
summary(lm(anno$mean~anno$Age))


#plot age vs mean size, vs full del distribution
pdf("Summary_supp_delsbyage.pdf",width=10,height=5)
par(mfrow = c(1, 2))
plot(anno$mean~anno$Age, xlab="Age in years",ylab="Mean deletion size",cex=1.5,pch=19,col="darkorange",main=paste0("Mean size, all deletions >100 bp"))
boxplot(Dels$len~Dels$Age,outline=FALSE)
boxplot(Dels$len~Dels$Age)

dev.off()


#moving forward will work with subset of data for representative plotting. setting seed to 
set.seed(1234)

## find deletions in subset of 30 k reads for consistency in plotting. Same code as above, but first subsetting meta alignments from subset of 30000 read names 
##Call deletions from reads with both primary and supplemental alignments 
#make empty df to append dels to
Dels1=data.frame('SuppName'=NA,'SuppPos'=NA,'SuppWidth'=NA,'SuppCig'=NA,'flag'=NA,'strand'=NA,'qwidth'=NA,'SuppEnd'=NA,'PrimPos'=NA,'PrimWidth'=NA,'PrimName'=NA,'PrimCig'=NA,'PrimEnd'=NA,'Sample'=NA)

#loop through reads
for (j in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[j],".rot.bam")
    reads1=scanBam(BamFile(bam))
    meta=reads1[[1]]
    meta1=data.frame(qname=meta$qname,pos=meta$pos,cigar=meta$cigar,flag=meta$flag,qwidth=meta$qwidth,strand=meta$strand)
    names=unique(meta1$qname)
    names1=names[sample(1:length(names),30000,replace=FALSE)]
    meta=meta1[which(meta1$qname%in%names1),]
    
#index reads which have supplemental alingment on chrM, pull out position and cigar from supp alignment , use cigar to calculate end
    idx=which(meta$flag%in%c(2048,2064))
    Dels=data.frame(SuppName=meta$qname[idx],SuppPos=meta$pos[idx],
                SuppWidth=cigarWidthAlongReferenceSpace(meta$cigar[idx]),SuppCig=meta$cigar[idx],flag=meta$flag[idx],strand=meta$strand[idx],qwidth=meta$qwidth[idx])
    Dels$SuppEnd=Dels$SuppPos+Dels$SuppWidth

#empty columns for stuff we will calculate   
    Dels$PrimPos=NA
    Dels$PrimWidth=NA
    Dels$PrimName=NA
    Dels$PrimCig=NA
    Dels$qwidth=NA

    for (i in 1:nrow(Dels)){
#index to match primary alignment with the supplemental
        ind=which(meta$qname==Dels$SuppName[i]& meta$flag%in%c(0,16))
#pull out primary alignment name, position and cigar. use cigar to calculate end
        Dels$PrimName[i]=meta$qname[ind]
        Dels$PrimPos[i]=meta$pos[ind]
        Dels$PrimWidth[i]=cigarWidthAlongReferenceSpace(meta$cigar[ind])
        Dels$PrimCig[i]=meta$cigar[ind]
        Dels$PrimEnd[i]=Dels$PrimPos[i]+Dels$PrimWidth[i]
#pull out length of query 
        Dels$qwidth[i]=meta$qwidth[ind]
}
#double check that primary and supplemental were from same read
identical(Dels$SuppName,Dels$PrimName)

#get rid of weird ones with NA locations
#Dels=Dels[-which(is.na(Dels$DelStart)),]
    
#get rid of Dels which start/end at cut site, these are likely fake    
#Dels=Dels[-c(which(Dels$SuppEnd>16568 & Dels$PrimPos<50), which(Dels$PrimEnd>16568 & Dels$SuppPos<50)),]

#add which sample the deletion came from, bind to larger df 
Dels$Sample=anno$Sample[j]
Dels1=rbind(Dels1,Dels)
    
#get rid of empty row at beginning    
Dels1=Dels1[-1,]
}
   
Dels=Dels1

#get rid of extra long reads, these are likely concatemers 
Dels=Dels[-which(Dels$qwidth>17000),]

#want to pull out locations on query sequence to look at arrangement of alignments, make blank lists for this 
sstart=rep(NA,nrow(Dels))
ssend=rep(NA,nrow(Dels))
#get query locations supp alignments 
for (i in 1:nrow(Dels)){
if(unlist(explodeCigarOps(Dels$SuppCig[i]))[1]%in%c("H","S")){
    sstart[i]=unlist(explodeCigarOpLengths(Dels$SuppCig[i]))[1]
}
else{
    sstart[i]=1
}
    ssend[i]=sstart[i]+cigarWidthAlongQuerySpace(Dels$SuppCig[i],after.soft.clipping=TRUE)
    }

#get query locations for primary alignments
pstart=rep(NA,nrow(Dels))
pend=rep(NA,nrow(Dels))
for (i in 1:nrow(Dels)){
if(unlist(explodeCigarOps(Dels$PrimCig[i]))[1]%in%c("H","S")){
    pstart[i]=unlist(explodeCigarOpLengths(Dels$PrimCig[i]))[1]
}
else{
    pstart[i]=1
}
    pend[i]=pstart[i]+cigarWidthAlongQuerySpace(Dels$PrimCig[i],after.soft.clipping=TRUE)
    }

#now want to look at order of alignments on query and whether alignments overlap on query or on reference to filter out overlapping ones
overlap=rep(NA,nrow(Dels))
overlap1=rep(NA,nrow(Dels))
order=rep("Prim",nrow(Dels))

for (i in 1:nrow(Dels)){
hits1=IRanges(start=sstart[i],end=ssend[i])
hits2=IRanges(start=pstart[i],end=pend[i])
if (length(findOverlaps(hits1,hits2))==0){
    overlap[i]=0    
}
    else{
    overlap[i]=width(intersect(hits1,hits2))}
    
if (sstart[i]<pstart[i]){
    order[i]="Supp"
}    
}

#look at whether alignments overlap on reference 
for (i in 1:nrow(Dels)){
hits1=IRanges(start=Dels$SuppPos[i],end=Dels$SuppEnd[i])
hits2=IRanges(start=Dels$PrimPos[i],end=Dels$PrimEnd[i])
if (length(findOverlaps(hits1,hits2))==0){
    overlap1[i]=0    
}
    else{
    overlap1[i]=width(intersect(hits1,hits2))}
}    
    
#calculate distance between alignments, start and end is based on order determined above
dist=NA
dist[which(order=="Supp")]=pstart[which(order=="Supp")]-ssend[which(order=="Supp")]

dist[which(order=="Prim")]=sstart[which(order=="Prim")]-pend[which(order=="Prim")]

#make new data frame with just pertinent info for deletions 
new=data.frame(name=Dels$SuppName,strand=Dels$strand,qwidth=Dels$qwidth,SuppPos=Dels$SuppPos,SuppWidth=Dels$SuppWidth,SuppEnd=Dels$SuppEnd,PrimPos=Dels$PrimPos,
               PrimWidth=Dels$PrimWidth,PrimEnd=Dels$PrimEnd,sstart,ssend,pstart,pend,overlap,overlap1,order,dist,Sample=Dels$Sample)    

#filter out potential dels which: overlap too much, are too far apart on query
new=new[which(new$overlap<50 & new$overlap1<50 & abs(new$dist)<300),]
#filter out potential dels for which one alignment is shorter than 200 bp
new=new[which(new$SuppWidth>200 & new$PrimWidth>200),]
    
#now that we have filtered list, calculate deletion location as space between alignments based on order determined above
new$DelStart=new$PrimEnd
new$DelEnd=new$SuppPos
new$DelStart[which(new$order=="Supp")]=new$SuppEnd[which(new$order=="Supp")]
new$DelEnd[which(new$order=="Supp")]=new$PrimPos[which(new$order=="Supp")]
new$newwidth=new$DelEnd-new$DelStart

#for those which alignments split near circle end (so start is much greater than end), correct for circularity
new$newwidth[which((new$DelStart-new$DelEnd)>50)]=new$DelEnd[which((new$DelStart-new$DelEnd)>50)]-(new$DelStart[which((new$DelStart-new$DelEnd)>50)]-16569)
    
#save both original list and filtered list    
write.csv(Dels1,"Deletions_supplemental_base_30k.csv")
write.csv(new,"Deletions_supplemental_filt_30k.csv")

##First make all plots with just supplemental read deletions 

Dels=read.csv("Deletions_supplemental_filt_30k.csv")
Dels$len=Dels$newwidth
summary(factor(Dels$Sample))
summary(factor(Dels$Sample[which(Dels$len>2000)]))

Dels$start=Dels$DelStart
Dels$end=Dels$DelEnd

# Prep to plot deletions from 30k subset
#first unrotate coordinates (cause aligned to rotated genome)
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#make gRanges 
gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$Sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

gr_Dels


#identify deletions that are the common deletion based on start and length, flag these as different color for plotting

Dels$Sample=factor(Dels$Sample,levels=anno$Sample)
gr_Dels$col=rep("black",length(gr_Dels))
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
gr_Dels$col[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="red"
gr_Dels=gr_Dels[which(gr_Dels$size>2000)]

Dels$len[which(Dels$cat=="common")]
cols <- c("grey33","red")

#plot using plot function defined at top
myplots <- lapply(levels(factor(Dels$Sample)), plot_data, data = Dels)

#plot deletions for 30k subset
arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)

pdf("supp_dels_30k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)
dev.off()

#now combine with primary read dels from 30k subset for plot
prim=read.csv("deletions_100bp_30k.csv")
prim=prim[,-c(1,3)]
Dels=read.csv("Deletions_supplemental_filt_30k.csv")
Dels$start=Dels$DelStart
Dels$end=Dels$DelEnd
Dels$len=Dels$newwidth
Dels$sample=Dels$Sample
Dels1=Dels[,c("sample","start","len","end")]

dels=rbind(prim,Dels1)


dim(Dels)
dim(prim)
dim(dels)

Dels=dels
summary(factor(Dels$sample))
summary(factor(Dels$sample[which(Dels$len>2000)]))

# Prep to plot deletions from 30k subset with combined primary and supplemental
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

dim(Dels)

gr_Dels=gr_Dels[which(gr_Dels$size>2000)]

#identify deletions that are the common deletion based on start and length, flag these as different color for plotting

Dels$sample=factor(Dels$sample,levels=anno$Sample)
gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size>4800 & gr_Dels$size<5100)]="red"
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]
cols <- c("grey33","red")

#plot using plot function defined at top
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)


pdf("combined_dels_30k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)
dev.off()

arrangeGrobByParsingLegend(myplots,nrow = 3, ncol = 5, legend.idx = 1)


sessionInfo()





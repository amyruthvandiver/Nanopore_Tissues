#load packages 
library(GenomicAlignments)
library(GenomicRanges)
library(gridExtra)
library(data.table)
library(ggbio)

#make annotation file for plotting
locs=read.csv("/Users/amyvandiver/Box/Nanopore/mitogenes.csv")
gr= GRanges(seqnames=rep("M",nrow(locs)), IRanges(start=locs$Start,end=locs$End),strand=rep('*',nrow(locs)),name=locs$Name,type=locs$Type,complex=locs$Complex)
gr$names1=gr$name
gr$names1[which(gr$type=="tRNA")]=NA
gr$names1[which(gr$type=="D-Loop")]=NA
gr$names1[which(gr$name=="ATP8")]=NA
gr$names1[which(gr$name=="ATP6")]="ATP8-ATP6"

seqlengths(gr)<-("M"=max(end(gr)))

setwd("/Users/amyvandiver/Box/Nanopore/Wanagat_Muscle/")

#make plot function
plot_data= function (data, sample) {
    age=anno$Age[which(anno$Name==sample)]
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

#annotation for placenta file
anno=read.csv("anno_other.csv")
anno$age=0
anno$Reads=c(34376,124911)

#define mode function to get summary stats
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#get summary statistics for placenta samples 
anno$chrMcount=rep(NA,nrow(anno))
anno$chrMmean=rep(NA,nrow(anno))
anno$chrMmedian=rep(NA,nrow(anno))
anno$chrMmode=rep(NA,nrow(anno))
anno$large=rep(NA,nrow(anno))

for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    param=ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("qual","flag"))
    dat=readGAlignments(bam,use.names=TRUE,param=param)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    anno$chrMcount[i]=length(datM)
    anno$chrMmean[i]=mean(qwidth(datM))
    anno$chrMmedian[i]=median(qwidth(datM))
    anno$chrMmode[i]=getmode(qwidth(datM))
    anno$large[i]=length(which(qwidth(datM)>15000))/length(datM)
}


anno$Sample=c("placenta1","placenta2")
write.csv(anno,"Anno_sum_placenta.csv")

anno


anno$chrMcount/anno$Reads

###Identify deletions in primary alignments 
#make empty data frame to add deletions to
dels=data.frame(sample=NA,read=NA,start=NA,len=NA,end=NA)

#loop through samples 
for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    dat=readGAlignments(bam,use.names=TRUE)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    
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
        del1=data.frame(sample=anno$Name[i],read=names(datM)[j],start=start,len=len,end=end)
        dels=rbind(dels,del1)    

    }else{
        if (length(which(tab$letter=="D"&tab$number>100))>1){
            idx=which(tab$letter=="D"&tab$number>100)
        #first get start, length, end for first deletion IDed in this read    
                p=tab[1:idx[1]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[1]]
                end=start+len
                del1=data.frame(sample=anno$Name[i],read=names(datM)[j],start=start,len=len,end=end)
        #now look at subsequent deletion, see if it overlaps   
            for (k in 2:length(idx)){
                p=tab[1:idx[2]-1,]
                start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
                len=tab$number[idx[k]]
                end1=start+len
                if(start-del1$end<301){
                    del1$len=end1-del1$start
                    del1$end=end1
                    dels=rbind(dels,del1)
                }
                else{
                rbind(dels,del1)
                del1=data.frame(sample=anno$Name[i],read=names(datM)[j],start=start,len=len,end=end1)
                dels=rbind(dels,del1)   
               }
                }
         #now look at subsequent deletion, see if it overlaps   
            if (length(which(tab$letter=="D"&tab$number>100))>2){
            for (k in 3:length(idx)){
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
                del1=data.frame(sample=anno$Name[i],read=names(datM)[j],start=start,len=len,end=end1)
                dels=rbind(dels,del1)   
               }                     
        }
      }     
    }

}            
}
        }

dels=dels[-1,]
write.csv(dels,paste0("deletions_placenta_100bp.csv"))

#look at deletions identified 
dels=read.csv("deletions_placenta_100bp.csv")
dels[which(dels$len>2000),]

Dels=dels

#unrotate coordinates for plotting 
Dels$len=Dels$end-Dels$start
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#Dels$end=Dels$start+Dels$len

gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

#make gRanges for plotting 
gr_Dels=gr_Dels[which(gr_Dels$size>2000 & gr_Dels$size <14000)]

Dels$sample=factor(Dels$sample,levels=anno$Sample)

gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size<5000 & gr_Dels$size>4000 )]="red"
gr_Dels[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size<5000 & gr_Dels$size>4000)]$size
cols <- c("grey33","red")


myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)

pdf("Placenta_prim_dels.pdf",width=12,height=5)
arrangeGrobByParsingLegend(myplots,nrow = 1, ncol = 2, legend.idx = 1)
dev.off()

arrangeGrobByParsingLegend(myplots,nrow = 1, ncol = 2, legend.idx = 1)


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
write.csv(Dels1,"Plac_Deletions_supplemental_base.csv")
write.csv(new,"Plac_Deletions_supplemental_filt.csv")

new=read.csv("Plac_Deletions_supplemental_filt.csv")


Dels=read.csv("deletions_placenta_100bp.csv")

#look at deletions identified 

dels=data.frame(start=c(Dels$start,new$DelStart),end=c(Dels$end,new$DelEnd),len=c(Dels$len,new$newwidth),Sample=c(Dels$sample,new$Sample))
dels[which(dels$len>2000),]

#unrotate coordinates for plotting 
Dels=dels
Dels$len=Dels$end-Dels$start
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#Dels$end=Dels$start+Dels$len

gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$Sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

#make gRanges for plotting 
gr_Dels=gr_Dels[which(gr_Dels$size>2000 & gr_Dels$size <14000)]

Dels$Sample=factor(Dels$Sample,levels=anno$Sample)

gr_Dels$col=rep("black",length(gr_Dels))
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]
cols <- c("grey33","red")


myplots <- lapply(levels(factor(Dels$Sample)), plot_data, data = Dels)

#plot
myplots


#make table of frequencies 
for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$Sample==anno$Sample[i] & Dels$len>2000))
}

anno

write.csv(anno,"placenta_freqs.csv")

sessionInfo()



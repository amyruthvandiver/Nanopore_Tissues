#load packages 
library(GenomicAlignments)
library(GenomicRanges)
library(gridExtra)
library(data.table)
library(ggbio)

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

#define plot function
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

#load annotation of brain samples 
anno=read.csv("anno_brain1.csv")
anno$Run=c("HB3_12092021","HB8_12092021","HB6_02022022","HB10_02022022","HB1_02082022","HB2_02082022")
anno$Reads=c(140548,235138,81411,161059,71466,53079)


#define mode function to get summary stats
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#loop through samples and get summary stats about chrM alignments 
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

anno


write.csv(anno,"Anno_sum_brain.csv")

anno=read.csv("Anno_sum_brain.csv")
anno

##Call Dels using "primary read" approach
#create empty df to append Deletions to
dels=data.frame(sample=NA,read=NA,start=NA,len=NA,end=NA)

#loop through samples, extract chrM alignments
for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    dat=readGAlignments(bam,use.names=TRUE)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    
    st=rep(NA,length(datM))
    len=rep(NA,length(datM))
#for each chrM alignment,convert cigar into table    
    for (j in 1:length(datM)){
        tab=data.frame(letter=unlist(explodeCigarOps(cigar(datM)[j])),
                   number=as.numeric(unlist(explodeCigarOpLengths(cigar(datM)[j]))))
#Identify deletions >100 bp, need to separately process reads with just 1 deletion vs more
#for those with 1 del, extract the location and length
        if (length(which(tab$letter=="D"&tab$number>100))==1){
        idx=which(tab$letter=="D"&tab$number>100)
        p=tab[1:idx-1,]
        start=start(datM)[j]+sum((p$number[which(p$letter%in%c("M","D"))]))
        len=tab$number[idx[1]]
        end=start+len
        del1=data.frame(sample=anno$Name[i],read=names(datM)[j],start=start,len=len,end=end)
        dels=rbind(dels,del1)    

    }else{
            
#for those with >1 del over 100 bp, look at each separately             
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
write.csv(dels,paste0("deletions_brain_100bp.csv"))

set.seed(1234)

##as above, call deletions using "primary" method, but in random subset of 12k reads 

dels=data.frame(sample=NA,read=NA,start=NA,len=NA,end=NA)

for (i in 1:nrow(anno)){
    bam=paste0("bams/",anno$Run[i],".rot.bam")
    dat=readGAlignments(bam,use.names=TRUE)
    datM=dat[which(seqnames(dat)=="chrM_rot")]
    datM=datM[sample(1:length(datM),12000,replace=FALSE),]
    
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
write.csv(dels,paste0("deletions_brain_12k_100bp.csv"))

#for subset deletions, unrotate coordinates to be on standard ref genomes
Dels=read.csv("deletions_brain_12k_100bp.csv")
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#create gRanges for plotting 
gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$len)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

#subset Granges to size
gr_Dels=gr_Dels[which(gr_Dels$size>2000)]

#flag deletions consistent with common deletion based on start and length 
gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8460 & start(gr_Dels) <8480 & gr_Dels$size<5000)]="red"
gr_Dels[which(start(gr_Dels)>8460 & start(gr_Dels) <8480 & gr_Dels$size<5000)]$size
cols <- c("grey33","red")

library(ggbio)

#order samples in annotation for plotting order
Dels$sample=factor(Dels$sample,levels=anno$Name[order(anno$Age)])

#plot using function above 
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)

pdf("brain_primaryread_dels_12k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 2, ncol = 3, legend.idx = 1)
dev.off()

write.csv(Dels,"deletions_100bp_12k_brain_unrot.csv")

anno$Sample=anno$Name

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
write.csv(Dels1,"Brain_deletions_supplemental_base.csv")
write.csv(new,"Brain_deletions_supplemental_filt.csv")

#unrotate coordinates for deletions for plotting and downstream
Dels$DelStart[which(Dels$DelStart > 15022)]=Dels$DelStart[which(Dels$DelStart > 15022)]-15022
Dels$DelStart[which(Dels$DelStart < 15022)]=Dels$DelStart[which(Dels$DelStart < 15022)]+1547

Dels$DelEnd[which(Dels$DelEnd > 15022)]=Dels$DelEnd[which(Dels$DelEnd > 15022)]-15022
Dels$DelEnd[which(Dels$DelEnd < 15022)]=Dels$DelEnd[which(Dels$DelEnd < 15022)]+1547

Dels$end=Dels$DelEnd

Dels$start=Dels$DelStart


set.seed(1234)

## find deletions in subset of 12k reads for consistency in plotting. Same code as above, but first subsetting meta alignments from subset of 30000 read names 
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
    names1=names[sample(1:length(names),12000,replace=FALSE)]
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
write.csv(Dels1,"Brain_Deletions_supplemental_base_12k.csv")
write.csv(new,"Brain_Deletions_supplemental_filt_12k.csv")

Dels=read.csv("Brain_Deletions_supplemental_filt_12k.csv")
Dels$len=Dels$newwidth

Dels$start=Dels$DelStart
Dels$end=Dels$DelEnd

Dels$sample=Dels$Sample

# Prep to plot deletions from 30k subset
#first unrotate coordinates (cause aligned to rotated genome)
Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

#make gRanges 
gr_Dels= GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$start,end=Dels$start),size=Dels$DelWidth)
gr_Dels$to.gr=GRanges(seqnames=rep("M",nrow(Dels)), IRanges(start=Dels$end,end=Dels$end))
gr_Dels$sample=Dels$sample
seqlengths(gr_Dels)<-("M"=16569)
seqlengths(gr_Dels)
seqlevels(gr_Dels)
seqlengths(gr_Dels$to.gr)<-("M"=max(end(gr)))

gr_Dels

gr_Dels=gr_Dels[which(Dels$len>2000)]
#identify deletions that are the common deletion based on start and length, flag these as different color for plotting

Dels$sample=factor(Dels$sample,levels=c("HB8", "HB6","HB10","HB1","HB2","HB3"))
gr_Dels$col=rep("black",length(gr_Dels))
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]
cols <- c("grey33","red")

#plot using plot function defined at top
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)

arrangeGrobByParsingLegend(myplots,nrow = 2, ncol = 3, legend.idx = 1)

pdf("supp_dels_brain_12k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 2, ncol = 3, legend.idx = 1)
dev.off()

#now combine with primary read dels from 30k subset for plot
prim=read.csv("deletions_100bp_12k_brain.csv")
prim=prim[,c("sample","start","len","end")]
Dels=read.csv("Brain_Deletions_supplemental_filt_12k.csv")
Dels$start=Dels$DelStart
Dels$len=Dels$newwidth
Dels$end=Dels$DelEnd
Dels$sample=Dels$Sample
Dels1=Dels[,c("sample","start","len","end")]

dels=rbind(prim,Dels1)


Dels=dels
summary(factor(Dels$sample))

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

gr_Dels=gr_Dels[which(Dels$len>2000)]
#identify deletions that are the common deletion based on start and length, flag these as different color for plotting

Dels$sample=factor(Dels$sample,levels=c("HB8", "HB6","HB10","HB1","HB2","HB3"))
gr_Dels$col=rep("black",length(gr_Dels))
gr_Dels$col[which(start(gr_Dels)>8400 & start(gr_Dels) <8500 & gr_Dels$size>4800 & gr_Dels$size<5100)]="red"
Dels$cat=NA
Dels$cat[which(Dels$start>8400 & Dels$start <8500 & Dels$len>4800 & Dels$len<5100)]="common"
Dels$len[which(Dels$cat=="common")]
cols <- c("grey33","red")

#plot using plot function defined at top
myplots <- lapply(levels(factor(Dels$sample)), plot_data, data = Dels)

arrangeGrobByParsingLegend(myplots,nrow = 2, ncol = 3, legend.idx = 1)

pdf("combined_dels_brain_12k.pdf",width=15,height=10)
arrangeGrobByParsingLegend(myplots,nrow = 2, ncol = 3, legend.idx = 1)
dev.off()

write.csv(Dels,"Brain_combined_dels_filt_unrot_12k.csv")



#make combined file for all
prim=read.csv("deletions_brain_100bp.csv")
supp=read.csv("Brain_deletions_supplemental_filt.csv")


prim=prim[,c("sample","start","len","end")]
supp$sample=supp$Sample
supp$len=supp$newwidth
supp$start=supp$DelStart
supp$end=supp$DelEnd
supp=supp[,c("sample","start","len","end")]

Dels=rbind(prim,supp)

#loop through samples to count number of dels >2kb per chrm read
anno$Dels=NA

for (i in 1:nrow(anno)){
anno$Dels[i]=length(which(Dels$sample==anno$Name[i] & Dels$len>2000))+1
}

anno$freq=anno$Dels/anno$chrMcount


#significance of difference of 2kb dels 
t.test(anno$freq~anno$Age,alternative="less")

summary(Dels[which(Dels$len>100),])

write.csv(anno,"Brain_anno_freqs.csv")

#unrotat coordinants for combined file 

Dels$start[which(Dels$start > 15022)]=Dels$start[which(Dels$start > 15022)]-15022
Dels$start[which(Dels$start < 15022)]=Dels$start[which(Dels$start < 15022)]+1547

Dels$end[which(Dels$end > 15022)]=Dels$end[which(Dels$end > 15022)]-15022
Dels$end[which(Dels$end < 15022)]=Dels$end[which(Dels$end < 15022)]+1547

Dels=Dels[which(Dels$len>2000),]

write.csv(Dels,"Brain_combined_Dels_unrot.csv")

sessionInfo()



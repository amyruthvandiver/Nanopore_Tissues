#load all packages 
library(GenomicAlignments)
library(GenomicRanges)
library(gridExtra)
library(data.table)

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

#load all Deletion lists, the ones where coordinates are unrotated 
prim=read.csv("prim_unrot.csv")
supp=read.csv("supplemental_new_unrot.csv")
blast=read.csv("all_blast_unrot.csv")
brain=read.csv("Brain_combined_Dels_unrot.csv")
combined=read.csv("Muscle_combined_unrot.csv")

combined$end=combined$start+combined$len

#filter out short dels 
primbreaks=c(prim$start[which(prim$len>2000)],prim$end[which(prim$len>2000)])
suppbreaks=c(supp$DelStart[which(supp$newwidth>2000)],supp$DelEnd[which(supp$newwidth>2000)])
blastbreaks=c(blast$start[which(blast$len>2000)],blast$end[which(blast$len>2000)])
brainbreaks=c(brain$start[which(brain$len>2000)],brain$end[which(brain$len>2000)])
combinedbreaks=c(combined$start[which(combined$len>2000)],combined$end[which(combined$len>2000)])

#load published deletion breakpoints
other=read.csv("dels_mitomap.csv")
other1=read.csv("mitobreak.csv")

#format those lists to get break points 
mapbreaks=as.numeric(c(unlist(tstrsplit(other$Junction,":")[1]),unlist(tstrsplit(other$Junction,":")[2])))
breakbreaks=c(other1[,2],other1[,3])

#make big table of break point location and type 
breaks=c(primbreaks,suppbreaks,blastbreaks,brainbreaks,combinedbreaks)
type=c(rep("prim",length(primbreaks)),rep("supp",length(suppbreaks)),
       rep("blast",length(blastbreaks)),rep("brain",length(brainbreaks)),rep("combined",length(combinedbreaks)))
       
breaks=c(breaks,mapbreaks,breakbreaks)

type=c(type,rep("map",length(mapbreaks)),rep("break",length(breakbreaks)))

locs$col=as.numeric(as.factor(locs$Type))+1
type=as.factor(type)

#density plots of each del type individually 
#pdf("Muscle_prim_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))
par(fig=c(0.1,0.9,0.45,0.9))
plot(density(breaks[which(type=="prim")],bw=100),main="Primary Read Deletions",col="dodgerblue3",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
par(fig=c(0.1,0.9,0.35,0.45), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.25,0.35), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
#dev.off()

#pdf("Muscle_supp_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))
par(fig=c(0.1,0.9,0.45,0.9))
plot(density(breaks[which(type=="supp")],bw=100),main="Supplemental Read Deletions",col="palegreen4",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
par(fig=c(0.1,0.9,0.35,0.45), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.25,0.35), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
#dev.off()

#pdf("Muscle_blast_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))
par(fig=c(0.1,0.9,0.45,0.9))
plot(density(breaks[which(type=="blast")],bw=100),main="Blast Read Deletions",col="darkorange3",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)

par(fig=c(0.1,0.9,0.35,0.45), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.25,0.35), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
#dev.off()


#pdf("Muscle_combined_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))
par(fig=c(0.1,0.9,0.45,0.9))
plot(density(breaks[which(type=="combined")],bw=100),main="Combined Deletions",col="darkred",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)

par(fig=c(0.1,0.9,0.35,0.45), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.25,0.35), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
#dev.off()

#plot prim, combined, blast and published dels together
pdf("Muscle_fig_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))

par(fig=c(0.1,0.9,0.75,0.9))

plot(density(breaks[which(type=="prim")],bw=100),main="Primary Read Deletions",col="dodgerblue3",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
rect(8370, -0.05, 8570, 0.01,col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, -0.05, 13550, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.6,0.75),new=TRUE)
plot(density(breaks[which(type=="combined")],bw=100),main="Combined Deletions",col="darkred",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
rect(8370, -0.05, 8570, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, -0.05, 13550, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.45,0.6),new=TRUE)
plot(density(breaks[which(type=="blast")],bw=100),main="Blast Deletions",col="darkorange3",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
rect(8370, -0.05, 8570, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, -0.05, 13550, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.35,0.45), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))
abline(v=1547,col="red",lty=2)
rect(8370, -0.05, 8570, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, -0.05, 13550, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.25,0.35), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))
abline(v=1547,col="red",lty=2)
rect(8370, -0.05, 8570, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, -0.05, 13550, 0.01, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
dev.off()

pdf("Brain_blast_densities.pdf",width=10,height=10)
par(cex=0.8, mai=c(0.05,0.05,0.1,0.1))
par(fig=c(0.1,0.9,0.645,0.9))
plot(density(breaks[which(type=="brain")],bw=100),main="Brain Deletions",col="darkorchid3",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
rect(8370, 0, 8570, 0.0002, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, 0, 13550, 0.0002, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(cex=0.8, mai=c(0.05,0.05,0.1,0.1), new=TRUE)
par(fig=c(0.1,0.9,0.39,0.645))
plot(density(breaks[which(type=="combined")],bw=100),main="Muscle Deletions ",col="darkred",xlim=c(0,16569),lwd=2.5)
abline(v=1547,col="red",lty=2)
rect(8370, 0, 8570, 0.0002, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)
rect(13350, 0, 13550, 0.0002, col="lightsalmon",border="lightsalmon",lty=3, xpd=FALSE)

par(fig=c(0.1,0.9,0.32,0.39), new=TRUE)
plot(density(breaks[which(type=="map")],bw=100),main="MitoMap Deletions",xlim=c(0,16569))
abline(v=1547,col="red",lty=2)

par(fig=c(0.1,0.9,0.25,0.32), new=TRUE)
plot(density(breaks[which(type=="break")],bw=100),main="MitoBreak Deletions",xlim=c(0,16569))
abline(v=1547,col="red",lty=2)

par(fig=c(0.1,0.9,0.2,0.25), new=TRUE)
plot(c(0,16569),c(0,50),type="n",ylab="",xlab="Location on ChrM",yaxt="n")
for (i in 1:(nrow(locs))){
rect(locs$Start[i],0,locs$End[i],50,col=locs$col[i])
}
abline(v=1547,col="red",lty=2)
dev.off()




library("sm")

set.seed(123)

#parametric testing using sm package to compare distributions, compare prim dels to others 

sm.density.compare(breaks[which(type%in%c("prim","map"))],type[which(type%in%c("prim","map"))],model="equal",h=50)
sm.density.compare(breaks[which(type%in%c("prim","break"))],type[which(type%in%c("prim","break"))],model="equal",h=50)
sm.density.compare(breaks[which(type%in%c("prim","combined"))],type[which(type%in%c("prim","combined"))],model="equal",h=50)


#parametric testing using sm package to compare distributions, compare supp dels to others 

sm.density.compare(breaks[which(type%in%c("combined","map"))],type[which(type%in%c("combined","map"))],model="equal",h=50)
sm.density.compare(breaks[which(type%in%c("combined","break"))],type[which(type%in%c("combined","break"))],model="equal",h=50)
sm.density.compare(breaks[which(type%in%c("combined","blast"))],type[which(type%in%c("combined","blast"))],model="equal",h=50)

sm.density.compare(breaks[which(type%in%c("combined","brain"))],type[which(type%in%c("combined","brain"))],model="equal",h=50)

head(brain)
head(combined)
combined$sample=combined$Sample

#make bed file of all deletion locations 
combined$chrom=rep("chrM",nrow(combined))
combined$type=rep("Muscle",nrow(combined))
bed=combined[,c('chrom','start','end','sample','type')]

brain$chrom=rep("chrM",nrow(brain))
brain$type=rep("Brain",nrow(brain))
bed2=brain[,c('chrom','start','end','sample','type')]


bedall=rbind(bed,bed2)

head(combined)

write.table(bedall,"Deletion_locs.bed",quote=F, sep="\t", row.names=F, col.names=F)

sessionInfo()



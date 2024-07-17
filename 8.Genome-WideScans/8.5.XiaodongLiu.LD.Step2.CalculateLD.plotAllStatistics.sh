library(scales)
source("~/Software/twisst/plot_twisst.R")



#Use larger oma for outer margin area
#Use mar for individual plots

## fst
#par(mai=c(0.2,0.4,0.2,0.4))
region="NC_030809.1.chr2"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
fst <- read.table("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/NC_030821.1.chr14.Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv",header=T)
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]

#peak<-c(fst$midbp[ which(fst$fst>0.61)][c(1,5)])
peak <- c(134700000, 135900000)
plot(NULL, xlim=c(start.pos,max(fst$midbp)), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,max(twisst_data_smooth$pos[[1]])), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,max(common$WinCenter)), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D

plot(NULL, xlim=c(start.pos,max(common$WinCenter)), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD

plot(NULL, xlim=c(start.pos,max(ld.output$midpos)), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, 140000000,by=20000000),labels=paste0(seq(start.pos/1000000, 140, by = 20), c(rep("",7)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

dev.off()

## chr16
region="NC_030823.1.chr16"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 80000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]

#peak<-c(fst$midbp[ which(fst$fst>0.60)][c(1,2)])
#peak2 <- c(fst$midbp[ which(fst$fst>0.60)][c(3,6)])
peak <- c(41400000, 41500000)
peak2<- c(43200000, 43300000)
peak3 <- c(49300000, 49600000)
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak3[1],0,peak3[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak3[1],0,peak3[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak3[1],0,peak3[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],-3,peak2[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak3[1],-3,peak3[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",4)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0.2,peak2[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak3[1],0.2,peak3[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr5
region="NC_030812.1.chr5"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 120000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]
#peak<-c(fst$midbp[ which(fst$fst>0.60)][c(1,2)])
#peak2 <- c(fst$midbp[ which(fst$fst>0.60)][c(3,6)])
peak <- c(55400000, 55500000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",6)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr6
region="NC_030813.1.chr6"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 120000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(400000, 1000000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",6)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr7
region="NC_030814.1.chr7"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 110000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(56200000, 56300000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",5)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr13
region="NC_030820.1.chr13"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 85000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(22400000, 22500000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",4)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr14
region="NC_030821.1.chr14"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 100000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(26500000, 26600000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=20000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 20), c(rep("",5)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr17
region="NC_030824.1.chr17"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 70000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(35800000, 36000000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=10000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 10), c(rep("",7)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr18
region="NC_030825.1.chr18"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 70000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(38300000, 39300000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=10000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 10), c(rep("",7)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr22
region="NC_030829.1.chr22"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 60000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(17000000, 17100000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)

box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=10000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 10), c(rep("",6)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()
dev.off()

## chr25
region="NC_030832.1.chr25"
region2 = gsub(".chr.+","", region,perl=T)


pdf(paste0(region,".pdf"))
#bitmap(paste0(region,".png"),h=4,w=3,res=300)
par(mfrow=c(5,1),mai = c(0.3, 0.6, 0.1, 0.3))
start.pos = 0
end.pos = 45000000
fst <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/",region,".Common.Defassa.Goatanc.fst.win100K.ABOVE10k.csv"),header=T)
fst$midbp[which.max(fst$fst)]
fst$midbp[which.max(fst$midbp)]

peak <- c(1100000, 1700000)
peak2 <- c(3200000, 3300000)

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.7), ylab=expression(F[ST]), xlab="",axes=F)
lines(fst$midbp,fst$fst)
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.7, by = 0.7))
abline(h=0.6122543,col="red",lty=2)
rect(peak[1],0,peak[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],0.7*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

## twisst
weights_file <- paste0("/home/users/xiaodong/Documents/Project/African1kg/waterbuck/","combined.",region2,"_100k.weights.mod.csv")

windows_data_file = paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/1.Fst_TWISST/","clean_regions_withNumSites_100k_",region2,".txt")

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=windows_data_file)
twisst_data_smooth <-smooth.twisst(twisst_data, span=0.02)

#plot.twisst(twisst_data_smooth,  ylim = c(0,2))
#par(mai=c(0.2,0.4,0.3,0.4))
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,1.5), ylab="TWISST", xlab="",axes=F)
colors=c(1:5)
sapply(1:5, function(i) lines(twisst_data_smooth$pos[[1]],twisst_data_smooth$weights[[1]][,i],col=colors[i]) )
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 2, by = 0.5))
rect(peak[1],0,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


## pi and tajima'sd
common <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)
defassa <- read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/2.PairwisePi_TajimaD/Defassa.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".thetas.win100k.thetasWindow.gz.pestPG"),header=T)

common$pi <- common$tP / common$nSites
defassa$pi <- defassa$tP / defassa$nSites

# plot pi

plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0,0.008), ylab="Pairwise Diversity", xlab="",axes=F)
lines(common$WinCenter, common$pi,col="#BEBADA")
lines(defassa$WinCenter, defassa$pi,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(0, 0.008, by = 0.002))
rect(peak[1],0,peak[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0,peak2[2],0.008*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()


# plot Tajima's D
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(-3,1.5), ylab="Tajima's D", xlab="",axes=F)
lines(common$WinCenter, common$Tajima,col="#BEBADA")
lines(defassa$WinCenter, defassa$Tajima,col="#FDB462")
axis(side=1,labels=FALSE)
axis(4, at = seq(-3, 1.5, by = 1.5))
rect(peak[1],-3,peak[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],-3,peak2[2],1.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()



## LD
library(data.table)
ld.1 <-read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Common/Common.Luangwa.Matetsi.Nairobi.Samburu.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.1$pos1 = as.numeric(gsub(".*:","",ld.1$V1))
ld.1$pos2 = as.numeric(gsub(".*:","",ld.1$V2))


window.size= 100000

breaks = seq(0, max(ld.1$pos1, ld.1$pos2), by = window.size)
ld.1$interval1 = findInterval(ld.1$pos1, breaks)
ld.1$interval2 = findInterval(ld.1$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.1$interval1,ld.1$interval2)) {

    index = which(ld.1$interval1 ==i & ld.1$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.1$V7[index],na.rm=T))
}

ld.output <- data.frame(midpos = positions, ld = ld.means)


ld.2 = read.table(paste0("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/3.ngsLD/Defassa_Down41inds/Defassa_Down41inds.Samole.Maswa.KVNP.QENP.Ugalla.Kafue.GoatRef.",region2,".0.01.41.ld"),header=F)

ld.2$pos1 = as.numeric(gsub(".*:","",ld.2$V1))
ld.2$pos2 = as.numeric(gsub(".*:","",ld.2$V2))

breaks = seq(0, max(ld.2$pos1, ld.2$pos2), by = window.size)
ld.2$interval1 = findInterval(ld.2$pos1, breaks)
ld.2$interval2 = findInterval(ld.2$pos2, breaks)

positions = c()
ld.means = c()
for (i in 1:max(ld.2$interval1,ld.2$interval2)) {

    index = which(ld.2$interval1 ==i & ld.2$interval2 == i)
    positions <- append(positions, (i-1)*window.size+window.size/2)
    ld.means <- append(ld.means, mean(ld.2$V7[index],na.rm=T))
}

ld.output2 <- data.frame(midpos = positions, ld = ld.means)

# plot LD                                                                                                                                                                                                                                            
plot(NULL, xlim=c(start.pos,end.pos), ylim=c(0.2,0.5), ylab=expression(r^2), xlab="",axes=F)
lines(ld.output$midpos, ld.output$ld,col="#BEBADA")
lines(ld.output2$midpos, ld.output2$ld,col="#FDB462")
axis(side=1, at = seq(start.pos, end.pos,by=10000000),labels=paste0(seq(start.pos/1000000, end.pos/1000000, by = 10), c(rep("",4)," Mb")))
axis(4, at = seq(0.2, 0.5, by = .15))
rect(peak[1],0.2,peak[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
rect(peak2[1],0.2,peak2[2],0.5*0.99, col = alpha('gray', 0.4),border=FALSE,lty=NULL)
box()

dev.off()

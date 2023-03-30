source("/home/xi/Desktop/scripts/Postdoc-scripts/Waterbuck/visFuns.R")
popall<-read.table("/home/xi/Desktop/NGSadmix.119inds/no_outgroup_crashed_duplicate_within_samples_Waterbuck.species.pops.txt",as.is=T, header=T)
outdir <- "/home/xi/Desktop/NGSadmix.119inds/"
waterbuckCols <- c ("Luangwa"="purple", "Matetsi"="red", "Nairobi"="dark green", "Samburu"="grey","Kafue"="dark blue","KVNP"="light blue","Maswa"="brown","QENP"="light green","Samole"="wheat", "Ugalla"="yellow")
popord <- c("Luangwa","Matetsi","Nairobi","Samburu",
            "Kafue","KVNP","Maswa","QENP","Samole","Ugalla")
pop <- popall$Pops
##### for each K
  par(mar=c(1,10,1,10) +0.1, xpd=F, mfrow=c(10,1), oma=c(3.0,0,0,0)) 
  waterbuckCols <- c ("Luangwa"="purple", "Matetsi"="red", "Nairobi"="dark green", "Samburu"="grey","Kafue"="dark blue","KVNP"="light blue","Maswa"="brown","QENP"="light green","Samole"="wheat", "Ugalla"="yellow")
  refpops2 <- c("Matetsi","KVNP") ###k=2
  fq <- list.files(path=paste0(indir, 2),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops2, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops2],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 2), cex = 2, xpd=NA)
  
  refpops3 <- c("Matetsi","Samburu","KVNP") ###k=3
  fq <- list.files(path=paste0(indir, 3),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops3, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops3],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 3), cex = 2, xpd=NA)
  
  refpops4 <- c("Matetsi","Samburu","KVNP", "Samole") ###k=4
  fq <- list.files(path=paste0(indir, 4),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops4, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops4],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 4), cex = 2, xpd=NA)
  
  refpops5 <- c("Matetsi","Samburu","KVNP", "Samole","Kafue") ###k=5
  fq <- list.files(path=paste0(indir, 5),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops5, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops5],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 5), cex = 2, xpd=NA)
  
  refpops6 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla") ###k=6
  fq <- list.files(path=paste0(indir, 6),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops6, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops6],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 6), cex = 2, xpd=NA)
  
  refpops7 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla","QENP") ###k=7
  fq <- list.files(path=paste0(indir, 7),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops7, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops7],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 7), cex = 2, xpd=NA)
  
  
  refpops8 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla","QENP","Nairobi") ###k=8
  fq <- list.files(path=paste0(indir, 8),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops8, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops8],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 8), cex = 2, xpd=NA)
  
  
  refpops9 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla","QENP","Nairobi","Maswa") ###k=9
  fq <- list.files(path=paste0(indir, 9),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops9, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops9],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=127, y=0.5, labels=paste("K =", 9), cex = 2, xpd=NA)
  
  refpops10 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla","QENP","Nairobi","Maswa","Luangwa") ###k=10
  fq <- list.files(path=paste0(indir, 10),recursive = TRUE, pattern="qopt_conv", full.names=T)
  q <- as.matrix(read.table(fq))
  ord <- orderInds(pop=pop, popord=popord)
  plotAdmix(q[, orderK(q, refpops=refpops10, pop=pop)], ord=ord,pop=pop, colorpal=waterbuckCols[refpops10],
            cex.lab=1.0, rotatelab=30,padj=0.02, main="")
  text(x=129, y=0.5, labels=paste("K =", 10), cex = 2, xpd=NA)
  
library(ggplot2)
library(grid)
library(gridExtra)
library(gridBase)
library(ggdendro)
library(lattice)
library(ape)
library('dendextend')
library(ggrepel)

setwd('/home/xi/Desktop/Fst_NJtree_PCAngsd_Heter/PCAngsd')

###  bar plot only based on 5 clusters from PC1 for both population and species, Chromosome 6 as example
info_all <- read.table("info.all", header = T)
ID <- read.table("108.ID.txt", header= T)
ID$id  <- 1:nrow(ID)
info_noOrder <- merge (ID, info_all , by='ID')
pop <- info_noOrder[order(info_noOrder$id), ] ####keep order with ID
mat <- as.matrix(read.table("chr7.Pcangsd.e9.cov", header = F)) ###matrix
e <- eigen(mat)

#pop
png(filename = "png_PC1/Chr6.pop.Fst.png", width = 800,height=800,)
head(d<-e$vectors[,1])
head(pop)
tab <-table(cut(d,5),pop$pop)
coll <- c("Samole" = "wheat","KVNP" = "light blue","QENP" = "light green","Maswa" = "brown","Ugalla" = "yellow3","Kafue" = "dark blue", "Samburu" = "grey","Nairobi" = "dark green","Matetsi" = "red","Luangwa" = "purple")
tab <-tab[,names(coll)]
df1 <- as.matrix(tab)
barplot(t(tab),col=coll,
        xlab = "Bins", 
        ylab = "Counts",
        cex.axis=1.5,
        cex.lab=1.5,
        main="Chr6:400000-1000000")
legend("topright", 
       legend = c("Samole","KVNP", "QENP","Maswa","Ugalla", "Kafue", "Samburu", "Nairobi","Matetsi", "Luangwa"), 
       fill = c("wheat","light blue", "light green","brown","yellow","dark blue", "grey", "dark green","red","purple"))
dev.off()
#species
png(filename = "png_PC1/Chr6.species.Fst.png", width = 800,height=800,)
tab <-table(cut(d,5),pop$species)
coll <- c("common" = "#BEBADA","defassa" = "#FDB462")
tab <-tab[,names(coll)]
barplot(t(tab),col=coll,
        xlab = "Bins", 
        ylab = "Counts",
        cex.axis=1.5,
        cex.lab=1.5,
        main="Chr6:400000-1000000")
legend("topright", 
       legend = c("Common","Defassa"), 
       fill = c("#BEBADA","#FDB462"))
dev.off()


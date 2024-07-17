library(scico)
library(gplots)
library("RColorBrewer")
library(tidyverse)

t<- read.table("/home/users/xi/Waterbuck_project/0.ColleaguesContribution/Xiaodong/right.tree.no.common.defassa.Fixed.name.Goat.108inds.8outgroups.Bohor_reedbuck.557.list.abbababa.txt",sep="\t",header=T,stringsAsFactors = FALSE)

t$pvalue = 2*pnorm(-abs(t$Z))
t$bonferroni_p = p.adjust(t$pvalue,method = "bonferroni")


mylist=strsplit(sort(unique(t$H1_Pops.H2_Pops.H3_Pops)),"-")[1:35]
cols = unique(unlist(lapply(mylist, '[', 3)))
rows = unique(unlist(lapply(mylist, function(x) paste(x[1], x[2], sep="-"))))

D.matrix = matrix(NA,ncol=length(cols),nrow= length(rows))

rownames(D.matrix) <- rows
#rownames(D.matrix) = c(cct,ddt)
colnames(D.matrix) = cols #7


for ( i in 1:ncol(D.matrix) ) { # by column
    p3 = colnames(D.matrix)[i]
    for ( j in 1:nrow(D.matrix) ) { # by row
        p12 = rownames(D.matrix)[j]
        p1 = strsplit(p12,"-")[[1]][1]
        p2 = strsplit(p12,"-")[[1]][2]
        index = which( t$H3_Pops == p3 & t$H1_Pops == p1 & t$H2_Pops == p2) 
        if( length(index) <1) {
            D.matrix[j,i] = NA
        } else {
            D.matrix[j,i] = mean(t$Dstat[index],na.rm=T)
        }
    }
}

Dmatrix = D.matrix[,c(3,4,7)]
#### two all-kob tests with the highest Dstat values do not make sense, 
#Dmatrix[11,3] = NA
#Dmatrix[12,4] =NA

heatmap.2(Dmatrix,dendrogram='none',Rowv=F,
          Colv = F,col=pal,tracecol=NA ,offsetRow =-50,colRow = 'black',colCol = 'black',key.title = "",lwid = c(1,2),
          keysize=1,cexRow = 1.5,cexCol = 1.5,colsep=1:ncol(Dmatrix),rowsep=1:nrow(Dmatrix),symbreaks = FALSE) 


### p value

P.matrix = matrix(NA,ncol=length(cols),nrow= length(rows))

rownames(P.matrix) <- rows
colnames(P.matrix) = cols


for ( i in 1:ncol(P.matrix) ) { # by column
    p3 = colnames(P.matrix)[i]
    for ( j in 1:nrow(P.matrix) ) { # by row
        p12 = rownames(P.matrix)[j]
        p1 = strsplit(p12,"-")[[1]][1]
        p2 = strsplit(p12,"-")[[1]][2]
        index = which( t$H3_Pops == p3 & t$H1_Pops == p1 & t$H2_Pops == p2) 
        if( length(index) <1) {
            P.matrix[j,i] = NA
        } else {
            pvalue = median(t$bonferroni_p[index],na.rm=T)
            if (pvalue <=0.05) {
                P.matrix[j,i] <- "*"
            } else {
                P.matrix[j,i] <- "n.s."
            }
        }
    }
}


Pmatrix = P.matrix[,c(3,4,7)]
#### two all-kob tests with the highest Dstat values do not make sense,

#Pmatrix[11,3] = NA
#Pmatrix[12,4] =NA


heatmap.2(Dmatrix,dendrogram='none',Rowv=F,
          Colv = F,col=pal,tracecol=NA ,offsetRow =-50,colRow = 'black',colCol = 'black',key.title = "",lwid = c(1,2),
          keysize=1,cexRow = 0.8,cexCol = 0.8,colsep=1:ncol(D.matrix),rowsep=1:nrow(Dmatrix),symbreaks = FALSE,cellnote=Pmatrix)

heatmap.2(P.matrix,dendrogram='none',Rowv=F,
          Colv = F,col=pal,tracecol=NA ,offsetRow =-50,colRow = 'black',colCol = 'black',key.title = "",lwid = c(1,2),
          keysize=1,cexRow = 1.5,cexCol = 1.5,colsep=1:ncol(P.matrix),rowsep=1:nrow(P.matrix),symbreaks = FALSE) 


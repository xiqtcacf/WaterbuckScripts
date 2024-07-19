source("/home/genis/github/apoh/apohFuns.R")

colpal <- c("purple", "red", "darkgreen", "grey", "darkblue", "lightblue", "brown", "lightgreen", "wheat", "yellow")

waterbuckCols <- c ("Luangwa"="purple", "Matetsi"="red", "Nairobi"="dark green", "Samburu"="grey","Kafue"="dark blue","KVNP"="light blue","Maswa"="brown","QENP"="light green","Samole"="wheat", "Ugalla"="yellow")
refpops10 <- c("Matetsi","Samburu","KVNP","Samole","Kafue","Ugalla","QENP","Nairobi","Maswa","Luangwa") ###k=10

colord <- c("KVNP", "Luangwa", "QENP", "Matetsi", "Samole", "Kafue", "Maswa", "Samburu", "Ugalla", "Nairobi")
colpal <- waterbuckCols[colord]



f <- "/home/genis/africa1kg/waterbuck/recent_hybrids/run2/tt_2boot.anccoef"

parentalQ <- read_ancestries_boot(f)[[1]][[1]]
#parentalQ <- parentalQ[, -c(1, ncol(parentalQ)-1, ncol(parentalQ))]

k <- parentalQ[[1]][[1]]


ids <- c("Samburu-415", "Luangwa-2535", "QENP-1161", "Samburu-393", "Nairobi-427", "Ugalla-2994", "Ugalla-2996", "Samburu-392", "Samburu-398", "Nairobi-429", "Samburu-408")


outdir <- "/home/genis/africa1kg/waterbuck/recent_hybrids/run2/plots"

for(i in 1:length(parentalQ)){

    qq <- parentalQ[[i]]
    pedigrees <- getAllSortedPedigrees(qq)
    outpng <- paste0(outdir, "/pedigree_", ids[i], ".png") 
    outpdf <- paste0(outdir, "/pedigree_", ids[i], ".pdf")

    bitmap(outpng, h=6, w=6, res=300)
    plotPedigree(pedigrees[[1]], title=paste0("Most compatible recent admixture pedigree\n", ids[i]), colpal=colpal)
    dev.off()

    outpng <- paste0(outdir, "/pairedanc_", ids[i], ".png")    

    pdf(outpdf, h=6, w=6)
    plotPedigree(pedigrees[[1]], title=paste0("Most compatible recent admixture pedigree\n", ids[i]), colpal=colpal)
    dev.off()

    outpng <- paste0(outdir, "/pairedanc_", ids[i], ".png")
    bitmap(outpng, h=4, w=6, res=300)
   # par(mar=c(5.1,6,4.1,2.1))
    plotOrderedEstimates(parentalQ=qq, pedigrees=pedigrees[1], main_title=paste0("Paired ancestry proportions\n", ids[i]), colpal=colpal)
    dev.off()

        outpdf <- paste0(outdir, "/pairedanc_", ids[i], ".pdf")
    pdf(outpdf, h=8, w=12)
   # par(mar=c(5.1,6,4.1,2.1))
    plotOrderedEstimates(parentalQ=qq, pedigrees=pedigrees[1], main_title=paste0("Paired ancestry proportions\n", ids[i]), colpal=colpal)
    dev.off()


}

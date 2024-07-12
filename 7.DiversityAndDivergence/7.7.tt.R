#!/usr/bin/env R

# Rscript tt.R <path to unfolded 2dsfs> <prefix for output files> <assumed mutation rate> <assumed generation time>

# Parameters used
#insfs <- "2dsfs/Matetsi1615-Samole7605_unfolded.2dsfs"
#mu <- "1.43e-8"
#g <- "7.1"

# Load data and parameters
args <- commandArgs(trailingOnly=T)

insfs <- args[1]
outprefix <- args[2]
mu <- as.numeric(args[3])
g <- as.numeric(args[4])

dat <- read.table(insfs,stringsAsFactors=F)
sfs <- dat$V2

# Reformat sfs 
names(sfs) <- sapply(lapply(strsplit(dat$V1,"_"), function(x) lapply(strsplit(x, split="/"), function(x) sum(as.integer(x)))), paste, collapse="")

# Apply formulae from Equation 3 in Sjodin et al. (2021)
res <- c(
    alpha1 = unname((2 * sfs["11"]) / (2 * sfs["21"] + sfs["11"])),
    alpha2 = unname((2 * sfs["11"]) / (2 * sfs["12"] + sfs["11"])),
    theta = unname(3/8 * (((2*sfs["21"] + sfs["11"]) * (2 * sfs["12"] + sfs["11"])) / (sfs["11"])) / sum(sfs)),
    t1 = unname((sfs["10"]/2 + sfs["20"] - (((2*sfs["21"] + sfs["11"]) * (6*sfs["12"] + sfs["11"])) / (8*sfs["11"]))) / sum(sfs)),
    t2 = unname((sfs["01"]/2 + sfs["02"] - (((6*sfs["21"] + sfs["11"]) * (2*sfs["12"] + sfs["11"])) / (8*sfs["11"]))) / sum(sfs)),
    v1 = unname(((sfs["10"] + sfs["12"]) / 2 - sfs["11"] * ((2*sfs["20"] - sfs["12"]) / (2*sfs["21"] - sfs["11"]))) / sum(sfs)),
    v2 = unname(((sfs["01"] + sfs["21"]) / 2 - sfs["11"] * ((2*sfs["02"] - sfs["21"]) / (2*sfs["12"] - sfs["11"]))) / sum(sfs))
)


# Write out results
outfile <- paste0(outprefix, ".tt.params.res")
write.table(x = t(t(res)), file = outfile,col.names=F, row.names=T,quote=F)

# Scale time in years
res_scaled<- c(
    T1 = unname(res["t1"] * g / mu),
    T2 = unname(res["t2"] * g / mu),
    Na = unname(res["theta"] / mu),
    `mu (assumed)` = mu,
    `g (assumed)` = g
)

outfile <- paste0(outprefix, ".tt.split.years")
write.table(x = t(t(res_scaled)), file = outfile ,col.names=F, row.names=T,quote=F)

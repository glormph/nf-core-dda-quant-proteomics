#!/usr/bin/env Rscript

library(DEqMS)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)
psms = args[1]
filtered_features = args[2]
setname = args[3]

psms = read.table(psms, header=T, sep="\t", comment.char="", quote="")
colnames(psms)[2] = "Feature"
features = read.table(filtered_features, header=T, sep="\t", comment.char="", quote="")

# this is to ensure we get values for each l
psms[psms == 0] <- NA
psms = na.omit(psms)
lastcol = dim(psms)[2]
psms[,3:lastcol] = log2(psms[,3:lastcol])
psm.counts = as.data.frame(table(psms$Feature))
countout = data.frame(feat=psm.counts$Var1, matrix(psm.counts$Freq, nrow(psm.counts), 1))
write.table(countout, 'psmcounts', sep='\t', quote=F, row.names=F, col.names=F)
rownames(psm.counts) = psm.counts$Var1

# Normalize and filter features against the passed filtered_features table
proteins.nm = medianSweeping(psms, group_col=2)
proteins.nm = proteins.nm[rownames(proteins.nm) %in% features[,1], ]
featcol = data.frame(feats=rownames(proteins.nm))
write.table(cbind(featcol, proteins.nm), 'normalized_feats', sep='\t', quote=F, row.names=F)

## DEqMS
#cond = c('a', 'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b')
#sampleTable <- data.frame(row.names=colnames(psms)[3:lastcol], cond=as.factor(cond))
#print(sampleTable) # REMOVE
#feats.mx = as.matrix(proteins.nm)
#design = model.matrix(~cond, sampleTable)
#dim(feats.mx) # REMOVE
#print(head(feats.mx)) # REMOVE
#
#fit1 <- eBayes(lmFit(feats.mx, design))
#fit1$count <- psm.counts[rownames(fit1$coefficients), 2]  # add an attribute containing PSM/peptide count for each gene
#print(fit1)

# FIXME does not work, mailed YZ
#fit2 = spectraCounteBayes(fit1)

## This is dev branch and not out yet
##VarianceBoxplot(fit2,n=30, main=paste("Set", setname), xlab="PSM count") # this function only available in latest devel branch

#DEqMS.results = outputResult(fit2, coef_col=3)
#write.table(DEqMS.results, "DEqMS.analysis.out.txt", quote=F,sep="\t",row.names = F)
#head(DEqMS.results,n=5)




# FIXME filter proteins 1%FDR list

#!/usr/bin/env Rscript

library(ggplot2)
library(forcats)
library(reshape2)
args = commandArgs(trailingOnly=TRUE)
nrsets = as.numeric(args[1])
feattype = args[2]
peptable = args[3]
make_normtable = FALSE
if(length(args) == 4) {
 make_normtable = TRUE
 normtable = args[4]
}
feats = read.table("feats", header=T, sep="\t", comment.char = "", quote = "")

featcol = list(peptides='Peptide.sequence', proteins='Protein.ID', genes='Gene.ID', assoc='Gene.Name')[[feattype]]

# nrpsms
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  nrpsms = melt(feats, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$variable)
  summary_psms = aggregate(value~Set, nrpsms, median)
  colnames(summary_psms) = c('Set', paste('no_psm_', feattype, sep=''))
  nrpsms = aggregate(value~get(featcol)+Set, nrpsms, max)
  nrpsms = transform(nrpsms, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
  png('nrpsms')
  print(ggplot(nrpsms, aes(y=value, x=setrank)) +
    geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
    theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
    scale_x_reverse())
    dev.off()
}


# featyield
# summary table: 
  # PROT/GENE:
  # DONE OK median # unique peptides / protein 
  # DONE OK median # unique peptides / protein (gene)
  # DONE OK # proteins
  # DONE OK # proteins (gene)
  # median # PSMs for quant used
  # PEP:
  # DONE and check # unique peptides -> unique for protein?
  # PSM
  # # PSMs total
  # Make prot/gene optional!

qcols = colnames(feats)[grep('_q.value', colnames(feats))]
overlap = na.exclude(feats[qcols])
overlap = dim(overlap[apply(overlap, 1, function(x) any(x<0.01)),])[1]
if (feattype == 'peptides') {
  am_prots = melt(feats, id.vars=c(featcol, "Protein.s."), measure.vars=qcols)
  am_prots$nrprots = lengths(regmatches(am_prots$Protein.s., gregexpr(';', am_prots$Protein.s.))) + 1
  # aggregate feats and remove col with ; where?
} else {
  am_prots = melt(feats, id.vars=featcol, measure.vars=qcols)
  pepcols = colnames(feats)[grep('Amount.Unique', colnames(feats))]
  pepprots = melt(feats, id.vars=featcol, measure.vars=pepcols)
  pepprots$variable = sub('_Amount.Unique.peptides', '', pepprots$variable)
  pepmed = aggregate(value~variable, pepprots, median)
  colnames(pepmed) = c('Set', paste('no_pep_', feattype, sep=''))
}
am_prots = am_prots[!is.na(am_prots$value) | am_prots$value < 0.01,]
am_prots$Set = sub('_q.value', '', am_prots$variable)
png('featyield', height=(nrsets + 2) * 72)
if (feattype == 'peptides') {
  totalunique = length(unique(subset(am_prots, nrprots == 1)$Peptide.sequence))
  unipepprotnr = aggregate(nrprots ~ Set, subset(am_prots, nrprots == 1), length)
  am_prots = aggregate(Peptide.sequence ~ Set, am_prots, length)
  am_prots = merge(am_prots, unipepprotnr)
  names(am_prots) = c('Set', 'All', 'Non-shared (unique)')
  write.table(am_prots[,c(1,3)], 'summary.txt', row.names=F, quote=F, sep='\t')
  am_prots = melt(am_prots)
  colnames(am_prots)[3] = 'accession'
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), legend.text=element_text(size=20), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=30)) +
    geom_bar(aes(fct_rev(Set), y=accession, fill=variable), stat='identity', position='dodge') +
    geom_text(data=subset(am_prots, variable=='All'), aes(fct_rev(Set), accession/2, label=accession), colour="white", size=7, nudge_x=-0.25) + 
    geom_text(data=subset(am_prots, variable=='Non-shared (unique)'), aes(fct_rev(Set), accession/2, label=accession), colour="white", size=7, nudge_x=+0.25) + 
    ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal uniques: ', totalunique)))
} else {
  if (length(grep('plex', names(feats)))) {
    # if isobaric, then show summary table of feats 1%FDR AND quant
    tmtcols = colnames(feats)[setdiff(grep('plex', colnames(feats)), grep('quanted', colnames(feats)))]
    not_fullna = feats[rowSums(is.na(feats[,tmtcols])) != length(tmtcols),]
    sum_prots = melt(not_fullna, id.vars=featcol, measure.vars=qcols)
    sum_prots = sum_prots[!is.na(sum_prots$value) | sum_prots$value < 0.01,]
    sum_prots$Set = sub('_q.value', '', sum_prots$variable)
    sum_prots = aggregate(get(featcol) ~ Set, sum_prots, length)
    summary = merge(pepmed, sum_prots, by='Set', all.y=T)
    colnames(summary)[ncol(summary)] = paste('nr_', feattype, '_q', sep='')
  } else {
    # just nr of proteins 1% FDR, no quant
    summary = merge(pepmed, am_prots, by='Set', all.y=T)
    colnames(summary)[ncol(summary)] = paste('nr_', feattype, sep='')
  }
  summary = merge(summary, summary_psms, by='Set', all.y=T)
  summary[is.na(summary)] = 0
  write.table(summary, 'summary.txt', row.names=F, quote=F, sep='\t')
  am_prots = aggregate(get(featcol) ~ Set, am_prots, length)
  colnames(am_prots)[2] = 'accession'
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), plot.title=element_text(size=30)) +
    geom_bar(aes(fct_rev(Set), y=accession), stat='identity') +
    geom_text(aes(fct_rev(Set), accession/2, label=accession), colour="white", size=10) + ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal identified: ', nrow(feats))))
}
dev.off()

# coverage
if (feattype == 'proteins') {
  covmed = median(feats$Coverage)
  png('coverage')
  plot = ggplot(feats) + geom_histogram(aes(Coverage), bins=50) + theme_bw()
  plot = plot + geom_label(x=max(ggplot_build(plot)$data[[1]]$x)*0.75, y=max(ggplot_build(plot)$data[[1]]$count)*0.75, label=sprintf('Median: %.3f', covmed))
  print(plot)
  dev.off()
}

#isobaric
if (length(grep('plex', names(feats)))) {
  tmtcols = colnames(feats)[setdiff(grep('plex', colnames(feats)), grep('quanted', colnames(feats)))]
  overlap = na.exclude(feats[c(tmtcols, qcols)])
  overlap = dim(overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),])[1]
  tmt = melt(feats, id.vars=featcol, measure.vars = tmtcols)
  tmt$Set = sub('_[a-z0-9]*plex.*', '', tmt$variable)
  tmt$variable = sub('.*_[a-z].*[0-9]*plex_', '', tmt$variable)
  outplot = ggplot(na.exclude(tmt)) + geom_boxplot(aes(fct_rev(Set), value, fill=fct_rev(variable)), position=position_dodge(width=1)) +
    coord_flip() + ylab('Fold change') + xlab('Channels') + theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), plot.title=element_text(size=30) ) + 
    theme(legend.text=element_text(size=20), legend.position="top", legend.title=element_blank()) +
    ggtitle(paste('Overlap with values in \nall ', length(tmtcols), 'channels: ', overlap))
  if (min(na.exclude(tmt$value)) >= 0) { outplot = outplot + scale_y_log10() }
  png('isobaric', height=(3 * nrsets + 1) * 72)
  print(outplot)
  dev.off()
}

#nrpsmsoverlapping
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[grep('plex', colnames(feats))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$variable)
  if (feattype == 'peptides') {
    nrpsms = aggregate(value~Peptide.sequence+Set, nrpsms, max)
  } else {
    nrpsms = aggregate(value~get(featcol)+Set, nrpsms, max)
  }
  nrpsms = transform(nrpsms, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
  png('nrpsmsoverlapping')
  print(ggplot(nrpsms, aes(y=value, x=setrank)) +
    geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
    theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
    scale_x_reverse())
    dev.off()
}


# percentage_onepsm
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[grep('quanted_psm_count', colnames(feats))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[grep('plex', colnames(feats))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_quanted_psm_count', '', nrpsms$variable) # this is enough for DEqMS pipeline
  # but need also following line for normal pipe bc nr-psms values are reported per channel
  nrpsms$Set = sub('_[a-z0-9]*plex.*', '', nrpsms$Set)
  feats_in_set = aggregate(value~Set, data=nrpsms, length) 
  feats_in_set$percent_single = aggregate(value~Set, data=nrpsms, function(x) length(grep('[^01]', x)))$value / feats_in_set$value * 100
  png('percentage_onepsm')
  print(ggplot(feats_in_set, aes(Set, percent_single)) +
    geom_col(aes(fill=Set)) + theme_bw() + ylab('% of identifications') +
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) )
  dev.off()
}

# ranked step plot MS1 peptide per protein
if (feattype != 'peptides') {
  peps = read.table(peptable, header=T, sep='\t', comment.char='', quote='')
  ms1qcols = grep('MS1.area', colnames(peps))
  nrpep_set = melt(peps, id.vars=c("Protein.s."), measure.vars=ms1qcols, na.rm=T)
  if (dim(nrpep_set)[1] != 0) {
    nrpep_set$Set = sub('_MS1.area.*', '', nrpep_set$variable)
    nrpep_set = aggregate(variable~Protein.s.+Set, nrpep_set, length) 
    nrpep_set = transform(nrpep_set, setrank=ave(variable, Set, FUN = function(x) rank(x, ties.method = "random")))
    png('ms1nrpeps')
    print(ggplot(nrpep_set, aes(y=variable, x=setrank)) +
      geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# peptides with MS1') +
      theme_bw() + 
      theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
      scale_x_reverse())
    dev.off()
  }
}

# precursorarea
if (length(grep('area', names(feats)))) {
    parea = melt(feats, id.vars=featcol, measure.vars = colnames(feats)[grep('area', colnames(feats))])
    parea$Set = sub('_MS1.*', '', parea$variable)
    png('precursorarea', height=(nrsets + 1) * 72)
    print(ggplot(parea) + 
      geom_boxplot(aes(fct_rev(Set), value)) + scale_y_log10() + coord_flip() + ylab("Intensity") + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank()))
    dev.off()
}


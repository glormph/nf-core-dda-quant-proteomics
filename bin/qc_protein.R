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

# featyield
qcols = colnames(feats)[grep('_q.value', colnames(feats))]
overlap = na.exclude(feats[qcols])
overlap = dim(overlap[apply(overlap, 1, function(x) any(x<0.01)),])[1]
if (feattype == 'peptides') {
  featcol = 'Peptide.sequence'
  am_prots = melt(feats, id.vars=c(featcol, "Protein.s."), measure.vars=qcols)
  am_prots$nrprots = lengths(regmatches(am_prots$Protein.s., gregexpr(';', am_prots$Protein.s.))) + 1
} else {
  featcol = 'Protein.accession'  
  am_prots = melt(feats, id.vars=featcol, measure.vars=qcols)
}
am_prots = am_prots[!is.na(am_prots$value),]
am_prots = am_prots[am_prots$value < 0.01,]
am_prots$Set = sub('_q.value', '', am_prots$variable)
png('featyield', height=(nrsets + 2) * 72)
if (feattype == 'peptides') {
  pepprotnr = aggregate(nrprots ~ Set, subset(am_prots, nrprots == 1), length)
  am_prots = aggregate(Peptide.sequence ~ Set, am_prots, length)
  am_prots = merge(am_prots, pepprotnr)
  names(am_prots) = c('Set', 'All', 'Non-shared (unique)')
  am_prots = melt(am_prots)
  colnames(am_prots)[3] = 'accession'
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), legend.text=element_text(size=20), legend.title=element_blank(), legend.position='top', plot.title=element_text(size=30)) +
    geom_bar(aes(fct_rev(Set), y=accession, fill=variable), stat='identity', position='dodge') +
    geom_text(data=subset(am_prots, variable=='All'), aes(fct_rev(Set), accession/2, label=accession), colour="white", size=7, nudge_x=-0.25) + ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal identified: ', dim(feats)[1])))
} else {
  am_prots = aggregate(Protein.accession ~ Set, am_prots, length)
  colnames(am_prots)[2] = 'accession'
  print(ggplot(am_prots) +
    coord_flip() + ylab('# identified') + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank(), plot.title=element_text(size=30)) +
    geom_bar(aes(fct_rev(Set), y=accession), stat='identity') +
    geom_text(aes(fct_rev(Set), accession/2, label=accession), colour="white", size=10) + ggtitle(paste('Overlap for all sets: ', overlap, '\nTotal identified: ', dim(feats)[1])))
}
dev.off()

# precursorarea
if (length(grep('area', names(feats)))) {
    parea = melt(feats, id.vars=featcol, measure.vars = colnames(feats)[grep('area', colnames(feats))])
    print(head(parea))
    parea$Set = sub('_MS1.*', '', parea$variable)
    png('precursorarea', height=(nrsets + 1) * 72)
    print(ggplot(parea) + 
      geom_boxplot(aes(fct_rev(Set), value)) + scale_y_log10() + coord_flip() + ylab("Intensity") + theme_bw() + theme(axis.title=element_text(size=30), axis.text=element_text(size=20), axis.title.y=element_blank()))
    dev.off()
}

# coverage
if (feattype == 'proteins') {
  covmed = median(feats$Coverage)
  png('coverage')
  plot = ggplot(feats) + geom_histogram(aes(Coverage), bins=50) + theme_bw()
  plot + geom_label(x=max(ggplot_build(plot)$data[[1]]$x)*0.75, y=max(ggplot_build(plot)$data[[1]]$count)*0.75, label=sprintf('Median: %.3f', covmed))
  print(plot)
  dev.off()
}

#isobaric
if (length(grep('plex', names(feats)))) {
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[intersect(grep('quanted', colnames(feats), invert=T), grep('plex', colnames(feats)))]
  overlap = na.exclude(feats[c(tmtcols, qcols)])
  overlap = dim(overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),])[1]
  tmt = melt(feats, id.vars=featcol, measure.vars = tmtcols)
  tmt$Set = sub('_[a-z0-9]*plex.*', '', tmt$variable)
  tmt$variable = sub('.*_[a-z].*[0-9]*plex_', '', tmt$variable)
  png('isobaric', height=(3 * nrsets + 1) * 72)
  print(ggplot(na.exclude(tmt)) + geom_boxplot(aes(fct_rev(Set), value, fill=fct_rev(variable)), position=position_dodge(width=1)) +
    scale_y_log10() + coord_flip() + ylab('Fold change') + xlab('Channels') + theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), plot.title=element_text(size=30) ) + 
    theme(legend.text=element_text(size=20), legend.position="top", legend.title=element_blank()) +
    ggtitle(paste('Overlap with values in \nall ', length(tmtcols), 'channels: ', overlap)))
    dev.off()
}


# nrpsms
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[intersect(grep('quanted', colnames(feats)), grep('plex', colnames(feats)))]
  nrpsms = melt(feats, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_[a-z0-9]*plex.*', '', nrpsms$variable)
  if (feattype == 'peptides') {
    nrpsms = aggregate(value~Peptide.sequence+Set, nrpsms, max)
  } else {
    nrpsms = aggregate(value~Protein.accession+Set, nrpsms, max)
  }
  nrpsms = transform(nrpsms, setrank=ave(value, Set, FUN = function(x) rank(x, ties.method = "random")))
  png('nrpsms')
  print(ggplot(nrpsms, aes(y=value, x=setrank)) +
    geom_step(aes(color=Set), size=2) + scale_y_log10() + xlab('Rank') + ylab('# PSMs quanted') +
    theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) +
    scale_x_reverse())
    dev.off()
}

#nrpsmsoverlapping
if (length(grep('plex', names(feats)))) {
  nrpsmscols = colnames(feats)[intersect(grep('quanted', colnames(feats)), grep('plex', colnames(feats)))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[intersect(grep('quanted', colnames(feats), invert=T), grep('plex', colnames(feats)))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_[a-z0-9]*plex.*', '', nrpsms$variable)
  if (feattype == 'peptides') {
    nrpsms = aggregate(value~Peptide.sequence+Set, nrpsms, max)
  } else {
    nrpsms = aggregate(value~Protein.accession+Set, nrpsms, max)
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
  nrpsmscols = colnames(feats)[intersect(grep('quanted', colnames(feats)), grep('plex', colnames(feats)))]
  qcols = colnames(feats)[grep('_q.value', colnames(feats))]
  tmtcols = colnames(feats)[intersect(grep('quanted', colnames(feats), invert=T), grep('plex', colnames(feats)))]
  overlap = na.exclude(feats[c(featcol, tmtcols, qcols, nrpsmscols)])
  overlap = overlap[apply(overlap[qcols], 1, function(x) any(x<0.01)),]
  nrpsms = melt(overlap, id.vars=featcol, measure.vars = nrpsmscols)
  nrpsms$Set = sub('_[a-z0-9]*plex.*', '', nrpsms$variable)
  feats_in_set = aggregate(value~Set, data=nrpsms, length) 
  feats_in_set$percent_single = aggregate(value~Set, data=nrpsms, function(x) length(grep('[^01]', x)))$value / feats_in_set$value * 100
  png('percentage_onepsm')
  print(ggplot(feats_in_set, aes(Set, percent_single)) +
    geom_col(aes(fill=Set)) + theme_bw() + ylab('% of identifications') +
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20), legend.position="top", legend.text=element_text(size=20), legend.title=element_blank()) )
  dev.off()
}

# normfac
if (make_normtable && length(grep('plex', names(feats)))) {
  norms = read.table(normtable, header=F, sep='\t', comment.char='', quote='')
  colnames(norms) = c('variable', 'Set', 'value')
  norms$variable = sub('[a-z].*[0-9]*plex_', '', norms$variable)
  png('normfac', height=(3 * nrsets + 1) * 72)
  print(ggplot(norms, aes(fct_rev(Set), value, group=variable)) + 
    geom_col(aes(fill=variable), position=position_dodge(width=1)) +
    geom_text(aes(y=0.1, label=value), position=position_dodge(width=1), colour="white", size=7, hjust=0) +
    coord_flip() + ylab('Normalizing factor') + xlab('Channels') + theme_bw() + 
    theme(axis.title=element_text(size=30), axis.text=element_text(size=20)) + 
    theme(legend.text=element_text(size=20), legend.position="top", legend.title=element_blank()))
  dev.off()
}

# ranked step plot MS1 peptide per protein
if (feattype != 'peptides') {
  peps = read.table(peptable, header=T, sep='\t', comment.char='', quote='')
  ms1qcols = grep('MS1.area', colnames(peps))
  nrpep_set = melt(peps, id.vars=c("Protein.s."), measure.vars=ms1qcols, na.rm=T)
  print(dim(nrpep_set), length(nrpep_set))
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

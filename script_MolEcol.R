# Author : Arnaud Meng

# Loading lib
library(igraph)
library(ggplot2)
library(parallel)
library(gridExtra)

# Custom operator for "NOT IN"
`%ni%` <- Negate(`%in%`)

# Import vertices & edges files
v = read.table("matches-wo-selfhits.filtered_pident60_qcov80_scov80.blast.vertices", h=T, sep=";", quote=NULL)
e = read.table("matches-wo-selfhits.filtered_pident60_qcov80_scov80.blast.edges", h=T, sep=";")

# Import taxonomy info
info = read.table("SUPP/infoDino4-good.csv", h=T, sep=";", na.strings=c("NA","NaN", " ", ""))
colnames(info)=c("phylum",
                 "class",
                 "order",
                 "family",
                 "genus",
                 "specie",
                 "prefix",
                 "strain",
                 "MMETSP",
                 "chloroplast",
                 "phototrophy",
                 "mixotrophy",
                 "harmful.for.human",
                 "symbiont",
                 "kleptoplastidy",
                 "ichyotoxin",
                 "parasite",
                 "dmsp",
                 "thecate",
                 "cyst.formation",
                 "transcriptome.quality")

# Merging tables vertices(v) & info (info)
v.info = merge(x=v, y=info, by="prefix", all.x=TRUE)
v.info = v.info[!duplicated(v.info[,"name"]),]
v.info = v.info[,c("name",
                   "prefix",
                   "phylum",
                   "class",
                   "order",
                   "family",
                   "genus",
                   "specie",
                   "strain",
                   "MMETSP",
                   "chloroplast",
                   "phototrophy",
                   "mixotrophy",
                   "parasite",
                   "dmsp",
                   "kleptoplastidy",
                   "symbiont",
                   "harmful.for.human",
                   "ichyotoxin",
                   "thecate",
                   "cyst.formation",
                   "transcriptome.quality")]

# Import annotation table
a = read.table("SUPP/all.sorted.filled.tsv", h=F, sep=";")
colnames(a)=c("name",               # col 0 : is the id of the input sequence.
              "checksum",           # col 1 : is the crc64 (checksum) of the protein sequence (supposed to be unique).
              "seqlength",          # col 2 : is the length of the sequence (in AA).
              "anaysis_method",     # col 3 : is the anaysis method launched.
              "db_entry",           # col 4 : is the database members entry for this match.
              "db_desc",            # col 5 : is the database member description for the entry.
              "domain_start",       # col 6 : is the start of the domain match.
              "domain_end",         # col 7 : is the end of the domain match.
              "evalue",             # col 8 : is the evalue of the match (reported by member database anayling method).
              "matching_statut",    # col 9 : is the status of the match (T: true, M: marginal).
              "date",               # col 10 : is the date of the run.
              "interpro_entry",     # col 11 : is the corresponding InterPro entry (if iprlookup requested by the user).
              "interpro_desc",      # col 12 : is the description of the InterPro entry.
              "go")                 # col 13 : is the GO (gene ontology) description for the InterPro entry.

# Merging tables : vertices(v.info) & annotations (a)
c = merge(x=v.info, y=a, by=c("name"), all.x=TRUE)

# Sorting : 
c = c[ order(c$name, c$go, c$interpro_desc, c$db_desc), ]

# Deleting duplicated annotations
c.uniq = c[!duplicated(c$name),]

# If needed : saving tables
# write.table(x=c.uniq, file="c.uniq.tsv", sep="\t", col.names=T, row.names=F, quote=F)

################################################################################
#
# Adding SLIM terms : CF script "go2slim.R"
#
# or load : my.slim = read.csv("my.slim.txt", h=T, sep=";")
#
################################################################################

my.slim = read.csv("SUPP/my.slim.txt", h=T, sep=";")

# Merging tables : vertices(c.uniq) & slim (my.slim)
c.uniq.slim = merge(x=c.uniq, y=my.slim, by=c("name", "go"), all.x=TRUE)
c.uniq.slim = c.uniq.slim[,c("name",
                             "prefix",
                             "phylum",
                             "class",
                             "order",
                             "family",
                             "genus",
                             "specie",
                             "strain",
                             "MMETSP",
                             "chloroplast",
                             "phototrophy",
                             "mixotrophy",
                             "harmful.for.human",
                             "symbiont",
                             "kleptoplastidy",
                             "ichyotoxin",
                             "parasite",
                             "dmsp",
                             "thecate",
                             "cyst.formation",
                             "transcriptome.quality",
                             "checksum",
                             "seqlength",
                             "anaysis_method",
                             "db_entry",
                             "db_desc",
                             "domain_start",
                             "domain_end",
                             "evalue",
                             "matching_statut",
                             "date",
                             "interpro_entry",
                             "interpro_desc",
                             "go",
                             "MFslim",
                             "BPslim",
                             "CCslim")]

# graph creation
g = graph.data.frame(e, directed=F, vertices=c.uniq.slim)
# IGRAPH UN-- 1283775 6213111 --

# Removing vertices from "bad quality" transcriptomes
transcriptome.bad = which(V(g)$transcriptome.quality=="bad")
g.good.only = delete.vertices(g, transcriptome.bad)
# IGRAPH UN-- 1275911 6142013 --

# Components (CCs) with >= 2 vertices
comp.good.only = decompose.graph(g.good.only, min.vertices = 2)
saveRDS(comp.good.only, file="comp.good.only.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
comp.good.only = readRDS("comp.good.only.rds")

# CCs sizes
comp.sizes = unlist(lapply(comp.good.only, function(x) vcount(x)))

# Saving
saveRDS(v, file="v.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(e, file="e.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(e.id60, file="e.id60.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(a, file="a.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(info, file="info.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(v.info, file="v.info.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(c, file="c.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(c.uniq, file="c.uniq.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(c.uniq.slim, file="c.uniq.slim.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(g, file="g.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(g.good.only, file="g.good.only.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

################################################################################
# ----------------------- ESTIMATE CCs ANNOTATIONS ----------------------------#
################################################################################

# Function to retrieve the list of CCs with "n" different annotations 
# at annotation level "golevel"
cc.n.annotation = function(graph, n, golevel) {
  
  # print(" -------------------------------------------------- new ")
  
  # get and unlist go term for the CC
  goterm = as.character(V(graph)$go)
  goterm = strsplit(goterm, split='|', fixed=T)
  goterm = unlist(goterm)
  
  # to find which goterm do not exist in goslim.g
  # for (go in goterm) { if (!is.na(go) & go %ni% V(goslim.g)$name) { print(paste(go, "NOT EXISTING in goslim.g")) } }
  
  # get slim term corresponding to GO id within 'goterm' list and also corresponding to 'golevel'
  slimterm = V(goslim.g)[which(V(goslim.g)$name %in% goterm & V(goslim.g)$Ontology==golevel)]$Term
  
  if (length(slimterm)==n) { return(graph) }
  
}

# Function to retrieve the list of CCs with more than "n" annotations
cc.more.annotation = function(graph, n, golevel) {
  
  # get and unlist go term for the CC
  goterm = as.character(V(graph)$go)
  goterm = strsplit(goterm, split='|', fixed=T)
  goterm = unlist(goterm)
  
  # get slim term corresponding to GO id within 'goterm' list and also corresponding to 'golevel'
  slimterm = V(goslim.g)[which(V(goslim.g)$name %in% goterm & V(goslim.g)$Ontology==golevel)]$Term
  
  if (length(slimterm)>n) { return(graph) }
  
}

# Function to retrieve the list of CCs "full NA" (fully unannotated)
cc.na.annotation = function(graph) {
  
  # get and unlist go term for the CC
  annotation.table =  as.data.frame(table(V(graph)$interpro_desc, useNA=c("ifany")))
  if (dim(annotation.table)[1]==1) {
    if (is.na(annotation.table$Var1)==T) {
      return(graph)
    }
  }
}

# Function to display categorized CCs
cc.new.annotation = function(graph, golevel) {
  
  print(paste("---------------------------------------------", vcount(graph)))
  
  # get and unlist go term for the CC
  goterm = as.character(V(graph)$go)
  goterm = strsplit(goterm, split='|', fixed=T)
  goterm = unlist(goterm)
  
  table = as.data.frame(table(goterm, useNA="ifany"), stringsAsFactors=F)
  table[,"term"] = NA
  
  # ccna   
  if (dim(table)[1]==1) { if(is.na(table$goterm)==T) { print("YES -- ccna") } }
  
  else {
    
    slimterm = V(goslim.g)[which(V(goslim.g)$name %in% goterm & V(goslim.g)$Ontology==golevel)]$Term
    
    # unfound goID 
    if (identical(slimterm, character(0))==TRUE) { print("GOid not found") }
    
    else {
      
      # cc1
      if (length(slimterm)==1) { print("YES ---- cc1") }
      # cc2 
      if (length(slimterm)==2) { print("YES ------ cc2") }
      # cc3 
      if (length(slimterm)==3) { print("YES -------- cc3") }
      # cc3+ 
      if (length(slimterm)>3) { print("YES ---------- cc3+") }
      
    }
  }
}

test = lapply(comp.good.only[1:50], cc.new.annotation, golevel="BP")
test = test[!sapply(test, is.null)]

# id60
cc.with.1.annotation.BP = mclapply(comp.good.only, cc.n.annotation, n=1, golevel="BP", mc.cores=4)
cc.with.2.annotation.BP = mclapply(comp.good.only, cc.n.annotation, n=2, golevel="BP", mc.cores=4)
cc.with.3.annotation.BP = mclapply(comp.good.only, cc.n.annotation, n=3, golevel="BP", mc.cores=4)
cc.with.na.annotation = mclapply(comp.good.only, cc.na.annotation, mc.cores=4)
cc.with.more3.annotation.BP = mclapply(comp.good.only, cc.more.annotation, n=3, golevel="BP", mc.cores=4)
cc.with.1.annotation.BP = cc.with.1.annotation.BP[!sapply(cc.with.1.annotation.BP, is.null)]
cc.with.2.annotation.BP = cc.with.2.annotation.BP[!sapply(cc.with.2.annotation.BP, is.null)]
cc.with.3.annotation.BP = cc.with.3.annotation.BP[!sapply(cc.with.3.annotation.BP, is.null)]
cc.with.na.annotation = cc.with.na.annotation[!sapply(cc.with.na.annotation, is.null)]
cc.with.more3.annotation.BP = cc.with.more3.annotation.BP[!sapply(cc.with.more3.annotation.BP, is.null)]
saveRDS(cc.with.1.annotation.BP, file="cc.with.1.annotation.BP.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(cc.with.2.annotation.BP, file="cc.with.2.annotation.BP.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(cc.with.3.annotation.BP, file="cc.with.3.annotation.BP.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(cc.with.na.annotation, file="cc.with.na.annotation.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

################################################################################
# ----------------------- CORE/ACCESSORY/PAN ----------------------------------#
################################################################################

# Get the number of sequences per "prefix"
cc.prefix.list = as.data.frame(table(unlist(sapply(comp.good.only, function(x) V(x)$prefix))))
colnames(cc.prefix.list) = c("prefix", "nbseq")

# Sort by the number of sequences (descending)
cc.prefix.list = cc.prefix.list[order(cc.prefix.list$nbseq, decreasing=T),]

# Barplot : number of sequences per transcriptome 
pdf("nbseq.per.transcriptome.pdf")
par(mar=c(10, 4.1, 4.1, 2.1))
barplot(cc.prefix.list$nbseq, ylab="number of sequences", names=cc.prefix.list$prefix, las=2, cex.axis=0.6, cex.names=0.6)
abline(h=9000, col="red")
dev.off()

# Keeping transcriptomes with more than 9000 sequences
cc.prefix.list.9000 = cc.prefix.list[which(cc.prefix.list$nbseq>9000),]

list.vec = list()           # transcriptomes vector
core.size = vector()        # core vector
pan.size = vector()         # pan vector
acc.size = vector()         # accessory vector

count = length(cc.prefix.list.9000$prefix)

for (transcriptome in cc.prefix.list.9000$prefix) {
  
  list.vec = unlist(c(list.vec, transcriptome))
  print(paste("remaining : ", count))
  
  # 1 transcriptome
  if (length(list.vec)==1) { 
    core.size = NA       
    pan.size = NA
    acc.size = NA
  }
  # at least 2 transcriptomes
  else if (length(list.vec)>=2) {
    v.spec = c.uniq.slim[which(c.uniq.slim$prefix %in% list.vec),]
    v.spec.names = unique(v.spec$name)
    e.spec = e[which(e$from %in% v.spec.names & e$to %in% v.spec.names),]
    g.spec = graph.data.frame(e.spec, directed=F, vertices=v.spec)
    comp.spec = decompose.graph(g.spec, min.vertices=2)
    comp.prefix.list = lapply(comp.spec, function(x) return(unique(V(x)$prefix)))
    comp.size = unlist(lapply(comp.prefix.list, function(x) length(x)))
    
    core.size = c(core.size, length(which(comp.size==length(list.vec))))
    acc.size = c(acc.size, length(which(comp.size==1)))
    pan.size = c(pan.size, length(comp.spec))
    
    print(core.size)
    print(acc.size)
    print(pan.size)
    
  }
  count = count - 1
}

# Retrieve CCs : core / acc / pan
fct = function (g, n) { if ( length(unique(V(g)$prefix))==n ) { return (g) }}

core.comp = lapply(comp.spec, fct, n=length(list.vec))
core.comp = core.comp[!sapply(core.comp, is.null)]
acc.comp = lapply(comp.spec, fct, n=1)
acc.comp = acc.comp[!sapply(acc.comp, is.null)]
pan.comp = comp.spec

saveRDS(core.comp, file="core.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(acc.comp, file="acc.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(pan.comp, file="pan.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(core.size, file="core.size.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(acc.size, file="acc.size.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(pan.size, file="pan.size.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

# Tables creation
core.pan.table = cbind(seq(1,length(core.size),1), core.size, acc.size, pan.size)
colnames(core.pan.table) = c("nbtranscriptome","core","accesory", "pan")
core.pan.table = as.data.frame(core.pan.table)
saveRDS(core.pan.table, file="core.pan.table.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

# Plot
pdf("core.pan.genome_CC_NEW-prefix.pdf")
lp = ggplot(core.pan.table, aes(nbtranscriptome))
lp = lp + geom_line(aes(y = core, colour = "Core")) + geom_point(aes(y = core, colour = "Core"))
lp = lp + geom_line(aes(y = accesory, colour = "Accesory")) + geom_point(aes(y = accesory, colour = "Accesory"))
lp = lp + geom_line(aes(y = pan, colour = "Pan")) + geom_point(aes(y = pan, colour = "Pan"))
lp = lp + geom_text(aes(y = core, label=core), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + geom_text(aes(y = accesory, label=accesory), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + geom_text(aes(y = pan, label=pan), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + xlab("number of transcriptomes")
lp = lp + ylab("number of connected components")
lp = lp + theme(legend.title=element_blank())
lp
dev.off()

pdf("core.pan.genome_CC_NEW-prefix_LOG.pdf")
lp = ggplot(core.pan.table, aes(nbtranscriptome))
lp = lp + geom_line(aes(y = core.log, colour = "Core")) + geom_point(aes(y = core.log, colour = "Core"))
lp = lp + geom_line(aes(y = accesory.log, colour = "Accesory")) + geom_point(aes(y = accesory.log, colour = "Accesory"))
lp = lp + geom_line(aes(y = pan.log, colour = "Pan")) + geom_point(aes(y = pan.log, colour = "Pan"))
lp = lp + geom_text(aes(y = core.log, label=core), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + geom_text(aes(y = accesory.log, label=accesory), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + geom_text(aes(y = pan.log, label=pan), vjust=2, hjust=1.5, size=1.5, angle=90)
lp = lp + xlab("number of transcriptomes")
lp = lp + ylab("number of connected components")
lp = lp + theme(legend.title=element_blank())
lp
dev.off()

# Predicting model core/accessory/pan of dinoflagellate proteome
# Purpose : try to predict how this proteome size could evolve 
# source : http://stackoverflow.com/questions/15535877/extrapolate-in-r-for-a-time-series-data

core.pan.table = readRDS("../RDS/core.pan.table.rds")
core.pan.table = core.pan.table[-1,]

x = core.pan.table$nbtranscriptome
y = core.pan.table$core

# Noise creation
noise = rnorm(length(x), mean=10, sd=80)
noisy.y = y + noise
plot(x,noisy.y,col='deepskyblue4',xlab='q',main='Observed data')
lines(x,y,col='firebrick1',lwd=3)

# According to my guess, the possible function could be : 
# y = N / x 
# with y the number of core genes, N the number of genes for the first point, 
# and x the number of transcriptome
# the more transcriptomes are involved, the more we could reach a "plateau" 
# that means : we reach a maximum number of core genes existing in the nature

dat = cbind(x, y, noisy.y)
colnames(dat) = c("nbtrans","nbcore","noisy_y")
dat = as.data.frame(dat)
dat$pred1 = predict(nls(nbcore~a/nbtrans, data=dat))
dat$pred2 = predict(nls(noisy_y~a/nbtrans, data=dat))

# pvalue of the a parameter : 
mod = nls(noisy_y~a/nbtrans, data=dat)
summary(mod)
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)    
# a   9008.3      318.1   28.32   <2e-16 ***

# correlation : prediction vs data 
cor(dat$nbcore, dat$pred1)
# [1] 0.976478

## extrapolate based on model
pred = data.frame(nbtrans=2:50)

mod = nls(nbcore~a/nbtrans, data=dat)
mod_noisy = nls(noisy_y~a/nbtrans, data=dat)

mod_noisy = nls(noisy_y~a/(nbtrans+b), data=dat, start=c(a=1, b=2))

pred$nbcore = predict(mod, newdata=pred)
pred$nbcore_noisy = predict(mod_noisy, newdata=pred)

p1 =    ggplot(dat, aes(x=nbtrans, y=noisy_y)) +
  geom_point(color="darkgrey") +
  geom_line(aes(y=nbcore, color="nbcore"), color="black" ,size=1) +
  geom_point(aes(y=nbcore, color="nbcore"), color="black") +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  geom_line(aes(y=pred1, colour="data_model"), size=1) +
  geom_line(data=pred, aes(x=nbtrans, y=nbcore_noisy, colour="prediction"), color="blue", linetype=2)
print(p1)

pdf("extrapolate_core_300.pdf")
print(p1)
dev.off()

write.table(pred, file="extrapolate_core.txt",row.names=F, quote=F, sep="\t")

# function to create random string of length = 12
MHmakeRandomString <- function(n=1, lenght=12) {
  # initialize vector
  randomString <- c(1:n)                 
  for (i in 1:n) { randomString[i] <- paste(sample(c(0:9, letters, LETTERS), lenght, replace=TRUE), collapse="") }
  return(randomString)
}

# function to extract sequences from each components of a list 
extract.cc.sequence = function(cc) {
  
  cc.nodes.index = ""
  
  # create component index
  index_name = MHmakeRandomString()
  index_nodes = vcount(cc)
  index_edges = ecount(cc)
  cc.nodes.index = as.character(paste(index_nodes, "_", index_edges, "_", index_name, sep=""))
  
  # recupere les noms des sequences
  # cc.nodes.names = as.data.frame(V(core.comp[[1]])$name)
  cc.nodes.names = as.data.frame(V(cc)$name)
  
  write.table(cc.nodes.names, 
              file=paste(cc.nodes.index, ".txt", sep=""), 
              append=FALSE, 
              sep="\n",
              row.names=FALSE,
              col.names=FALSE,
              quote=FALSE)
  
  return(cc.nodes.names)
}

# Retrieving IDs for core CCs in .txt files
lapply(core.comp, extract.cc.sequence)

# SHELL
#
# Concatenate ID files 
# cat *.txt > core.seq
#
# Extracting core sequence based on their IDs
# perl extractFromFasta.pl ALL.clean.pep list core.seq > core.fasta
# 

core.blastp.nr = read.table("core.vs.nr.diamond.1e-25.blastp", h=F)
colnames(core.blastp.nr) = c("name","nr.subject","nr.id_perc","nr.align_len","nr.mismatches","nr.gap_openings","nr.qstart","nr.qend","nr.sstart","nr.send","nr.evalue","nr.bitscore")
comp.good.only_ccindex = readRDS("../../RDS/comp.good.only_ccindex.rds")
core.comp = readRDS("../../RDS/core.core.rds")
core.comp.seq = as.data.frame(unlist(lapply(core.comp, function(x) return(V(x)$name))))
colnames(core.comp.seq) = "name"
c.uniq.slim = readRDS("../../RDS/c.uniq.slim.rds")

a = merge(x=core.comp.seq, y=comp.good.only_ccindex, by="name", all.x=T)
b = merge(x=a, y=c.uniq.slim, by="name", all.x=T)
c = merge(x=b, y=core.blastp.nr, by="name", all.x=T)
c = c[!duplicated(c[,"name"]),]


cclist = list()
for(cc in c$cc_index){
  seq = c[which(c$cc_index==cc),]
  a = as.numeric(dim(seq[which(!is.na(seq$nr.subject)),])[1])
  u = as.numeric(dim(seq[which(is.na(seq$nr.subject)),])[1])
  
  if(a!=0 & u!=0) {
    prop = a/u
    if (prop>=1.0) { cclist = c(cclist,cc) }
  }
}

cclist = as.data.frame(unique(unlist(cclist)))
colnames(cclist) = "cc_index"
res = merge(x=cclist, y=a, by="cc_index", all.x=T)
res2 = merge(x=res, y=core.blastp.nr, by="name", all.x=T)
res2 = res2[!duplicated(res2[,"name"]),]
res2_NA = res2[which(is.na(res2$nr.subject)),]
res2_NA = res2_NA[!duplicated(res2_NA[,"name"]),]

# Comparisons to nr / sp / BUSCO
core.blastp.nr = read.table("core.vs.nr.diamond.1e-25.blastp", h=F)
core.blastp.sp = read.table("core.vs.sp.diamond.1e-25.blastp", h=F)
core.blastp.busco = read.table("core.vs.busco.busco.txt", fill=T)
colnames(core.blastp.nr) = c("name","nr.subject","nr.id_perc","nr.align_len","nr.mismatches","nr.gap_openings","nr.qstart","nr.qend","nr.sstart","nr.send","nr.evalue","nr.bitscore")
colnames(core.blastp.sp) = c("name","sp.subject","sp.id_perc","sp.align_len","sp.mismatches","sp.gap_openings","sp.qstart","sp.qend","sp.sstart","sp.send","sp.evalue","sp.bitscore")
colnames(core.blastp.busco) = c("busco_group","status","scaffold","bitscore","length") 

core.blastp.nr.seq = as.data.frame(unique(core.blastp.nr$name), stringsAsFactors=FALSE)
core.blastp.sp.seq = as.data.frame(unique(core.blastp.sp$name), stringsAsFactors=FALSE)
core.blastp.busco.seq = as.data.frame(unique(core.blastp.busco$scaffold), stringsAsFactors=FALSE)
colnames(core.blastp.nr.seq) = "name"
colnames(core.blastp.sp.seq) = "name"
colnames(core.blastp.busco.seq) = "name"

# function to number sequences to each components of a list 
number.cc.sequence = function(cc) {
  
  # create component index of length N
  n=1
  lenght=12
  # initialize vector 
  randomString = c(1:n)                 
  for (i in 1:n) { randomString[i] = paste(sample(c(0:9, letters, LETTERS), lenght, replace=TRUE), collapse="") }
  
  seq.table = as.data.frame(V(cc)$name)
  seq.table$cc_index = rep(randomString, dim(seq.table)[1])
  colnames(seq.table) = c("name","cc_index")
  
  return(seq.table)
}

# Retrieving CORE sequences with their CCs association ID
core.comp = readRDS("RDS/core.comp.rds")
core.seq = mclapply(core.comp, number.cc.sequence, mc.cores=8)
core.seq = do.call("rbind", core.seq)

# Merging tables
core.comp.nr.blastp = merge(x=core.blastp.nr.seq, y=core.seq, all.x=TRUE)
core.comp.sp.blastp = merge(x=core.blastp.sp.seq, y=core.seq, all.x=TRUE)
core.comp.busco.blastp = merge(x=core.blastp.busco.seq, y=core.seq, all.x=TRUE)

matched.seq = rbind(core.comp.nr.blastp, core.comp.sp.blastp, core.comp.busco.blastp)
unmatched.seq = core.seq[which(core.seq$name %ni% matched.seq$name),]

# How many core full NA CCs ?
core.full.na = lapply(core.comp, cc.na.annotation)
core.full.na = core.full.na[!sapply(core.full.na, is.null)]
core.full.na.seq = lapply(core.full.na, number.cc.sequence)
core.full.na.seq = do.call("rbind", core.full.na.seq)
core.no.hit = merge(x=unmatched.seq, y=core.full.na.seq, by="name")

# Get the CCs IDs
unique(core.no.hit$cc_index.y)
# [1] "RXJHlIhcw9Cm" "5tQiTiLhfS48" "AAe3pgnSyLlU" "DyqXnkvAU9x4" "ZXNkylzG7wZZ"
# [6] "9Gz7Y7C3Imgd" "radACsCnwOri" "Ik38rtFb8FBM" "kVCzD7jEKtA6" "tVXQSIcySEnl"
# [11] "V6vuicCBoxfx" "1MDzFe4QWHhn" "48FcxWMbaivn" "8toCsU0zFozp" "C7Utibrxnf8y"
# [16] "6QcuqfoUXL6O"

nacore1 = core.full.na.seq[which(core.full.na.seq$cc_index=="RXJHlIhcw9Cm"),]
nacore2 = core.full.na.seq[which(core.full.na.seq$cc_index=="5tQiTiLhfS48"),]
nacore3 = core.full.na.seq[which(core.full.na.seq$cc_index=="AAe3pgnSyLlU"),]
nacore4 = core.full.na.seq[which(core.full.na.seq$cc_index=="DyqXnkvAU9x4"),]
nacore5 = core.full.na.seq[which(core.full.na.seq$cc_index=="ZXNkylzG7wZZ"),]
nacore6 = core.full.na.seq[which(core.full.na.seq$cc_index=="9Gz7Y7C3Imgd"),]
nacore7 = core.full.na.seq[which(core.full.na.seq$cc_index=="radACsCnwOri"),]
nacore8 = core.full.na.seq[which(core.full.na.seq$cc_index=="Ik38rtFb8FBM"),]
nacore9 = core.full.na.seq[which(core.full.na.seq$cc_index=="kVCzD7jEKtA6"),]
nacore10 = core.full.na.seq[which(core.full.na.seq$cc_index=="tVXQSIcySEnl"),]
nacore11 = core.full.na.seq[which(core.full.na.seq$cc_index=="V6vuicCBoxfx"),]
nacore12 = core.full.na.seq[which(core.full.na.seq$cc_index=="1MDzFe4QWHhn"),]
nacore13 = core.full.na.seq[which(core.full.na.seq$cc_index=="48FcxWMbaivn"),]
nacore14 = core.full.na.seq[which(core.full.na.seq$cc_index=="8toCsU0zFozp"),]
nacore15 = core.full.na.seq[which(core.full.na.seq$cc_index=="C7Utibrxnf8y"),]
nacore16 = core.full.na.seq[which(core.full.na.seq$cc_index=="6QcuqfoUXL6O"),]

write(as.character(nacore1$name), file="nacore1.txt")
write(as.character(nacore2$name), file="nacore2.txt")
write(as.character(nacore3$name), file="nacore3.txt")
write(as.character(nacore4$name), file="nacore4.txt")
write(as.character(nacore5$name), file="nacore5.txt")
write(as.character(nacore6$name), file="nacore6.txt")
write(as.character(nacore7$name), file="nacore7.txt")
write(as.character(nacore8$name), file="nacore8.txt")
write(as.character(nacore9$name), file="nacore9.txt")
write(as.character(nacore10$name), file="nacore10.txt")
write(as.character(nacore11$name), file="nacore11.txt")
write(as.character(nacore12$name), file="nacore12.txt")
write(as.character(nacore13$name), file="nacore13.txt")
write(as.character(nacore14$name), file="nacore14.txt")
write(as.character(nacore15$name), file="nacore15.txt")
write(as.character(nacore16$name), file="nacore16.txt")

# Core CCs fonctions summary (INTERPRO - GOslim)
core.interpro = merge(x=core.seq, y=c.uniq.slim, all.x=T)

core.interpro.MFslim = as.data.frame(table(core.interpro$MFslim)) 
core.interpro.BPslim = as.data.frame(table(core.interpro$BPslim)) 
core.interpro.CCslim = as.data.frame(table(core.interpro$CCslim))
core.interpro.MFslim =  core.interpro.MFslim[order(core.interpro.MFslim$Freq, decreasing=T),]
core.interpro.BPslim =  core.interpro.BPslim[order(core.interpro.BPslim$Freq, decreasing=T),]
core.interpro.CCslim =  core.interpro.CCslim[order(core.interpro.CCslim$Freq, decreasing=T),]
core.interpro.MFslim = cbind(head(core.interpro.MFslim, 10), rep("MF",10))
core.interpro.BPslim = cbind(head(core.interpro.BPslim, 10), rep("BP",10))
core.interpro.CCslim = cbind(head(core.interpro.CCslim, 10), rep("CC",10))
colnames(core.interpro.MFslim) = c("annotation","Freq","level")
colnames(core.interpro.BPslim) = c("annotation","Freq","level")
colnames(core.interpro.CCslim) = c("annotation","Freq","level")

core.interpro.go = rbind(core.interpro.MFslim, core.interpro.BPslim, core.interpro.CCslim)
core.interpro.go$annotation = factor(core.interpro.go$annotation, levels=core.interpro.go$annotation)

pdf("core_interpro_annotation_top10.pdf")
ggplot(core.interpro.go, aes(x=annotation, y=Freq, fill=level)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_x_discrete(labels=core.interpro.go$annotation) + 
  theme(axis.text=element_text(size=8)) +
  coord_flip()
dev.off()

# STATISTIC : is there some under-represented core CCs ?
mat = lapply(core.comp, function(x) as.data.frame(as.numeric(table(V(x)$prefix))))
mat = lapply(mat, function(x) t(x))
div.shannon.prefix = unlist(lapply(mat, function(x) round(diversity(x, "shannon"),3)))
hmax =  log(43) # maximum number of species
div.pielou.prefix = div.shannon.prefix/hmax

################################################################################
# ------------------------- FUNCTIONAL TRAITS ---------------------------------#
################################################################################

kleptoplastidy.comp = function(graph, param) { 
  res = as.data.frame(table(V(graph)$kleptoplastidy, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

mixotrophy.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$mixotrophy, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

symbionte.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$symbiont, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

ichyotoxin.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$ichyotoxin, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

parasite.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$parasite, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

dmsp.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$dmsp, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

thecate.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$thecate, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

cyst.forming.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$cyst.formation, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

phototrophy.comp = function(graph, param) {
  res = as.data.frame(table(V(graph)$phototrophy, useNA="ifany"))
  if (dim(res)[1]==1) { if (is.na(res$Var1)==F & res$Var1==param) { return(graph) }}
}

chloroplast.comp = function(graph) {
  res = as.data.frame(table(V(graph)$chloroplast, useNA="ifany"))
  if (dim(res)[1]>=1) {
    if ("R" %ni% res$Var1 & NA %ni% res$Var1) {
      return(graph)
    }
  }
}

harmful4humans.comp = function(graph) {
  res = as.data.frame(table(V(graph)$harmful.for.human, useNA="ifany"))
  if (dim(res)[1]>=1) {
    if ("no" %ni% res$Var1 & NA %ni% res$Var1) {
      return(graph)
    }
  }
}

# Full CCs defintion
klepto.comp = lapply(comp.good.only, kleptoplastidy.comp, param="yes")
klepto.comp = klepto.comp[!sapply(klepto.comp, is.null)]

mixo.comp = lapply(comp.good.only, mixotrophy.comp, param="yes")
mixo.comp = mixo.comp[!sapply(mixo.comp, is.null)]

symb.comp = lapply(comp.good.only, symbionte.comp, param="yes")
symb.comp = symb.comp[!sapply(symb.comp, is.null)]

h4h.comp = sapply(comp.good.only, harmful4humans.comp)
h4h.comp = h4h.comp[!sapply(h4h.comp, is.null)]

ichyo.comp = lapply(comp.good.only, ichyotoxin.comp, param="yes")
ichyo.comp = ichyo.comp[!sapply(ichyo.comp, is.null)]

parasite.comp = lapply(comp.good.only, parasite.comp, param="yes")
parasite.comp = parasite.comp[!sapply(parasite.comp, is.null)]

dmsp.comp = lapply(comp.good.only, dmsp.comp, param="yes")
dmsp.comp = dmsp.comp[!sapply(dmsp.comp, is.null)]

thecate.comp = lapply(comp.good.only, thecate.comp, param="yes")
thecate.comp = thecate.comp[!sapply(thecate.comp, is.null)]

cyst.forming.comp = lapply(comp.good.only, cyst.forming.comp, param="yes")
cyst.forming.comp = cyst.forming.comp[!sapply(cyst.forming.comp, is.null)]

chloro.comp = lapply(comp.good.only, chloroplast.comp)
chloro.comp = chloro.comp[!sapply(chloro.comp, is.null)]

phototrophy.comp = lapply(comp.good.only, phototrophy.comp, param="yes")
phototrophy.comp = phototrophy.comp[!sapply(phototrophy.comp, is.null)]

saveRDS(klepto.comp, file="klepto.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(h4h.comp, file="h4h.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(mixo.comp, file="mixo.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(symb.comp, file="symb.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(ichyo.comp, file="ichyo.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(parasite.comp, file="parasite.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(dmsp.comp, file="dmsp.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(thecate.comp, file="thecate.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(cyst.forming.comp, file="cyst.forming.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)
saveRDS(chloro.comp, file="chloro.comp.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

# Extracting h4h CCs
h4h.seq = unlist(lapply(h4h.comp, function(x) return(V(x)$name)))
write(h4h.seq, file="h4h.seq")

# Extracting non-h4h CCs
all.seq = unlist(V(g.goog.only)$name)
all.wo.h4h.seq = all.seq[which(all.seq %ni% h4h.seq)]
write(all.wo.h4h.seq, file="all.wo.h4h.seq")

# How many transcriptomes per functional trait ?
info = readRDS("../RDS/info.rds")
symb.pre = length(unique(info[which(info$symbiont=="yes" & info$transcriptome.quality=="good"),"prefix"]))
chlo.pre = length(unique(info[which(info$chloroplast!="R"  & info$transcriptome.quality=="good"),"prefix"]))
mixo.pre = length(unique(info[which(info$mixo=="yes" & info$transcriptome.quality=="good"),"prefix"]))
h4h.pre = length(unique(info[which(info$harmful.for.human!="no" & info$transcriptome.quality=="good"),"prefix"]))
klep.pre = length(unique(info[which(info$kleptoplastidy=="yes" & info$transcriptome.quality=="good"),"prefix"]))
ichy.pre = length(unique(info[which(info$ichyotoxin=="yes" & info$transcriptome.quality=="good"),"prefix"]))
para.pre = length(unique(info[which(info$parasite=="yes" & info$transcriptome.quality=="good"),"prefix"]))
dmsp.pre = length(unique(info[which(info$dmsp=="yes" & info$transcriptome.quality=="good"),"prefix"]))
thec.pre = length(unique(info[which(info$thecate=="yes" & info$transcriptome.quality=="good"),"prefix"]))
cyst.pre = length(unique(info[which(info$cyst.formation=="yes" & info$transcriptome.quality=="good"),"prefix"]))

# max symbiont transcriptomes = 12
# max chloroplast transcriptomes = 43
# max mixotrophy transcriptomes = 18
# max harmful4humans transcriptomes = 14
# max kleptoplastidy transcriptomes = 2
# max ichyotoxin transcriptomes = 3
# max parasite transcriptomes = 1
# max dmsp transcriptomes = 25
# max thecate transcriptomes = 28 
# max cyst transcriptomes = 16

# Seeking full NA CCs for each functional trait

chloro.na.cc = mclapply(chloro.comp, cc.na.annotation, mc.cores=4)
chloro.na.cc = chloro.na.cc[!sapply(chloro.na.cc, is.null)]
saveRDS(chloro.na.cc, file="chloro.na.cc.rds", ascii=FALSE, version=NULL, compress=TRUE, refhook=NULL)

# CC chloro - full NA = 220,081
# CC symbiont - full NA = 64,203
# CC thecate - full NA = 111,538
# CC DMSP - full NA = 99,906
# CC mixotrophy - full NA = 40,042
# CC h4h - full NA = 31,496
# CC klepto - full NA = 5,948
# CC ichyo - full NA = 3,404
# CC cyst - full NA = 54,796
# CC parasite - full NA = 640

# Looking for CCs having 1/2/3... TRANSCRIPTOMES for each functional trait
symb.comp = readRDS("../RDS/symb.comp.rds")
chloro.comp = readRDS("../RDS/chloro.comp.rds")
mixo.comp = readRDS("../RDS/mixo.comp.rds")
h4h.comp = readRDS("../RDS/h4h.comp.rds")
klepto.comp = readRDS("../RDS/klepto.comp.rds")
ichyo.comp = readRDS("../RDS/ichyo.comp.rds")
para.comp = readRDS("../RDS/para.comp.rds")
dmsp.comp = readRDS("../RDS/dmsp.comp.rds")
thecate.comp = readRDS("../RDS/thecate.comp.rds")
cyst.comp = readRDS("../RDS/cyst.comp.rds")

# Function to retrieve the list  of CCs with "n" different transcriptomes 
cc.n.transcriptome = function(graph, n) {
  
  # get, unlist and count the number of prefix found in the graph 
  pre = as.character(V(graph)$prefix)
  pre = unlist(pre)
  pre = unique(pre)
  nbpre = length(pre)
  
  if (nbpre==n) { return(graph) }
}

# symb CCs (DO THE SAME THING FOR OTHER FUNCTIONAL TRAITS)
symb.comp.12.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 12), mc.cores=8) ; symb.comp.12.prefix = symb.comp.12.prefix[!sapply(symb.comp.12.prefix, is.null)] ; saveRDS(symb.comp.12.prefix, file="symb.comp.12.prefix.rds")
symb.comp.11.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 11), mc.cores=8) ; symb.comp.11.prefix = symb.comp.11.prefix[!sapply(symb.comp.11.prefix, is.null)] ; saveRDS(symb.comp.11.prefix, file="symb.comp.11.prefix.rds")
symb.comp.10.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 10), mc.cores=8) ; symb.comp.10.prefix = symb.comp.10.prefix[!sapply(symb.comp.10.prefix, is.null)] ; saveRDS(symb.comp.10.prefix, file="symb.comp.10.prefix.rds")
symb.comp.9.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 9), mc.cores=8) ; symb.comp.9.prefix = symb.comp.9.prefix[!sapply(symb.comp.9.prefix, is.null)] ; saveRDS(symb.comp.9.prefix, file="symb.comp.9.prefix.rds")
symb.comp.8.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 8), mc.cores=8) ; symb.comp.8.prefix = symb.comp.8.prefix[!sapply(symb.comp.8.prefix, is.null)] ; saveRDS(symb.comp.8.prefix, file="symb.comp.8.prefix.rds")
symb.comp.7.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 7), mc.cores=8) ; symb.comp.7.prefix = symb.comp.7.prefix[!sapply(symb.comp.7.prefix, is.null)] ; saveRDS(symb.comp.7.prefix, file="symb.comp.7.prefix.rds")
symb.comp.6.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 6), mc.cores=8) ; symb.comp.6.prefix = symb.comp.6.prefix[!sapply(symb.comp.6.prefix, is.null)] ; saveRDS(symb.comp.6.prefix, file="symb.comp.6.prefix.rds")
symb.comp.5.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 5), mc.cores=8) ; symb.comp.5.prefix = symb.comp.5.prefix[!sapply(symb.comp.5.prefix, is.null)] ; saveRDS(symb.comp.5.prefix, file="symb.comp.5.prefix.rds")
symb.comp.4.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 4), mc.cores=8) ; symb.comp.4.prefix = symb.comp.4.prefix[!sapply(symb.comp.4.prefix, is.null)] ; saveRDS(symb.comp.4.prefix, file="symb.comp.4.prefix.rds")
symb.comp.3.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 3), mc.cores=8) ; symb.comp.3.prefix = symb.comp.3.prefix[!sapply(symb.comp.3.prefix, is.null)] ; saveRDS(symb.comp.3.prefix, file="symb.comp.3.prefix.rds")
symb.comp.2.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 2), mc.cores=8) ; symb.comp.2.prefix = symb.comp.2.prefix[!sapply(symb.comp.2.prefix, is.null)] ; saveRDS(symb.comp.2.prefix, file="symb.comp.2.prefix.rds")
symb.comp.1.prefix = mclapply(symb.comp, function(x) cc.n.transcriptome(x, 1), mc.cores=8) ; symb.comp.1.prefix = symb.comp.1.prefix[!sapply(symb.comp.1.prefix, is.null)] ; saveRDS(symb.comp.1.prefix, file="symb.comp.1.prefix.rds")

table = cbind(seq(1,12,1),c(length(symb.comp.1.prefix),
                            length(symb.comp.2.prefix),
                            length(symb.comp.3.prefix),
                            length(symb.comp.4.prefix),
                            length(symb.comp.5.prefix),
                            length(symb.comp.6.prefix),
                            length(symb.comp.7.prefix),
                            length(symb.comp.8.prefix),
                            length(symb.comp.9.prefix),
                            length(symb.comp.10.prefix),
                            length(symb.comp.11.prefix),
                            length(symb.comp.12.prefix)))
colnames(table) = c("nb_prefix","nb_comp")
write.table(table, file="symb.table.txt", row.names=F, quote=F)

# Plot
tab = read.table("symb.table.txt", h=T)

pdf("symb_bp.pdf")
bp = barplot(tab$nb_comp, 
             names.arg=tab$nb_prefix,
             ylim=c(0,max(tab$nb_comp+2000)))
text(x=bp, y=tab$nb_comp, label=tab$nb_comp, pos=3, cex=0.8)
dev.off()

pdf("symb_rose.pdf")
bp = ggplot(tab, aes(factor(nb_prefix), log(nb_comp)))
bp = bp + geom_bar(stat="identity", width=0.8, color="black", fill="#999999")
bp = bp + scale_y_continuous(breaks=0:10)
bp = bp + coord_polar()
bp = bp + labs(x = "", y = "")
bp = bp + geom_text(aes(label=paste("(",nb_comp,")",sep=""), y=log(nb_comp)+3), size=2.2)
bp = bp + theme(panel.grid.minor = element_blank())
bp
dev.off()

# Function to draw the graph
plot.my.graph = function(net) {
  
  par(mfrow=c(3, 3))
  
  # assign the "genus" attribute as the vertex color
  vpref = unique(V(net)$genus)
  vcol  = rainbow(length(vpref))[rank(vpref)]
  color.vertex = data.frame(v1=unique(vpref), v2=unique(vcol))
  colnames(color.vertex) = c("genus", "color")
  
  # assign color to each species
  V(net)$color=V(net)$genus
  for (i in 1:dim(color.vertex)[1]) { V(net)$color=gsub(color.vertex[i,1], color.vertex[i,2], V(net)$color) }
  
  # plot the network
  network_plot = plot(net,
                      vertex.label = "",
                      vertex.label.color='black',
                      vertex.label.cex=0.4,
                      vertex.frame.color='white',
                      vertex.size=15,
                      edge.color='grey',
                      edge.curved=FALSE,
                      layout=layout.lgl,
                      main=paste("Component : ",vcount(net),"v",ecount(net),"e", sep=""))
  
  # plot the genus distribution
  genus.distr = as.data.frame(table(V(net)$genus))
  colnames(genus.distr) = c("genus", "Freq")
  genus.distr = merge(x=genus.distr, y=color.vertex, all.x=TRUE, by="genus")
  rownames(genus.distr) = genus.distr$genus
  bp = barplot(genus.distr$Freq, 
               names.arg=genus.distr$genus,
               cex.names=0.8,
               las=2,
               col=as.character(genus.distr$color),
               main="Number of vertices per genus",
               ylab="Number of vertices")
  
  # plot the toppology indexes
  # assortativity
  assortativity = round(assortativity.degree(net, directed=F),3)
  # calculate diversity at genus level
  mat = as.data.frame(table(V(net)$genus))
  colnames(mat) = c("genus", "Freq")
  rownames(mat) = mat$genus
  mat$genus = NULL
  mat = t(mat)
  div.shannon.genus = round(diversity(mat, "shannon"),3)
  div.simpson.genus = round(diversity(mat, "simpson"),3)
  divplot =  plot(1, type="n", axes=F, xlab="", ylab="")
  text(x = 1, 
       y = 1,
       cex = 0.8,
       col = "black",
       paste("Assortativity = ", assortativity, "\n",
             "Shannon diversity (at genus level) = ", div.shannon.genus, "\n",
             "Simpson diversity (at genus level) = ", div.simpson.genus))
  
  # annotation
  MFslim = as.data.frame(table(V(net)$MFslim, useNA='ifany')) 
  BPslim = as.data.frame(table(V(net)$BPslim, useNA='ifany')) 
  CCslim = as.data.frame(table(V(net)$CCslim, useNA='ifany'))
  MFslim =  MFslim[order(MFslim$Freq, decreasing=T),]
  BPslim =  BPslim[order(BPslim$Freq, decreasing=T),]
  CCslim =  CCslim[order(CCslim$Freq, decreasing=T),]
  
  colnames(MFslim) = c("annotation","Freq")
  colnames(BPslim) = c("annotation","Freq")
  colnames(CCslim) = c("annotation","Freq")
  
  # MFslim[is.na(MFslim$annotation),"annotation"] = "na"
  
  barplot(MFslim$Freq, names.arg=MFslim$annotation, cex.names=0.8, las=2, main="MF GOslim", ylab="Number of vertices")
  barplot(BPslim$Freq, names.arg=BPslim$annotation, cex.names=0.8, las=2, main="BP GOslim", ylab="Number of vertices")
  barplot(CCslim$Freq, names.arg=CCslim$annotation, cex.names=0.8, las=2, main="CC GOslim", ylab="Number of vertices")
  
}

symb.comp.8.prefix = readRDS("symb.comp.8.prefix.rds")

pdf("symb_comp8.pdf")
lapply(symb.comp.8.prefix, plot.my.graph)
dev.off()

# top 10 annotations for each trait
symb.annot.mf = as.data.frame(table(unlist(lapply(symb.comp, function(x) return(V(x)$MFslim)))))
symb.annot.bp = as.data.frame(table(unlist(lapply(symb.comp, function(x) return(V(x)$BPslim)))))
symb.annot.cc = as.data.frame(table(unlist(lapply(symb.comp, function(x) return(V(x)$CCslim)))))
symb.annot.mf =  symb.annot.mf[order(symb.annot.mf$Freq, decreasing=T),]
symb.annot.bp =  symb.annot.bp[order(symb.annot.bp$Freq, decreasing=T),]
symb.annot.cc =  symb.annot.cc[order(symb.annot.cc$Freq, decreasing=T),]
symb.annot.mf = cbind(head(symb.annot.mf, 10), rep("MF",10))
symb.annot.bp = cbind(head(symb.annot.bp, 10), rep("BP",10))
symb.annot.cc = cbind(head(symb.annot.cc, 10), rep("CC",10))
colnames(symb.annot.mf) = c("annotation","Freq","level")
colnames(symb.annot.bp) = c("annotation","Freq","level")
colnames(symb.annot.cc) = c("annotation","Freq","level")
symb.go = rbind(symb.annot.mf, symb.annot.bp, symb.annot.cc)
symb.go$annotation = factor(symb.go$annotation, levels=symb.go$annotation)

pdf("symb_interpro_annotation_top10.pdf")
ggplot(symb.go, aes(x=annotation, y=Freq, fill=level)) + 
  geom_bar(stat="identity", position="dodge") + 
  scale_x_discrete(labels=symb.go$annotation) + 
  ylab("Frequency") +
  xlab("Functional annotations") +
  ggtitle("Top 10 functional annotations per level in symbiotic CC") +
  theme(axis.text=element_text(size=8)) +
  coord_flip() + 
  annotate("text", x=28, y=10000, label= "NA sequences at MF level = 172557") +
  annotate("text", x=27, y=10000, label= "NA sequences at BP level = 205382") + 
  annotate("text", x=26, y=10000, label= "NA sequences at CC level = 221346")
dev.off()
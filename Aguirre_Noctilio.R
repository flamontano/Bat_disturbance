
##########################################################################
### ANALYSIS FOR MANUSCRIPT ON BOLIVIAN BAT ASSEMBLAGES - OCTOBER 2015
##########################################################################

#Flavia Montano

### II for Noctilionoidae (superfamily only)

### basics 

setwd("E:/2. MANUSCRIPTS/2. BAT ASSEMBLAGES_Aguirre et al/analyses")
load("E:/2. MANUSCRIPTS/2. BAT ASSEMBLAGES_Aguirre et al/analyses/analyses2015.RData")
#packages needed

library(ecodist)
library(vegan)
library(picante)
library(cluster)
library(ape)
library(geiger)

###   0. measuring the effort, S and number of individuals 

effort<-read.csv("effort.csv", header=T)
plot(effort)

library(Hmisc)
effort.<-as.matrix(effort)
rcorr(effort., type="pearson")


###  1. obtaining data

databats<-read.csv("databats.csv", header=T, na.strings=c(""))

# preparing taxonomic data from databats

noctilio.bats<-subset(databats, family!="Vespertilionidae")
noctilio.bats<-subset(noctilio.bats, family!="Molossidae")
noctilio.bats<-subset(noctilio.bats, family!="Emballonuridae")
noctilio.bats<-subset(noctilio.bats, family!="Thyropteridae")
#write.csv(noctilio.bats, file = "noctilio_bats.csv")

#creating all data.frames for noctilio

noctilio.taxo.abundance<-noctilio.bats[,-(2:4)]#final dataframe
noctilio.taxo.abundance <- data.frame(noctilio.taxo.abundance[,-1], row.names=noctilio.taxo.abundance[,1])#placing row names

noctilio.taxo.presabs<-data.frame(noctilio.taxo.abundance)
noctilio.taxo.presabs[noctilio.taxo.presabs>0]<-1
noctilio.taxo.presabs #final data frame

noctilio.taxo.pcent.sites<-prop.table(as.matrix(noctilio.taxo.abundance), margin=2)*100
colSums(noctilio.taxo.pcent.sites)# prueba de que las columnas suman 100
noctilio.taxo.pcent.sites<-as.data.frame(noctilio.taxo.pcent.sites) #final data frame

#final data frames (will need transpose)
noctilio.taxo.abundance
noctilio.taxo.presabs
noctilio.taxo.pcent.sites

# preparing functional data from databats

guilds.noctilio<-aggregate(.~guild, data=noctilio.bats, FUN=sum)

## for sigma plot curves
noctilio.guild.pres<-noctilio.bats
noctilio.guild.pres[noctilio.guild.pres>0]<-1
noctilio.guild.pres<-aggregate(.~guild, data=noctilio.guild.pres, FUN=sum)
write.csv(noctilio.guild.pres, file = "noctilio_guilds.csv")

noctilio.func.k.abundance<-guilds.noctilio[,-(2:4)]#final dataframe
noctilio.func.k.abundance<- data.frame(noctilio.func.k.abundance[,-1], row.names=noctilio.func.k.abundance[,1])#placing row names
noctilio.func.k.abundance # final

noctilio.func.k.presabs <-data.frame(noctilio.func.k.abundance)
noctilio.func.k.presabs[noctilio.func.k.presabs>0]<-1
noctilio.func.k.presabs #final

noctilio.func.k.abund.pcent.site<-data.frame(noctilio.func.k.abundance)
noctilio.func.k.abund.pcent.site<-prop.table(as.matrix(noctilio.func.k.abundance), margin=2)*100
colSums(noctilio.func.k.abund.pcent.site)# prueba de que las columnas suman 100
noctilio.func.k.abund.pcent.site<-as.data.frame(noctilio.func.k.abund.pcent.site) #final data frame
noctilio.func.k.abund.pcent.site#final data

#final data frames (will need transpose)
noctilio.func.k.abundance # final
noctilio.func.k.presabs #final
noctilio.func.k.abund.pcent.site#final data

### preparing phylogenetic data

library(picante)
library(geiger)
library(ape)
noctilio.taxo.abundance

#con geiger
noctilio.tree<-treedata(bolivia.tree, noctilio.taxo.abundance, warnings=F)
plot(noctilio.tree$phy, cex=0.5)
axisPhylo(cex=0.9)
windows()
plot.phylo(noctilio.tree$phy, type = "fan", show.tip.label = TRUE, cex=0.7, show.node.label=T, edge.color= "orange", edge.width=1, tip.color="black")
noctilio.tree<-noctilio.tree$phy


#final tree 
noctilio.tree

### 2. distance measures

#taxonomic distance matrices based on bray curtis dissimilarity
noc.dm.bc.taxo.pres.abs <-bcdist(t(noctilio.taxo.presabs))
noc.dm.bc.taxo.abundance <-bcdist(t(noctilio.taxo.abundance))
noc.dm.bc.taxo.abundance.pcensite <-bcdist(t(noctilio.taxo.pcent.sites))

#functional distance matrices based on bray curtis dissimilarity
noc.dm.bc.func.k.pres.abs <-bcdist(t(noctilio.func.k.presabs))
noc.dm.bc.func.k.abundance  <-bcdist(t(noctilio.func.k.abundance))
noc.dm.bc.func.k.abund.pcent.site <-bcdist(t(noctilio.func.k.abund.pcent.site))

#phylogenetic distance matrices based on sorensens similarity (then 1- to get distances)

x<-as.data.frame(t(noctilio.taxo.presabs))
noc.dm.sor.phylo.pres.abs<-GUniFrac(x,noctilio.tree, alpha=c(0,0.5,1))
noc.dm.sor.phylo.pres.abs<-noc.dm.sor.phylo.pres.abs$unifracs[,,2]#extracting matrix for alpha=0.5
noc.dm.sor.phylo.pres.abs<-as.matrix(noc.dm.sor.phylo.pres.abs)
noc.dm.sor.phylo.pres.abs<-as.dist(noc.dm.sor.phylo.pres.abs,diag=FALSE,upper=FALSE)#getting half matrix


x<-as.data.frame(t(noctilio.taxo.abundance))
noc.dm.sor.phylo.abundance<-GUniFrac(x,noctilio.tree, alpha=c(0,0.5,1))
noc.dm.sor.phylo.abundance<-noc.dm.sor.phylo.abundance$unifracs[,,2]#extracting matrix for alpha=0.5
noc.dm.sor.phylo.abundance<-as.matrix(noc.dm.sor.phylo.abundance)
noc.dm.sor.phylo.abundance<-as.dist(noc.dm.sor.phylo.abundance,diag=FALSE,upper=FALSE)#getting half matrix

x<-as.data.frame(t(noctilio.taxo.pcent.sites))
noc.dm.sor.phylo.abund.pcent.site<-GUniFrac(x,noctilio.tree, alpha=c(0,0.5,1))
noc.dm.sor.phylo.abund.pcent.site<-noc.dm.sor.phylo.abund.pcent.site$unifracs[,,2]#extracting matrix for alpha=0.5
noc.dm.sor.phylo.abund.pcent.site<-as.matrix(noc.dm.sor.phylo.abund.pcent.site)
noc.dm.sor.phylo.abund.pcent.site<-as.dist(noc.dm.sor.phylo.abund.pcent.site,diag=FALSE,upper=FALSE)#getting half matrix

### 3. Cluster analyses 

# taxo
noctilio.cluster.taxo.presabs<-agnes(noc.dm.bc.taxo.pres.abs, method="average")
summary(noctilio.cluster.taxo.presabs)
plot(noctilio.cluster.taxo.presabs)
windows()
pltree(noctilio.cluster.taxo.presabs)

noctilio.cluster.taxo.abundance<-agnes(noc.dm.bc.taxo.abundance, method="average")
summary(noctilio.cluster.taxo.abundance)
#plot(noctilio.cluster.taxo.abundance)

noctilio.cluster.taxo.ab.pcent.site<-agnes(noc.dm.bc.taxo.abundance.pcensite, method="average")
#plot(noctilio.cluster.taxo.ab.pcent.site)
pltree(noctilio.cluster.taxo.ab.pcent.site)


# func

noctilio.cluster.func.k.presabs<-agnes(noc.dm.bc.func.k.pres.abs, method="average")
#plot(noctilio.cluster.func.k.presabs)
noctilio.cluster.func.k.abundance<-agnes(noc.dm.bc.func.k.abundance, method="average")
#plot(noctilio.cluster.func.k.abundance)
noctilio.cluster.func.k.ab.pcent.site<-agnes(noc.dm.bc.func.k.abund.pcent.site, method="average")
#plot(noctilio.cluster.func.k.ab.pcent.site)

# phylo

noctilio.cluster.phylo.presabs<-agnes(noc.dm.sor.phylo.pres.abs, method="average")
#plot(noctilio.cluster.phylo.presabs)
noctilio.cluster.phylo.abundance<-agnes(noc.dm.sor.phylo.abundance, method="average")
#plot(noctilio.cluster.phylo.abundance)
noctilio.cluster.phylo.ab.pcent.site<-agnes(noc.dm.sor.phylo.abund.pcent.site, method="average")
#plot(noctilio.cluster.phylo.ab.pcent.site)

#Results of clusters

# agglomerative coefficients 

ag1<-noctilio.cluster.taxo.presabs$ac#conglomerate
ag2<-noctilio.cluster.taxo.abundance$ac
ag3<-noctilio.cluster.taxo.ab.pcent.site$ac
ag4<-noctilio.cluster.func.k.presabs$ac
ag5<-noctilio.cluster.func.k.abundance$ac
ag6<-noctilio.cluster.func.k.ab.pcent.site$ac
ag7<-noctilio.cluster.phylo.presabs$ac
ag8<-noctilio.cluster.phylo.abundance$ac
ag9<-noctilio.cluster.phylo.ab.pcent.site$ac

agc<-round(c(ag1,ag2,ag3,ag4,ag5,ag6,ag7,ag8,ag9),4)

# cophenetic correlation coefficient

cc1<-cor(noc.dm.bc.taxo.pres.abs, cophenetic(noctilio.cluster.taxo.presabs))
cc2<-cor(noc.dm.bc.taxo.abundance, cophenetic(noctilio.cluster.taxo.abundance))
cc3<-cor(noc.dm.bc.taxo.abundance.pcensite, cophenetic(noctilio.cluster.taxo.ab.pcent.site))
cc4<-cor(noc.dm.bc.func.k.pres.abs, cophenetic(noctilio.cluster.func.k.presabs))
cc5<-cor(noc.dm.bc.func.k.abundance, cophenetic(noctilio.cluster.func.k.abundance))
cc6<-cor(noc.dm.bc.func.k.abund.pcent.site, cophenetic(noctilio.cluster.func.k.ab.pcent.site))
cc7<-cor(noc.dm.sor.phylo.pres.abs, cophenetic(noctilio.cluster.phylo.presabs))
cc8<-cor(noc.dm.sor.phylo.abundance, cophenetic(noctilio.cluster.phylo.abundance))
cc9<-cor(noc.dm.sor.phylo.abund.pcent.site, cophenetic(noctilio.cluster.phylo.ab.pcent.site))
cophCor<-round(c(cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9),4)

# putting cluster results in one table

clusters<-c("Taxon.presabs","Taxon.abund","Taxon.scaled", "Func.presabs", "Func.abund", "Func.scaled","Phylo.presabs", "Phylo.abun","Phylo.scaled")
resultsTable.noc<-cbind(clusters,cophCor,agc)
write.csv(resultsTable.noc, file = "clusterresults_noc.csv")


### plotting restults

windows()
par(mfrow=c(1,3))
pltree(noctilio.cluster.taxo.presabs, main="presence/absence", )
pltree(noctilio.cluster.taxo.abundance, main="abundance")
pltree(noctilio.cluster.taxo.ab.pcent.site, main="scaled abundance")

windows()
par(mfrow=c(1,3))
pltree(noctilio.cluster.func.k.presabs, main="presence/absence", )
pltree(noctilio.cluster.func.k.abundance, main="abundance")
pltree(noctilio.cluster.func.k.ab.pcent.site, main="scaled abundance")
dev.off()

windows()
par(mfrow=c(1,3))
pltree(noctilio.cluster.phylo.presabs, main="presence/absence", )
pltree(noctilio.cluster.phylo.abundance, main="abundance")
pltree(noctilio.cluster.phylo.ab.pcent.site, main="scaled abundance")
dev.off()

# plot con labels (presence absence)

noctilio.cluster.func.k.presabs$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas medio", 
                                                      "Pie de Monte (Vargas)", "Yungas alto", "Amazonica-Ichilo", 
                                                      "Pie de Monte (Terán)", "Pie de Monte (Flores)", "Cerrado", "Madre de Dios")
noctilio.cluster.taxo.presabs$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas medio", 
                                                     "Pie de Monte (Vargas)", "Amazonica-Ichilo", 
                                                     "Pie de Monte (Flores)", "Cerrado", "Pie de Monte (Terán)", "Madre de Dios","Yungas alto")
noctilio.cluster.phylo.presabs$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas medio", 
                                                      "Pie de Monte (Vargas)", "Amazonica-Ichilo", 
                                                      "Pie de Monte (Terán)","Pie de Monte (Flores)", "Cerrado", "Madre de Dios","Yungas alto")
windows()
par(mfrow=c(1,3))
pltree(noctilio.cluster.func.k.presabs, main="Functional Diversity")
pltree(noctilio.cluster.taxo.presabs, main="Taxonomic Diversity")
pltree(noctilio.cluster.phylo.presabs, main="Phylogenetic diversity")
dev.off()


# plot con labels (relative abundance)

noctilio.cluster.func.k.ab.pcent.site$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas Medio", 
                                                       "Pie de Monte - Vargas", "Yungas Alto", "Amazonica-Ichilo", 
                                                       "Pie de Monte - Terán", "Pie de Monte - Flores", "Cerrado", "Madre de Dios")
noctilio.cluster.taxo.ab.pcent.site$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas Medio", 
                                                     "Pie de Monte - Vargas", "Amazonica-Ichilo", 
                                                     "Pie de Monte - Flores", "Cerrado", "Pie de Monte - Terán", "Madre de Dios","Yungas Alto")
noctilio.cluster.phylo.ab.pcent.site$order.lab<-cbind("Sabana Inundable", "Chaco", "Yungas Medio", 
                                                      "Pie de Monte - Vargas", "Amazonica-Ichilo", 
                                                      "Pie de Monte - Terán","Pie de Monte - Flores", "Cerrado", "Madre de Dios","Yungas Alto")
windows()
par(mfrow=c(1,3))
pltree(noctilio.cluster.func.k.ab.pcent.site, main="Functional Diversity")
pltree(noctilio.cluster.taxo.ab.pcent.site, main="Taxonomic Diversity")
pltree(noctilio.cluster.phylo.ab.pcent.site, main="Phylogenetic diversity")
dev.off()

### 4. CORRELATION-PEARSON & MANTEL TESTS 

#presence-absence


mantel.taxofunc.pa<-mantel(noc.dm.bc.taxo.pres.abs, noc.dm.bc.func.k.pres.abs, permutations=9999)
mantel.taxophylo.pa<-mantel(noc.dm.bc.taxo.pres.abs, noc.dm.sor.phylo.pres.abs, permutations=9999)
mantel.funcphylo.pa<-mantel(noc.dm.bc.func.k.pres.abs, noc.dm.sor.phylo.pres.abs, permutations=9999)

partialmantel.taxofunc.pa<-mantel.partial(noc.dm.bc.taxo.pres.abs, noc.dm.bc.func.k.pres.abs, noc.dm.sor.phylo.pres.abs, permutations=9999)
partialmantel.taxophylo.pa<-mantel.partial(noc.dm.bc.taxo.pres.abs, noc.dm.sor.phylo.pres.abs, noc.dm.bc.func.k.pres.abs, permutations=9999)
partialmantel.funcphylo.pa<-mantel.partial(noc.dm.sor.phylo.pres.abs, noc.dm.bc.func.k.pres.abs, noc.dm.bc.taxo.pres.abs, permutations=9999)

#abundance
mantel.taxofunc.ab<-mantel(noc.dm.bc.taxo.abundance, noc.dm.bc.func.k.abundance, permutations=9999)
mantel.taxophylo.ab<-mantel(noc.dm.bc.taxo.abundance, noc.dm.sor.phylo.abundance, permutations=9999)
mantel.funcphylo.ab<-mantel(noc.dm.bc.func.k.abundance, noc.dm.sor.phylo.abundance, permutations=9999)

partialmantel.taxofunc.ab<-mantel.partial(noc.dm.bc.taxo.abundance, noc.dm.bc.func.k.abundance, noc.dm.sor.phylo.abundance, permutations=9999)
partialmantel.taxophylo.ab<-mantel.partial(noc.dm.bc.taxo.abundance, noc.dm.sor.phylo.abundance,noc.dm.bc.func.k.pres.abs, permutations=9999)
partialmantel.funcphylo.ab<-mantel.partial(noc.dm.bc.func.k.abundance, noc.dm.sor.phylo.abundance, noc.dm.bc.taxo.abundance, permutations=9999)


#abundance.percent.per.site
mantel.taxofunc.per<-mantel(noc.dm.bc.taxo.abundance.pcensite, noc.dm.bc.func.k.abund.pcent.site, permutations=9999)
mantel.taxophylo.per<-mantel(noc.dm.bc.taxo.abundance.pcensite, noc.dm.sor.phylo.abund.pcent.site, permutations=9999)
mantel.funcphylo.per<-mantel(noc.dm.bc.func.k.abund.pcent.site, noc.dm.sor.phylo.abund.pcent.site, permutations=9999)

partialmantel.taxofunc.per<-mantel.partial(noc.dm.bc.taxo.abundance.pcensite, noc.dm.bc.func.k.abund.pcent.site, noc.dm.sor.phylo.abund.pcent.site, permutations=9999)
partialmantel.taxophylo.per<-mantel.partial(noc.dm.bc.taxo.abundance.pcensite, noc.dm.sor.phylo.abund.pcent.site,noc.dm.bc.func.k.abund.pcent.site, permutations=9999)
partialmantel.funcphylo.per<-mantel.partial(noc.dm.bc.func.k.abund.pcent.site, noc.dm.sor.phylo.abund.pcent.site, noc.dm.bc.taxo.abundance.pcensite, permutations=9999)


#Results of Mantel

#presence-absence

m1<-mantel.taxofunc.pa$statistic
m2<-mantel.taxophylo.pa$statistic
m3<-mantel.funcphylo.pa$statistic

pm1<-partialmantel.taxofunc.pa$statistic
pm2<-partialmantel.taxophylo.pa$statistic
pm3<-partialmantel.funcphylo.pa$statistic

#abundance
m4<-mantel.taxofunc.ab$statistic
m5<-mantel.taxophylo.ab$statistic
m6<-mantel.funcphylo.ab$statistic

pm4<-partialmantel.taxofunc.ab$statistic
pm5<-partialmantel.taxophylo.ab$statistic
pm6<-partialmantel.funcphylo.ab$statistic

#abundance.percent.per.site
m7<-mantel.taxofunc.per$statistic
m8<-mantel.taxophylo.per$statistic
m9<-mantel.funcphylo.per$statistic

pm7<-partialmantel.taxofunc.per$statistic
pm8<-partialmantel.taxophylo.per$statistic
pm9<-partialmantel.funcphylo.per$statistic

mantel<-round(c(m1,m2,m3,m4,m5,m6,m7,m8,m9),4)
parcial.mantel<-round(c(pm1,pm2,pm3,pm4,pm5,pm6,pm7,pm8,pm9),4)



#presence-absence

sm1<-mantel.taxofunc.pa$signif
sm2<-mantel.taxophylo.pa$signif
sm3<-mantel.funcphylo.pa$signif

spm1<-partialmantel.taxofunc.pa$signif
spm2<-partialmantel.taxophylo.pa$signif
spm3<-partialmantel.funcphylo.pa$signif

#abundance
sm4<-mantel.taxofunc.ab$signif
sm5<-mantel.taxophylo.ab$signif
sm6<-mantel.funcphylo.ab$signif

spm4<-partialmantel.taxofunc.ab$signif
spm5<-partialmantel.taxophylo.ab$signif
spm6<-partialmantel.funcphylo.ab$signif

#abundance.percent.per.site
sm7<-mantel.taxofunc.per$signif
sm8<-mantel.taxophylo.per$signif
sm9<-mantel.funcphylo.per$signif

spm7<-partialmantel.taxofunc.per$signif
spm8<-partialmantel.taxophylo.per$signif
spm9<-partialmantel.funcphylo.per$signif

signif.mantel<-round(c(sm1,sm2,sm3,sm4,sm5,sm6,sm7,sm8,sm9),6)
signif.parcial.mantel<-round(c(spm1,spm2,spm3,spm4,spm5,spm6,spm7,spm8,spm9),6)

#Table of mantel and partial mantel restuls

mantels<-c("Taxon-Funcc.presabs", "Taxon-Phylo.presabs", "func-Phylo.presabs", "Taxon-Funcc.ab", "Taxon-Phylo.abs", "func-Phylo.ab","Taxon-Funcc.scaled", "Taxon-Phylo.scaled", "func-Phylo.scaled")
partial.mantels<-c("Taxon-Funcc.presabs", "Taxon-Phylo.presabs", "func-Phylo.presabs", "Taxon-Funcc.ab", "Taxon-Phylo.abs", "func-Phylo.ab","Taxon-Funcc.scaled", "Taxon-Phylo.scaled", "func-Phylo.scaled")

Tablemantel.noc<-cbind(mantels,mantel,signif.mantel)
Table.partialmantel.noc<-cbind(partial.mantels,parcial.mantel,signif.parcial.mantel)

write.csv(Tablemantel.noc, file = "mantelresults_noc.csv")
write.csv(Table.partialmantel.noc, file = "partial_mantelresults_noc.csv")

##########
## ALTERNATIVE. using MRM

## exploring distances

func.max<-max(noc.dm.bc.func.k.abund.pcent.site)
func.min<-min(noc.dm.bc.func.k.abund.pcent.site)
func.mean<-mean(noc.dm.bc.func.k.abund.pcent.site)
func.sdev<-sd(noc.dm.bc.func.k.abund.pcent.site)

func.max
func.min
func.mean
func.sdev

## taxonomic

taxo.max<-max(noc.dm.bc.taxo.abundance.pcensite)
taxo.min<-min(noc.dm.bc.taxo.abundance.pcensite)
taxo.mean<-mean(noc.dm.bc.taxo.abundance.pcensite)
taxo.sdev<-sd(noc.dm.bc.taxo.abundance.pcensite)

taxo.max
taxo.min
taxo.mean
taxo.sdev

## phylo

phylo.max<-max(noc.dm.sor.phylo.abund.pcent.site)
phylo.min<-min(noc.dm.sor.phylo.abund.pcent.site)
phylo.mean<-mean(noc.dm.sor.phylo.abund.pcent.site)
phylo.sdev<-sd(noc.dm.sor.phylo.abund.pcent.site)

phylo.max
phylo.min
phylo.mean
phylo.sdev

library(ecodist)

#presence-absence

mrm.func.pa<-MRM(noc.dm.bc.func.k.pres.abs~noc.dm.bc.taxo.pres.abs+ noc.dm.sor.phylo.pres.abs, nperm=9999)

#abundance
mrm.func.ab<-MRM(noc.dm.bc.func.k.abundance~noc.dm.bc.taxo.abundance+ noc.dm.sor.phylo.abundance, nperm=9999)

#abundance.percent.per.site
mrm.func.per<-MRM(noc.dm.bc.func.k.abund.pcent.site~noc.dm.bc.taxo.abundance.pcensite+noc.dm.sor.phylo.abund.pcent.site, nperm=9999)


### partitioning the relative contribution of each variable
## with hierarchical partitioning algorithm (Chevan and Shuterland)

library(hier.part)

##running all possible models models ## only for per.site
## noc.dm.bc.taxo.abundance.pcensite =  Variable 1, 
## noc.dm.sor.phylo.abund.pcent.site = variable 2

mrm.func.per3<-MRM(noc.dm.bc.func.k.abund.pcent.site~noc.dm.bc.taxo.abundance.pcensite+noc.dm.sor.phylo.abund.pcent.site, nperm=9999)
mrm.func.per2<-MRM(noc.dm.bc.func.k.abund.pcent.site~noc.dm.sor.phylo.abund.pcent.site, nperm=9999)
mrm.func.per1<-MRM(noc.dm.bc.func.k.abund.pcent.site~noc.dm.bc.taxo.abundance.pcensite, nperm=9999)


#extracting R squared from each model ##


R3<-mrm.func.per3$r.squared[1:1]
R2<-mrm.func.per2$r.squared[1:1]
R1<-mrm.func.per1$r.squared[1:1]
R0<-0
gof<-matrix(nrow=4, ncol=1)
gof[,1]<-as.numeric(c(R0,R1,R2,R3))
gof

part.dist<- partition(gof, pcan = 2)
part.dist


#### otros analisis... 

simpson.noctilio<-read.csv("guilds_simpson.csv", header=T, row.names=1)
simpson<-diversity(simpson.noctilio, "inv")


save.image("E:/2. MANUSCRIPTS/2. BAT ASSEMBLAGES_Aguirre et al/analyses/analyses2015.RData")



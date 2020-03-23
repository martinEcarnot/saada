library(rnirs)
library(dplyr)
library(sampling)
library(ggplot2)

source("~/Documents/INRA/R/pre.R")
source("~/Documents/INRA/R/orthEf2.R")

if (T) {
load("~/Documents/INRA/melanges/2019/saada2019_sp_Brimrose")

mono=droplevels(sp[grep("mono",sp$Culture),])
# mono=droplevels(mono[mono$Trait=="S",])
ngeno=nlevels(mono$Genotype1)

# Pré
p=rbind(list('red',c(10,1,1)),list('snv',''),list('detr',''))
# p=rbind(list('red',c(10,1,1)),list('snv',''),list('sder',c(1,3,9)))
mono$xp=pre(mono$x,p)
comp=15
# 
# Porj orth sur l effet rep a partir des 29 premiers geno
mo29=droplevels(mono[mono$Genotype1 %in% levels(mono$Genotype1)[1:30],])
kW=orthEf2(mo29,"xp",c(4),5)  

i=40
j=41
g1=levels(mono$Genotype1)[i]
g2=levels(mono$Genotype1)[j]

ideux=mono$Genotype1 %in% c(g1,g2) #& mono$traitement=="T-S"
mono2=droplevels(mono[ideux,])

  ngeno2=nlevels(mono2$Genotype1)
  
  # mono2$xp = mono2$xp - (mono2$xp %*% kW %*% t(kW))
  
  Xr=mono2[mono2$Rep == 1,]
  Xu=mono2[mono2$Rep == 2,]
  
  # ?quilibrage des effectifs
  t=table(Xu$Genotype1)
  gmin=which.min(t)
  gmax=which.max(t)
  imax=which(Xu$Genotype1==names(gmax))
  iout=imax[sample(t[gmax],abs(t[1]-t[2]))]
  if (length(iout)!=0) {Xu=Xu[-iout,]}
  
  t=table(Xr$Genotype1)
  gmin=which.min(t)
  gmax=which.max(t)
  imax=which(Xr$Genotype1==names(gmax))
  iout=imax[sample(t[gmax],abs(t[1]-t[2]))]
  if (length(iout)!=0) {Xr=Xr[-iout,]}
  
  ngenoXr=nlevels(Xr$Genotype1)
  segm=segmcvkfold(n = nrow(Xr$xp), K = 2, nrep = 3)
  
  # fm <- fitcv(Xr$xp, Xr$Genotype1,fun = plsda,  ncomp = comp,  segm = segm, da="dalm")
  # z <- err(fm, ~ ncomp)
  # minw=selncomp.wold(z, nam="errp",plot=F,typ = "smooth")$sel
  # if (is.na(minw)) {minw=which.min(z$errp)}
  
  # fm <- plsda(Xr$xp,Xr$Genotype1,Xu$xp,Xu$Genotype1, ncomp=comp, da=dadis, diss="mahalanobis")
  fm <- plsda(Xr$xp,Xr$Genotype1,Xu$xp,Xu$Genotype1, ncomp=comp, da=dalm)
}
  

df=as.data.frame(rbind(fm$fm$Tr,fm$fm$Tu))
df=cbind(df,rbind(Xr,Xu))
ggplot(data = df,aes(x=comp1,y=comp2, col=Genotype1,size=Rep)) +geom_point() + stat_ellipse(level=0.6)


# correction de la moyenne du jeu de cal
mu=colMeans(fm$fm$Tu)
fm$fm$Tuc=fm$fm$Tu-mu

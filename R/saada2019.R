# Spectres SAADA 2019
# Essai de discri de genotypes

library(rnirs)
library(dplyr)
library(sampling)
library(e1071)
# source('Script_R_2020/SIGNE_maha0.R')
source("~/Documents/INRA/R/pre.R")
source("~/Documents/INRA/melanges/2019/segmFact.R")
source('~/Documents/INRA/melanges/2019/saada2019_prepa.R')

sp=saada2019_prepa()
mono=droplevels(sp[grep("monogénotype",sp$culture),])

# Pré
# p=rbind(list('red',c(1300, 1, 2)),list('snv',''),list('sder',c(2,3,11)))
p=rbind(list('red',c(100,1,2)),list('snv',''),list('sder',c(2,3,41)))
mono$xp=pre(mono$x,p)
comp=12
rpca=pca(mono$xp,ncomp=comp)
plotxy(rpca$Tr[, c(1,2)],group=mono$Rep)

stop()
# Discri
n=5
nbgeno=round(seq(2,10,length.out = n))
zt=list()
for (i in 1:length(nbgeno)) {
  print(i)
  genook=levels(mono$Genotype1)[sample(nbgeno[i])]
  monor=droplevels(mono[mono$Genotype1 %in% genook,])
  ngeno=nlevels(monor$Genotype1)
  segm=segmFact(monor, var=list("Genotype1","plante"),nel=list(ngeno,rep(4,ngeno)), nrep=1 )
  fm <- fitcv(monor$xp, monor$Genotype1,fun = plsda,  ncomp = 20,  segm = segm, print = F, da=dalm)
  z <- err(fm, ~ ncomp)
  zt[[i]]=z
}


z <- err(fm, ~ ncomp)
plotmse(z, nam = "errp")
min(z$errp)

# Avec 2 geno
r=list()
c=1
for (i in 1) { #:ngeno) {
  for (j in 1:ngeno) {
    print(i)
    g1=levels(mono$Genotype1)[i]
    g2=levels(mono$Genotype1)[j]
    ideux=mono$Genotype1 %in% c(g1,g2)
    mono2=droplevels(mono[ideux,])
    ngeno2=nlevels(mono2$Genotype1)
    segm=segmFact(mono2, var=list("Genotype1","plante"),nel=list(ngeno2,rep(4,ngeno2)), nrep=2 )
    fm <- fitcv(mono2$xp, mono2$Genotype1,fun = plsdalm,  ncomp = 10,  segm = segm)
    r[[c]]=fm
    # z <- err(fm, ~ ncomp)
    # plot(z$errp)
    c=c+1
  }
}

z <- err(fm, ~ ncomp)
# plotmse(z, nam = "errp")
# i6=fm$y$ncomp==5
# table(fm$y$x1[i6],fm$fit$x1[i6])

load(file = "~/Documents/INRA/R/sorties_saada2019")
rdf=data.frame(minerr=single(length(r)))
for (i in 1:length(r)) {
  print(i)
  z <- err(r[[i]], ~ ncomp)
  mini=which(z$errp == min(z$errp))[1]
  minw=selncomp.wold(z, nam="errp")$sel
  rdf$minerr[i]=z[minw, 4]
  rdf$mincmp[i]=mini
  rdf$minw[i]=minw
}

fm4=lapply(fm,function (x) {x[x$ncomp==4,]})


# SVM


ngeno=nlevels(monor$Genotype1)

Xr=monor[-segm[[1]][[1]],]
Xu=monor[segm[[1]][[1]],]
r=svm(Xr$xp,Xr$Genotype1)
table(predict(r,Xu$xp),Xu$Genotype1)





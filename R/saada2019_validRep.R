# On calibre sur 1 rep, on valide sur l'autre

library(rnirs)
library(dplyr)
library(sampling)

# source("~/Documents/INRA/R/pre.R")
# source("~/Documents/INRA/R/orthEf.R")
# source("~/Documents/INRA/melanges/2019/segmFact.R")
# source('~/Documents/INRA/melanges/2019/saada2019_prepa.R')
source("C:/Users/robot/Documents/Martin/melanges/pre.R")
source("C:/Users/robot/Documents/Martin/melanges/orthEf.R")
source("C:/Users/robot/Documents/Martin/melanges/2019/segmFact.R")
# source('C:/Users/robot/Documents/Martin/melanges/saada2019_prepa.R')

# sp=saada2019_prepa()
# sp$code.champ=as.factor(sp$code.champ)

load("C:/Users/robot/Documents/Martin/melanges/2019/saada2019_sp")
mono=droplevels(sp[grep("mono",sp$Culture),])
# mono=droplevels(mono[mono$Trait=="S",])
ngeno=nlevels(mono$Genotype1)

# Pré
# p=rbind(list('red',c(400,1,1)),list('snv',''),list('sder',c(2,2,65)))
p=rbind(list('red',c(400,1,1)),list('snv',''),list('detr',''))
# p=rbind(list('red',c(400,1,1)),list('detr',''))
mono$xp=pre(mono$x,p)
comp=15

# Porj orth sur l effet rep a partir des 29 premiers geno
mo29=droplevels(mono[mono$Genotype1 %in% levels(mono$Genotype1)[1:30],])

# Avec 2 geno
c=1
zt=list()
mt=list()
nam=list()
erVCV=data.frame(g1=character(),g2=character(),SECV=double(),SEV=double(),SEVmin=double(),SEV2=double(),VLcv=integer(),VLminV=integer(),SEV5=integer()) 
for (i in round(seq(31,ngeno, 1))) {  # 1:ngeno) {
  print(i)
  for (j in round(seq(30,ngeno-1, 1))) {  # 1:ngeno
    g1=levels(mono$Genotype1)[i]
    g2=levels(mono$Genotype1)[j]
    ideux=mono$Genotype1 %in% c(g1,g2) #& mono$traitement=="T-S"
    mono2=droplevels(mono[ideux,])
    
    mono2=orthEf(mo29,mono2,"xp",c(4),6)  # mono2=orthEf(mono2,mono2,"xp",c(1,5,9),1)
    mono2$xp=mono2$xort
    
    if (nlevels(mono2$Genotype1)!=2 | nlevels(mono2$Rep)!=2 | length(unique((mono2$code.champ)))<4) {
      sprintf('levels Geno = %d ;   levels rep = \n',nlevels(mono2$Genotype1),nlevels(mono2$rep))
      } else {
    ngeno2=nlevels(mono2$Genotype1)
    # segm <- vector("list", length = 2)
    # segm[[1]][[1]]=which(mono2$rep == "R1")
    # segm[[1]][[2]]=which(mono2$rep == "R2")

    # On fait une Crossval pour avoir le meilleur nVL, qu'on utilise pour faire une vraie valid R1 - R2
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
    segm=segmFact(Xr, var=list("Genotype1","plante"),nel=list(ngenoXr,rep(4,ngenoXr)), nrep=3 )
    fm <- fitcv(Xr$xp, Xr$Genotype1,fun = plsda,  ncomp = comp,  segm = segm, da="dadis", diss="mahalanobis")
    # fm <- fitcv(Xr$xp, Xr$Genotype1,fun = plsda,  ncomp = comp,  segm = segm, da="dalm")
    z <- err(fm, ~ ncomp)
    minw=selncomp.wold(z, nam="errp",plot=F,typ = "smooth")$sel
    if (is.na(minw)) {minw=which.min(z$errp)}

    fm <- plsda(Xr$xp,Xr$Genotype1,Xu$xp,Xu$Genotype1, ncomp=comp, da=dadis, diss="mahalanobis")
    # fm <- plsda(Xr$xp,Xr$Genotype1,Xu$xp,Xu$Genotype1, ncomp=comp, da=dalm)
    # r[[c]]=fm
    zV <- err(fm, ~ ncomp)
    # print(zV)
    erVCV1=data.frame(g1=g1,g2=g2,SECV=z[minw,"errp"],SEV=zV[minw,"errp"],SEVmin=min(zV$errp),SEV2=zV[2,"errp"], VLcv=minw, VLminV=which.min(zV$errp),SEC5=zV[5,"errp"])
    erVCV=rbind(erVCV,erVCV1)
    zt[[c]]=cbind(z[,4],zV[,4])
    mt[[c]]=c(minw,which.min(zV$errp))
    nam[[c]]=c(g1,g2)
    c=c+1
    }
  }
}
# erVCV=erVCV[-1,]
  print(colMeans(erVCV[,3:9]))
stop()

rdf=data.frame(minerr=single(length(r)-1), mincmp=single(length(r)-1))
plerr=matrix(,ncol=15)
res=data.frame()
for (i in 2:length(r)) {
  print(i)
  z <- err(r[[i]], ~ ncomp)
  plot(z$errp)
  mini=which(z$errp == min(z$errp))[1]
  # minw=selncomp.wold(z, nam="errp")$sel
  rdf$minerr[i]=z[mini, 4]
  rdf$mincmp[i]=mini
  plerr=rbind(plerr,z$errp)
  fm=r[[i]][1:3]
  fm2=lapply(fm,function (x) {x[x$ncomp==2,]})
  res=rbind(res,fm2$r[,c("rownam","y1")])
  # rdf$minw[i]=minw
}

stop()
pm=pmatch(res$rownam,rownames(mono), duplicates.ok = T)
res=cbind(res,mono[pm,])



rpca=pca(mono2$xp,ncomp=comp)
plotxy(rpca$Tr[, c(1,2)],group=mono2$Genotype1)



rt = matrix(unlist(lapply(r, "err", , "~ncomp")),ncol=20,byrow = T)
plotsp(r)


# LOO sur 2 geno, 1 rep, 1 traitement
mono2r1=droplevels(mono2[mono2$traitement=="T-S" & mono2$rep=="R1",])
segm=segmcvkfold(nrow(mono2r1),K=nrow(mono2r1))
fm <- fitcv(mono2r1$xp, mono2r1$Genotype1,fun = plsdalm,  ncomp = 20,  segm = segm, print = F)
z <- err(fm, ~ ncomp)
plot(z$errp)
fm4=lapply(fm,function (x) {x[x$ncomp==4,]})
table(fm4$fit$x1,fm4$y$x1)




# Sensibilité au nombre de Genotypes
# n=5
# nbgeno=round(seq(2,10,length.out = n))
# zt=list()
# for (i in 1:length(nbgeno)) {
#   print(i)
#   genook=levels(mono$Genotype1)[sample(nbgeno[i])]
#   monor=droplevels(mono[mono$Genotype1 %in% genook,])
#   ngeno=nlevels(monor$Genotype1)
#   segm <- vector("list", length = 1)
#   segm[[1]][[1]]=which(monor$rep == "R1")
#   segm[[1]][[2]]=which(monor$rep == "R2")
#   fm <- fitcv(monor$xp, monor$Genotype1,fun = plsdalm,  ncomp = 20,  segm = segm, print = T)
#   z <- err(fm, ~ ncomp)
#   zt[[i]]=z
# }



# Attribue une rep (1 ou 2) au traitement


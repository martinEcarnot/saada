segmFact <- function (dat,var=list(...),nel=list(...),nrep=1) {

  G=levels(dat$Genotype1)
  CC=levels(dat$code.champ)
  
  datg1=droplevels(dat[dat$Genotype1==G[1],])
  datg2=droplevels(dat[dat$Genotype1==G[2],])

  
    # On créé les groupes de CV en s'assurant que tous les niveaux d'un facteurs sont présents dans un groupe
  segm <- vector("list", length = nrep)

  for(i in 1:nrep) {)
    seggp <- vector("list", length = ngp)

    for(j in 1:ngp) {

      
      m=mstage(datdef, stage=list("cluster","stratified"),varnames=var,size=nel, method=list("srswor","srswor"))
      # seggp[[j]]=m$`2`$ID_unit  # Abandonné, car, un fois deflatté, les indices en sont plus ceux du df initial
      datval=getdata(datdef,m)[[2]]
      idval=which(rownames(dat)  %in%  rownames(datval))
      seggp[[j]]=idval
      datdef=datdef[-which(rownames(datdef)  %in%  rownames(datval)),]
    }
    segm[[i]] <- seggp
  }
  return(segm)
}

segmFact <- function (dat,var=list(...),nel=list(...),nrep=1) {

  # On créé les groupes de CV en s'assurant que tous les niveaux d'un facteurs sont présents dans un groupe
  segm <- vector("list", length = nrep)

  for(i in 1:nrep) {
    #  Nombre de groupe  = (effectif du + petit gpe de strate 1) / nbre d'element de ce stage demandé par groupe de strate n
    # ngp=floor(mean(table(dat[,var[[1]]]))/unique(nel[[2]]))
    ngp=floor(min(table(dat[,var[[1]]]))/unique(nel[[2]]))
    seggp <- vector("list", length = ngp)
    datdef=dat
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

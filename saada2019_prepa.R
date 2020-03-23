saada2019_prepa <- function () {

d='/home/ecarnot/Documents/INRA/melanges/2019'
# d="/run/user/1000/gvfs/smb-share:server=stocka2,share=agap-ble/Ble/PerkinElmer/frederic-2019/SAADA2019_epis/"
# sp=PE_read_dir(d)    # save(sp,file=file.path(d,'saada2019_spectres.txt'))
## sp=read.table(file.path(d,'saada2019_spectres.txt'), header=T, row.names=1, sep=';', dec='.')    # Plus utilisé
load(file.path(d,'saada2019_spectres.txt'))


# Enleve ce qui n'est pas de SAADA
sp=sp[-seq(1,29),]
sp=sp[grep("es",rownames(sp), invert=T),]

md=read.table(file.path(d,'Copie de SAADA_2019-Cahier_TOUT.csv'),header=T,  sep=';')[,c(1,2,4:7,10)]

# On duplique les infos de md sur sp
colnames(sp)=1e7/as.numeric(gsub("X","",colnames(sp)))
sp=data.frame(x=I(as.matrix(sp)))
class(sp$x)="matrix"
rownames(sp)[3988]=gsub("-","-8",rownames(sp)[3988])
rownames(sp)[660]=gsub("7.","7",rownames(sp)[660])
rownames(sp)[3686]=gsub("73","7",rownames(sp)[3686])

sp=cbind(code.champ=as.numeric(sub("\\-.*", "", rownames(sp))), plante=sub(".*-", "", rownames(sp)),sp)

sp$plante=as.factor(trimws(gsub("\\(1\\)","",sp$plante)))

sp=left_join(sp,md,by ="code.champ")
sp$Rep=as.factor(sp$Rep)

return(sp)

}

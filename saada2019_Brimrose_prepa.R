saada2019_Brimrose_prepa <- function () {

d='/home/ecarnot/Documents/INRA/melanges/2019'
# d="/run/user/1000/gvfs/smb-share:server=stocka2,share=agap-ble/Ble/PerkinElmer/frederic-2019/SAADA2019_epis/"
# sp=PE_read_dir(d)    # save(sp,file=file.path(d,'saada2019_spectres.txt'))
sp=read.table(file.path(d,'saada2019_Brimrose.txt'), header=T,  sep=';', dec='.')    # Plus utilisé


md=read.table(file.path(d,'Copie de SAADA_2019-Cahier_TOUT.csv'),header=T,  sep=';')[,c(1,2,4:7,10)]

# On duplique les infos de md sur sp
nam=gsub("grains_proteines_etalons-","",sp$X.1)
nam=gsub("grains_proteines_etalons_","",nam)
colnames(sp)=as.numeric(gsub("X","",colnames(sp)))
sp=data.frame(x=I(as.matrix(sp[,3:ncol(sp)])))
class(sp$x)="matrix"
sp$code.champ=as.numeric(gsub("-.*","",nam))

sp=left_join(sp,md,by ="code.champ")

return(sp)

}

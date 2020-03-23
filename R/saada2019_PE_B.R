saada2019_PE_B = function(zt_mt_B) {
  
library(ggplot2)

# load('zt_mt_B')
# load('zt_mt_PE')

ncpl=length(zt_mt_B[[1]])
# errPE=matrix(nrow=ncpl, ncol=8)  #matrix(nrow=ncpl,colnames=c('erCV','erVw','erVm'))
errB=matrix(nrow=ncpl, ncol=8)  #matrix(nrow=ncpl,colnames=c('erCV','erVw','erVm'))
for (i in 1:ncpl) {
  errB1=zt_mt_B[[1]][[i]]
  errB[i,1]=errB1[zt_mt_B[[2]][[i]][1],1]
  errB[i,2]=errB1[zt_mt_B[[2]][[i]][1],2]
  errB[i,3]=min(errB1[,2])
  errB[i,4]=zt_mt_B[[2]][[i]][1]
  errB[i,5]=zt_mt_B[[2]][[i]][2]
  errB[i,6]=errB1[6,2]
  errB[i,7]=abs(errB[i,4]- errB[i,5])
  
  # errPE1=zt_mt_PE[[1]][[i]]
  # errPE[i,1]=errPE1[zt_mt_PE[[2]][[i]][1],1]
  # errPE[i,2]=errPE1[zt_mt_PE[[2]][[i]][1],2]
  # errPE[i,3]=min(errPE1[,2])
  # errPE[i,4]=zt_mt_PE[[2]][[i]][1]
  # errPE[i,5]=zt_mt_PE[[2]][[i]][2]
  # errPE[i,6]=errPE1[6,2]
  # errPE[i,7]=abs(errPE[i,4]- errPE[i,5])
  
}

# plot(errPE[,1],errPE[,3], col=errPE[,5])
df=as.data.frame(cbind(errPE,errB))

return(df)
}

# pl= ggplot(data = df,aes(x=V1,y=V2))+ geom_point(size=3, alpha=1)
# pl= ggplot(data = df,aes(x=V1,y=V2,color =V4))+ geom_point(size=3, alpha=1)
# # pl= ggplot(data = df,aes(x=V6,y=V14,color =V4))+ geom_point(size=3, alpha=1)
# 
# # histogrammes
# pl=ggplot(df, aes(x=V10)) +geom_histogram(color="darkblue", fill="lightblue", bins=15)  + theme(text = element_text(size=20))




  
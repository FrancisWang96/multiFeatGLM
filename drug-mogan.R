library(rcdk)
mols<-load.molecules(c("/home/gulab-01/桌面/Structure2D_CID_31703.sdf","/home/gulab-01/桌面/Structure2D_CID_2907.sdf","/home/gulab-01/桌面/Structure2D_CID_36314.sdf","/home/gulab-01/桌面/Structure2D_CID_41867.sdf","/home/gulab-01/桌面/Structure2D_CID_84691.sdf"))
fps<-lapply(mols,get.fingerprint,type="extended",size=256,fp.mode='bit')
fp.sim<-fingerprint::fp.sim.matrix(fps,method="tanimoto")
fp.dist<-1-fp.sim
cls<-hclust(as.dist(fp.dist))

num<-fps[[1]]@bits
dt<-data.frame(mol='fpt',fp=1:256)
dt$fp=0
dt$mol<-paste0(dt$mol,1:256)
dt[num,'fp']<-1
rownames(dt)<-dt$mol
dt$mol<-NULL
names(dt)[1]<-paste0('molecule',i)
yaowu1<-rep(1,time=45) %*%t(dt)

num<-fps[[2]]@bits
dt<-data.frame(mol='fpt',fp=1:256)
dt$fp=0
dt$mol<-paste0(dt$mol,1:256)
dt[num,'fp']<-1
rownames(dt)<-dt$mol
dt$mol<-NULL
names(dt)[1]<-paste0('molecule',i)
yaowu2<-rep(1,time=42) %*%t(dt)

num<-fps[[3]]@bits
dt<-data.frame(mol='fpt',fp=1:256)
dt$fp=0
dt$mol<-paste0(dt$mol,1:256)
dt[num,'fp']<-1
rownames(dt)<-dt$mol
dt$mol<-NULL
names(dt)[1]<-paste0('molecule',i)
yaowu3<-rep(1,time=46) %*%t(dt)

num<-fps[[4]]@bits
dt<-data.frame(mol='fpt',fp=1:256)
dt$fp=0
dt$mol<-paste0(dt$mol,1:256)
dt[num,'fp']<-1
rownames(dt)<-dt$mol
dt$mol<-NULL
names(dt)[1]<-paste0('molecule',i)
yaowu4<-rep(1,time=41) %*%t(dt)

num<-fps[[5]]@bits
dt<-data.frame(mol='fpt',fp=1:256)
dt$fp=0
dt$mol<-paste0(dt$mol,1:256)
dt[num,'fp']<-1
rownames(dt)<-dt$mol
dt$mol<-NULL
names(dt)[1]<-paste0('molecule',i)
yaowu5<-rep(1,time=46) %*%t(dt)

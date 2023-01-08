library(stats)
library(h2o)
h2o.init()
combine.feature1<-read.csv('/home/gulab-01/×ÀÃæ/mrmr600_combine_doc.csv',header=T)
combine.feature2<-read.csv('/home/d/GDSCH2O/mrmr600_combine_cyc.csv',header=T)
combine.feature3<-read.csv('/home/d/mrmr600_combine_pac.csv',header=T)
combine.feature4<-read.csv('/home/d/mrmr600_combine_epi.csv',header=T)
combine.feature5<-read.csv('/home/d/mrmr600_combine_cis.csv',header=T)
print(combine.feature2$train_label)
#dim(combine.feature1)
combine.feature1<-cbind(rna.data,cnv.data,mut.data,yaowu1,train_label)
combine.feature2<-cbind(rna.data,cnv.data,mut.data,yaowu2,train_label)
combine.feature3<-cbind(rna.data,cnv.data,mut.data,yaowu3,train_label)
combine.feature4<-cbind(rna.data,cnv.data,mut.data,yaowu4,train_label)
combine.feature5<-cbind(rna.data,cnv.data,mut.data,yaowu5,train_label)
#combine.feature1$train_label
combine.feature1<-cbind(combine.feature1,train_label)
combine.feature2<-cbind(combine.feature2,train_label)
combine.feature3<-cbind(combine.feature3,train_label)
combine.feature4<-cbind(combine.feature4,train_label)
combine.feature5<-cbind(combine.feature5,train_label)
combine.feature<-rbind.fill(combine.feature1,combine.feature2,combine.feature3,combine.feature4,combine.feature5)
combine.feature[is.na(combine.feature)]=0

combine.feature<-read.csv('/home/gulab-01/×ÀÃæ/mrmr500_combine_cis.csv',header=T)
yaowu<-rbind(yaowu4,yaowu1,yaowu2,yaowu5,yaowu3)
combine.feature<-cbind(combine.feature,yaowu)

fea<-combine.feature
fea$train_label<-cellLine.class$IC50
data<-fea[sample(nrow(fea)),]
featuress<-colnames(fea)
rtrain<-0
rsquared<-0
rtrmse<-0
rtemse<-0
#xulie<-c(1:4)
pre<-data.frame()
#pre<-pre[-1,0]
folds<-cut(seq(1,nrow(data)),breaks = 10,labels = FALSE)
for(i in 1:10){
  #h2o.init()
  testIndexes<-which(folds==i,arr.ind = TRUE)
  testdata<-data[testIndexes,]
  #testx<-testdata[,-"train_label"]
  traindata<-data[-testIndexes,]
  #trainx<-traindata[,-"train_label"]
  print(i)
  #label<-traindata$label
  traindata<-as.h2o(traindata)
  #h2model=h2o.loadModel(path="/home/d/mrmr500_com_yaowu_xgb_ga_h2o/DeepLearning_model_R_1655714017812_10")
  h2model=h2o.deeplearning(x=featuress,y="train_label",training_frame =traindata,activation="TanhWithDropout",hidden=c(1000,800,500,100),l1=0.1,l2=0.1,epochs=10,rate=0.000005,distribution ="gamma",stopping_metric="RMSE")
  
  rtrain<-rtrain+h2o.r2(h2model)
  rtrmse<-rtrmse+h2o.rmse(h2model)
  #print(rtrain)
  summary(h2model)
  testdata<-as.h2o(testdata)
  pre_mu<-h2o.predict(h2model,newdata=testdata)
  p2<-h2o.performance(h2model,newdata = testdata)
  pre_sig<-p2@metrics[["mean_residual_deviance"]]
  pre_mu<-as.matrix(pre_mu)
  #pre_sig<-as.matrix(pre_sig)
  #pre_sig<-predict(model,what='sigma',type='response',newdata=testdata)
  #pre_q<-predict(gammodel,what='nu',type='response',newdata=testdata)
  ig_range<-gammaqujian(0.975,pre_mu,pre_sig)
  #pre_mu<-ig_range$Mu
  #ig_range<-expqujian(0.975,pre_mu)
  #Q.stats(resid = resid(gammodel),xvar=testdata[,-"train_label"],n.inter=10)
  testdata<-as.matrix(testdata) 
  pre1=cbind(pre_mu,ig_range$Down,ig_range$Up,testdata[,"train_label"],pre_sig)
  pre<-rbind(pre,pre1)
  #testdata<-as.matrix(testdata)
  da<-corr.test(pre_mu,testdata[,"train_label"])
  train_label<-as.data.frame(testdata[,"train_label"])
  pre_mu<-as.data.frame(pre_mu)
  rtemse<-rtemse+rmse(pre_mu,train_label)
  rsquared<-(da$r)^2+rsquared
  print((da$r)^2)
  #print(rsquared)
  print(da$p)
}

names(pre)[1:5]<-c('pre_mu','V3','V2','V4','pre_sig')
cp<-cpi(pre,220)
cp
da<-corr.test(pre$pre_mu,pre$V4)
rsquared<-(da$r)^2
rsquared
pre_mu<-log(pre$pre_mu)
log_y<-log(pre$V4)
rmse<-rmse(pre_mu,log_y)
rmse
h2o.saveModel(h2model,path="/home/gulab-01/æ¡Œé¢/mrmr1500_com_ga_all_h2o",force=TRUE)
h2o.clearLog()
h2o.shutdown()
#ggsave(p1,filename="mrmr800_com_cyc_ga_h2o_pre.png",width=12,height=9)
write.csv (pre,file="//home/gulab-01/æ¡Œé¢/mrmr1500_com_all_h2o_pre.csv",row.names = F)

library(randomForest)
library(xgboost)
library(Boruta)
library(Matrix)
library(psych)
library(gamlss)
library(plyr)

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
dim(combine.feature)

combine.feature[,"train_label"]<-log10(combine.feature[,"train_label"])
#combine.feature<-na.omit(combine.feature)
#df<-sparse.model.matrix(train_label~.-1,data=combine.feature)
df<-combine.feature[,-which(names(combine.feature)%in%c("train_label"))]

#df<-rna.data
#df<-cis_ccle
train_label<-combine.feature[,"train_label"]
df<-data.matrix(df)
#x<-df
df<-Matrix(df,sparse = T)
#train<-list(data=dfx,label=df[,"train_label"])
xgbx<-xgb.DMatrix(data=df,label=train_label)
#xgbx<-xgb.DMatrix(data=df,label=train_label)
xgb<-xgboost(data=xgbx,max_depth=6,eta=0.4,nrounds =200)
i<-xgb.importance(model=xgb)
newdata<-combine.feature[,i$Feature]
train_label<-combine.feature$train_label
newdata<-cbind(newdata,train_label)

bor<-Boruta(train_label~.,data=newdata,pvalue=0.01,doTrace=2,maxRuns=500,getImp=getImpLegacyRfZ)
z<-getSelectedAttributes(bor)
#z[-1]

fea<-newdata[,z[-1]]
#fea<-subset(fea,select = -c(PPP1R2P6_rna))
l<-newdata['train_label']
fea<-cbind(fea,l)
#fea<-combine.feature
#fea$train_label<-cellLine.class$IC50
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
  testIndexes<-which(folds==i,arr.ind = TRUE)
  testdata<-data[testIndexes,]
  #testx<-testdata[,-"train_label"]
  traindata<-data[-testIndexes,]
  #trainx<-traindata[,-"train_label"]
  print(i)
  #label<-traindata$label
  gammodel=gamlss(train_label~.,data=traindata,family =LO)
  #print(Rsq(gammodeGG
  rtrain<-rtrain+Rsq(gammodel)
  #plot(gammodel)
  summary(gammodel)
  
  pre_mu<-predict(gammodel,newdata=testdata,type="response")
  pre_sig<-predict(gammodel,what='sigma',type='response',newdata=testdata)
  #pre_nu<-predict(gammodel,what='nu',type='response',newdata=testdata)
  #pre_q<-predict(gammodel,what='nu',type='response',newdata=testdata)
  ig_range<-loqujian(0.975,pre_mu,pre_sig)
  #pre_mu<-ig_range$Mu
  #ig_range<-expqujian(0.975,pre_mu)
  #Q.stats(resid = resid(gammodel),xvar=testdata[,-"train_label"],n.inter=10)
  pre1=cbind(pre_mu,ig_range$Down,ig_range$Up,testdata$train_label)
  pre<-rbind(pre,pre1)
  da<-corr.test(pre_mu,testdata$train_label)
  rsquared<-(da$r)^2+rsquared
  print((da$r)^2)
  #print(rsquared)
  print(da$p)
}
#cp<-cpi(pre,46)
#cp
da<-corr.test(pre$pre_mu,pre$V4)
rsquared<-(da$r)^2
rsquared
rmse<-rmse(pre$pre_mu,pre$V4)
rmse
write.csv (pre,"/home/gulab-01/×ÀÃæ/mrmr500_com_cis_h2o_pre.csv",row.names = F)

##################AIC test####################

mEXP<-gamlss(train_label~.,data=fea,family = EXP)
mGA<-gamlss(train_label~.,data=fea,family = GA)
mIG<-gamlss(train_label~.,data=fea,family = IG)
mLOGNO<-gamlss(train_label~.,data=fea,family = LOGNO)
mLOGNO2<-gamlss(train_label~.,data=fea,family = LOGNO2)
mGG<-gamlss(train_label~.,data=fea,family = GG)
GAIC(mEXP,mGA,mIG,mLOGNO,mLOGNO2,mGG,k=2)

mNO<-gamlss(train_label~.,data=fea,family = NO)
mGU<-gamlss(train_label~.,data=fea,family = GU)
mLO<-gamlss(train_label~.,data=fea,family = LO)
GAIC(mNO,mGU,mLO,k=2)
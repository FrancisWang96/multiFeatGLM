library(mRMRe)
data<-read.csv("/home/gulab-01/桌面/com_cis.csv",header = T)
data2<-read.csv('/home/base/com_cyc.csv',header=T)
data3<-read.csv('/home/base/com_pac.csv',header=T)
data4<-read.csv('/home/base/com_epi.csv',header=T)
data5<-read.csv('/home/base/com_cis.csv',header=T)
#dim(combine.feature1)
#data1<-cbind(data1,yaowu1)
#data2<-cbind(data2,yaowu2)
#data3<-cbind(data3,yaowu3)
#data4<-cbind(data4,yaowu4)
#data5<-cbind(data5,yaowu5)
d<-data1$train_label
is.na(d)
data<-rbind.fill(data1,data2,data3,data4,data5)
tr<-data$train_label
is.na(tr)
data[is.na(data)]=0
dim(data)
feature_num=ncol(data)-1
train_feature<-data[,0:feature_num]
train_label=data[,ncol(data)]

#train_feature<-cbind(rna.data,cnv.data,mut.data)
mrmr_feature<-train_feature
mrmr_feature$y<-train_label

#筛选的数据要加上Y值
target_indices = which(names(mrmr_feature)=='y')

for (m in which(sapply(mrmr_feature, class)!="numeric")){
  mrmr_feature[,m]=as.numeric(unlist(mrmr_feature[,m]))
} 

#转化成mRMR.data的形式才能被使用
Data <- mRMR.data(data = data.frame(mrmr_feature))
#data就是数据，target_indices就是Y（label）值，也就是需要对比相关度的列
#feature_count设置选择特征数
mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, feature_count = 1000, solution_count = 1)
#获得mrmr选择后的特征索引
#获取筛选出来的特征的列，包含在mrmr@filters中，mrmr@filters[原特征个数]这个list里
index=mrmr@filters[[as.character(mrmr@target_indices)]]
#index
#新数据提取
new_data <- nrow(data):ncol(index)

new_data = data[,index]
dim(new_data)
new_data_0 = cbind(new_data,train_label)

#写入新数据  自行更改最后保存路径
write.csv (new_data_0,"/home/d/mrmr600_combine_cis.csv",row.names = F)

library(Boruta)
library(caret)
library(randomForest)
library(ROCR)
library(dplyr)
library(h2o)
library(openxlsx)
h2o.init()

######  training data prepare   ######
dataSource <- "GDSCbrca"
impRate <- 0.5
correlationRate <- 0.8
seeds <- 777
############   cellLine data  ############
cellLine.class <- read.xlsx("/home/gulab-01/桌面/Cisplatin_new.xlsx",sheet = 1)
cellLine.class <- cellLine.class[order(rownames(cellLine.class),decreasing = F),]

############## expression data ###########
rna.data <- read.csv("/home/gulab-01/桌面/GDSC_rna_cis.csv",row.names = 1,check=F)
colnames(rna.data) <- paste(colnames(rna.data),"rna",sep = "_")
colnames(rna.data) <- colnames(rna.data)
rna.data <- rna.data[order(rownames(rna.data),decreasing = F),]

features <- colnames(rna.data)
train.data.h2o <- as.h2o(rna.data)
anomaly_model <- h2o.deeplearning(x = features, training_frame = train.data.h2o, 
                                  activation = "Tanh",   autoencoder = TRUE,   hidden = 100,   
                                  epochs = 100,variable_importances=TRUE,seed=seeds,standardize=TRUE)
rna.import.matrix <- anomaly_model@model$variable_importances
import.feature <- rna.import.matrix$variable[1:(nrow(rna.import.matrix)*impRate)]
rna.data <- rna.data[,which(colnames(rna.data) %in% import.feature)]

omic.name <- "rna"
variable.import <- anomaly_model@model$variable_importances
out.variable.file <- paste(dataSource,"_",omic.name,"_autoVarImport.txt",sep = "")
write.table(variable.import,file = out.variable.file,row.names = F,col.names = T,sep = "\t",quote = F)

# let's take the third hidden layer
train_features <- as.data.frame(h2o.deepfeatures(anomaly_model, train.data.h2o, layer = 1)) 
feature.file.name <- paste(dataSource,"_",omic.name,"_autoLayer1.txt",sep = "")
write.table(train_features,file = feature.file.name,row.names = F,col.names = T,sep = "\t",quote = F)



###############  cnv  data   ##############
cnv.data <- read.csv("/home/d/GDSCH2O/cnv_neww_cyc.csv",row.names = 1,check=F)
colnames(cnv.data) <- paste(colnames(cnv.data),"cnv",sep = "_")
cnv.data <- cnv.data[order(rownames(cnv.data),decreasing = F),]
features <- colnames(cnv.data)
train.data.h2o <- as.h2o(cnv.data)
anomaly_model <- h2o.deeplearning(x = features,   training_frame = train.data.h2o, 
                                  activation = "Tanh",   autoencoder = TRUE,   hidden = 100,   
                                  epochs = 100,variable_importances=TRUE,seed=seeds,standardize=TRUE)
cnv.import.matrix <- anomaly_model@model$variable_importances
import.feature <- cnv.import.matrix$variable[1:(nrow(cnv.import.matrix)*impRate)]
cnv.data <- cnv.data[,which(colnames(cnv.data) %in% import.feature)]

omic.name <- "cnv"
variable.import <- anomaly_model@model$variable_importances
out.variable.file <- paste(dataSource,"_",omic.name,"_autoVarImport.txt",sep = "")
write.table(variable.import,file = out.variable.file,row.names = F,col.names = T,sep = "\t",quote = F)

# let's take the third hidden layer
train_features <- as.data.frame(h2o.deepfeatures(anomaly_model, train.data.h2o, layer = 1)) 
feature.file.name <- paste(dataSource,"_",omic.name,"_autoLayer1.txt",sep = "")
write.table(train_features,file = feature.file.name,row.names = F,col.names = T,sep = "\t",quote = F)


##############   mutation   data   ##################
mut.data <- read.csv("/home/d/GDSCH2O/mut_new_cyc.csv",row.names = 1,check=F)
colnames(mut.data) <- paste(colnames(mut.data),"mut",sep = "_")
mut.data <- mut.data[order(rownames(mut.data),decreasing = F),]

#############  correlation detection   ########
combine.feature <- cbind(rna.data,cnv.data)
cor.matrix <- cor(combine.feature)
highlyCorrelated <- findCorrelation(cor.matrix, cutoff=correlationRate)

if(length(highlyCorrelated) > 0){
  combine.feature <- combine.feature[,-highlyCorrelated]
}
combine.feature <- cbind(combine.feature,mut.data)
combine.feature <- combine.feature[order(rownames(combine.feature),decreasing = F),]
combine.feature$train_label<-cellLine.class[,2]
write.csv (combine.feature,"/home/d/com_cyc.csv",row.names = FALSE)

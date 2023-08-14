library(glmnet)
library(Matrix)
library(foreach)
library(pROC)
library(PRROC)

left_anchor_data <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\Example_left_anchor_data.csv",header = T,row.names = 1)
right_anchor_data <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\Example_right_anchor_data.csv",header = T,row.names = 1)
label <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\Example_label.csv",header = T,row.names = 1)
ppi <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\TF_PPI.csv",header = F)
cross_tissue_expression <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\TF_cross_tissue_expression.csv",header = T,row.names = 1)
TF_expression <- read.csv("D:\\个人\\Project2\\data\\GitHub\\data\\TF_gene_expression.csv",header = T,row.names = 1)
Label <- label$Label

#Set cutoff 
cutoff <- 0.8
#Construct the adjacency matrix of TF PPI network
TF <- colnames(cross_tissue_expression)
TF_number <- ncol(cross_tissue_expression)
ppi <- as.matrix(ppi)
ppi_matrix <- matrix(0,TF_number,TF_number)
rownames(ppi_matrix) <- TF
colnames(ppi_matrix) <- TF
for(i in 1:TF_number){
  match_index1 <- which(TF[i]==ppi[,1])
  match_index_vector1 <- ppi[match_index1,2]
  match_index2 <- which(TF[i]==ppi[,2])
  match_index_vector2 <- ppi[match_index2,1]
  match_index_vector <- c(match_index_vector1,match_index_vector2)
  match_index_final <- c()
  for(j in match_index_vector){
    match_index3 <- which(j==TF)
    match_index_final <- c(match_index_final,match_index3)
  }
  ppi_matrix[match_index_final,i] <- 1
}
#Generating TF cross cell type co-expression matrix
co_expression <- cor(cross_tissue_expression)
co_expression[is.na(co_expression)] <- 0
#Combine PPI and co-expression network 
PPIandCoExPression <- ppi_matrix* co_expression
PPIandCoExPression[which(PPIandCoExPression <= cutoff)]<-0
contral_expression <- which(TF_expression[,1]<1)
PPIandCoExPression[,contral_expression]<-0
PPIandCoExPression[contral_expression,]<-0
#Generate genomic feature TF complex activity
TF_pair <- c()
n_row <- nrow(PPIandCoExPression)
for (i in 1:nrow(PPIandCoExPression)) {
  index_pair <- which(PPIandCoExPression[i:n_row,i]>0)
  if (length(index_pair)>0){
    w_vec <- PPIandCoExPression[index_pair+(i-1),i]
    index_vec <- cbind(i,index_pair+(i-1),w_vec)
    TF_pair<-rbind(TF_pair,index_vec)
  }
}
combin_matrix <- t(TF_pair)
TF_complex_data <- matrix(NA,nrow(left_anchor_data),ncol(combin_matrix))
TF_complex_name <- c()
TF_name <- colnames(cross_tissue_expression)
for(i in 1:ncol(combin_matrix)){
  feature_vec <- combin_matrix[,i]
  data_sub = (left_anchor_data[,feature_vec[1]]*right_anchor_data[,feature_vec[2]]+left_anchor_data[,feature_vec[2]]*right_anchor_data[,feature_vec[1]])
  sub_name = TF_name[feature_vec]
  merge_name<- paste(sub_name[1],"_",sub_name[2],sep = "")
  TF_complex_data[,i] <- data_sub
  TF_complex_name[i] <- merge_name
}
colnames(TF_complex_data)<-TF_complex_name
tfc_data <- as.data.frame(TF_complex_data)
#Extract the features corresponding to the candidate TF network
TF_weight <- colSums(PPIandCoExPression)
candidateTF_index <- which(TF_weight > 0)
candidateTF <- TF[candidateTF_index]
candidateTF_weight <- TF_weight[candidateTF_index]
left_anchor_data_candidate <- data.frame(left_anchor_data[,candidateTF_index],left_anchor_data[,(ncol(left_anchor_data)-3):ncol(left_anchor_data)])
right_anchor_data_candidate <- data.frame(right_anchor_data[,candidateTF_index],right_anchor_data[,(ncol(right_anchor_data)-3):ncol(right_anchor_data)])
non_TF_weight <- rep(0,4)
non_TF_feature <- c("openness","CTCF","SMC3","RAD21")
anchor_weight <- c(candidateTF_weight,non_TF_weight)
TF_complex_weight <- TF_pair[,3]
model_input_data <- cbind(left_anchor_data_candidate,right_anchor_data_candidate,tfc_data,Label)
feature_weight <- c(anchor_weight,anchor_weight,TF_complex_weight)

#SpecLoop Model
set.seed(2)
require(caret)
folds <- createFolds(y=model_input_data$Label,k=5)
fold_predict <- c()
true_value <- c()
feature_data <- ncol(model_input_data) -1
for(i in 1:5){
  print(i)
  fold_test <- model_input_data[folds[[i]],]
  fold_train <- model_input_data[-folds[[i]],]
  y1<-fold_test$Label
  y2<-fold_train$Label
  x1<-fold_test[,1:feature_data]
  x2<-fold_train[,1:feature_data]
  x1<-as.matrix(x1)
  x2<-as.matrix(x2)
  fold_pre <-cv.glmnet(x2,y2,family = "binomial",penalty.factor=feature_weight)
  fold_predict1 <- predict(fold_pre,type='response',x1)
  fold_predict1<- as.numeric(fold_predict1)
  true_value1 <- y1
  fold_predict <- append(fold_predict,fold_predict1)
  true_value <- append(true_value,true_value1)
}
result_model <- data.frame(true_value,fold_predict)
#ROC curve
modelroc <- roc(result_model$true_value,result_model$fold_predict)
auroc<-auc(modelroc)
auroc<-round(auroc,4)
par(pty="s")
plot(modelroc, col="red",main = "SpecLoop_ROCcurve",auc.polygon=TRUE,auc.polygon.col="white",legacy.axes=TRUE,
     grid.col=c("black", "black"),xaxs="i",yaxs="i",font.lab=2,cex.main=2,cex.lab=2,lwd=3)
auc_name<-c("AUC=")
auroc<-c(auroc)
auroc.<-paste(auc_name,auroc)
legend("bottomright",legend = auroc.,col=c("red"),lwd = 2,cex=1.8)
#PR curve
fg = result_model$fold_predict[which(result_model$true_value==1)]
bg = result_model$fold_predict[which(result_model$true_value==0)]
pr <-pr.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
aupr <- round(pr$auc.integral,4)
par(pty="s")
plot(pr, col="blue",font.lab=2,cex.main=2,cex.lab=2,main ="SpecLoop_PRcurve" ,auc.main = F,xaxs="i",yaxs="i")
aupr1 <-c(aupr)
aupr.<-paste(auc_name,aupr)
legend("bottomleft",legend = aupr.,col=c("blue"),lwd = 2,cex=1.8)
fit<-fold_pre
coefficients <- coef(fit,s=fit$lambda.min)
feature_coefficents <- coefficients[,1]
Active.Index<-which(feature_coefficents!=0)
output_feature <- names(Active.Index[-1])
output_feature <- sub("_(la|ra)$", "", output_feature)
output_feature <- unique(output_feature)
SpecLoop_output <- setdiff(output_feature,non_TF_feature)

print(SpecLoop_output)

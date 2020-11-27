library(data.table)
library(h2o)
library(dplyr)
library(ROCR)
library(ggplot2)

#Random Forest Code
#LJC is the generated file from community assignment topic model weights and other metadata.
dim(LJC)
results <- list()
ImportanceA<-list()
perfs<-perf<-list()

for (i in 1:1000) {
  #MakeTrainingSetVSTestSet
  sid<-sample(nrow(LJC),0.7*nrow(LJC))
  train<-LJC[sid,]
  test<- LJC[-sid,]

  #test<-dplyr::sample_n(train, 10)
  #No. of rows and columns in train
  dim(train)
  str(train)

  #No. of rows and columns in test
  dim(test)
  str(train)

  #Load onto H20
  localH2O <- h2o.init(nthreads = -1)
  h2o.init()
  #data to h2o cluster
  train.h2o <- as.h2o(train)
  test.h2o <- as.h2o(test)
  #check column index number
  colnames(train.h2o)

  #Dependant Variable(Labor_Fever)
  y.dep <- 5

  #Independant Varibale
  x.indep <-c(1:4,6:21)

  #Do random Forest
  system.time(
    rforest.model <- h2o.randomForest(y=y.dep, x=x.indep, training_frame = train.h2o,
                                      ntrees = 1000, mtries = 3, max_depth = 20,
                                      seed = 1122))

  h2o.performance(rforest.model)
  #check variable importance
  Importance<- h2o.varimp(rforest.model)
  #WriteTable
  write.csv(Importance, file = "/RF_PredictionOfFebrileStatus.csv", row.names = FALSE)

  #making predictions on unseen data
  rf_test<-system.time(
    predict.rforest<- as.data.frame(h2o.predict(rforest.model, test.h2o)))
  Labor_Fever_rf <- data.frame( orig_cat=test$Labor_Fever, predicted_impact=predict.rforest$predict)

  h2o.performance(rforest.model,test.h2o)

  write.csv(Labor_Fever_rf, file = "/RF_PredictionOfFebrileStatus_test.csv", row.names = FALSE)


  # without the loop (Febrile)
  q <- data.frame(
    predicted = as.integer(predict.rforest$predict == "Febrile"),
    labels = as.integer(test$Labor_Fever == "Febrile"),
    p = predict.rforest[, "Febrile"]
  )
  pred <- prediction(q$p, q$labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  plot(perf)
  #Afebrile
  q <- data.frame(
    predicted = as.integer(predict.rforest$predict == "Afebrile"),
    labels = as.integer(test$Labor_Fever == "Afebrile"),
    p = predict.rforest[, "Afebrile"]
  )
  pred <- prediction(q$p, q$labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  plot(perf)
  results[[i]] <- auc
  ImportanceA[[i]]<-Importance
  perfs[[i]]<-perf
}
###
boxplot(`2`~ Labor_Fever,data=LJC)
boxplot(`1`~ Labor_Fever,data=LJC)
boxplot(MAGE~ Labor_Fever,data=LJC)
###PlotAUCIntervak

AUC<-unlist(results, recursive = FALSE)
summary(AUC)
plot(density(AUC),main="AUC of 1000 RF Models for Labor Prediction")


abline(v = median(AUC), col="red", lwd=2, lty=1)
abline(v = max(AUC), col="blue", lwd=2, lty=2)
abline(v = min(AUC), col="darkgreen", lwd=2, lty=2)

#
localH2O <- h2o.init(nthreads = -1)
h2o.init()
results <- list()
ImportanceA<-list()
perfs<-list()

for (i in 1:1000) {
  #MakeTrainingSetVSTestSet
  sid<-sample(nrow(LJC),0.7*nrow(LJC))
  train<-LJC[sid,]
  test<- LJC[-sid,]

  #test<-dplyr::sample_n(train, 10)
  #No. of rows and columns in train
  dim(train)
  str(train)

  #No. of rows and columns in test
  dim(test)
  str(train)

  #Load onto H20
  localH2O <- h2o.init(nthreads = -1)
  h2o.init()
  #data to h2o cluster
  train.h2o <- as.h2o(train)
  test.h2o <- as.h2o(test)
  #check column index number
  colnames(train.h2o)

  #Dependant Variable(Labor_Fever)
  y.dep <- 5

  #Independant Varibale
  x.indep <-c(1:4,6:21)

  #Do random Forest
  system.time(
    rforest.model <- h2o.randomForest(y=y.dep, x=x.indep, training_frame = train.h2o,
                                      ntrees = 1000, mtries = 3, max_depth = 20,
                                      seed = 1122))

  h2o.performance(rforest.model)
  #check variable importance
  Importance<- h2o.varimp(rforest.model)
  #WriteTable
  write.csv(Importance, file = "/RF_PredictionOfFebrileStatus.csv", row.names = FALSE)

  #making predictions on unseen data
  rf_test<-system.time(
    predict.rforest<- as.data.frame(h2o.predict(rforest.model, test.h2o)))
  Labor_Fever_rf <- data.frame( orig_cat=test$Labor_Fever, predicted_impact=predict.rforest$predict)

  h2o.performance(rforest.model,test.h2o)

  write.csv(Labor_Fever_rf, file = "/RF_PredictionOfFebrileStatus_test.csv", row.names = FALSE)

  # without the loop (Febrile)
  q <- data.frame(
    predicted = as.integer(predict.rforest$predict == "Febrile"),
    labels = as.integer(test$Labor_Fever == "Febrile"),
    p = predict.rforest[, "Febrile"]
  )
  pred <- prediction(q$p, q$labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  plot(perf)
  #Afebrile
  q <- data.frame(
    predicted = as.integer(predict.rforest$predict == "Afebrile"),
    labels = as.integer(test$Labor_Fever == "Afebrile"),
    p = predict.rforest[, "Afebrile"]
  )
  pred <- prediction(q$p, q$labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")@y.values[[1]]
  plot(perf)
  results[[i]] <- auc
  ImportanceA[[i]]<-Importance
  perfs[[i]]<-perf
}

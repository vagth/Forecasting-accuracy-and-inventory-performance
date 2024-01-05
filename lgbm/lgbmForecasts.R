## Path ----
setwd("/home/lgbm")

## Libraries ----
library(lubridate)
library(data.table)
library(stringr)
library(lightgbm)
library(methods)
library(Matrix)
library(rpart)
library(tsutils)
library(forecast)
library(dplyr)                                                 
library(plyr)                                                   
library(readr) 

## Read data ----
calendar <- head(fread("calendar.csv", stringsAsFactors = F),1941)
M5 <- fread("sales_train_evaluation.csv", stringsAsFactors = F)
M5$store_id = M5$dept_id = M5$cat_id <- NULL
M5<-aggregate(.~item_id+state_id, M5, sum)
save.image("rolling data/Dataset.Rdata")

load("rolling data/Dataset.Rdata")

## Engineer Features ----
include_id = exclude_id <- c()
time_fe <- c()
for (Sid in unique(M5$state_id)){
  
  t0 <- Sys.time()
  sub_set <- M5[M5$state_id==Sid,]
  row.names(sub_set) <- NULL
  dataset_for_GBM_model <- NULL
  
  for (tsid in 1:nrow(sub_set)){
    
    #Sales
    sales <- as.numeric(sub_set[tsid,3:ncol(sub_set)])
    
    #Calendar features
    dates <- as.Date(calendar$date)
    day_id <- as.factor(day(dates)) 
    month_id <- as.factor(month(dates))
    weekday_id <- as.factor(wday(dates)) 
    week_id <- as.factor(week(dates))
    
    #SNAP and special days
    if (sub_set$state_id[1] == "CA"){
      snap <- as.factor(calendar$snap_CA)
    }else if (sub_set$state_id[1] == "TX"){
      snap <- as.factor(calendar$snap_TX)
    }else{
      snap <- as.factor(calendar$snap_WI)
    }
    event_name <- calendar$event_name_1
    event_type <- calendar$event_type_1
    event_name[is.na(event_name)] <- "NoEvent"
    event_type[is.na(event_type)] <- "NoEvent"
    event_name <- as.factor(event_name)
    event_type <- as.factor(event_type)
    
    #This is the main data set for series tsid
    temp <- data.frame(sales, day_id, month_id, weekday_id, week_id, snap, event_name, event_type)
    temp$weekend <- 0 #Flag weekends
    temp[(temp$weekday_id==7)|(temp$weekday_id==1),]$weekend <- 1
    
    #Include product_id
    temp$prod_id <- sub_set$item_id[tsid]
    temp$state_id <- sub_set$state_id[tsid]
    
    #FLag day id
    temp$origin <- c(1:nrow(temp))
    
    #Drop days for which sales has not started yet
    start_p <- as.numeric(row.names(temp[temp$sales>0,]))[1]
    temp <- temp[start_p:nrow(temp),]
    
    temp$lag7 = temp$lag14 = temp$lag28 <- NA 
    temp$M_7 = temp$M_14 = temp$M_30 = temp$M_60 = temp$M_180 <- NA 
    temp$S_7 = temp$S_14 = temp$S_30 = temp$S_60 = temp$S_180 <- NA 
    if (nrow(temp)>=730){
      
      for (i in 210:nrow(temp)){
        temp$lag7[i] <- temp$sales[i-30]
        temp$lag14[i] <- temp$sales[i-(30+7)]
        temp$lag28[i] <- temp$sales[i-(30+7+7)]
        
        temp$M_7[i] <- mean(temp$sales[(i-30):(i-30-7+1)])
        temp$S_7[i] <- sd(temp$sales[(i-30):(i-30-7+1)])
        
        temp$M_14[i] <- mean(temp$sales[(i-30):(i-30-14+1)])
        temp$S_14[i] <- sd(temp$sales[(i-30):(i-30-14+1)])
        
        temp$M_30[i] <- mean(temp$sales[(i-30):(i-30-30+1)])
        temp$S_30[i] <- sd(temp$sales[(i-30):(i-30-30+1)])
        
        temp$M_60[i] <- mean(temp$sales[(i-30):(i-30-60+1)])
        temp$S_60[i] <- sd(temp$sales[(i-30):(i-30-60+1)])
        
        temp$M_180[i] <- mean(temp$sales[(i-30):(i-30-180+1)])
        temp$S_180[i] <- sd(temp$sales[(i-30):(i-30-180+1)])
      }
      temp <- na.omit(temp)
      
      dataset_for_GBM_model <- rbind(dataset_for_GBM_model, temp)
      include_id <- c(include_id, paste0(sub_set$item_id[tsid], "_", sub_set$state_id[tsid])) 
      
    }else{
      
      exclude_id <- c(exclude_id, paste0(sub_set$item_id[tsid], "_", sub_set$state_id[tsid]))
      
    }
    
  }
  
  t1 <- Sys.time()
  time_fe <- c(time_fe, as.numeric(difftime(t1,t0,units="secs")))
  
  write.csv(dataset_for_GBM_model, paste0("rolling data/",Sid,".csv"), row.names = F)
  rm(dataset_for_GBM_model)
  
  print(c(Sid,length(include_id)+length(exclude_id)))
  
}
rm(Sid, i, t0, t1, tsid, start_p, temp, sub_set, calendar)
rm(dates, sales, day_id, month_id, weekday_id, week_id, snap, event_name, event_type)

#Load data and create complete set
dataset_for_GBM_model_all <- list.files(path = "rolling data",     
                                        pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                                            
  bind_rows                                                        

save.image("rolling data/FE.Rdata")


## Train LightGBM models and forecast ----
origin <- 1940

while (origin>=(1941-393-21)) {
  
  load("rolling data/FE.Rdata")
  
  rm(M5,exclude_id,include_id)
  
  print(paste("Origin",origin))
  print(Sys.time())
  
  fh <- 30
  
  # (1) State models
  # Train model
  
  train_set <- dataset_for_GBM_model_all[(dataset_for_GBM_model_all$origin<=origin)&(dataset_for_GBM_model_all$origin>=(origin-365*2.25)),]
  test_set <- dataset_for_GBM_model_all[(dataset_for_GBM_model_all$origin>origin)&(dataset_for_GBM_model_all$origin<=(origin+fh)),]
  complete_set <- rbind(train_set, test_set)
  complete_set$origin <- NULL
  
  x_train <- head(sparse.model.matrix(sales ~ .-1, data = complete_set), nrow(train_set))
  x_test <- tail(sparse.model.matrix(sales ~ .-1, data = complete_set), nrow(test_set))
  y_train = train_set$sales
  y_test = test_set$sales
  
  rm(train_set,complete_set)
  
  bst_tweedie <- lightgbm(
    data = x_train
    , label = y_train
    , objective = "tweedie"
    , tweedie_variance_power = 1.1
    , subsample = 0.5
    , subsample_freq = 1
    , learning_rate = 0.015
    , num_leaves = 2^8-1
    , min_data_in_leaf = 2^8-1
    , feature_fraction = 0.5
    , max_bin = 100
    , n_estimators = 1000
    , verbose = -1
    , boost_from_average = FALSE
  )
  
  
  #Forecast
  pred_tweedie <- predict(bst_tweedie, x_test)
  rm(bst_tweedie)
  
  #Store predictions
  eval_frame <- data.frame(test_set$prod_id, test_set$state_id, pred_tweedie)
  colnames(eval_frame) <- c("item_id", "state_id","LightGBM_S")
  if (nrow(eval_frame)<(8613*fh)){
    eval_frame$fh <- rep(c(1:(1941-origin)),nrow(eval_frame)/(1941-origin))
  }else{
    eval_frame$fh <- rep(c(1:fh), nrow(eval_frame)/fh)
  }
  
  rm(test_set, x_train, x_test, y_train, y_test)
  rm(pred_tweedie)
  
  #Order based on horizon
  eval_frame <- eval_frame[order(eval_frame$state_id, eval_frame$item_id, eval_frame$fh),]
  row.names(eval_frame) <- NULL
  
  #Save forecasts and errors
  write.csv(eval_frame, paste0("rolling results/Forecasts_",origin,".csv"), row.names = F)
  
  origin <- origin - 1
  
  rm(list=setdiff(ls(), c("origin")))
  gc()
}

rm(list=setdiff(ls(), c("origin")))
.rs.restartR()

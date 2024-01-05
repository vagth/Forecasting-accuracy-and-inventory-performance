path<-"/home/Paper"
setwd(path)

#### Libraries ####
library(zoo)
library(ggplot2)
library(dplyr)
library(plyr)
library(readr)
library(vroom)
library(tidyr)
library(stringr)
library(lubridate)
library(forecast)
library(reshape2)
library(tsfeatures)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(readr)
library(moments)
library(tseries)
library(fracdiff)
library(Matrix)
library(MLmetrics)
library(smooth)

#### Functions ####
intervals<-function(x){
  y<-c()
  k<-1
  counter<-0
  for (tmp in (1:length(x))){
    if(x[tmp]==0){
      counter<-counter+1
    }else{
      k<-k+1
      y[k]<-counter
      counter<-1
    }
  }
  y<-y[y>0]
  y[is.na(y)]<-1
  return(y)
}
demand<-function(x){
  y<-x[x!=0]
  return(y)
}
statistics<-function(data,tsid){
  input <- data[tsid,]
  id<-input$id
  input <- as.numeric(input)
  if (max(input,na.rm = T)>0){
    input<-input[min(which(input!=0)):length(input)] #remove leading zeros
  }else{
    input<-0
  }
  
  
  lngth <- length(input)
  D <- demand(input)
  ADI <- mean(intervals(input))
  CV2 <- (sd(D)/mean(D))^2
  Min <- min(input)
  Low25 <- as.numeric(quantile(input,0.25))
  Mean <- mean(input)
  Median <- median(input)
  Up25 <- as.numeric(quantile(input,0.75))
  Max <- max(input)
  pz <- length(input[input==0])/lngth
  
  if (!is.na(CV2) & !is.na(ADI)){
    if (ADI > (4/3)){
      if (CV2 > 0.5){
        Type <- "Lumpy"
      }else{
        Type <- "Intermittent"
      }
    }else{
      if (CV2 > 0.5){
        Type <- "Erratic"
      }else{
        Type <- "Smooth"
      }
    }
  }else{
    Type <- NA
  }
  
  matrix_s <- data.frame(id, lngth, ADI, CV2, pz, Type, Min, Low25, Mean, Median, Up25, Max)
  
  return(matrix_s)
}
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

Naive <- function(x, h, type="simple"){
  frcst <- rep(tail(x,1), h)
  if (type=="seasonal"){
    frcst <- head(rep(as.numeric(tail(x,7)), h), h) 
  }
  return(frcst)
}
SexpS <- function(x, h){
  a <- optim(c(0), SES, x=x, h=1, job="train", lower = 0.1, upper = 0.3, method = "L-BFGS-B")$par
  y <- SES(a=a, x=x, h=1, job="forecast")$mean
  forecast <- rep(as.numeric(y), h)
  return(forecast)
}
SES <- function(a, x, h, job){
  y <- c()  
  y[1] <- x[1] #initialization
  
  for (t in 1:(length(x))){
    y[t+1] <- a*x[t]+(1-a)*y[t]
  }
  
  fitted <- head(y,(length(y)-1))
  forecast <- rep(tail(y,1),h)
  if (job=="train"){
    return(mean((fitted - x)^2))
  }else if (job=="fit"){
    return(fitted)
  }else{
    return(list(fitted=fitted,mean=forecast))
  }
}
Croston <- function(x, h, type){
  if (type=="classic"){
    mult <- 1 
    a1 = a2 <- 0.1 
  }else if (type=="optimized"){
    mult <- 1 
    a1 <- optim(c(0), SES, x=demand(x), h=1, job="train", lower = 0.1, upper = 0.3, method = "L-BFGS-B")$par
    a2 <- optim(c(0), SES, x=intervals(x), h=1, job="train", lower = 0.1, upper = 0.3, method = "L-BFGS-B")$par
  }else if (type=="sba"){
    mult <- 0.95
    a1 = a2 <- 0.1
  }
  yd <- SES(a=a1, x=demand(x), h=1, job="forecast")$mean
  yi <- SES(a=a2, x=intervals(x), h=1, job="forecast")$mean
  forecast <- rep(as.numeric(yd/yi), h)*mult
  return(forecast)
}
TSB <- function(x, h){
  n <- length(x)
  p <- as.numeric(x != 0)
  z <- x[x != 0]
  
  a <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.8) 
  b <- c(0.01,0.02,0.03,0.05,0.1,0.2,0.3)
  MSE <- c() ; forecast <- NULL
  for (atemp in a){
    for (btemp in b){
      zfit <- vector("numeric", length(x))
      pfit <- vector("numeric", length(x))
      zfit[1] <- z[1] ; pfit[1] <- p[1]
      
      for (i in 2:n) {
        pfit[i] <- pfit[i-1] + atemp*(p[i]-pfit[i-1])
        if (p[i] == 0) {
          zfit[i] <- zfit[i-1]
        }else {
          zfit[i] <- zfit[i-1] + btemp*(x[i]-zfit[i-1])
        }
      }
      yfit <- pfit * zfit
      forecast[length(forecast)+1] <- list(rep(yfit[n], h))
      yfit <- c(NA, head(yfit, n-1))
      MSE <- c(MSE, mean((yfit-x)^2, na.rm = T) )
    }
  }
  return(forecast[[which.min(MSE)]])
}
ADIDA <- function(x, h){
  al <- round(mean(intervals(x)),0) #mean inter-demand interval
  #Aggregated series (AS)
  AS <- as.numeric(na.omit(as.numeric(rollapply(tail(x, (length(x) %/% al)*al), al, FUN=sum, by = al))))
  forecast <- rep(SexpS(AS, 1)/al, h)
  return(forecast)
}
iMAPA <- function(x, h){
  mal <- round(mean(intervals(x)),0)
  frc <- NULL
  for (al in 1:mal){
    frc <- rbind(frc, rep(SexpS(as.numeric(na.omit(as.numeric(rollapply(tail(x, (length(x) %/% al)*al), al, FUN=sum, by = al)))), 1)/al, h))
  }
  forecast <- colMeans(frc)
  return(forecast)
}
MA <- function(x, h){
  mse <- c()
  for (k in 2:14){
    y <- rep(NA, k)
    for (i in (k+1):length(x)){
      y <- c(y, mean(x[(i-k):(i-1)]))
    }
    mse <- c(mse, mean((y-x)^2, na.rm = T))
  }
  k <- which.min(mse)+1
  forecast <- rep(mean(as.numeric(tail(x, k))), h)
  return(forecast)
}

SS<-function(SL, sd, L){
  if( any(SL <= 0 | SL >= 1) ) stop('SL not between 0 and 1')
  z<- qnorm(SL, 0, 1)
  rst <- z*sd*sqrt(L)
  rs <-   round(rst, 2)
  return(rs)
}
RS<-function(rsid,initial_stock, actual_demand, FT, L, SL, S=NULL, method){
  if( any(is.null(initial_stock) | is.null(actual_demand) | is.null(FT) | is.null(L) | is.null(SL)) ) stop('Variable missing')
  fh<-FT+L
  periods<-length(actual_demand)-1
  actual_demand<-head(actual_demand,periods)
  
  known_periods<-periods-393 #365
  known_demand<-head(actual_demand,known_periods)
  
  orders=sent=backorders=frst_diffs<-rep(0,known_periods)
  forecasts<-rep(0,known_periods+1)
  s<-rep(initial_stock,(known_periods))
  start<-TRUE
  count_ft<-0
  
  actual_demand_formetrics=pinball_frc=pinball_act<-c()
  aek_actual_demand=aek_forecasts<-NULL
  for (periodid in c((known_periods+1):periods)){
    sent<-c(sent,min(actual_demand[periodid], s[periodid-1]))
    
    if (periodid %in% round(seq(periods,1,by=-(FT)))){
      count_ft<-count_ft+1
      known_demand<-head(actual_demand,periodid)
      ss<-SS(SL, sd(known_demand), (L+FT))
      if (is.null(S) || S==0){
        
        if (method=="Croston"){
          forecast<-try((Croston(known_demand,fh,"classic")),silent = T)
        }else if (method=="optCroston"){
          forecast<-try((Croston(known_demand,fh,"optimized")),silent = T)
        }else if (method=="SBA"){
          forecast<-try((Croston(known_demand,fh,"sba")),silent = T)
        }else if (method=="TSB"){
          forecast<-try((TSB(tail(known_demand,56),fh)),silent = T)
        }else if (method=="ADIDA"){
          forecast<-try((ADIDA(known_demand,fh)),silent = T)
        }else if (method=="iMAPA"){
          forecast<-try((iMAPA(known_demand,fh)),silent = T)
        }else if (method=="SES"){
          forecast<-try((as.numeric(SexpS(known_demand,fh))),silent = T)
        }else if (method=="ETS"){
          if (start | (count_ft*FT %% 28)==0){
            model<-es(ts(tail(known_demand,56), frequency = 7), h=fh)
            model<-substr(model$model,5,7)
            forecast<-try((as.numeric(es(ts(tail(known_demand,56), frequency = 7), model = model, h=fh)$forecast)),silent = T)
            start<-FALSE
          }else{
            forecast<-try((as.numeric(es(ts(tail(known_demand,56), frequency = 7), model = model, h=fh)$forecast)),silent = T)
          }
        }else if (method=="Naive"){
          forecast<-try((Naive(known_demand,fh)),silent = T)
        }else if (method=="SNaive"){
          forecast<-try((Naive(known_demand,fh,"seasonal")),silent = T)
        }else if (method=="MA"){
          forecast<-try((MA(known_demand,fh)),silent = T)
        }else if (method=="ARIMA"){
          forecast<-try((as.numeric(forecast::forecast(auto.arima(ts(tail(known_demand,56), frequency = 7)), h=fh)$mean)),silent = T)
        }else if (method=="Comb"){
          if (start | (count_ft*FT %% 28)==0){
            model<-es(ts(tail(known_demand,56), frequency = 7), h=fh)
            model<-substr(model$model,5,7)
            forecast<-try((as.numeric(es(ts(tail(known_demand,56), frequency = 7), model = model, h=fh)$forecast)),silent = T)
            start<-FALSE
          }else{
            forecast<-try((as.numeric(es(ts(tail(known_demand,56), frequency = 7), model = model, h=fh)$forecast)),silent = T)
          }
          forecast2<-try((as.numeric(forecast::forecast(auto.arima(ts(tail(known_demand,56), frequency = 7)), h=fh)$mean)),silent = T)
          forecast<-(forecast+forecast2)/2
        }else if (method=="LGBM"){
          forecast<-try(fread(file=paste0("/home/lgbm/rolling results/Forecasts_",periodid,".csv")),silent = T)
          forecast$id<-paste0(forecast$item_id,"_",forecast$state_id)
          forecast<-forecast[forecast$id==rsid,]
          forecast<-forecast[order(forecast$fh),]
          forecast<-head(forecast$LightGBM_S,fh)
        }
        
        if (class(forecast)=="try-error"){forecast<-0}
        if (any(forecast<0)){forecast[forecast<0]<-0}
        Q<-sum(forecast)+ss
        
        pinball_frc<-c(pinball_frc, Q)
        pinball_act<-c(pinball_act, sum(actual_demand[(periodid+1):(periodid+fh)]) )
      }else{
        Q<-S
      }
      forecasts<-c(forecasts,head(forecast,fh)) #rep(forecast,FT))
      actual_demand_formetrics<-c(actual_demand_formetrics, actual_demand[(periodid+1):(periodid+fh)])
      frst_diffs<-c(frst_diffs,rep(mean(diff(known_demand)^2),FT))
      orders<-c(orders, round(max(Q-s[periodid-1]+sent[periodid]-sum(orders[max(periodid-L,1):(periodid-1)]),0)))
    }else{
      orders<-c(orders, 0)
      if (count_ft==0){
        forecasts<-c(forecasts, 0)
      }
    }
    
    
    if (periodid <= L){
      s<-c(s, s[periodid-1]-sent[periodid])
    }else{
      s<-c(s, s[periodid-1]-sent[periodid]+orders[periodid-L])
    }
  }
  
  actual_demand_formetrics<-actual_demand_formetrics[!is.na(actual_demand_formetrics)]
  forecasts<-forecasts[(known_periods+2):length(forecasts)]
  forecasts<-head(forecasts,length(actual_demand_formetrics))
  index<-c(1:(periods-365))
  
  inv<-s[-index] ; ord<-orders[-index]
  
  #Pinboll loss
  pinball_act<-as.numeric(na.omit(pinball_act))
  pinball_frc<-head(pinball_frc,length(pinball_act))
  temp_sum<-c()
  for (tpnbl in c(1:length(pinball_act))){
    if (pinball_act[tpnbl]>=pinball_frc[tpnbl]){
      temp_sum <- c(temp_sum, (pinball_act[tpnbl] - pinball_frc[tpnbl])*as.numeric(SL))
    }else{
      temp_sum <- c(temp_sum, (pinball_frc[tpnbl] - pinball_act[tpnbl])*(1-as.numeric(SL)))
    }
  }
  scaler<-mean((frst_diffs[-index]))
  PLd <- mean(temp_sum/scaler)
  PL <- mean(temp_sum)/scaler
  
  return(c(
    FT,
    L,
    SL,
    sum(sent[-index])/sum(actual_demand[-index]),
    mean(sent[-index][actual_demand[-index]>0]/actual_demand[-index][actual_demand[-index]>0]),
    mean(s[-index]),
    max(s[-index]),
    sum(sum(actual_demand[-index])-sum(sent[-index])),
    length(inv[inv==0]),
    mean(inv)/mean(actual_demand[-index]),
    sum(ord),
    length(ord[ord>0]),
    sum(actual_demand_formetrics),
    sum(forecasts),
    mean(actual_demand_formetrics-forecasts),
    mean(abs(actual_demand_formetrics-forecasts)),
    sqrt(mean((actual_demand_formetrics-forecasts)^2)),
    sqrt(mean((actual_demand_formetrics-forecasts)^2)/mean((frst_diffs[-index]))),
    
    sqrt(mean((actual_demand_formetrics-forecasts)^2)/(mean(actual_demand_formetrics)^2)),
    
    mean(abs(actual_demand_formetrics-forecasts))/mean(abs(frst_diffs[-index])),
    PLd,
    PL
  ))
}

#### Read Data ####
data<-read.csv("sales_train_evaluation.csv")
data$id<-NULL
data$store_id<-NULL
data<-aggregate(.~item_id+dept_id+cat_id+state_id,data,sum)
data$id<-paste0(data$item_id,"_",data$state_id)
data<-data[moveme(names(data),"id first")]

data$item_id=data$dept_id=data$cat_id=data$state_id<-NULL

#### Read Setups ####
stats<-read.csv("stats.txt")
stats<-stats[stats$lngth>=(393+21),]
excluded<-read.csv("excluded_ids.csv") ; excluded<-excluded$x
stats<-stats[!(stats$id %in% excluded),]

data<-data[data$id %in% unique(stats$id),]

#### Start Experiment ####
m_names<-c("LGBM","Naive","SNaive","MA","SES","Croston","TSB","SBA","ETS","ARIMA","ADIDA","iMAPA")

for (sel_method in m_names){
  results<-NULL
  for (tmp in c(1:nrow(data))){
    print(c(sel_method,tmp))
    inf<-c(1)
    infos<-data[tmp,inf]
    d_ts<-as.numeric(data[tmp,-inf])
    
    d_ts<-ts(d_ts,frequency = 7)
    res<-NULL
    
    for (t_t in c(3,7,14,21)){
      for (serl in c(0.9, 0.95, 0.99)){
        for (l in c(3, 9)){
          res<-rbind(res,data.frame(infos,sel_method,t(RS(infos,0,d_ts,t_t,l,serl,S=NULL,sel_method))))
        }
      }
    }
    
    colnames(res)<-c("id","method","FT","L","target_service_level","fill_rate","service_level","mean_inventory","max_inventory","lost_sales","out_of_stock_days",
                     "inventory_days","sum_orders","count_orders","sum_demand","sum_forecasts","ME","MAE","RMSE","RMSSE","RMSsE","FD2","PLd","PL")
    results<-rbind(results,res)
    gc()
  }
  gc()
  
  write.csv(results,paste0("rolling results/results_",sel_method,".csv"),row.names = F)
}
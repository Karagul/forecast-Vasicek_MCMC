############################
##### MCMC testing 2  ######
############################

#######################
#### Preprocessing  ###
#######################

#all packages needed
requirements <- c("lattice",'MASS','ggplot2','grid','snow','pbapply',"parallel","reshape","RColorBrewer",
                  'WriteXLS')
#install packages if you have not
for(requirement in requirements){if( !(requirement %in% installed.packages())) install.packages(requirement)}
#load all required packages
lapply(requirements, require, character.only=T);

#######################################################################################
#### for testing purpose, simulate observed short rate and zero coupon bond yields ####
#######################################################################################


### set up the directory and load data
data =read.csv('Data/Monthly Curve Data USD.csv',header=TRUE)


#######################################
#### Metropolis-Hastings algorithm ####
#######################################

B_function = function(delta_t,a,b){
  B = (1-exp(-a*(delta_t)))/a
  return(B)
}

A_function = function(delta_t,a,b){
  A= (b-sigma_r^2/(2*a^2))*(B_function(delta_t,a,b)-delta_t)- (sigma_r^2)*B_function(delta_t,a,b)^2/(4*a)
  return(A)
}

run_MCMC<- function(startvalue,iterations,n_out_of_sample=0,train_start='3/31/1998'
                    ,train_end='12/29/2017',data=data){
  
  data = data[1:(nrow(data)-n_out_of_sample),]
  data = data[which(data$Date==train_start):which(data$Date==train_end),]
  
  obs_r = data$X1M
  y = data[,c('X1Y','X2Y','X5Y','X10Y','X15Y','X20Y','X30Y')]
  y_delta_t = c(1,2,5,10,15,20,30)     
  r = obs_r[2:(length(obs_r))]
  
  ### X, sigma_r, sigma_y and Sigma matrix
  X = t(rbind(rep(1,length(r)),obs_r[1:(length(obs_r)-1)]))
  sigma_r = sd(diff(obs_r))
  #sigma_r = sd((obs_r))
  sigma_y = cov(y)
  Sigma = sigma_r^2*solve(t(X)%*%X)
  
  #### MLE estimates of alpha (alpha0 and alpha1)
  alpha_hat = solve(t(X)%*%X)%*%t(X)%*%r
  
  #############################################
  ### 3 sub functions will be used in MCMC ####
  #############################################
  
  #### function A and B
  B_function = function(delta_t,a,b){
    B = (1-exp(-a*(delta_t)))/a
    return(B)
  }
  
  A_function = function(delta_t,a,b){
    A= (b-sigma_r^2/(2*a^2))*(B_function(delta_t,a,b)-delta_t)- (sigma_r^2)*B_function(delta_t,a,b)^2/(4*a)
    return(A)
  }
  
  ##### function f(y|alpha) used in step 3 of MCMC simulation
  ##### -1/2*(y-beta0-beta0*r)T * (sigma^y)^(-1) * (y-beta0-beta0*r)
  
  f_y = function(alpha0,alpha1,n){
    ## change
    a = -12*log(alpha1)
    b = alpha0/(1-alpha1)
    beta_0 = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)/(-delta_t)}))
    beta_1 = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)/delta_t}))
    f = (-0.5*as.matrix((y[2:nrow(y),][n,])-beta_0-t(beta_1 %*% t(r[n]))) %*% solve(sigma_y) %*% t(as.matrix((y[2:nrow(y),][n,])-beta_0-t(beta_1 %*% t(r[n])))))
    return(f)
  }
  
  ##### initialize lists to store the values of alpha0 and alpha1
  alpha0_dict = c()
  alpha1_dict = c()
  accept_prob_dict = c()
  a_dict = c()
  b_dict = c()
  
  ##### 1.set an initial vector of values (alpha(0)), j= 0 
  alpha_0 = startvalue
  alpha_j = alpha_0
  
  ##### MCMC loop
  for (k in 1:iterations){
    set.seed(k)
    ##### 2.generate (alpha (j+1)) from pi(alpha(j))
    ##### alpha_j from previous iteration
    pi_Mu = solve(solve(Sigma) + t(X)%*%X) %*% (solve(Sigma)%*%alpha_j+t(X)%*%r)
    pi_Sigma = solve(solve(Sigma)+ (sigma_r)^(-2)*t(X)%*%X)
    ##### simulate alpha0 and alpha1 from bivariate normal distribution
    alpha_j_1 = mvrnorm(n=1,mu=pi_Mu,Sigma=pi_Sigma)
    
    #### alpha1 should always be positive, this ensures it is a positive number
    while(alpha_j_1[2]<0){
      alpha_j_1 = mvrnorm(n=1,mu=pi_Mu,Sigma=pi_Sigma)
    }
    
    ##### 3.calculate the proposal acceptance probability alpha 
    ##### Sample number from uniform (0,1) distribution, if greater than acceptance probability, then reject, otherwise accept.
    
    #product_j_1 = 0
    #product_j = 0
    #for(i in 1:length(r)){
    #  product_j = product_j + f_y(alpha_j[1],alpha_j[2],i)
    #  product_j_1 = product_j_1 + f_y(alpha_j_1[1],alpha_j_1[2],i)
    #}
    product_j = sum(unlist(lapply(1:length(r),function(i){f_y(alpha_j[1],alpha_j[2],i)})))
    product_j_1 = sum(unlist(lapply(1:length(r),function(i){f_y(alpha_j_1[1],alpha_j_1[2],i)})))
    
    accept_prob = exp(product_j_1-product_j)
    accept_prob_dict = c(accept_prob_dict,accept_prob)
    
    #### determine if accept or not, compared with a random simulated value from uniform distribution (0,1)
    if(accept_prob > runif(1,0,1) & !is.na(accept_prob)){
      alpha0_dict = c(alpha0_dict,alpha_j_1[1])
      alpha1_dict = c(alpha1_dict,alpha_j_1[2])
      alpha_j = alpha_j_1
      #### change
      a_dict = c(a_dict,-log(alpha_j_1[2])*12)
      b_dict = c(b_dict,alpha_j_1[1]/(1-alpha_j_1[2]))
    }else{
      alpha0_dict = c(alpha0_dict,alpha_j[1])
      alpha1_dict = c(alpha1_dict,alpha_j[2])
      alpha_j = alpha_j
      a_dict = c(a_dict,-log(alpha_j[2])*12)
      b_dict = c(b_dict,alpha_j[1]/(1-alpha_j[2]))
    }
    print(paste('Progress: ',k*100/iterations,'%'))  
  } ###  end of iteration
  #return(list(alpha0_dict,alpha1_dict,accept_prob_dict))
  return(list(alpha0_dict,alpha1_dict,accept_prob_dict,a_dict,b_dict))
}### end of MCMC function

#################################################################
#### function to perform in-sample and out-of-sample testing ####
#################################################################

Vasciek_prediction = function(a_vec,b_vec,train_start='12/31/2008',train_end='12/29/2017'
                              ,data,type = 'In-Sample'){
  pdf(paste(type,'_calibration_1.pdf',sep=''))
  for(i in 1:length(a_vec)){
    #### set up
    plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
        xaxt='n', yaxt='n', xlab='', ylab='')
    text(1,4,paste('Parameter a:',a_vec[i]), pos=4)
    text(1,3,paste('Parameter b:',b_vec[i]), pos=4)
    text(1,2,paste('Data starting Date:', train_start), pos=4)
    text(1,1,paste('Data ending Date:', train_end), pos=4)
    points(rep(1,4),1:4, pch=15)
    a = a_vec[i]
    b = b_vec[i]
    alpha0= b*(1-exp(-a/12))
    alpha1= exp(-a/12)
  #data_test = data[(nrow(data)-n_out_of_sample+1):nrow(data),]
    data_train = data[which(data$Date==train_start):which(data$Date==train_end),]
  
    obs_r = data_train$X1M
    y = data_train[,c('X1Y','X2Y','X5Y','X10Y','X15Y','X20Y','X30Y')]
    y_delta_t = c(1,2,5,10,15,20,30)     
    r = obs_r[2:(length(obs_r))]
  
  ### X, sigma_r, sigma_y and Sigma matrix
    X = t(rbind(rep(1,length(r)),obs_r[1:(length(obs_r)-1)]))
    sigma_r = sd(diff(obs_r))
    sigma_y = cov(y)
    Sigma = sigma_r^2*solve(t(X)%*%X)
    
    B_function = function(delta_t,a,b){
      B = (1-exp(-a*(delta_t)))/a
      return(B)
    }
    
    A_function = function(delta_t,a,b){
      A= (b-sigma_r^2/(2*a^2))*(B_function(delta_t,a,b)-delta_t)- (sigma_r^2)*B_function(delta_t,a,b)^2/(4*a)
      return(A)
    }
  
    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)
    clusterExport(cl, c("alpha0","alpha1","a","b","rnorm","obs_r","sigma_r","n_out_of_sample","data_train",
                        "A_function","B_function","y_delta_t","type","n_pred"), envir=environment())
    predict_in_sample_list <- pblapply(1:n_iter, cl = cl, function(k){
      predict_r_dict = c()
      predict_y_dict = c()
      if(type !='In-Sample'){
        predict_r_dict[1] = data_train$X1M[nrow(data_train)]
        project_length = n_pred
      }else{
        project_length = length(obs_r)
      }
      for(j in 2:project_length){
        #set.seed(j)
        epilson_r = rnorm(1,mean=0,sd=1)
        if(type == 'In-Sample'){
          predict_r = alpha0 + alpha1*obs_r[j-1] + sigma_r*epilson_r*sqrt(12)}else{
            predict_r = alpha0 + alpha1*predict_r_dict[j-1] + sigma_r*epilson_r*sqrt(12)
          }
        ### add into the list
        predict_r_dict = c(predict_r_dict,predict_r)
      }
      return((predict_r_dict))
    })
    in_sample_r=do.call(rbind,predict_in_sample_list)
    in_sample_r= melt(t(in_sample_r))
    names(in_sample_r) = c('Date','Num_iteration','Short_rate')
  
    A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
    B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
  
    if(type == 'In-Sample'){  
      in_sample_r_97_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.975)}))
      in_sample_r_02_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.025)}))
      in_sample_r_mean = unlist(lapply(1:(nrow(y)-1),function(x){mean(unlist(lapply(predict_in_sample_list, `[[`, x)))}))
  
      in_sample_y_97_5 = do.call(rbind,lapply(in_sample_r_97_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
      in_sample_y_02_5 = do.call(rbind,lapply(in_sample_r_02_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
      in_sample_y_mean = do.call(rbind,lapply(in_sample_r_mean, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
      colnames(in_sample_y_mean) = names(y)
      yield_curve_est_mean = melt(in_sample_y_mean)
      yield_curve_est_mean$cat = 'Estimated (Mean)'
  
      colnames(in_sample_y_02_5) = names(y)
      yield_curve_est_02_5 = melt(in_sample_y_02_5)
      yield_curve_est_02_5$cat = 'Estimated (2.5%)'
  
      colnames(in_sample_y_97_5) = names(y)
      yield_curve_est_97_5 = melt(in_sample_y_97_5)
      yield_curve_est_97_5$cat = 'Estimated (97.5%)'
  
      #### predicted yield curve with actual y
      in_sample_r_actual = do.call(rbind,lapply(obs_r[-1],function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
      colnames(in_sample_r_actual) = names(y)
      yield_curve_est_actual_r = melt(in_sample_r_actual)
      yield_curve_est_actual_r$cat = 'Estimated Yield Curve with Actual Short Rate'
  
      yield_curve_est = rbind(yield_curve_est_mean,yield_curve_est_02_5,yield_curve_est_97_5,yield_curve_est_actual_r)
      yield_curve_est = cbind(data_train$Date[-1],yield_curve_est[,-1])
      names(yield_curve_est) = c('Month','Maturity','Value','cat')
  
      ### actual yield curve
      yield_curve_actual = melt(data_train[-1,c('Date',names(y))],1)
      names(yield_curve_actual) = c('Month','Maturity','Value')
      yield_curve_actual$cat = 'Actual'
  
      ### estimated yield curve and actual yield curve
      yield_curve_all = rbind(yield_curve_actual,yield_curve_est)
      yield_curve_all$Month = as.character(yield_curve_all$Month)
      yield_curve_all$Maturity = as.numeric(gsub('X|Y','',yield_curve_all$Maturity))
  
      date_list= unique(yield_curve_all$Month)
      for(i in 1:ceiling(length(date_list)/4)){
        ind1 = i*4 - 3 
        ind2 = min(i*4,length(date_list))
        p = xyplot(Value~Maturity|Month, group=cat,data= yield_curve_all[yield_curve_all$Month %in% date_list[ind1:ind2],]
                  ,ylab='Yield Curves', xlab='Maturity'
                  ,main = paste('Yield Curves 1 Year to 30 Year Maturity (a:',a,', b:',b,')',by='')
                  #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
                  ,auto.key = list(x=0.7,y=0.85, columns=1, cex=0.3,border = FALSE, lines = TRUE, points=F)
                  ,type=c('g', 'l')
                  ,par.settings = list(superpose.line = list(
                    lty = c(1,2,2,1,1)
                    ,lwd=c(1, rep(1, 5)))
                    ,superpose.symbol = list(col= brewer.pal(10, 'RdYlGn')[2:10]))
        )
      print(p)
      }
    } else{
      in_sample_r_mean = unlist(lapply(1:(n_pred),function(x){mean(unlist(lapply(predict_in_sample_list, `[[`, x)))}))
      in_sample_y_mean = do.call(rbind,lapply(in_sample_r_mean, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
    
      colnames(in_sample_y_mean) = names(y)
      yield_curve_est_mean = melt(in_sample_y_mean)
      names(yield_curve_est_mean) = c('Month','Maturity','Value')
      yield_curve_all = yield_curve_est_mean

      yield_curve_all$Month = as.numeric(yield_curve_all$Month)
      yield_curve_all$Maturity = as.numeric(gsub('X|Y','',yield_curve_all$Maturity))
    
      p=xyplot(Value~Month, group=Maturity,data= yield_curve_all
              ,ylab='Yield Curves', xlab='Months'
              ,main = paste('Yield Curves 1 Year to 30 Year Maturity (monthly mean from',n_iter,'Simulations)')
              #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
              ,auto.key = list(space='right', columns=1, border = FALSE, lines = TRUE, points=F)
              ,type=c('g', 'l')
              ,par.settings = list(superpose.line = list(
                lwd=c(3, rep(1.5, 5)))
                ,superpose.symbol = list(col= brewer.pal(6, 'RdYlGn')))
      )
      print(p)
      }
    }
    dev.off()
    return(list(predict_in_sample_list,A_delta_t,B_delta_t,y_delta_t,y))
}

###################################################
###### calibration with shorter time periods ######
###################################################

### arguments set up
n_out_of_sample = 0
train_start= '12/31/2008'
train_end = '12/29/2017'
n_MCMC = 2000
n_iter = 10000
n_pred = 0

### inital values of alpha
alpha_initial = c(0.0008298707,0.9917013)   ### starting from a = 0.1 and b = 0.1

#result = run_MCMC(alpha_initial,n_MCMC,n_out_of_sample,train_start)
result = run_MCMC(alpha_initial,n_MCMC,n_out_of_sample,train_start,train_end=train_end,data)


#a_vec = c(0.4283286)             ### control the speed
#b_vec = c(0.02207140)          ### control the mean

#a_vec = c(0.3743440)             ### control the speed
#b_vec = c(0.02230472)          ### control the mean

a_vec = c(0.3677407)             ### control the speed
b_vec = c(0.02338157)          ### control the mean

output=Vasciek_prediction(a_vec,b_vec,train_start='12/31/2008',train_end='12/29/2017'
                          ,data,type = 'In-Sample')

pdf('In_of_sample_calibration_2009.pdf')
for(i in 1:length(a_vec)){
  #### set up
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(1,4,paste('Parameter a:',a_vec[i]), pos=4)
  text(1,3,paste('Parameter b:',b_vec[i]), pos=4)
  text(1,2,paste('Data starting Date:', train_start), pos=4)
  text(1,1,paste('Data ending Date:', data$Date[nrow(data)]), pos=4)
  points(rep(1,4),1:4, pch=15)
  a = a_vec[i]
  b = b_vec[i]
  alpha0= b*(1-exp(-a/12))
  alpha1= exp(-a/12)
  #data_test = data[(nrow(data)-n_out_of_sample+1):nrow(data),]
  data_train = data[which(data$Date==train_start):nrow(data),]
  
  obs_r = data_train$X1M
  y = data_train[,c('X1Y','X2Y','X5Y','X10Y','X15Y','X20Y','X30Y')]
  y_delta_t = c(1,2,5,10,15,20,30)     
  r = obs_r[2:(length(obs_r))]
  
  ### X, sigma_r, sigma_y and Sigma matrix
  X = t(rbind(rep(1,length(r)),obs_r[1:(length(obs_r)-1)]))
  sigma_r = sd(diff(obs_r))
  sigma_y = cov(y)
  Sigma = sigma_r^2*solve(t(X)%*%X)
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("alpha0","alpha1","a","b","rnorm","obs_r","sigma_r","n_out_of_sample","data_train",
                      "A_function","B_function","y_delta_t"), envir=environment())
  predict_in_sample_list <- pblapply(1:n_iter, cl = cl, function(k){
    predict_r_dict = c()
    predict_y_dict = c()
    for(j in 2:length(obs_r)){
      #set.seed(j)
      epilson_r = rnorm(1,mean=0,sd=1)
      predict_r = alpha0 + alpha1*obs_r[j-1] + sigma_r*epilson_r
      #A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
      #B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
      #predict_y = (B_delta_t*predict_r - A_delta_t)/y_delta_t
      ### add into the list
      predict_r_dict = c(predict_r_dict,predict_r)
      #predict_y_dict = rbind(predict_y_dict,predict_y)
    }
    return((predict_r_dict))
  })
  in_sample_r=do.call(rbind,predict_in_sample_list)
  in_sample_r= melt(t(in_sample_r))
  names(in_sample_r) = c('Date','Num_iteration','Short_rate')
  
  A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
  B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
  
  in_sample_r_97_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.975)}))
  in_sample_r_02_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.025)}))
  in_sample_r_mean = unlist(lapply(1:(nrow(y)-1),function(x){mean(unlist(lapply(predict_in_sample_list, `[[`, x)))}))
  
  in_sample_y_97_5 = do.call(rbind,lapply(in_sample_r_97_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  in_sample_y_02_5 = do.call(rbind,lapply(in_sample_r_02_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  in_sample_y_mean = do.call(rbind,lapply(in_sample_r_mean, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
  colnames(in_sample_y_mean) = names(y)
  yield_curve_est_mean = melt(in_sample_y_mean)
  yield_curve_est_mean$cat = 'Estimated (Mean)'
  
  colnames(in_sample_y_02_5) = names(y)
  yield_curve_est_02_5 = melt(in_sample_y_02_5)
  yield_curve_est_02_5$cat = 'Estimated (2.5%)'
  
  colnames(in_sample_y_97_5) = names(y)
  yield_curve_est_97_5 = melt(in_sample_y_97_5)
  yield_curve_est_97_5$cat = 'Estimated (97.5%)'
  
  #### predicted yield curve with actual y
  in_sample_r_actual = do.call(rbind,lapply(obs_r[-1],function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
  colnames(in_sample_r_actual) = names(y)
  yield_curve_est_actual_r = melt(in_sample_r_actual)
  yield_curve_est_actual_r$cat = 'Estimated Yield Curve with Actual Short Rate'
  
  yield_curve_est = rbind(yield_curve_est_mean,yield_curve_est_02_5,yield_curve_est_97_5,yield_curve_est_actual_r)
  yield_curve_est = cbind(data_train$Date[-1],yield_curve_est[,-1])
  names(yield_curve_est) = c('Month','Maturity','Value','cat')
  
  ### actual yield curve
  yield_curve_actual = melt(data_train[-1,c('Date',names(y))],1)
  names(yield_curve_actual) = c('Month','Maturity','Value')
  yield_curve_actual$cat = 'Actual'
  
  ### estimated yield curve and actual yield curve
  yield_curve_all = rbind(yield_curve_actual,yield_curve_est)
  yield_curve_all$Month = as.character(yield_curve_all$Month)
  yield_curve_all$Maturity = as.numeric(gsub('X|Y','',yield_curve_all$Maturity))
  
  date_list= unique(yield_curve_all$Month)
  for(i in 1:ceiling(length(date_list)/4)){
    ind1 = i*4 - 3 
    ind2 = min(i*4,length(date_list))
    p = xyplot(Value~Maturity|Month, group=cat,data= yield_curve_all[yield_curve_all$Month %in% date_list[ind1:ind2],]
               ,ylab='Yield Curves', xlab='Maturity'
               ,main = paste('Yield Curves 1 Year to 30 Year Maturity (a:',a,', b:',b,')',by='')
               #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
               ,auto.key = list(x=0.7,y=0.85, columns=1, cex=0.3,border = FALSE, lines = TRUE, points=F)
               ,type=c('g', 'l')
               ,par.settings = list(superpose.line = list(
                 lty = c(1,2,2,1,1)
                 ,lwd=c(1, rep(1, 5)))
                 ,superpose.symbol = list(col= brewer.pal(10, 'RdYlGn')[2:10]))
    )
    print(p)
  }
}
dev.off()

###########################
#####  Stress Testing #####
###########################

### arguments set up
n_out_of_sample = 0
train_start= '12/31/2008'
n_MCMC = 2000
n_iter = 10000
shock = - 0.01

data_ST = cbind(data$Date,data[,-1]+shock)
names(data_ST) = names(data)

### inital values of alpha
alpha_initial = c(0.0008298707,0.9917013)   ### starting from a = 0.1 and b = 0.1

result = run_MCMC(alpha_initial,n_MCMC,n_out_of_sample,train_start,data=data_ST)

a_vec = c(0.3705106)             ### estimated parameter from shocking 100bps down
b_vec = c(0.01223165)          ### estimated parameter from shocking 100bps down

output=Vasciek_prediction(a_vec,b_vec,train_start='12/31/2008',train_end='12/29/2017'
                          ,data_ST,type = 'In-Sample')


pdf('Stress_Testing.pdf')
for(i in 1:length(a_vec)){
  #### set up
  plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
       xaxt='n', yaxt='n', xlab='', ylab='')
  text(1,5,paste('Parameter a:',a_vec[i]), pos=4)
  text(1,4,paste('Parameter b:',b_vec[i]), pos=4)
  text(1,3,paste('Data starting Date:', train_start), pos=4)
  text(1,2,paste('Data ending Date:', data_ST$Date[nrow(data)]), pos=4)
  text(1,1,paste('Shock:', shock), pos=4)
  points(rep(1,5),1:5, pch=15)
  a = a_vec[i]
  b = b_vec[i]
  alpha0= b*(1-exp(-a/12))
  alpha1= exp(-a/12)
  #data_test = data[(nrow(data)-n_out_of_sample+1):nrow(data),]
  data_train = data_ST[which(data_ST$Date==train_start):nrow(data_ST),]
  
  obs_r = data_train$X1M
  y = data_train[,c('X1Y','X2Y','X5Y','X10Y','X15Y','X20Y','X30Y')]
  y_delta_t = c(1,2,5,10,15,20,30)     
  r = obs_r[2:(length(obs_r))]
  
  ### X, sigma_r, sigma_y and Sigma matrix
  X = t(rbind(rep(1,length(r)),obs_r[1:(length(obs_r)-1)]))
  sigma_r = sd(diff(obs_r))
  sigma_y = cov(y)
  Sigma = sigma_r^2*solve(t(X)%*%X)
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("alpha0","alpha1","a","b","rnorm","obs_r","sigma_r","n_out_of_sample","data_train",
                      "A_function","B_function","y_delta_t"), envir=environment())
  predict_in_sample_list <- pblapply(1:n_iter, cl = cl, function(k){
    predict_r_dict = c()
    predict_y_dict = c()
    for(j in 2:length(obs_r)){
      #set.seed(j)
      epilson_r = rnorm(1,mean=0,sd=1)
      predict_r = alpha0 + alpha1*obs_r[j-1] + sigma_r*epilson_r
      #A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
      #B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
      #predict_y = (B_delta_t*predict_r - A_delta_t)/y_delta_t
      ### add into the list
      predict_r_dict = c(predict_r_dict,predict_r)
      #predict_y_dict = rbind(predict_y_dict,predict_y)
    }
    return((predict_r_dict))
  })
  in_sample_r=do.call(rbind,predict_in_sample_list)
  in_sample_r= melt(t(in_sample_r))
  names(in_sample_r) = c('Date','Num_iteration','Short_rate')
  
  A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
  B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
  
  in_sample_r_97_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.975)}))
  in_sample_r_02_5 = unlist(lapply(1:(nrow(y)-1),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),0.025)}))
  in_sample_r_mean = unlist(lapply(1:(nrow(y)-1),function(x){mean(unlist(lapply(predict_in_sample_list, `[[`, x)))}))
  
  in_sample_y_97_5 = do.call(rbind,lapply(in_sample_r_97_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  in_sample_y_02_5 = do.call(rbind,lapply(in_sample_r_02_5, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  in_sample_y_mean = do.call(rbind,lapply(in_sample_r_mean, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
  colnames(in_sample_y_mean) = names(y)
  yield_curve_est_mean = melt(in_sample_y_mean)
  yield_curve_est_mean$cat = 'Estimated (Mean)'
  
  colnames(in_sample_y_02_5) = names(y)
  yield_curve_est_02_5 = melt(in_sample_y_02_5)
  yield_curve_est_02_5$cat = 'Estimated (2.5%)'
  
  colnames(in_sample_y_97_5) = names(y)
  yield_curve_est_97_5 = melt(in_sample_y_97_5)
  yield_curve_est_97_5$cat = 'Estimated (97.5%)'
  
  #### predicted yield curve with actual y
  in_sample_r_actual = do.call(rbind,lapply(obs_r[-1],function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  
  colnames(in_sample_r_actual) = names(y)
  yield_curve_est_actual_r = melt(in_sample_r_actual)
  yield_curve_est_actual_r$cat = 'Estimated Yield Curve with Actual Short Rate'
  
  yield_curve_est = rbind(yield_curve_est_mean,yield_curve_est_02_5,yield_curve_est_97_5,yield_curve_est_actual_r)
  yield_curve_est = cbind(data_train$Date[-1],yield_curve_est[,-1])
  names(yield_curve_est) = c('Month','Maturity','Value','cat')
  
  ### actual yield curve
  yield_curve_actual = melt(data_train[-1,c('Date',names(y))],1)
  names(yield_curve_actual) = c('Month','Maturity','Value')
  yield_curve_actual$cat = 'Actual'
  
  ### estimated yield curve and actual yield curve
  yield_curve_all = rbind(yield_curve_actual,yield_curve_est)
  yield_curve_all$Month = as.character(yield_curve_all$Month)
  yield_curve_all$Maturity = as.numeric(gsub('X|Y','',yield_curve_all$Maturity))
  
  date_list= unique(yield_curve_all$Month)
  for(i in 1:ceiling(length(date_list)/4)){
    ind1 = i*4 - 3 
    ind2 = min(i*4,length(date_list))
    p = xyplot(Value~Maturity|Month, group=cat,data= yield_curve_all[yield_curve_all$Month %in% date_list[ind1:ind2],]
               ,ylab='Yield Curves', xlab='Maturity'
               ,main = paste('Yield Curves 1 Year to 30 Year Maturity (a:',a,', b:',b,')',by='')
               #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
               ,auto.key = list(x=0.7,y=0.85, columns=1, cex=0.3,border = FALSE, lines = TRUE, points=F)
               ,type=c('g', 'l')
               ,par.settings = list(superpose.line = list(
                 lty = c(1,2,2,1,1)
                 ,lwd=c(1, rep(1, 5)))
                 ,superpose.symbol = list(col= brewer.pal(10, 'RdYlGn')[2:10]))
    )
    print(p)
  }
}
dev.off()


########################################
###### Out of sample simulation ########
########################################

### arguments set up
n_out_of_sample = 0
train_start= '12/31/2008'
train_end = '12/29/2017'
n_MCMC = 2000
n_iter = 10000
n_pred = 25*12

### inital values of alpha
alpha_initial = c(0.0008298707,0.9917013)   ### starting from a = 0.1 and b = 0.1

result = run_MCMC(alpha_initial,n_MCMC,n_out_of_sample,train_start,train_end=train_end,data)

a_vec = c(0.3677407)             ### estimated parameter from shocking 100bps down
b_vec = c(0.02338157)          ### estimated parameter from shocking 100bps down


output=Vasciek_prediction(a_vec,b_vec,train_start='12/31/2008',train_end='12/29/2017'
                              ,data,type = 'Out-Of-Sample')

predict_in_sample_list = output[[1]]
A_delta_t              = output[[2]]
B_delta_t              = output[[3]]
y_delta_t              = output[[4]]
y                      = output[[5]]


#### calculate the quantile
cal_quantile = function(quantile,maturity){
  simulated_r = unlist(lapply(1:(n_pred),function(x){quantile(unlist(lapply(predict_in_sample_list, `[[`, x)),quantile)}))
  simulated_y = do.call(rbind,lapply(simulated_r, function(r){(B_delta_t*r - A_delta_t)/y_delta_t}))
  simulated_y2 = simulated_y[c(1*12,5*12,10*12,15*12,20*12,25*12),]
  names(simulated_y2) = names(y)
  output = simulated_y2[,which(names(simulated_y2)==maturity)]
  return(output)
}

table_quantile = function(maturity='X1Y'){
  quantile_list = c(0,0.001,0.005,0.01,0.05,0.1,0.25,0.5,0.75
                    ,0.9,0.95,0.99,0.995,0.999,1)
  table = do.call(rbind
                  ,lapply(quantile_list,function(x){cal_quantile(x,maturity)}))
  table = as.data.frame(table)
  table = cbind(paste(quantile_list*100,'%',sep=''),table)
  names(table) = c(' ','1','5','10','15','20','25')
  return(table)
}

#### Libor Scenario Statistics
Libor_1Y=table_quantile('X1Y')
Libor_2Y=table_quantile('X2Y')
Libor_5Y=table_quantile('X5Y')
Libor_10Y=table_quantile('X10Y')
Libor_15Y=table_quantile('X15Y')
Libor_20Y=table_quantile('X20Y')
Libor_30Y=table_quantile('X30Y')


WriteXLS(c("Libor_1Y",
           "Libor_2Y",
           "Libor_5Y",
           "Libor_10Y",
           "Libor_15Y",
           "Libor_20Y",
           "Libor_30Y"),
         "Libor_quantile_table.xlsx")

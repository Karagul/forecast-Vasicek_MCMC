#all packages needed
requirements <- c("lattice",'MASS','ggplot2','grid','snow','pbapply',"parallel","reshape","RColorBrewer")
#install packages if you have not
for(requirement in requirements){if( !(requirement %in% installed.packages())) install.packages(requirement)}
#load all required packages
lapply(requirements, require, character.only=T);

#######################################################################################
#### for testing purpose, simulate observed short rate and zero coupon bond yields ####
#######################################################################################

### set up the directory and load data
setwd('C:/Users/rstrepparava/Documents/Projects_2018/TIAA')
#data = read.csv('Libor 1M and zero rates data for Jia.csv',header=TRUE)
data =read.csv('Monthly Curve Data USD.csv',header=TRUE)

### only use most recent 5 years data
#data = data[(nrow(data)-120):nrow(data),]

### short rate and y 
# obs_r = data$Libor.1M
# y = data[,c('X3M','X6M','X1Y','X2Y','X5Y','X10Y','X15Y','X20Y','X30Y')]
# y_delta_t = c(0.25,0.5,1,2,5,10,15,20,30)     
# r = obs_r[2:(length(obs_r))]

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
  A= (b-sigma_r^2/(2*a^2))*(B_function(delta_t,a,b)-delta_t)- (sigma_r^2)*B_function(delta_t,a,b)^2/4*a
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

#######################################
#### Metropolis-Hastings algorithm ####
#######################################

run_MCMC<- function(startvalue,iterations){
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


##################
#### Testing #####
##################

### inital values of alpha
#alpha_initial = c(0.0001198097,0.9968304)
alpha_initial = c(0.002438529,0.9512294)   ### starting from a = 0.05 and b = 0.05

result = run_MCMC(alpha_initial,2000)

save(result,file='result_4000.RData')

### histogram of the simulated alpha0 and alpha1
pdf('MCMC_Vasicek.pdf')
plot1 <- 
  qplot(x=1:length(result[[4]]),y=result[[4]]
        ,main = 'Convergence of the parameter a'
        ,xlab = '# of iteration'
        ,ylab = 'a'
        ,geom=c("line","point")
        ,fill=I('red')
        ,size=I(0.4)
        ,col=I('blue'))
plot2 <- 
  qplot(x=1:length(result[[5]]),y=result[[5]]
        ,main = 'Convergence of the parameter b'
        ,xlab = '# of iteration'
        ,ylab = 'b'
        ,geom=c("line","point")
        ,fill=I('red')
        ,size=I(0.4)
        ,col=I('blue'))
plot1
plot2
#grid.arrange(plot2, plot3, ncol=2)


#############################################################
#### Prediction for short rate r(t) and different yields ####
#############################################################

#alpha0= 6.059986e-05
#alpha1= 0.9987300

#alpha0 = 0.001175000   #### 1000 iterations
#alpha1 = 0.9674624   #### 1000 iterations

alpha0 = 0.001109794   #### 4000 iterations
alpha1 = 0.9685025   #### 1000 iterations

#alpha0 = 0.0007219632
#alpha1 = 0.9660325

a = -log(alpha1)*12
b = alpha0/(1-alpha1)

#### one step ahead in-sample prediction
n_iter = 1000
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("alpha0","alpha1","a","b","rnorm","obs_r","sigma_r",
                    "A_function","B_function","y_delta_t"), envir=environment())
predict_list <- pblapply(1:n_iter, cl = cl, function(k) {
  predict_r_dict = c()
  predict_y_dict = c()
  for(i in 2:length(obs_r)){
    set.seed(i)
    epilson_r = rnorm(1,mean=0,sd=1)
    predict_r = alpha0 + alpha1*obs_r[i-1] + sigma_r*epilson_r
    A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
    B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
    predict_y = (B_delta_t*predict_r - A_delta_t)/y_delta_t
    ### add into the list
    predict_r_dict = c(predict_r_dict,predict_r)
    predict_y_dict = rbind(predict_y_dict,predict_y)
  }
  return(list(predict_r_dict,predict_y_dict))
})

#### extract the 2nd element
predict_list_temp=lapply(predict_list, `[[`, 2)
predict_list_avg = predict_list_temp[[1]]
for(i in 2:length(predict_list_temp)){
  predict_list_avg = predict_list_avg + predict_list_temp[[i]]
}

#### estimated yield curve
est_yield_curve = predict_list_avg/length(predict_list_temp)

#####################################################################
##### Yield curve estimated by yield curve fitting (In-sample) ######
######################################################################

plotPrep = function(date){
  ind = which(data$Date==date)
  RealYield = y[ind,]
  EstYield = est_yield_curve[ind-1,]
  return(cbind(c(RealYield,EstYield),c(rep('real yield curve',length(RealYield)),rep('estimated yield curve',length(RealYield)))))
}

plotYieldCurve = function(date){
        dataPlot = data.frame(
                          Years = rep(y_delta_t,2),
                          Rate =  unlist(plotPrep(date)[,1]),
                          Category = unlist(plotPrep(date)[,2])
                    )

        p4 <- xyplot(Rate~Years, group=Category,data= dataPlot
                      ,ylab='Rate', xlab='Years'
                      ,main = paste('Yield curve estimated by yield curve fitting on',date)
                      ,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
                      ,type=c('g', 'l')
                      ,par.settings = list(superpose.line = list(
                                      lwd=c(3, rep(1.5, 5)))
                      ,superpose.symbol = list(col= 'blue')))
print(p4)
}

plotYieldCurve('2/28/2018')
plotYieldCurve('9/30/2013')
plotYieldCurve('10/30/2009')
plotYieldCurve('1/31/2017')
plotYieldCurve('1/31/2018')

dev.off()

#######################################
##### Out-of-sample Prediction  #######
#######################################

#### arguments
n_iter = 1000
n_pred = 100
start_date = '4/30/2018'
#alpha0 = 5.802376e-05
#alpha1 = 0.9987713
a = -log(alpha1)*12
b = alpha0/(1-alpha1)

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterExport(cl, c("alpha0","alpha1","a","b","rnorm","obs_r","sigma_r","n_pred","start_date","data",
                    "A_function","B_function","y_delta_t"), envir=environment())
predict_out_sample_list <- pblapply(1:n_iter, cl = cl, function(k){
                                ind = which(data$Date==start_date)
                                start_r = obs_r[ind]
                                predict_r = alpha0 + alpha1*start_r + sigma_r*rnorm(1,mean=0,sd=1)
                                for(i in 2:n_pred){
                                  epilson_r = rnorm(1,mean=0,sd=1)
                                  predict_r = c(predict_r,alpha0 + alpha1*predict_r[i-1] + sigma_r*epilson_r)
                                }
                                A_delta_t = unlist(lapply(y_delta_t,function(delta_t){A_function(delta_t,a,b)}))
                                B_delta_t = unlist(lapply(y_delta_t,function(delta_t){B_function(delta_t,a,b)}))
                                predict_r_dict = predict_r
                                predict_y_dict = do.call(rbind,lapply(predict_r,function(x){return((B_delta_t*x - A_delta_t)/y_delta_t)}))
                                return(list(predict_r_dict,predict_y_dict))
                            })
out_of_sample_r=do.call(rbind,lapply(predict_out_sample_list, `[[`, 1))
out_of_sample_r= melt(t(out_of_sample_r))
names(out_of_sample_r) = c('Date','Num_iteration','Short_rate')

out_of_sample_y = lapply(predict_out_sample_list, `[[`, 2)
out_of_sample_y_sum = out_of_sample_y[[1]]
for(i in 2:length(out_of_sample_y)){
  out_of_sample_y_sum = out_of_sample_y_sum + out_of_sample_y[[i]]
}

#### out-0f-sample estimated yield curve
est_yield_curve_out_sample = out_of_sample_y_sum/length(out_of_sample_y)

#### all simulated out-of-sample short rate 
xyplot(Short_rate~Date, group=Num_iteration,data= out_of_sample_r
       ,ylab='Rate', xlab='Years',col='blue'
       ,main = paste('short rate out-of-sample prediction starting from ', start_date)
       #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
       ,type=c('g', 'l')
       ,par.settings = list(superpose.line = list(
         lwd=c(3, rep(1.5, 5)))
         ,superpose.symbol = list(col= 'blue')))

colnames(est_yield_curve_out_sample) = names(y)
est_yield_curve_future = melt(est_yield_curve_out_sample)
names(est_yield_curve_future) = c('Month','Maturity','Value')

est_yield_curve_future = est_yield_curve_future[order(match(est_yield_curve_future$Maturity,names(y))),]
est_yield_curve_future$Maturity = factor(est_yield_curve_future$Maturity,levels = names(y))

#### all simulated out-of-sample short rate 
xyplot(Value~Month, group=Maturity,data= est_yield_curve_future
       ,ylab='Yield Curves', xlab='Months'
       ,main = paste('Yield Curves 1 Year to 30 Year Maturity (monthly mean from',n_iter,'Simulations)')
       #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
       ,auto.key = list(space='right', columns=1, border = FALSE, lines = TRUE, points=F)
       ,type=c('g', 'l')
      ,par.settings = list(superpose.line = list(
        lwd=c(3, rep(1.5, 5)))
        ,superpose.symbol = list(col= brewer.pal(6, 'RdYlGn')))
       )


#### plot the first thirty out-of-sample curves to display Principal Components
N <- 30
index <- c()
for (k in seq(1,601,100))
  index <- c(index, seq(k, k+N-1))

first_thirty <- est_yield_curve_future[index, ]


xyplot(Value~Maturity, group=Month,data=first_thirty
       ,ylab='Yield Curves', xlab='Tenor'
       ,main = paste('First 30 Monthly OOS Yield Curves 1 Year to 30 Year Maturity')
       #,auto.key = list(space='top', columns=1, border = FALSE, lines = TRUE, points=F)
       ,auto.key = list(space='right', columns=1, border = FALSE, lines = TRUE, points=F)
       ,type=c('g', 'l')
       ,par.settings = list(superpose.line = list(
         lwd=c(3, rep(1.5, 5)))
         ,superpose.symbol = list(col= brewer.pal(6, 'RdYlGn')))
)


df   <- cast(first_thirty, Month ~ Maturity)
pca1 <- prcomp(df, scale. = TRUE)
# sqrt of eigenvalues
print(pca1$sdev)

# loadings
head(pca1$rotation)

# PCs (aka scores)
head(pca1$x)

# PCA plots

scores = as.data.frame(pca1$x)

# plot of observations
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of OOS Yield Curves")


# Circle of correlations

# function to create a circle
circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(df, pca1$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0, 0), y1 = c(0, 0, 0, 0, 0, 0, 0), 
                    x2 = correlations$PC1, y2 = correlations$PC2)

# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, 
                                                             colour = "gray65") + xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs", 
                                                                                                                           y = "pc2 axis") + ggtitle("Circle of correlations")



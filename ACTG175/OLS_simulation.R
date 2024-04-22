source("utils.R")
source("fun_model.R")
source("fun_data.R")

library(MASS)
library(dplyr)
library(rmutil)
library(epiDisplay)
library(ATE)
library(matrixcalc)
library(tmvtnorm)
library(corpcor)
library(prettyR)

B = 200
r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
n_r = length(r)
rep = 500


level = 0.05 #4.18改
n = 400
k = 2

######################################
model = 'Gaussian'

e = 200
theta_true = c(0,0,0,0,0,0,0,1)
k_par = 4
simulation(rep,B,level,e,r,theta_true,k_par,cov_matrix=diag(length(theta_true)),model)

##################################

OLS_simulation(rep,B,level,e,r,theta_true,k_par)


theta_true = c(0,0,1,0,1)
k_par = 2
OLS_simulation(rep,B,level,e,r,theta_true,k_par)

source("fun_data.R")
theta_true = rep(0,14)
k_par = 4
rep = 1
level = 0.1
e = 1000

group_method = 2
level = 0.22
e = 1000


source("fun_data.R")
group_method = 3
level = 0.05
e = 1000
#e_grid = seq(from=50, to=2000, by=100)
e_grid = c(1,5,30,50,100,500,1000,2000,10^4)
k_par = 4
OLS_grid_search(e_grid)

e_list = c(5, 500,5000)
for (i in (1:3)){
  e = e_list[i]
  result = OLS_application(rep,B,level,e,r,group_method,k_par)
  print(result)
  print("----------------------------------------------------")
}
OLS_application(rep,B,level,e,r,group_method,k_par)


source("fun_data.R")
k_par = 4
group_method = 3
OLS_application(rep,B,level,e,r,group_method,k_par)



group_method = 4
k_par = 8
level = 0.05
e = 5000
OLS_application(rep,B,level,e,r,group_method,k_par)

source("fun_logistic_regression.R")
group_method = 4
k_par = 8
level = 0.2
e = 5000
OLS_application_logistic(rep,B,level,e,r,group_method,k_par)
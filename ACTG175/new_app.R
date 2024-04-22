nCores = 16
myCluster <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(myCluster)

source("fun_model.R")
library(speff2trial)

library(foreach)
library(doParallel)
library(MASS)
library(dplyr)
library(rmutil)
library(epiDisplay)

del_var <- c("pidnum", "arms", "cd40", "cd420", "cd496", "r", "cd80", "cens", "days")
ACTG175 <- na.omit(ACTG175[,!colnames(ACTG175) %in% del_var])
X <- ACTG175[, !colnames(ACTG175) %in% c("cd820","zprior","treat")]
Y <- ACTG175[,"cd820"]

Ti<- ACTG175$treat

B = 200
r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
n_r = length(r)
rep = 500
level = 0.05 
e = 10000


##################### group by race & wtkg #####################
# n = dim(ACTG175)[1]
# n_subgroup <- 2 # race & wtkg
# k_n <- 2^n_subgroup
# k_par = k_n

# Z <- matrix(0,n,k_n)
# for (i in 1:n){ # subgrou = 2
#     indicator_2base = rep(0,n_subgroup) 
#     if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1 # dim1
#     if (X$wtkg[i] > median(X$wtkg)) indicator_2base[2] = indicator_2base[2]+1  # dim2
#     indicator_10base = indicator_2base[2]+indicator_2base[1]*2
#     Z[i,indicator_10base+1] = 1
# }
# colnames(Z) <- c("G1","G2","G3","G4") 

# DZ <- matrix(0,n,k_n)
# for (i in 1:n){
#     DZ[i,] = Ti[i]*Z[i,]
# }
# colnames(DZ) <- c("TG1","TG2","TG3","TG4") 

# X = cbind(Z,X) 
# D = cbind(DZ,X) 

# D <- D[, !colnames(D) %in% c("race")]
# D <- D[, !colnames(D) %in% c("G4")]
# D = as.matrix(D)

# star_time <- Sys.time() 
# tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
# best_subgroup = tmp$best_subgroup
# l_naive = tmp$nonpri_naive
# result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("dplyr","MASS","rmutil","epiDisplay","tmvtnorm","matrixcalc","corpcor")) %dopar% {
#     unname(fun_linear_regression_parallel(D,Y,B,level,e,r,BS_true,k_par))
# }
# parallel_result = rowMeans(result)
# final_result = c(parallel_result,l_naive,best_subgroup)
# final_result
# write.csv(result,file = "result/race_wtkg.csv")
# end_time <- Sys.time() 
# run_time <- end_time - star_time  ## 计算程序运行时间
# print("Time for the 100 epoch:")
# print(run_time)

##################### group by age & race #####################
# n = dim(ACTG175)[1]
# n_subgroup <- 2 
# k_n <- 2^n_subgroup
# k_par = k_n

# Z <- matrix(0,n,k_n)
# for (i in 1:n){ 
#     indicator_2base = rep(0,n_subgroup) 
#     if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1 
#     if (X$age[i] > median(X$age)) indicator_2base[2] = indicator_2base[2]+1  
#     indicator_10base = indicator_2base[2]+indicator_2base[1]*2
#     Z[i,indicator_10base+1] = 1
# }
# colnames(Z) <- c("G1","G2","G3","G4") 

# DZ <- matrix(0,n,k_n)
# for (i in 1:n){
#     DZ[i,] = Ti[i]*Z[i,]
# }
# colnames(DZ) <- c("TG1","TG2","TG3","TG4") 

# X = cbind(Z,X) 
# D = cbind(DZ,X) 

# D <- D[, !colnames(D) %in% c("race")]
# D = as.matrix(D)

# star_time <- Sys.time() 
# tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
# best_subgroup = tmp$best_subgroup
# l_naive = tmp$nonpri_naive
# result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("dplyr","MASS","rmutil","epiDisplay","tmvtnorm","matrixcalc","corpcor")) %dopar% {
#     unname(fun_linear_regression_parallel(D,Y,B,level,e,r,BS_true,k_par))
# }
# parallel_result = rowMeans(result)
# final_result = c(parallel_result,l_naive,best_subgroup)
# final_result

# write.csv(result,file = "result/race_age.csv")
# end_time <- Sys.time() 
# run_time <- end_time - star_time  ## 计算程序运行时间
# print("Time for the 100 epoch:")
# print(run_time)  # 52 sec


##################### group by age & race & wtkg #####################
n = dim(ACTG175)[1]
n_subgroup <- 3 
k_n <- 2^n_subgroup
k_par = k_n

print(colnames(X))
Z <- matrix(0,n,k_n)
for (i in 1:n){ 
    indicator_2base = rep(0,n_subgroup) 
    if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1 
    if (X$age[i] > median(X$age)) indicator_2base[2] = indicator_2base[2]+1  
    if (X$wtkg[i] > median(X$wtkg)) indicator_2base[3] = indicator_2base[3]+1  
    indicator_10base = indicator_2base[3]+indicator_2base[2]*2+indicator_2base[1]*4
    Z[i,indicator_10base+1] = 1
}
colnames(Z) <- c("G1","G2","G3","G4","G5","G6","G7","G8") 

DZ <- matrix(0,n,k_n)
for (i in 1:n){
    DZ[i,] = Ti[i]*Z[i,]
}
colnames(DZ) <- c("TG1","TG2","TG3","TG4","TG5","TG6","TG7","TG8") 

X = cbind(Z,X) 
# for (j in 1:dim(X)[2]){
#     X[,j] = (X[,j]-mean(X[,j]))/sd(X[,j])
# }

D = cbind(DZ,X) 
D <- D[, !colnames(D) %in% c("race")]

D = as.matrix(D)

# star_time <- Sys.time() 
# tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
# best_subgroup = tmp$best_subgroup
# l_naive = tmp$nonpri_naive
# result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("dplyr","MASS","rmutil","epiDisplay","tmvtnorm","matrixcalc","corpcor")) %dopar% {
#     unname(fun_linear_regression_parallel(D,Y,B,level,e,r,BS_true,k_par))
# }
# parallel_result = rowMeans(result)
# final_result = c(parallel_result,l_naive,best_subgroup)
# final_result

# write.csv(result,file = "result/race_age_wtkg.csv")
# end_time <- Sys.time() 
# run_time <- end_time - star_time  ## 计算程序运行时间
# print("Time for the 100 epoch:")
# print(run_time)
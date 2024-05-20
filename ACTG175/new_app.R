nCores = 8
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

set.seed(980616)

# group = c('symptom','hemo') #!
# group = c('age','race') #?
path_name = paste(group, collapse = '_')
log_file <- paste("newlog/", path_name, ".txt", sep = "")
sink(log_file, append = TRUE)
cat("\n\n\n")
start_time = Sys.time()
print("Start time is:")
print(start_time)
cat("with cross-validation \n")
# cat("with cd420 as response variable")
r = c(-Inf,1/30,1/15,1/10,1/5,0.5)
level = 0.2
e = 4000   #! update

cat("correction term:", r, "\n")
cat("confident level:", level, "\n")
cat("privacy budget:", e, "\n")


del_var <- c("pidnum", "arms", "cd40", "cd420", "cd496", "r", "cd80", "cens", "days")
# del_var <- c("pidnum", "arms", "cd40", "cd820", "cd496", "r", "cd80", "cens", "days")

ACTG175 <- na.omit(ACTG175[,!colnames(ACTG175) %in% del_var])
X <- ACTG175[, !colnames(ACTG175) %in% c("cd820","zprior","treat")]
Y <- ACTG175[,"cd820"]
# Y <- ACTG175[,"cd420"]
Ti<- ACTG175$treat

B = 200
n_r = length(r)
rep = 100


n = dim(ACTG175)[1]
n_subgroup <- length(group)
k_n <- 2^n_subgroup
k_par = k_n

Z <- matrix(0,n,k_n)
for (i in 1:n){ 
    indicator_2base = rep(0,n_subgroup) 
    # if (X$karnof[i] == 1) indicator_2base[1] = indicator_2base[1]+1 
    # if (X$wtkg[i] > median(X$wtkg)) indicator_2base[1] = indicator_2base[1]+1 
    # if (X$offtrt[i] == 0) indicator_2base[1] = indicator_2base[1]+1 #!
    # if (X$race[i] == 0) indicator_2base[1] = indicator_2base[1]+1  #?
    if (X$hemo[i] == 0) indicator_2base[1] = indicator_2base[1]+1 
    # if (X$preanti[i] > median(X$preanti)) indicator_2base[1] = indicator_2base[1]+1 


    # if (X$age[i] > median(X$age)) indicator_2base[2] = indicator_2base[2]+1   #?
    # if (X$age[i] > 26.3) indicator_2base[2] = indicator_2base[2]+1  
    if (X$symptom[i] == 0) indicator_2base[2] = indicator_2base[2]+1 #!
    # if (X$karnof[i] > median(X$karnof)) indicator_2base[2] = indicator_2base[2]+1
    indicator_10base = indicator_2base[2]+indicator_2base[1]*2
    Z[i,indicator_10base+1] = 1
}
colnames(Z) <- c("G1","G2","G3","G4") 

# cat("X$offtrt[i] == 0 \n") #!
cat("X$symptom[i] == 0 \n") #!

# cat("X$age[i] > median(X$age)", median(X$age), "\n") #?
# cat("X$race[i] == 0 \n") #?
# cat("X$26.3[i] > 40 \n")
cat("X$hemo[i] == 0 \n")
# cat("X$preanti[i] > median(X$preanti)", median(X$preanti), "\n") 
# cat("X$wtkg[i] > median(X$wtkg)", median(X$wtkg), "\n")
# cat("age setline: median", median(X$age), "\n")

DZ <- matrix(0,n,k_n)
for (i in 1:n){
    DZ[i,] = Ti[i]*Z[i,]
}
colnames(DZ) <- c("TG1","TG2","TG3","TG4") 

X = cbind(Z,X) 
D = cbind(DZ,X) 

D <- D[, !colnames(D) %in% c("symptom")]  #!
# D <- D[, !colnames(D) %in% c("offtrt")]  #!
# D <- D[, !colnames(D) %in% c("race")]  #?
D <- D[, !colnames(D) %in% c("hemo")]  #!revised
# D <- D[, !colnames(D) %in% c("G4")]  #!revised
# dim(D)[2]
# qr(D)$rank
# cat("remove symptom, offtrt \n") #!
cat("remove symptom,hemo \n")

# cat("remove race") #?

normalize_min_max <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
D$wtkg <- normalize_min_max(D$wtkg)
D$age <- normalize_min_max(D$age)
D$karnof <- normalize_min_max(D$karnof)
D$preanti <- normalize_min_max(D$preanti)
# print(D)
D = as.matrix(D)
tmp = fun_linear_regression(D,Y,B,level,e,r,BS_true,k_par)
best_subgroup = tmp$best_subgroup
l_naive = tmp$nonpri_naive
k_fold = 3

result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("dplyr","MASS","rmutil","epiDisplay","tmvtnorm","matrixcalc","corpcor")) %dopar% {
    unname(fun_linear_regression_parallel(D,Y,B,level,e,r,BS_true,k_par,k_fold))
}
parallel_result = rowMeans(result)
final_result = c(parallel_result,l_naive,best_subgroup)
# final_result
sprintf("%.3f", final_result)
final_result

# write.csv(result,file = paste("result/", path_name, ".csv", sep = ""))
end_time <- Sys.time() 
run_time <- end_time - start_time  ## 计算程序运行时间
print("Time for the 100 epoch:")
print(run_time)
sink()
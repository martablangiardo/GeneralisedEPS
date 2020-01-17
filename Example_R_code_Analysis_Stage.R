Packages <- c("copula","plyr","reshape2","sp","abind","R2WinBUGS","rgdal","spdep")
lapply(Packages, library, character.only = TRUE)


setwd("C:/Simul/MAR/Imput_Health") # set working directory
load("DataGen_MAR_EPS.RData")


######################
### ANALYSIS STAGE ###
######################

EPS.mis <- array(NA,dim=list(numareas,nsim))
for(j in 1:nsim){
  EPS.mis[IndMis==0,j] <- eps.mean[,j] # see Figure 2 of the paper
}

eps.mean.new <- array(0,dim=list(numareas,nsim))
for(j in 1:nsim){
  eps.mean.new[IndMis==0,j] <- eps.mean[,j]
}
eps.mean.new[,1]

eps.prec.new <- array(NA,dim=list(numareas,nsim))
for(j in 1:nsim){
  eps.prec.new[IndMis==0,j] <- 1/eps.var[,j]
}

for(j in 1:nsim){
  eps.prec.new[is.na(eps.prec.new[,j])]<- 1/mean(eps.var[,j])
}
eps.prec.new[,1]

# NOTE: for eps.mean and eps.prec we used a trick, as the loop goes from 1 to the total number of areas=625.
# In particular, for the eps.mean we set the values of the out-of-sample areas equal to 0, and for eps.prec we set the 
# values of the out-of-sample areas equal to inverse of the mean of the posterior variance. 
# Importantly, this does not affect the computation, as for the out-of-sample areas the
# value used at each iteration for the generalized EPS is predicted by the imputation model. 
# This is reached by the use of the step functions in the model.


### Generate terms for ICAR ###

London.map <- readOGR(dsn = ".", layer = "London_Ward_CityMerged")
coords.London <- coordinates(London.map)
ID <- row.names(coords.London)

London.nb <- knn2nb(knearneigh(coords.London, k = 4), row.names = ID)

dsts <- unlist(nbdists(London.nb, coords.London))
summary(dsts)
max.dst.L <- max(dsts)

London.D.nb <- dnearneigh(coords.London, d1 = 0, d2 = 1 * max.dst.L, row.names = ID) 
London.D.nb

numL <- sapply(London.D.nb, length)
adjL <- unlist(London.D.nb)
allparL <- listw2WB(nb2listw(London.D.nb, style="B"))
weightsL <- dput(allparL$weights, control=NULL)

### Generate inits ###

initsEpsMis1 <- array(NA,dim=list(numareas,nsim))
for(j in 1: nsim){
  initsEpsMis1[IndMis==1,] <- c(runif(309, min=min(eps.mean[,j]), max=max(eps.mean[,j])))
}

initsEpsMis1[,1]

set.seed(3570)
initsEpsMis2 <- array(NA,dim=list(numareas,nsim))
for(j in 1: nsim){
  initsEpsMis2[IndMis==1,] <- c(runif(309, min=min(eps.mean[,j]), max=max(eps.mean[,j])))
}

initsEpsMis2[,1]

### run model ###

setwd("C:/Simul/MAR/Imput_Health/BETA") # set directory to store the output

for (j in 1:nsim){
  print(j) 
  dataIM<-list("numareas"=numareas, "EPS.mis"=EPS.mis[,j], 
               "eps.mean"=eps.mean.new[,j],"eps.prec"=eps.prec.new[,j],
               "X"=X[,j],"C"=C[,j], "adjSp"=adjL, "numSp"=numL, 
               "weights"=weightsL, "E"=expected, "Y"=Y[,j],
               "ind.mis"=IndMis) 
 
  inits <- list(
    list(EPS.mis=initsEpsMis1[,j], 
         eta=0.1,  gamma1=0.5, gamma2=0.3, tau.eps=30, 
         theta=rep(-1,numareas), 
         beta0=0.1, beta1=0.1, beta2=0.1, beta3=0.1, 
         sigma=0.8, V=runif(625,-1, 1)),
    list(EPS.mis=initsEpsMis2[,j], 
         eta=0.2,  gamma1=0.3, gamma2=0.1, tau.eps=40,
         theta=rep(1,numareas), 
         beta0=0.3, beta1=0.2, beta2=0.3, beta3=0.3, 
         sigma=0.25, V=runif(625,-1, 1)))
  
  
  parameters=c("beta0", "beta1", "beta2","beta3","eta","gamma1","gamma2")
  
  modelsimIm <- bugs(dataIM, inits, parameters, model.file="Analysis_Stage.txt",
                     n.chains = 2, n.iter = 55000, n.burnin = 50000, n.thin = 5,
                     debug=FALSE, bugs.directory="C:/WinBUGS14",  DIC = FALSE,
                     working.directory =getwd())
  
  write.csv(modelsimIm$sims.matrix,paste0("para",j,".csv"))
}

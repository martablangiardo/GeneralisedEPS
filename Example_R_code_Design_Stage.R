Packages <- c("copula","plyr","reshape2","sp","abind","R2WinBUGS","rgdal","spdep")
lapply(Packages, library, character.only = TRUE)


setwd("C:/Simul/MAR/UpScaling_EPSgen") # set working directory
load("DataSimul_MAR.RData")


################
# DESIGN STAGE # 
################

table(IndMis) # 0=in-sample areas, 1=out-of-sample areas

K1 <- 3 # continuous individual-level variables
K2 <- 5 # total number of individual-level variables
numvars <- 5 

corr_mat <- matrix(c(1,0.7,0.3,0.2,0.2, 0.7,1,0.2,0.1,0.1, 0.3,0.2,1,0.4,0.4, 0.2,0.1,0.4,1,0.5, 0.2,0.1,0.4,0.5,1), nrow=numvars, ncol=numvars, byrow=TRUE)
B <- solve(corr_mat)

C.Nomis  <- array(NA,dim=list(numareasNOmis, nsim)) # subset of C (for in-sample areas)
for(j in 1:nsim){
  C.Nomis[,j] <- C[which(IndMis==0),j]
}

X.Nomis <- array(NA,dim=list(numareasNOmis, nsim)) # subset of X (for in-sample areas)
for(j in 1:nsim){
  X.Nomis[,j] <- X[which(IndMis==0),j]
}

dim(datm.Nomis) 
datm.Nomis.forBugs <- aperm(datm.Nomis,c(3,1,2,4)) 
dim(datm.Nomis.forBugs)

### Generate terms for Multivariate ICAR ###

SubLondon.map <-readOGR(dsn = ".", layer = "subsetLondon")
length(SubLondon.map)

# extract the coordinates and row names 
coords.SubLondon <- coordinates(SubLondon.map)
IDs <- row.names(coords.SubLondon)

SubLondon.nb <- knn2nb(knearneigh(coords.SubLondon, k=4), row.names = IDs)
dsts1 <- unlist(nbdists(SubLondon.nb, coords.SubLondon))
summary(dsts1)
max.dst <- max(dsts1)

# distnce-based neighbours

SubLondon3.nb <- dnearneigh(coords.SubLondon, d1 = 0, d2 = 1 * max.dst, row.names = IDs) 
numSub3 <- sapply(SubLondon3.nb, length)
adjSub3 <- unlist(SubLondon3.nb)
allpar3 <- listw2WB(nb2listw(SubLondon3.nb, style="B")) #listw2WB crates a list of weights for winbugs
weights3 <- dput(allpar3$weights, control=NULL)

### Generate inits ###

prec.xi <- matrix(c(1.44,-2.66,0.73,-0.54,-3.32, -2.66,6.75,-1.93,-0.40,7.16, 0.73,-1.93,0.63,-0.02,-1.74, -0.54,-0.40,-0.02,10.05,-1.22, -3.32,7.16,-1.74,-1.22,16.20), nrow=5, ncol=5, byrow=TRUE)
prec.xi2 <- matrix(c(1.83,-4.24,0.90,2.46,-3.47,-4.24,12.20,-2.65,-6.43,10.28,0.90,-2.65,0.71,0.88,-1.20,2.46,-6.43,0.88,14.77,-7.83,-3.47,10.28,-1.20,-7.83,20.57), nrow=5, ncol=5, byrow=TRUE)

matxi <- rep(0,1580) #1580=316(in-sample areas)*5(nb individual-level variables)
matxi <- matrix(matxi, nrow=numvars, ncol=numareasNOmis)
matxib <- rep(0,1580)
matxib <- matrix(matxib, nrow=numvars, ncol=numareasNOmis)
dc <- 0.1
dc1 <- 0.2
dm <- c(0.1, 0.1, 0.1, 0.1, 0.1) 
dm1 <- c(0.2, 0.2, 0.2, 0.3, 0.3)

### run model ###

setwd("C:/Simul/MAR/UpScaling_EPSgen/EPS") #set directory to store the output

for (j in 1:100){
  print(j)
  datamodel <- list("numareasNOmis"=numareasNOmis, "K1"=K1, "K2"=K2, "numvars"=numvars, 
                 "numind"=numind, "m"=round(datm.Nomis.forBugs[,,,j],5),
                 "adj"=adjSub3, "num"=numSub3, "weights"=weights3,
                 "B"=B, "X"=X.Nomis[,j],"C"=C.Nomis[,j])
  
  inits <- list(
    list(delta0=0.1,  nu=c(0.1, -0.1, 0.1, -0.1, 0.1), prec.xi=prec.xi, xi=matxi, 
         deltac=dc, deltam=dm, tauX=5, tau=c(10.06, 4.73, 2.01)),
    list(delta0=0.05,  nu=c(0.3, -0.3, 0.3, -0.3, 0.3),  prec.xi=prec.xi2, xi=matxib, 
         deltac=dc1, deltam=dm1, tauX=10, tau=c(8, 5, 4)))
  
  parameters <- c("eps","delta0","deltac","deltam","tauX")
  
  modelsim <- bugs(datamodel, inits, parameters, model.file="Design_Stage.txt",
                   n.chains = 2, n.iter = 50000, n.burnin = 45000, n.thin =5,
                   debug=FALSE, bugs.directory="C:/WinBUGS14", DIC=FALSE,
                   working.directory =getwd())
   
  write.csv(modelsim$sims.matrix,paste0("epspar",j,".csv"))
}

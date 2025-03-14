################################################################################
###################### Loading the packages ########################################
################################################################################
library(dplyr)
library(Matrix)
library(splines)
library(splines2)
library(optimParallel)
library(MASS)
library(mgcv)
library(scam)
library(scales)
library(invgamma)
library(pracma)
library(parallel)
library(purrr)
library("fastmatrix")
library("DescTools")
library(plyr)
library("LaplacesDemon")
library("pracma")
library(expm)
library(fBasics)
library("interp")
library(ggplot2)



###################################################################################
dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(dir)
data = read.csv("DataFinal.csv")

data$sexo = replace(data$sexo,data$sexo == "0","Female")
data$sexo = replace(data$sexo,data$sexo == "1","Male")

print(ggplot(data = data, aes(x = time, y = y,group  = id))+
        facet_grid(sexo ~ Age)+ 
        geom_line(lty= 2)+
        labs(y = "Mortality rate",x = "Time",cex=1.2)+ 
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=16),
                                                                     axis.title=element_text(size=16,face="bold"), strip.text = element_text(
                                                                       size = 16)))



###################### Preparing the data
############################################
data$sexo = factor(data$sexo)
data$Age = factor(data$Age)
data$yscale = (scale(data$y))
data$time2 <- rescale(data$time, to = c(0, 1))
datafinal = data
############################################
############################################

############################################
#### Creating the trasnformation functions


############################################
############################################
############################################
### Transformation in y
sspec <- s(yscale, k =10 , bs = "ps")
sspec$mono <- 1
sm2 <- smoothCon(sspec, data = data, scale.penalty=FALSE, absorb.cons = T)[[1]]

A1 = sm2$X


(sm2$deriv <- 1)
Ad1 <- PredictMat(sm2, data)
index1 = c(sm2$g.index)
K1 = sm2$S[[1]]


#### Product between id and time

Xre = model.matrix(~-1+as.factor(id),data = data)

data$Time2 = scale(data$time)
yscale = data$yscale
tensor_cat = function(Time2,q,Xc,data){
  sspec = s(Time2, k = q, bs = "cs")
  sm = smoothCon(sspec, data = data, scale.penalty = FALSE, absorb.cons = F)[[1]]
  q1 = q 
  q2 = ncol(Xc)
  n <- length(yscale)
  
  X <- matrix(0, n, q1 * q2)
  for (i in 1:n) {
    X[i, ] <- sm$X[i, ] %x% Xc[i,]
  }
  
  #D <- diag(q1 * q2)
  #D <- D[, -q2]
  #D1 <- t(diff(diag(q2)))
  #D[1:q2, 1:(q2 - 1)] <- D1
  #X <- X %*% D
  return(list(X = X, knots = sm$knots, K = sm$S[[1]]))
}

sm = tensor_cat(data$Time2,5,Xre,data)
A2 = sm$X
dim(A2)
index2 = rep(F, ncol(A2))

Ad2 = matrix(0, nrow(A2), ncol(A2))
create_bspline_matrix <- function(n) {
  # Inicializa uma matriz de zeros
  matrix <- matrix(0, n, n)
  
  # Preenche a diagonal principal com 2
  diag(matrix) <- 2
  
  # Preenche as diagonais superiores e inferiores com -1
  for (i in 1:(n-1)) {
    matrix[i, i+1] <- -1
    matrix[i+1, i] <- -1
  }
  
  # Ajusta os valores das extremidades conforme o exemplo dado
  matrix[1, 1] <- 1
  matrix[n, n] <- 1
  
  return(matrix)
}

# Criando a matriz 5x5
K_dim <- 5
K_matrix <- create_bspline_matrix(K_dim)
K_matrix = K_matrix%*%K_matrix
K12 = K_matrix%x%diag(ncol(Xre))
q = 6
K22 = (diag(1,q-1)%x%diag(10^(-6),ncol(Xre)))

######### Linear effects

X = model.matrix(~ 1  + Time2*sexo + Time2*Age+ Age+sexo, data = data)
head(X)
X = X[,c(-1,-2)]
#X = X[,-1]
Xre = model.matrix(~as.factor(id),data = data)
X = cbind(X)
Xi0 = matrix(0, ncol = ncol(X),nrow = nrow(X))



###### Aggregating all
A = cbind(A1,A2, X)
Ad = cbind(Ad1,Ad2,Xi0)

indext = c(index1,index2,rep(FALSE,ncol(X)))
K2 = matrix(0,ncol(K1),nrow(K1))
K3 = K12 + K22 

y = data$y

############################################
############################################
############################################
#################### Fitting the model
library(Rcpp)
library("RcppParallel")


sourceCpp("funcs.cpp")
source("Model.R")

system.time({m12 = bctm(A,Ad,X,K1,K2,K3,q,indext)})


####### Getting the results

post_multi = list()
points_t = m12[[1]]$sigma
maximos = m12[[2]]$maximos
mu_ini = m12[[2]]$mu_ini
sigma_ini = m12[[2]]$sigma_ini
lista = as.list(as.data.frame(t(points_t)))
weight = c(abs(diff(points_t[,2])[1])*abs(diff(points_t[,1],6)[1]))
unl = (sapply(m12, '[[', 'unn_log_post')) - mean(sapply(m12, `[[`, 'unn_log_post'))

#############Ecnontrar o maximo apenas uma vez e replicar os pontos
xis_c2 = function(index,max){
  if(index %in% which(indext == FALSE)){
    return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
  }else{
    return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
  }
}

post_alpha  = exp(unl)/(sum(exp(unl)*weight)) #Normalizing the hyperpetameter posterior
lpost = post_alpha
lpost2 <- cut((post_alpha), 10, label = FALSE)
colr <- rev(heat.colors(10))


### Plot of hyperparameter posterior
plot(points_t[,1] ~ points_t[,2],
     col = colr[lpost2], pch = 20)


################################################################
#######################################################################
################################ Getting the posterior marginals of parameters
pontos = points_t
lpost2_norm = lpost

for(j in 1:length(pontos[,1])){
  post_multi[[j]] = matrix(unlist(m12[[j]]$posti_cond),nrow = 31,ncol = length(indext))
}

marginal = function(index,pila){
  pila_marginal = vector()
  xis = xis_c2(index,maximos)
  mar = function(j){
    return(sum(sapply(lapply(1:(length(lista)), function(i) pila[[i]][,index]),'[[',j)*(lpost*weight)))
  }
  pila_marginal = unlist(map(1:length(xis), mar))
  return(list(pila_marginal = as.vector(pila_marginal), xis = xis))
}

marginals2 = lapply(1:length(mu_ini),function(ii) marginal(ii,post_multi))

####Marginals

for(i in 1:length(indext)){

  plot(marginal(i,post_multi)$xis,marginal(i,post_multi)$pila_marginal,type = "l", col = "blue")

}


################################################################
#######################################################################
################################ Simulating values from marginals of parameters

dens_est = function(x,norm,xis){
  if(x > max(xis) || x< min(xis)){
    return(0)
  }else{
    ind1 = which(abs(xis - x) == sort((abs(xis - x)))[1])
    ind2 = which(abs(xis - x) == sort((abs(xis - x)))[2])
    pesos = abs(x - xis[c(ind1,ind2)])/(sum(abs(x - xis[c(ind1,ind2)])))
    return(interp1(x = xis, y = norm,xi = x,method = "nearest"))}
}

MH = function(index){
  dens_norm = marginal(index,post_multi)
  norm = dens_norm$pila_marginal
  xis = dens_norm$xis
  x = rep(0,5000)
  x[1] = xis[norm == max(norm)]     #initialize; I've set arbitrarily set this to 3
  for(i in 2:5000){
    current_x = x[i-1]
    proposed_x = current_x + rnorm(1,mean=0,sd=1)
    #proposed_x = runif(1,current_x -0.2,current_x+0.2) 
    A = min(1,dens_est(proposed_x,norm,xis)/dens_est(current_x,norm,xis))
    #A = interp1(x = xis,y = norm, xi = proposed_x,method = "spline")/interp1(x = xis,y = norm, xi = current_x,method = "spline")
    if(runif(1)<A){
      x[i] = proposed_x       # accept move with probabily min(1,A)
    } else {
      x[i] = current_x        # otherwise "reject" move, and stay where we are
    }
  }
  
  return(x[seq(1500,5000,2)])
}



Beta = mclapply(seq(1,length(mu_ini),1),MH,mc.cores = 6)
Beta = matrix(unlist(Beta),ncol = length(mu_ini),nrow = length(Beta[[1]]))
dim(A)

gamma_f = function(beta,index){
  gammat = vector(length = length(index))
  gammat[index] = exp(beta[index])
  gammat[!index] = beta[!index]
  return(gammat)
}


Gamma = apply(Beta,1,gamma_f,index = indext)

library(coda)
## Table of estimates and HPD
library(xtable)
dim(A)
head(A)
xtable(cbind(apply(Gamma,1,mean)[145:dim(A)[2]],HPDinterval(as.mcmc(t(Gamma)))[145:dim(A)[2],]),digits = 4)
head(X)
apply(Gamma,1,mean)

###############################################################################
##################### Estimating the conditional densities and distributions
head(A)
sspec <- s(yscale, k = 10 , bs = "ps")
sspec$mono <- 1
ynew = seq(min(yscale)-5,max(yscale)+5,0.01)
ynew_o = ynew*(attr(yscale, "scaled:scale")) + attr(yscale,"scaled:center")
dfnew = as.data.frame(cbind(ynew))

names(dfnew) = c("yscale")
sm2 <- smoothCon(sspec, data = data, scale.penalty=FALSE, absorb.cons = T)[[1]]
B =PredictMat(sm2, dfnew)

(sm2$deriv <- 1)
Bl <- PredictMat(sm2, dfnew)
Xm = cbind(A2,X)
dim(Xm)

dens_estim= function(ind,B,Bl,Gamma,k){
  
  gamma = t(Gamma)
  
  
  h1 = B%*%(t(gamma[,1:(k - 1)])) 
  
  h1l = Bl%*%(t(gamma[,1:(k - 1)])) 
  
  h2 = Xm[ind,]%*%t(gamma[,(k):(ncol(gamma))])
  
  h1 = apply(h1,1,mean)
  
  h1l = apply(h1l,1,mean)
  
  h2 = mean(h2)
  
  
  f = dnorm(h1 + h2)*(h1l)
  
  plot(ynew,f,type = "l",col = "blue")
  return(f)
  
}


prob_estim = function(ind,Gamma,B,Bl,k){
  gamma = t(Gamma)
  
  h1 = B%*%(t(gamma[,1:(k - 1)])) 
  
  h1l = Bl%*%(t(gamma[,1:(k - 1)])) 
  
  h2 = Xm[ind,]%*%t(gamma[,(k):(ncol(gamma))])
  
  h1 = apply(h1,1,mean)
  
  h1l = apply(h1l,1,mean)
  
  h2 = mean(h2)
  
  Fa = pnorm(h1 + h2)
  
  Faord = sort(Fa)
  
  qi = (ynew_o[Faord > 0.025])[1]
  qs = (ynew_o[Faord > 0.975])[1]
  qm = (ynew_o[Faord > 0.5])[1]
  
  
  
  return(list(qi = qi,  qs = qs, qm = qm))
  
  
}

Estim = lapply(seq(1:(dim(data)[1])),prob_estim, Gamma = Gamma, B = B, Bl = Bl,10)

qm = sapply(Estim, "[[","qm")
qi = sapply(Estim, "[[","qi")
qs = sapply(Estim, "[[","qs")

datapred = as.data.frame(cbind(data$y,qm,qi,qs))
datafinal = as.data.frame(cbind(data,datapred))

datapred = datapred[order(datapred$qm),]

par(mar = c(5, 7, 4, 2)) 
plot(datapred$qm,type = "l", col =  "darkblue",lwd = 2,lty =2, cex.axis = 2.5,cex.lab = 2.5, ylab = "Mortality rate")
points(datapred$qs, type = "l",col = "darkblue",lwd = 2)
points(datapred$qi, type = "l",col = "darkblue",lwd = 2)
points(datapred$V1,pch = 20)
sum(datapred$V1 < datapred$qs & datapred$V1 > datapred$qi)/length(datapred$V1)

i= 10
#Maranhão

data1 = datafinal[(datafinal$id == i) & (datafinal$sexo == "Male"),]
data1$MHDI
print(ggplot(data = data1, aes(x = time, y = qm))+
        facet_wrap(~Age )+ 
        geom_line(col= "black")+ labs(x = "Time",cex=1.5)+ 
        labs(y = "Mortality rate",cex=1.5)+
        geom_point(data=data1,aes(x=time,y=y),pch=1, colour = "black", fill = "black") + 
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=17),
                                                                     axis.title=element_text(size=17,face="bold"),
                                                                     strip.text = element_text(size = 17, face = "bold"))+ 
                           geom_ribbon(aes(ymin = qi, ymax = qs),alpha = 0.2))



#Rio de Janeiro
i = 21
data1 = datafinal[(datafinal$id == i) & (datafinal$sexo == "Male"),]

print(ggplot(data = data1, aes(x = time, y = qm))+
        facet_wrap(~Age )+ 
        geom_line(col= "black")+ labs(x = "Time",cex=1.5)+ 
        labs(y = "Mortality rate",cex=1.5)+
        geom_point(data=data1,aes(x=time,y=y),pch=1, colour = "black", fill = "black") + 
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=16),
                                                                     axis.title=element_text(size=16,face="bold"),
      strip.text = element_text(size = 17, face = "bold"))+ 
  geom_ribbon(aes(ymin = qi, ymax = qs),alpha = 0.2))

################################################################################
############ Producing the graphs


estados = read.csv("DataFinal2.csv", header=TRUE, sep = ";")

estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Cear�"),"Ceará")
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Amap�"),"Amapá")
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Esp�rito Santo"),"Espírito Santo")     
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Goi�s"),"Goiás") 
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Maranh�o"),"Maranhão") 
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Paran�"),"Paraná") 
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Para�ba"),"Paraíba") 
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Par�"),"Pará") 
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Piau�"),"Piauí")
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "Rond�nia"),"Rondônia")
estados$Estado = replace(estados$Estado,which(estados$Estado ==  "S�o Paulo"),"São Paulo")


tensor_cat_deriv = function(Time2,q,Xc,data){
  sspec = s(Time2, k = q, bs = "cs")
  sm = smoothCon(sspec, data = data, scale.penalty = FALSE, absorb.cons = F)[[1]]
  sm$deriv = 1
  q1 = q 
  q2 = ncol(Xc)
  n <- length(yscale)
  dB <- PredictMat(sm, data)  # Derivada da base de splines
  
  X <- matrix(0, n, q1 * q2)
  for (i in 1:n) {
    X[i, ] <- dB[i, ] %x% Xc[i,]
  }
  
  #D <- diag(q1 * q2)
  #D <- D[, -q2]
  #D1 <- t(diff(diag(q2)))
  #D[1:q2, 1:(q2 - 1)] <- D1
  #X <- X %*% D
  return(list(X = X, knots = sm$knots, K = sm$S[[1]]))
}
data

smd = tensor_cat_deriv(data$Time2,5,Xre,data)
A2d = smd$X
A2d
datay = as.data.frame(ynew)

sspec <- s(yscale, k =10 , bs = "ps")
sspec$mono <- 1
names(datay) = "yscale"
sm2 <- smoothCon(sspec, data = datay, scale.penalty=FALSE, absorb.cons = T)[[1]]
(sm2$deriv <- 1)
A1d <- PredictMat(sm2, data)


gamma = apply(Gamma,1,mean)
data$idf = as.factor(data$id)
head(estados)
estados$Estado = factor(estados$Estado, levels = unique(estados$Estado))
Xtime = data.frame(cbind(-A2%*%t(gamma)[,(dim(A1)[2]+1):(dim(A1)[2]+ 27+ 27*4)],data$time, data$yscale,estados$Estado))
#Xtime = data.frame(cbind(-A2%*%t(gamma)[,(dim(A1)[2]+1):(dim(A1)[2]+ 27*4)],data$idf,data$time, data$yscale,data$id))
names(Xtime)= c("ytime","Time","yscale","Estados")
#names(Xtime)= c("ytime","id","Time","yscale")

Xtime$ytime= as.numeric(Xtime$ytime)
Xtime$Time = as.numeric(Xtime$Time)
#Xtime$ytime = (A2d%*%t(gamma)[,(dim(A1)[2]+1):(dim(A1)[2]+ 27+ 27*4)])/(A1d%*%t(gamma)[,1:9])
Xtime$Estados = as.factor(Xtime$Estados)

Xtime$Estados= factor(estados$Estado, levels = unique(estados$Estado))
library(ggplot2)

head(Xtime)
Xtime$id
print(ggplot(data = Xtime, aes(x = Time, y = ytime))+
        facet_wrap(~Estados )+
        geom_line(col= "black")+ labs(x = "Time", y = "",cex=1.2)+
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=12),
                                                                     axis.title=element_text(size=14,face="bold")))




dens_estim= function(ind,B,Bl,Gamma,k){
  
  gamma = t(Gamma)
  
  
  h1 = B%*%(t(gamma[,1:(k - 1)])) 
  
  h1l = Bl%*%(t(gamma[,1:(k - 1)])) 
  
  h2 = Xm[ind,]%*%t(gamma[,(k):(ncol(gamma))])
  
  h1 = apply(h1,1,mean)
  
  h1l = apply(h1l,1,mean)
  
  h2 = mean(h2)
  
  
  f = dnorm(h1 + h2)*(h1l)*(1/attr(data$yscale,"scaled:scale"))
  
  return(f)
  
  
}

####### Comparing different times
inds = data[(data$time == 1) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 4),]$X
estados[inds,]
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,550))

inds = data[(data$time == 4) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 4),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lty = 5,lwd = 2)

inds = data[(data$time == 8) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 4),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l",lty = 9,lwd= 2)
legend("topright", lty = c(1,5,9), legend = c("Time = 1", "Time = 4", "Time = 8"), cex = 1.4)

#########
inds = data[(data$time == 1) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 6),]$X
estados[inds,]
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,550))

inds = data[(data$time == 4) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 6),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lty = 5,lwd = 2)

inds = data[(data$time == 8) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 6),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l",lty = 9,lwd= 2)
legend("topright", lty = c(1,5,9), legend = c("Time = 1", "Time = 4", "Time = 8"), cex = 1.4)



####### Comparing different Sex

inds = data[(data$time == 1) & (data$Age3 == 1) & (data$sexo == "Female") & (data$id == 16),]$X
estados[inds,]#which state
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,200))

inds = data[(data$time == 1) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 16),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 5)
legend("topright", lty = c(1,5),lwd = c(2,2), legend = c("Female", "Male"), cex = 1.4)





##### Comparing different Ages
inds = data[(data$time == 4) & (data$Age1 == 0)  & (data$Age2 == 0)  & (data$Age3 == 0)  & (data$sexo == "Male") & (data$id == 18),]$X
estados[inds,]#which state
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,200))

inds = data[(data$time == 4) & (data$Age1 == 1) & (data$sexo == "Male") & (data$id == 18),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 5)

inds = data[(data$time == 4) & (data$Age2 == 1) & (data$sexo == "Male") & (data$id == 18),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 9)

inds = data[(data$time == 4) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 18),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 4)
legend("topright", lty = c(1,5,9,4), legend = c("50-59", "60-69", "70-79", ">80"), cex = 1.4)



####### Time*Age
inds = data[(data$time == 4) & (data$Age1 == 1) & (data$sexo == "Male") & (data$id == 15),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,300))

inds = data[(data$time == 4) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 15),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 5)

inds = data[(data$time == 7) & (data$Age1 == 1) & (data$sexo == "Male") & (data$id == 15),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 9,lwd = 2, xlim = c(0,200))

inds = data[(data$time == 7) & (data$Age3 == 1) & (data$sexo == "Male") & (data$id == 15),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 4,lwd = 2, xlim = c(0,200))
legend("topright", lty = c(1,5,9,4), legend = c("Age = 60-59 and Time = 4", "Age > 80 and Time = 4","Age = 60-59 and Time = 7", "Age > 80 and Time = 7"), cex = 1.4)


####### Time*Sex
inds = data[(data$time == 1) & (data$Age1 == 0)  & (data$Age2 == 0)  & (data$Age3 == 0)& (data$sexo == "Female") & (data$id == 4),]$X
estados[inds,]
f = dens_estim(inds,B,Bl,Gamma,10)
plot(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 1,lwd = 2, xlim = c(0,90))

inds = data[(data$time == 1) & (data$Age1 == 0)  & (data$Age2 == 0)  & (data$Age3 == 0)& (data$sexo == "Male") & (data$id == 4),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", lwd = 2, lty = 5)

inds = data[(data$time == 7) & (data$Age1 == 0)  & (data$Age2 == 0)  & (data$Age3 == 0) & (data$sexo == "Female") & (data$id == 4),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 9,lwd = 2, xlim = c(0,200))

inds = data[(data$time == 7) & (data$Age1 == 0)  & (data$Age2 == 0)  & (data$Age3 == 0)& (data$sexo == "Male") & (data$id == 4),]$X
f = dens_estim(inds,B,Bl,Gamma,10)
points(ynew_o,f, type= "l", xlab="Mortality rate", ylab = "", cex.axis = 1.6, cex.sub= 1.6, cex.lab = 1.6, lty = 4,lwd = 2, xlim = c(0,200))
legend("topright", lty = c(1,5,9,4), legend = c("Female and Time = 4", "Male and Time = 4","Female and Time = 7", "Male and Time = 7"), cex = 1.4)




################################################
######## Producing figures ...
ey = vector()
sd = vector()
skew = vector()


for(i in 1:dim(data)[1]){
  f = dens_estim(i,B,Bl,Gamma,10)
  ey[i] = sum(diff(ynew_o)[1]*f*ynew_o)
  sd[i] = sqrt(sum(diff(ynew_o)[1]*f*(ynew_o^2)) - (sum(diff(ynew_o)[1]*(f)*ynew_o))^2)
  skew[i] = sum(diff(ynew_o)[1]*f*(((ynew_o - ey[i])/sd[i])^3))
  print(i)
}

data$EY = ey
data$VY = sd
data$SKEW = skew
id_labels <- c("1" = "Acre", "20" = "Rio Grande do Sul")  # Replace with your actual IDs and labels

data_1 = data[data$id == 1 | data$id == 20,]
data_1$id = as.factor(data_1$id)
library(ggplot2)
print(ggplot(data = data_1, aes(x = time, y = VY,group = id,colour = id))+
        facet_grid(sexo ~ Age)+ 
        geom_line(lty= 4)+
        labs(y = "Standard Deviation",x = "Time", colour = "State",cex=1.2)+ 
        scale_colour_manual(values = c("1" = "black", "20" = "red"), labels = id_labels) + 
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=16),
                                                                     axis.title=element_text(size=16,face="bold"), strip.text = element_text(
                                                                       size = 16),legend.text = element_text(size = 17),legend.title = element_text(size = 19)))

print(ggplot(data = data_1, aes(x = time, y = SKEW,group = id,colour = id))+
        facet_grid(sexo ~ Age)+ 
        geom_line(lty= 4)+
        labs(y = "Skewness",x = "Time", colour = "State",cex=1.2)+ 
        scale_colour_manual(values = c("1" = "black", "20" = "red"), labels = id_labels) + 
        theme_bw() +theme(axis.line = element_line(colour = "black"),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())+ theme(axis.text=element_text(size=16),
                                                                     axis.title=element_text(size=16,face="bold"), strip.text = element_text(
                                                                       size = 16),legend.text = element_text(size = 17),legend.title = element_text(size = 19)))




####################
library(Rcpp)
library(splines2)
library(LaplacesDemon)
library(Matrix)
library(parallel)
library(dplyr)
library(ggplot2)
##################
dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(dir)
df = read.csv("datafinal.csv")
sourceCpp("funcs2_3.cpp")
dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(dir)
source("Allmodels.R")
source("func_aux.R")
sourceCpp("arg_min.cpp")
sourceCpp("gradiente.cpp")
############# Modelo 2 ################################################

### Transformacao em y
df_backup = df
df = df[1:floor(dim(df)[1]*0.9),]
#### Centrar e na centrar da a mesma coisa. Com ou sem intercepto tmb
df = as.data.frame(df)
A <- bernsteinPoly(df$case_count, degree = 20, intercept = FALSE,Boundary.knots = c(1,max(df$case_count)+10),data = df)
dim(A)
#umv = rep(1,dim(A)[1])
Ad = deriv(A)
#Xi0 = matrix(0,nrow(A), ncol(X))
q = ncol(A)
Sigma <- matrix(1,q,q)                             # Duplicate matrix
Sigma[upper.tri(Sigma)] <- 0          # Change lower triangular part

#Nao zerar
cmxy <- colMeans(A)
A <- sweep(A,2,cmxy)
A = A%*%Sigma
Ad = Ad%*%Sigma


index1 = c(rep(TRUE,q))
P <- diff(diag(q), differences = 1)
S <- crossprod(P)
P <- diff(diag(q), differences = 3)
S <- crossprod(P)
K1 = S
num_basis1 = q
Dp1 = diag(10^(-6),ncol(K1),ncol(K1))

dim(A)
dim(Ad)
K1 = K1 + Dp1
dim(K1)


############### Splines
head(df)
summary(df$dia)
Xyear <- bernsteinPoly(df$year, degree = 5, intercept = FALSE,Boundary.knots = c(1,7),data = df)
cmxYear <- colMeans(Xyear)
Xyear <- sweep(Xyear,2,cmxYear)
Xmes <- bernsteinPoly(df$mes, degree = 5, intercept = FALSE,Boundary.knots = c(1,12),data = df)
cmxMes <- colMeans(Xmes)
Xmes <- sweep(Xmes,2,cmxMes)
Xdia <- bernsteinPoly(df$dia, degree = 5, intercept = FALSE,Boundary.knots = c(1,7),data = df)
cmxDia <- colMeans(Xdia)
Xdia <- sweep(Xdia,2,cmxDia)


q = ncol(Xyear)

P <- diff(diag(q), differences = 3)
S <- crossprod(P)
K2 = S
num_basis1 = q

Dp1 = diag(10^(-6),ncol(K2),ncol(K2))

dim(A)
dim(Ad)
K2 = K2 + Dp1
#### fixed effects 

Xperiod = model.matrix(~ 1 +Period2, data = df)

X = cbind(Xyear,Xmes,Xdia,Xperiod)

dim(X)

dim(K2)
#X = X[,-1]
index2 = rep(F,dim(X)[2])


############################################
Xrei = matrix(0,nrow(X),ncol(X))
A = cbind(A,X)
Ad = cbind(Ad,Xrei)
indext = c(index1,index2)
length(indext)
dim(A)
dim(Ad)

K3 = K2
########## Fitting
system.time({m2 = model2(A,Ad,indext,K1,K2,K3,Xperiod,n)})

#####################################

beta = t(m2)

Quantil = matrix(0,dim(A)[1],9)

gamma2 = apply(beta,1,gamma_f,index = indext)
yseq = seq(1,max(df$case_count)+10,0.1)
Ay <- bernsteinPoly(yseq, degree = 20, intercept = FALSE,Boundary.knots = c(1,max(df$case_count)+10))
Ayd <- deriv(Ay)
Ay <- sweep(Ay,2,cmxy)


q = ncol(Ay)
Sigma <- matrix(1,q,q)                             # Duplicate matrix
Sigma[upper.tri(Sigma)] <- 0          # Change lower triangular part
#qi = vector()
#qs = vector()
#qm = vector()
Quantil = matrix(0,dim(A)[1],9)
ey = vector()
Ay = Ay%*%Sigma
Ayd = Ayd%*%Sigma

plot((Ay%*%gamma2[1:ncol(Ay),])[,1])
#gamma2 = gamma_f(mu_post_old, indext)
h1 = apply(Ay%*%gamma2[1:ncol(Ay),],1,mean)
h1l = apply(Ayd%*%gamma2[1:ncol(Ay),],1,mean)

gamma3 = apply(gamma2,1,mean)
qqnorm(qnorm(pnorm(A%*%gamma3)))
abline(a = 0,b = 1)  


for(j in 1:dim(X)[1]){
  #print(j)
  
  #j = sample(1:dim(A)[1],1)
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  

  h2 = mean(c(X[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + gamma[c(length(gamma))] + c(Bx2[j,]%*%gamma[(ncol(Apred) + 1):(length(gamma) - 1)]))
  Fy = pnorm(h1 + h2)
  fy = dnorm(h1 + h2)*h1l
  
  ey[j] = sum(fy*diff(yseq)[1]*yseq)
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + c(X[j,]*gamma[c(length(gamma))]))
  # usar 0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975
  
  
  Quantil[j,1] = (yseq[Fy > 0.025])[1]
  Quantil[j,2] = (yseq[Fy > 0.05])[1]
  Quantil[j,3] = (yseq[Fy > 0.1])[1]
  Quantil[j,4] =(yseq[Fy > 0.25])[1]
  Quantil[j,5] = (yseq[Fy > 0.5])[1]
  Quantil[j,6] =  (yseq[Fy > 0.75])[1]
  Quantil[j,7] = (yseq[Fy > 0.9])[1]
  Quantil[j,8] =  (yseq[Fy > 0.95])[1]
  Quantil[j,9] = (yseq[Fy > 0.975])[1]
  
  
  
  
}
#datapred_model1 = as.data.frame(cbind(df,qm,qi,qs,ey))
datapred_model1 = as.data.frame(cbind(df,Quantil,ey))
datapred = datapred_model1
head(datapred)
cr = 1 - sum(datapred$case_count > datapred$`9` | datapred$case_count < datapred$`1`)/(nrow(datapred))
cr

#df$predict2 = datapred$qm
mean(abs(ey - df$case_count))

quantiles = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
pred = (Quantil)



plotquantil = function(pred,df,quantiles){
indices = order(pred[,5])
pred = t(pred[order(pred[,5]),])
pred = c(pred)
#
i <- 1:nrow(df)
nd <- expand.grid(quantile = factor(quantiles), i = i)
#
df_sorted = df[indices,]
nd <- cbind(nd, df_sorted[nd$i,,drop = FALSE])
nd$predict <- pred
newdat <- nd
head(newdat)
#newdat$Subject <- factor(paste("Subject =", newdat$Subject))
#data$Subject <- factor(paste("Subject =", data$Subject))

## Usar esse gráfico
ggplot(newdat, aes(x = seq_along(predict), y = predict, color = quantile)) +
  geom_line(size = 1) +
  geom_point(aes(x = seq_along(predict), y = case_count), data = newdat, color = "black", alpha = 0.3,
             size = 0.75) +
  #geom_line(aes(x = Days, y = Reaction), data = data, color = "black", alpha = 0.3) +
  facet_wrap(~ month, nrow = 3) +
  ylab("Average reaction time (ms)") +
  labs(colour = "Estimated conditional quantiles") +
  colorspace::scale_color_discrete_diverging("Blue-Red2") +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.3,
                                        linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = rel(0.9)),
        legend.title = element_text(size = rel(0.9)),
        legend.position = "bottom",
        legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))
}

plotquantil(pred,df,quantiles)
############# New values #################

head(df)
df_teste = df_backup[floor(dim(df)[1] + 1):dim(df_backup)[1],]

summary(df$dia)
Xyear <- bernsteinPoly(df_teste$year, degree = 5, intercept = FALSE,Boundary.knots = c(1,7))
#cmx <- colMeans(Xyear)
Xyear <- sweep(Xyear,2,cmxYear)
Xmes <- bernsteinPoly(df_teste$mes, degree = 5, intercept = FALSE,Boundary.knots = c(1,12))
#cmx <- colMeans(Xmes)
Xmes <- sweep(Xmes,2,cmxMes)
Xdia <- bernsteinPoly(df_teste$dia, degree = 5, intercept = FALSE,Boundary.knots = c(1,7))
#cmx <- colMeans(Xdia)
Xdia <- sweep(Xdia,2,cmxDia)
Quantil = matrix(0,dim(Xdia)[1],9)
ey = vector()


Xperiod = model.matrix(~ 1 +Period2, data = df_teste)

Xteste = cbind(Xyear,Xmes,Xdia,Xperiod)

dim(X)

#X = X[,-1]
index2 = rep(F,dim(X)[2])

for(j in 1:dim(Xteste)[1]){
  print(j)
  
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  dim(Xteste)
  
  h2 = mean(c(Xteste[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  h2
  head(X)
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + gamma[c(length(gamma))] + c(Bx2[j,]%*%gamma[(ncol(Apred) + 1):(length(gamma) - 1)]))
  Fy = pnorm(h1 + h2)
  fy = dnorm(h1 + h2)*h1l
  
  ey[j] = sum(fy*diff(yseq)[1]*yseq)
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + c(X[j,]*gamma[c(length(gamma))]))
  
  Quantil[j,1] = (yseq[Fy > 0.025])[1]
  Quantil[j,2] = (yseq[Fy > 0.05])[1]
  Quantil[j,3] = (yseq[Fy > 0.1])[1]
  Quantil[j,4] =(yseq[Fy > 0.25])[1]
  Quantil[j,5] = (yseq[Fy > 0.5])[1]
  Quantil[j,6] =  (yseq[Fy > 0.75])[1]
  Quantil[j,7] = (yseq[Fy > 0.9])[1]
  Quantil[j,8] =  (yseq[Fy > 0.95])[1]
  Quantil[j,9] = (yseq[Fy > 0.975])[1]
  
}
mean(abs(ey - df_teste$case_count))

quantiles = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
pred = (Quantil)


plotquantil(pred,df_teste,quantiles)

indices = order(pred[,5])
pred = t(pred[order(pred[,5]),])
pred = c(pred)
#
i <- 1:nrow(df_teste)
nd <- expand.grid(quantile = factor(quantiles), i = i)
#
df_sorted = df_teste[indices,]
nd <- cbind(nd, df_sorted[nd$i,,drop = FALSE])
nd$predict <- pred
newdat <- nd
head(newdat)
#newdat$Subject <- factor(paste("Subject =", newdat$Subject))
#data$Subject <- factor(paste("Subject =", data$Subject))

## Usar esse gráfico
ggplot(newdat, aes(x = seq_along(predict), y = predict, color = quantile)) +
  geom_line(size = 1) +
  geom_point(aes(x = seq_along(predict), y = case_count), data = newdat, color = "black", alpha = 0.3,
             size = 0.75) +
  #geom_line(aes(x = Days, y = Reaction), data = data, color = "black", alpha = 0.3) +
  facet_wrap(~ month, nrow = 3) +
  ylab("Average reaction time (ms)") +
  labs(colour = "Estimated conditional quantiles") +
  colorspace::scale_color_discrete_diverging("Blue-Red2") +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.3,
                                        linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = rel(0.9)),
        legend.title = element_text(size = rel(0.9)),
        legend.position = "bottom",
        legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))


## Usar esse gráfico
datapred = as.data.frame(cbind(df_teste,Quantil,ey))
#datapred = datapred_model1
head(datapred)
#cr = 1 - sum(datapred$case_count > datapred$qs | datapred$case_count < datapred$qi)/(nrow(datapred))
#cr

#df$predict2 = datapred$qm
mean(abs(ey - df_teste$case_count))

quantiles = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
pred = t(Quantil)

#indices = order(pred[,5])
#pred = t(pred[order(pred[,5]),])
pred = c(pred)
#
i <- 1:nrow(df_teste)
nd <- expand.grid(quantile = factor(quantiles), i = i)
#

nd <- cbind(nd, df_teste[nd$i,,drop = FALSE])
nd$predict <- pred
newdat <- nd
head(newdat)
head(datapred)
head(newdat)
newdata_prev = newdat[1:700,]
ggplot(newdata_prev, aes(x = i, y = predict, color = quantile)) +
  geom_line(size = 1) +
  geom_line(aes(x = i, y = case_count), data = newdata_prev, color = "black", alpha = 0.9,
            size = 0.75) +
  #geom_line(aes(x = Days, y = Reaction), data = data, color = "black", alpha = 0.3) +
  #facet_wrap(~ month, nrow = 3) +
  ylab("Average reaction time (ms)") +
  labs(colour = "Estimated conditional quantiles") +
  colorspace::scale_color_discrete_diverging("Blue-Red2") +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.3,
                                        linetype = "dotted"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = rel(0.9)),
        legend.title = element_text(size = rel(0.9)),
        legend.position = "bottom",
        legend.key = element_rect(fill = "transparent", colour = "transparent")) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))

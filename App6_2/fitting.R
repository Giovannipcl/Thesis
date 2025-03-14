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
sourceCpp("funcs1.cpp")
source("Allmodels.R")
source("func_aux.R")
sourceCpp("arg_min.cpp")
sourceCpp("gradiente.cpp")

################## MOdelo 1##############

### Transformacao em y

head(df)

df = as.data.frame(df)
df_backup = df
df = df[1:floor(dim(df)[1]*0.9),]

A <- bernsteinPoly(df$case_count, degree = 20, intercept = FALSE,Boundary.knots = c(1,max(df$case_count)+10),data = df)
dim(A)
#umv = rep(1,dim(A)[1])
Ad = deriv(A)
#Xi0 = matrix(0,nrow(A), ncol(X))
q = ncol(A)
Sigma <- matrix(1,q,q)                             # Duplicate matrix
Sigma[upper.tri(Sigma)] <- 0          # Change lower triangular part

#Nao zerar
cmx <- colMeans(A)
A <- sweep(A,2,cmx)
A = A%*%Sigma
Ad = Ad%*%Sigma
q = q
index1 = c(rep(TRUE,q))
P <- diff(diag(q), differences = 1)
S <- crossprod(P)
P <- diff(diag(q), differences = 3)
S <- crossprod(P)
K1 = S
num_basis1 = q
n = dim(A)[1]
Dp1 = diag(10^(-6),ncol(K1),ncol(K1))
K1 = K1 + Dp1

#### fixed effects 

X = model.matrix(~ year + mes + dia + dia*Period2, data = df)
index2 = rep(F,dim(X)[2])

############################################
Xrei = matrix(0,nrow(X),ncol(X))
A = cbind(A,X)
Ad = cbind(Ad,Xrei)
indext = c(index1,index2)
length(indext)
dim(Ad)
n = dim(A)[1]
K2 = K1


##############################################
##### Model fitting
system.time({m1 = model1(A,Ad,indext,K1,X,n)})
##############################################


beta <- do.call(rbind, m1)
beta = t(beta)

yseq = seq(1,max(df$case_count)+10,0.1)
Ay <- bernsteinPoly(yseq, degree = 20, intercept = FALSE,Boundary.knots = c(1,max(df$case_count)+10))

Ay <- sweep(Ay,2,cmx)
Ayd <- deriv(Ay)

q = ncol(Ay)
Sigma <- matrix(1,q,q)                             # Duplicate matrix
Sigma[upper.tri(Sigma)] <- 0          # Change lower triangular part


qi = vector()
qs = vector()
qm = vector()
ey = vector()
Ay = Ay%*%Sigma
Ayd = Ayd%*%Sigma
length(qm)
Quantil = matrix(0,dim(A)[1],9)

gamma2 = apply(beta,1,gamma_f,index = indext)


h1 = apply(Ay%*%gamma2[1:ncol(Ay),],1,mean)
h1l = apply(Ayd%*%gamma2[1:ncol(Ay),],1,mean)

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
ey
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

############################################################################################################################
############################################################################################################################
###############################################Predict new values#############################################################################


df_backup
df_teste = df_backup[floor(dim(df)[1] + 1):dim(df_backup)[1],]


ey = vector()

Xpred =  model.matrix(~ year + mes + dia + dia*Period2, data = df_teste)

Quantil = matrix(0,dim(Xpred)[1],9)


for(j in 1:dim(Xpred)[1]){
  print(j)
  
  #j = sample(1:dim(A)[1],1)
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  
  
  h2 = mean(c(Xpred[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  
  
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

datapred = as.data.frame(cbind(df_teste,Quantil,ey))
#datapred = datapred_model1
head(datapred)
#cr = 1 - sum(datapred$case_count > datapred$qs | datapred$case_count < datapred$qi)/(nrow(datapred))
#cr

#df$predict2 = datapred$qm
mean(abs(ey - df_teste$case_count))

quantiles = c(0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975)
pred = (Quantil)

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

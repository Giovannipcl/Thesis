################ Loading the packages ###################
library("tramME")
library(lme4)
library(splines2)
library(LaplacesDemon)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv)
library(coda)
library(MASS)
#####################################
### Getting the data
data("sleepstudy")
sleepstudy$Subject
####################################
##### Loading the functions and .cpp files
dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
dir_models <- file.path(dir, "Models") 
setwd(dir_models)
source("func_auxiliar.R")
source("model.R")

dir_cpp <- file.path(dir, "cpp") 
setwd(dir_cpp)
sourceCpp("arg_min_Armad.cpp")
sourceCpp("gradiente_Armad.cpp")
sourceCpp("funcs.cpp")
sourceCpp("cpo.cpp")

###################################################
####### Starting the model

#### Transformation in y

sleepstudy$y = sleepstudy$Reaction
A <- bernsteinPoly(sleepstudy$Reaction, degree = 10, intercept = FALSE,Boundary.knots = c(min(sleepstudy$Reaction)-12,max(sleepstudy$Reaction)+12),data = sleepstudy)
Ad = deriv(A)
q = ncol(A)
Sigma <- matrix(1,q,q)                            
Sigma[upper.tri(Sigma)] <- 0          

#Nao zerar
#cmx <- colMeans(A)
#A <- sweep(A,2,cmx)
A = A%*%Sigma
Ad = Ad%*%Sigma
index1 = c(rep(TRUE,q))


##### Getting the matrix K1
P <- diff(diag(q), differences = 3)
S <- crossprod(P)
K1 = S

Dp1 = diag(10^(-6),ncol(K1),ncol(K1))
K1 = K1 + Dp1


#### Getting the fixed effects and random effects
sleepstudy$Dayss = scale(sleepstudy$Days)
Xre = model.matrix(~-1 + Subject,data = sleepstudy)
X = model.matrix(~ 1+Days, data = sleepstudy)
X = cbind(X,Xre,c(sleepstudy$Days)*Xre)
index2 = rep(F,dim(X)[2])
nf = length(levels(sleepstudy$Subject)) 
Xrei = matrix(0,nrow(X),ncol(X))   # matrix of first derivative



#### Final matrix
A = cbind(A,X)
Ad = cbind(Ad,Xrei)
indext = c(index1,index2)


################################
##############Fitting the model
system.time({m1 = trame_fit(A,Ad,indext,K1,nf)})
################################
################################

############# Extracting the samples from the marginal posterior distributions
beta = m1$beta
beta <- do.call(rbind, beta)
beta = t(beta)


##### Extracting the random effects
b0 = apply(beta,2,mean)[13:(13+17)]
b1 = apply(beta,2,mean)[(13+18):ncol(A)]
data <- data.frame(b0 = b0, b1 = b1)

ggplot(data, aes(x = b0, y = b1)) +
  geom_point(color = "black", size = 2) +  
  labs(
    x = expression(b[0]),
    y = expression(b[1])) +
  theme_minimal() +                       
  theme(
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),   
    axis.text.x = element_text(size = 14),    
    axis.text.y = element_text(size = 14),    
    plot.title = element_text(size = 18)      
  )


###################################################################################################
#######################################################################################################
###########################################################################

yseq = seq(min(sleepstudy$Reaction)-12,max(sleepstudy$Reaction)+12,0.051)

Ay<- bernsteinPoly(yseq, degree = 10, intercept = FALSE,Boundary.knots = c(min(sleepstudy$Reaction)-12,max(sleepstudy$Reaction)+12))
Ayd <- deriv(Ay)
q = ncol(Ay)
Sigma <- matrix(1,q,q)                             
Sigma[upper.tri(Sigma)] <- 0          
Ay = Ay%*%Sigma
Ayd = Ayd%*%Sigma
#sm$deriv <- 0
#datanew = as.data.frame(yseq)
#names(datanew) = "Reaction"
#head(datanew)
#Ay <- PredictMat(sm, datanew)
#sm$deriv = 1
#Ayd <- PredictMat(sm, datanew)
nd <- expand.grid(Days = seq(min(sleepstudy$Days), max(sleepstudy$Days),
                             length.out = 100),
                  Subject = unique(sleepstudy$Subject))
Xre2 = model.matrix(~-1 + Subject,data = nd)
X2 = model.matrix(~ 1 + Days, data = nd)
gamma2 = apply(beta,1,gamma_f,index = indext)
X2 = cbind(X2,Xre2,c(nd$Days)*Xre2)
Quant = matrix(0,9,dim(X2)[1])
h1 = apply(Ay%*%gamma2[1:ncol(Ay),],1,mean)


for(j in 1:dim(X2)[1]){
  
  
  #j = 1
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  
  

  
  h2 = mean(c(X2[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  dim(X2)
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + gamma[c(length(gamma))] + c(Bx2[j,]%*%gamma[(ncol(Apred) + 1):(length(gamma) - 1)]))
  Fy = pnorm(h1 + h2)
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + c(X[j,]*gamma[c(length(gamma))]))
  
  
  
  Quant[9,j] = (yseq[Fy > 0.975])[1]
  Quant[8,j] = (yseq[Fy > 0.95])[1]
  Quant[7,j] = (yseq[Fy > 0.9])[1]
  Quant[6,j] = (yseq[Fy > 0.75])[1]
  Quant[5,j] = (yseq[Fy > 0.5])[1]
  Quant[4,j] = (yseq[Fy > 0.25])[1]
  Quant[3,j] = (yseq[Fy > 0.1])[1]
  Quant[2,j] = (yseq[Fy > 0.05])[1]
  Quant[1,j] = (yseq[Fy > 0.025])[1]
  
}

pr <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
plot_cquant(sleepstudy, nd, Quant, pr)


##################################################################
#####################indiv 308 e 309 ####################################################
##################################################################


### Individual 308
nd <- expand.grid(Days = seq(min(sleepstudy$Days), max(sleepstudy$Days),
                             length.out = 10),
                  Subject = unique(sleepstudy$Subject))
nd
nd$Subject == 308

Xre2 = model.matrix(~-1 + Subject,data = nd)[ nd$Subject == 308,]
X2 = model.matrix(~ 1 + Days, data = nd)[nd$Subject == 308,]


X2 = cbind(X2,Xre2,c(nd$Days[nd$Subject == 308])*Xre2)
fy_den = matrix(0,length(yseq),dim(X2)[1])
h1 = apply(Ay%*%gamma2[1:ncol(Ay),],1,mean)
h1l = apply(Ayd%*%gamma2[1:ncol(Ay),],1,mean)

for(j in 1:dim(X2)[1]){
  print(j)
  
  #j = 1
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  
  
  
  
  h2 = mean(c(X2[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + gamma[c(length(gamma))] + c(Bx2[j,]%*%gamma[(ncol(Apred) + 1):(length(gamma) - 1)]))
  fy = dnorm(h1 + h2)*h1l
  fy_den[,j] = fy
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + c(X[j,]*gamma[c(length(gamma))]))
  
  
  
}

df_set1 <- as.data.frame(fy_den)
df_set1$reaction_time <- yseq
df_set1$set <- "308"



#######################################################################################################################################################
### Individual 309
nd <- expand.grid(Days = seq(min(sleepstudy$Days), max(sleepstudy$Days),
                             length.out = 10),
                  Subject = unique(sleepstudy$Subject))

nd$Subject == 309

Xre2 = model.matrix(~-1 + Subject,data = nd)[ nd$Subject == 309,]
X2 = model.matrix(~ 1 + Days, data = nd)[nd$Subject == 309,]
#X = model.matrix(~1+ Days*Subject, data = sleepstudy)

X2 = cbind(X2,Xre2,c(nd$Days[nd$Subject == 309])*Xre2)
fy_den = matrix(0,length(yseq),dim(X2)[1])
h1 = apply(Ay%*%gamma2[1:ncol(Ay),],1,mean)
h1l = apply(Ayd%*%gamma2[1:ncol(Ay),],1,mean)
#gamma = gamma_f(mu_post, indext)
for(j in 1:dim(X2)[1]){
  print(j)
  
  #j = 1
  
  #A <- sweep(A,2,cmx)
  
  #h1 = Ay%*%gamma[1:ncol(Ay)]
  
  
  
  
  h2 = mean(c(X2[j,]%*%gamma2[(ncol(Ay)+1):(ncol(A)),]))
  
  
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + gamma[c(length(gamma))] + c(Bx2[j,]%*%gamma[(ncol(Apred) + 1):(length(gamma) - 1)]))
  fy = dnorm(h1 + h2)*h1l
  fy_den[,j] = fy
  #Fy = pnorm(c(Apred%*%gamma[1:ncol(Apred)]) + c(X[j,]*gamma[c(length(gamma))]))
  
  
  
}


df_set2 <- as.data.frame(fy_den)
df_set2$reaction_time <- yseq
df_set2$set <- "309"



df_combined <- rbind(df_set1, df_set2)

df_long <- pivot_longer(df_combined, cols = -c(reaction_time, set), names_to = "group", values_to = "density")
blue_palette <- colorRampPalette(c("darkblue", "lightblue"))(10)
orange_palette <- colorRampPalette(c("darkorange", "lightyellow"))(10)


df_long$group = factor(df_long$group, levels = c("V1" ,"V2" ,"V3" ,"V4" ,"V5", "V6" ,"V7" ,"V8" ,"V9", "V10"))
#############################################################################################################################3

combined_palette <- c(blue_palette, orange_palette)
ggplot(df_long, aes(x = reaction_time, y = density, group = interaction(group, set), color = interaction(group, set))) +
  geom_line(size = 1) +  
  scale_color_manual(values = combined_palette, labels = c(seq(0,9,1), seq(0,9,1))) +  
  theme_minimal() +  
  labs(
    x = "Average reaction time (ms)",
    y = "Density"
  ) +
  guides(color = guide_legend(
    ncol = 2,
    title = c("308      309")
  )) +
  theme( axis.text = element_text(size = rel(1.2)),
         axis.title = element_text(size = rel(1.2)),
         legend.position = "right",  
         legend.justification = c("left", "center"),  
         legend.title = element_text(size = 15),  
         legend.text = element_text(size = 15),  
         panel.grid = element_blank()  
  ) 


#############################################################################

d <- expand.grid(Days = seq(min(sleepstudy$Days), max(sleepstudy$Days),
                            length.out = 10),
                 Subject = 308)


Xm = model.matrix(~ 1 + Days, data = d)



days_effects = (Xm[,c(1:2)]%*%apply(gamma2,1,mean)[c(11:12)])

ncol(A)
marginal_days = function(i){  
  
  nd <- expand.grid(Days = i,
                    Subject = unique(sleepstudy$Subject))
  
  Xre2 = model.matrix(~-1 + Subject,data = nd)
  X2 = cbind(Xre2,c(nd$Days)*Xre2)
  
  
  h1 = apply(sweep(Ay%*%gamma2[1:ncol(Ay),], 2, c(apply(days_effects[i] + X2%*%(gamma2[13:ncol(A),]),2,mean)) , "+"),1,mean)
  fy = dnorm(h1)*h1l
  return(fy)
}
h1 = lapply(seq(1:10),marginal_days)



df <- data.frame(
  valores = unlist(h1),                
  posicao = rep(1:length(h1), times = sapply(h1, length))  
)

df$x = rep(yseq,10)
dim(Ay)
length(h1)
df$posicao =factor(df$posicao)
names(df)[2] = c("Days")
head(df)
ggplot(df, aes(x = x, y = valores, group = Days, color = Days)) +
  geom_line(size = 1) +  
  scale_color_manual(values = blue_palette) + 
  theme_minimal() +  
  labs(
    x = "Average reaction time (ms)",
    y = "Density"
  ) +
  theme( axis.text = element_text(size = rel(1.2)),
         axis.title = element_text(size = rel(1.2)),
         legend.position = "right",  
         legend.justification = c("left", "center"),  
         legend.title = element_text(size = 15),
         legend.text = element_text(size = 15),  
         panel.grid = element_blank()  
  ) 

##########################################################

yseq = seq(min(sleepstudy$Reaction)-12,max(sleepstudy$Reaction)+12,0.051)
Aq <- bernsteinPoly(yseq, degree = 10, intercept = FALSE,Boundary.knots = c(min(sleepstudy$Reaction)-12,max(sleepstudy$Reaction)+12))
#Ay <- sweep(Ay,2,cmx)
Ayq <- deriv(Aq)

q = ncol(Aq)
Sigma <- matrix(1,q,q)                             # Duplicate matrix
Sigma[upper.tri(Sigma)] <- 0          # Change lower triangular part

Aq = Aq%*%Sigma
Ayq = Ayq%*%Sigma

###############################################

#sm$deriv <- 0
datanew = as.data.frame(yseq)
names(datanew) = "Reaction"
head(datanew)

#######################################################################
indbeta = 12
median_impact = apply((-beta[,indbeta]/t(Ayq%*%gamma2[1:ncol(Ayq),])),2,median)
#median_impact = apply(as.mcmc((beta[,12]/t(Ayq%*%gamma2[1:ncol(Ayq),]))),1,HPDinterval)


me = as.mcmc((-beta[,indbeta]/t(Ayq%*%gamma2[1:ncol(Ayq),])))
intervals=HPDinterval(me,prob = 0.9)
days_impact = as.data.frame(cbind(intervals,median_impact,yseq))
names(days_impact)
ggplot(days_impact, aes(x = yseq)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray", alpha = 0.5) +  # Sombra cinza
  geom_line(aes(y = median_impact), color = "darkblue", size = 1.5) +  # Linha da mediana azul
  labs(
    x = "Median of average reaction time median", 
    y = "Impact of Days") +
  theme_minimal()  +                       # Tema minimalista
  theme(
    axis.title.x = element_text(size = 16),  # Tamanho do rótulo do eixo x
    axis.title.y = element_text(size = 16),  # Tamanho do rótulo do eixo y
    axis.text.x = element_text(size = 14),    # Tamanho dos valores do eixo x
    axis.text.y = element_text(size = 14),    # Tamanho dos valores do eixo y
    plot.title = element_text(size = 18)      # Tamanho do título
  )

############# For a particular individual

#######################################################################
median_impact308 = apply((-(beta[,12] + beta[,31])/t(Ayq%*%gamma2[1:ncol(Ayq),])),2,mean)
median_impact309 = apply((-(beta[,12] + beta[,32])/t(Ayq%*%gamma2[1:ncol(Ayq),])),2,mean)


days_impact = as.data.frame(cbind(median_impact308,median_impact309,yseq))
names(days_impact)

# Transformar os dados no formato longo (long format)
dados_long <- gather(days_impact, key = "linha", value = "valor", median_impact308, median_impact309)

head(dados_long)
names(dados_long) = c("yseq", "Subject", "valor")
dados_long$Subject
dados_long[dados_long$Subject == "median_impact308",2] <- "308"
dados_long[dados_long$Subject == "median_impact309",2] <- "309"


ggplot(dados_long, aes(x = yseq, y = valor, color = Subject)) +
  geom_line(size = 1) +  # Adiciona as linhas
  labs(
    x = "Median of average reaction time median", 
    y = "Impact of Days") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))  +                       # Tema minimalista
  theme(
    axis.title.x = element_text(size = 16),  # Tamanho do rótulo do eixo x
    axis.title.y = element_text(size = 16),  # Tamanho do rótulo do eixo y
    axis.text.x = element_text(size = 14),    # Tamanho dos valores do eixo x
    axis.text.y = element_text(size = 14),    # Tamanho dos valores do eixo y
    plot.title = element_text(size = 18),      # Tamanho do título
    legend.title = element_text(size = 16),   # Tamanho do título da legenda
    legend.text = element_text(size = 14)     # Tamanho do texto da legenda
  )

##################### PIT

means = m1$means
mu = m1$mu
pontos2 = m1$pontos2
post_alpha = m1$post_alpha
weights = m1$weights


pit3 = function(k,ind,wY){
  
  #matriz_corrigida <- as.matrix((solve(calc_neg_hess_ff4(means[,k],pontos2[k,],A = A, Al = Ad, K1 = K1,K2=K2,K3=K3,X=Xf,index = indext))))
  matriz_corrigida = solve(calc_neg_hess_ff4(means[,k],pontos2[k,],A = A, Al = Ad, K1 = K1,index = indext,nf = nf))
  beta_cpo = mvrnorm(1000,mu[,k], matriz_corrigida)
  beta_cpo = apply(beta_cpo,1,gamma_f,indext)
  indic_true = yseq < sleepstudy$Reaction[ind]
  newyseq =  yseq[yseq < sleepstudy$Reaction[ind]]
  Aqn = Aq[indic_true,]
  Aqnd = Ayq[indic_true,]
  
  densidade_bctm = function(ind){
    (1/(dnorm(A[ind,]%*%beta_cpo)*(Ad[ind,]%*%beta_cpo)))
    
  }
  
  densidade_bctm2 = function(ind){
    mean(1/(dnorm(A[ind,]%*%beta_cpo)*(Ad[ind,]%*%beta_cpo)))
    
  }
  
  h1 = sweep(Aqn%*%beta_cpo[1:ncol(Aqn),],2,A[ind,(ncol(Aqn)+1):(ncol(A))]%*%beta_cpo[(ncol(Aqn)+1):(ncol(A)),],"+")
  fy = dnorm(h1)*(Aqnd%*%beta_cpo[1:ncol(Aqn),])
  
  pit =mean((apply(diff(newyseq)[1]*fy,2,sum))*densidade_bctm(ind))
  
  return(list(CPO =post_alpha[k]*densidade_bctm2(ind)*wY[k], PIT =  pit))
  
  
}



CPO_values_all = vector()
PIT_values_all = vector()
dim(A)

for(j in 1:dim(A)[1]){
  print(j)
  pit_ind_1 = (mclapply(1:(dim(pontos2)[1]), function(i) pit3(i,j, weights),mc.cores = 6))
  CPO_values <- sapply(pit_ind_1, function(x) x$CPO)
  PIT_values <- sapply(pit_ind_1, function(x) x$PIT)
  
  CPO_values_all[j] = 1/sum(CPO_values);
  PIT_values_all[j] = sum(PIT_values*post_alpha*weights)*CPO_values_all[j];
  print(PIT_values_all[j] )
  print(CPO_values_all[j])
}

# Calcular os quantis teóricos da uniforme
n <- length(PIT_values_all)
theoretical_quantiles <- qunif(ppoints(n), min = 0, max = 1)
sorted_pits = sort(PIT_values_all)
# Fazer o Q-Q plot
plot(theoretical_quantiles, sorted_pits, 
     main = "Q-Q Plot dos PITs", 
     xlab = "Quantis teóricos (Uniforme)", 
     ylab = "Quantis observados (PIT)", 
     pch = 19, col = "blue")

# Adicionar linha de referência
abline(0, 1, col = "red", lwd = 2)


plot(ecdf(PIT_values_all), main = "", 
     xlab = "PIT", ylab = "Empirical CDF", col = "black", lwd = 2,cex.axis =1.5)
curve(punif(x, 0, 1), add = TRUE, col = "red", lwd = 2) # Theoretical CDF
legend("bottomright", legend = c("Empirical CDF", "Theoretical CDF"), col = c("black", "red"), lwd = 2,cex=1.3)

hist(PIT_values_all, xlab = "PIT", breaks =10, cex.lab=1.5,cex.axis =1.5,main ="")
abline(h = 18, lwd = 2,lty = 2)
####################################

# KS
ks_teste <- ks.test(PIT_values_all, "punif", min = 0, max = 1)

print(ks_teste)

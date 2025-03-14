
library(LaplacesDemon)
library(parallel)
library(dplyr)
library(Rcpp)


model1 = function(A,Al,index,K1,X,n){

###################################################################

xin2 = calc_x0(c(1),rep(1,length(index)),A,Al,index,K1,X)


d = optim(par = c(0), fn = calc_lpost3,xin = xin2, A = A,Al = Al,index = index,K1 = K1,X = X,n =n, method = "L-BFGS-B",
                         control=list(fnscale=-1,maxit = 5000,trace = TRUE),hessian = TRUE)


min_values = d$par - 3*sqrt(diag(solve(-(d$hessian))))
max_values = d$par + 3*sqrt(diag(solve(-(d$hessian))))

result <- generate_uniform_grid(1, 21, min_values, max_values)
sigmas = c(result$points)$Var1

v.f.w = sapply(1:length(sigmas), function(i) calc_lpost_cpp(sigmas[i], xin2,A,Al,index,K1,X,n))

unl = v.f.w - mean(v.f.w)
pontos_x = do.call(cbind,result$ranges)
pontos_novos = exp(pontos_x)

weights = diff(pontos_novos)
weights = c(weights,weights[length(weights)])

post_alpha =exp(unl)/sum(exp(unl)*weights)
means = (sapply(1:length(sigmas), function(i) calc_x0(sigmas[i], xin2,A,Al,index,K1,X)))



sigmas_beta = sapply(1:length(sigmas), function(i) vari = diag(chol2inv(chol(calc_neg_hess_ff4(means[,i],sigmas[i],A = A, Al = Al, K1 = K1,index = index,X = X)))))

nsigmas = length(sigmas)

system.time({correções4 = mclapply(1:nsigmas, function(i) teste = optim(rep(0,length(index)),arg_min22, gr = gradiente_f,
                                                                      H = solve(calc_neg_hess_ff4(means[,i],sigmas[i],A = A, Al = Al, K1 = K1, index = index,X = X)),  A = A,
                                                                      mu_ini = c(means[,i]),Al = Al, K_p = K(sigmas[i],K1,X),index = index,
                                                                      method = "BFGS")$par,mc.cores = 6)})

correções4 <- do.call(cbind, correções4)
mu = means + correções4

mu_post = sapply(1:ncol(A), function(ii) weighted.mean(mu[ii,], w=post_alpha*weights))



mu_post_old = sapply(1:ncol(A), function(ii) weighted.mean(means[ii,], w=post_alpha*weights))


sigma_post = 
  sapply(1:ncol(A), function(ii) 
    weighted.mean((mu[ii,] - mu_post[ii])^2 + sigmas_beta[ii,], w=post_alpha*weights)
  ) %>% 
  sqrt

sigma_post_old = 
  sapply(1:ncol(A), function(ii) 
    weighted.mean((means[ii,] - mu_post_old[ii])^2 + sigmas_beta[ii,], w=post_alpha*weights)
  ) %>% 
  sqrt

beta_post = mclapply(1:dim(A)[2],function(i) MH(i, mu_post = mu_post, mu = mu, sigmas_beta = sigmas_beta, post_alpha = post_alpha, weights = weights),mc.cores = 6)
return(beta_post)
}


#################################################
########## Model 2 e 3

model2 = function(A,Al,index,K1,K2,K3,X,n){


  xin = rep(1,length(index))

xin2 = calc_x0(c(-1,-1,-1,-1),rep(1,length(index)),A,Al,index,K1,K2,K3,X)


n = dim(df)[1]
d = optim(par = c(0,0,0,0), fn = calc_lpost_cpp,xin = xin2, A = A,Al = Al,index = index,K1 = K1,K2 = K2,K3 = K3,X = X,n =n, method = "L-BFGS-B",
                        control=list(fnscale=-1,maxit = 5000,trace = TRUE),hessian = TRUE)


min_values = d$par - 2*sqrt(diag(solve(-(d$hessian))))
max_values = d$par + 2*sqrt(diag(solve(-(d$hessian))))

result <- generate_uniform_grid(4, 4, min_values, max_values)

# Visualiza os pontos e os pesos
pontos2 <- as.matrix(result$points)
weights <- c(result$weights)

v.f.w = sapply(1:(dim(pontos2)[1]), function(i) calc_lpost_cpp(pontos2[i,], xin2,A,Al,index,K1,K2,K3,X,n))

unl = v.f.w - mean(v.f.w)
post_alpha =exp(unl)/sum(exp(unl)*weights)

means = (sapply(1:dim(pontos2)[1], function(i) calc_x0(pontos2[i,], xin2,A,Al,index,K1,K2,K3,X)))


sigmas_beta = sapply(1:dim(pontos2)[1], function(i) vari = diag(chol2inv(chol(calc_neg_hess_ff4(means[,i],pontos2[i,],A = A, Al = Al, K1 = K1,K2 = K2,K3,index = index,X = X)))))
sigmas = pontos2

nsigmas = dim(sigmas)[1]

system.time({correções3 = sapply(1:nsigmas, function(i) teste = optim(rep(0,length(index)),arg_min22, gr = gradiente_f,
                                                                      H = solve(calc_neg_hess_ff4(means[,i],pontos2[i,],A = A, Al = Al, K1 = K1, K2 = K2, K3 = K3, index = index,X = X)),  A = A,
                                                                      mu_ini = c(means[,i]),Al = Al, K_p = K(pontos2[i,],K1,K2,K3,X),index = index,
                                                                      method = "BFGS")$par)})




#####################################
mu = means + correções3
mu_post = sapply(1:ncol(A), function(ii) weighted.mean(mu[ii,], w=post_alpha))
mu_post_old = sapply(1:ncol(A), function(ii) weighted.mean(means[ii,], w=post_alpha))

sigma_post = 
  sapply(1:ncol(A), function(ii) 
    weighted.mean((mu[ii,] - mu_post[ii])^2 + sigmas_beta[ii,], w=post_alpha)
  ) %>% 
  sqrt

sigma_post_old = 
  sapply(1:ncol(A), function(ii) 
    weighted.mean((means[ii,] - mu_post_old[ii])^2 + sigmas_beta[ii,], w=post_alpha)
  ) %>% 
  sqrt

beta_post = lapply(1:dim(A)[2],function(i) MH(i, mu_post = mu_post, mu = mu, sigmas_beta = sigmas_beta, post_alpha = post_alpha, weights = weights))
beta_post <- do.call(rbind, beta_post)
return(beta_post)

}


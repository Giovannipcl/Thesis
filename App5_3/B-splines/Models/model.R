


trame_fit = function(A,Al,index,K1,nf){

system.time({xin2 = calc_x0(c(0,0,0,0),rep(1,length(index)),A,Al,index,K1,nf)})
plot(xin2)
#########################################################
n = length(sleepstudy$y)
system.time({ d = optim(par = c(0,0,0,0), fn = calc_lpost_cpp,xin = xin2, A = A,Al = Al,index = index,K1 = K1,nf = nf,n =n, method = "L-BFGS-B",
                        control=list(fnscale=-1,maxit = 5000,trace = TRUE),hessian = TRUE)})


min_values = d$par - 3*sqrt(diag(solve(-(d$hessian))))
max_values = d$par + 3*sqrt(diag(solve(-(d$hessian))))

result <- generate_uniform_grid(4, 4, min_values, max_values)
# Visualiza os pontos e os pesos
points = as.matrix(result$points)

###########################################################
library(parallel)
library(dplyr)
pontos2 = (points)
nsigmas = dim(pontos2)[1]
sigmamod = d$par
print(paste("Moda:",sigmamod))
v.f.w = unlist(mclapply( as.data.frame(t(pontos2)), calc_lpost_cpp,xin = xin2 ,A = A, Al = Al, index = index, K1 = K1,n = n,nf = nf))
unl = v.f.w - mean(v.f.w)
pontos_x = do.call(cbind,result$ranges)
pontos_novos = exp(pontos_x)

weights = calculate_cell_weights(pontos_novos)


post_alpha =exp(unl)/sum(exp(unl)*weights)


means = apply(pontos2,1,calc_x0,x0 = xin2 ,A = A, Al = Al, index = index, K1 = K1,nf = nf, simplify = TRUE)
sigmas_beta = sapply(1:nsigmas, function(i) vari = diag(chol2inv(chol(calc_neg_hess_ff4(means[,i],c(pontos2[i,]),A = A, Al = Al, K1 = K1,index = index,nf = nf)))))



sigmas = pontos2

system.time({correções3 = sapply(1:nsigmas, function(i) teste = optim(rep(0,length(indext)),arg_min22, gr = gradiente_f,
                                                                      H = solve(calc_neg_hess_ff4(means[,i],pontos2[i,],A = A, Al = Al, K1 = K1,index = index,nf = nf)),  A = A,
                                                                      mu_ini = c(means[,i]),Al = Al, K_p = K(pontos2[i,],K1,nf),index = index,
                                                                      method = "BFGS")$par)})


mu = means + correções3

mu_post = sapply(1:ncol(A), function(ii) weighted.mean(mu[ii,], w=post_alpha*weights))



mu_post_old = sapply(1:ncol(A), function(ii) weighted.mean(means[ii,], w=post_alpha*weights))
library(dplyr)
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

beta_post = mclapply(1:dim(A)[2],function(i) MH(i, mu_post = mu_post, mu = mu, sigmas_beta = sigmas_beta, post_alpha = post_alpha, weights = weights),mc.cores = 1)

return(list(beta=beta_post, mu = mu, means = means, pontos2 = pontos2,post_alpha = post_alpha,weights = weights))
}


library(parallel)
library(purrr)
library("fastmatrix")
library("DescTools")
library(plyr)
library("LaplacesDemon")
library("optimParallel")
library("pracma")
library(expm)
#install.packages("fBasics")
library(fBasics)
library(coda)
library("interp")



bctm = function(A,Al,X,K1,K2,K3,q,index){
  
  n = length(y)
  Dp1 = diag(10^(-6),ncol(K1),ncol(K1))
  Dp2 = diag(10^(-6),ncol(K2),ncol(K2))
  Dp3 = diag(10^(-6),ncol(K3),ncol(K3))
  index = indext
  
  Cbeta = function(beta,index){
    cbetat = rep(1,length = length(index))
    cbetat[index] = exp(beta[index])
    return(diag(cbetat))
  }
  
  Cbetal = function(beta,index){
    cbetat = rep(0,length = length(index))
    cbetat[index] = exp(beta[index])
    return(diag(cbetat))
  }
  
  
  K1 = K1 + Dp1
  K2 = K2 + Dp2
  #K3 = K3 + Dp3
  #K = function(sigma,K1,K2){
  #  as.matrix(drop(bdiag((1/(sigma[1]))*K1 + (1/(sigma[2]))*K2,diag(c(rep(0.1,ncol(X))))))) 
  #}
  
  
  K = function(sigma,K1,K2,K3){
    as.matrix(drop(bdiag((exp(sigma[1]))*(K1 + K2),exp(sigma[2])*K3,diag(c(rep(10^-6,ncol(X))))))) 
  }
  
  #K = function(sigma,K1,K2,K3){
  #  as.matrix(drop(bdiag((exp(sigma[1]))*(K1 + K2),exp(sigma[2])*K3,diag(c(rep(10^{-6},ncol(X))))))) 
  #}
  
  
  calc_lprior = function(sigma){
    sum(dhalfcauchy(exp(sigma),scale = 25,log = TRUE))
  }
  
  #calc_lprior = function(sigma){
  #  sum(dnorm(exp(sigma),0,1,log = TRUE))
  #}
  
  #calc_lprior = function(sigma){
  #  sum(dinvgamma(exp(sigma),1,0.001,log = TRUE))
  #}
  
  
  gamma_f = function(beta,index){
    gammat = vector(length = length(index))
    gammat[index] = exp(beta[index])
    gammat[!index] = beta[!index]
    return(gammat)
  }
  
  calc_ljoint = function(A,beta, sigma,  Al,K1,K2,K3){
    gamma = gamma_f(beta,index)
    chol_Q = K(sigma,K1,K2,K3) %>% chol
    logdet_Q_half = chol_Q %>% diag %>% log %>% sum
    quad_form = crossprod(chol_Q %*% beta) %>% drop
    
    res = sum(dnorm(A%*%gamma,0,1,log = TRUE)) + sum(log(Al%*%gamma)) + logdet_Q_half - 0.5*quad_form + calc_lprior(sigma) 
    return(res)
    
  }
  
  calc_ff = function(beta,sigma){
    gamma = gamma_f(beta,index)
    if(all(Al%*%gamma > 0)){
      sum(dnorm(A%*%gamma,0,1,log = TRUE)) + sum(log(Al%*%gamma)) - 
        0.5 * (as.matrix(drop(beta%*% K(sigma,K1,K2,K3) %*% beta))) + sum(beta[index])
    }else{-10000000}
  }
  
  
  calc_grad_ff = function(beta, sigma){                                                             
    gamma = gamma_f(beta,index)
    cbetat = rep(0,length = length(index))
    cbetat[index] = rep(1, sum(index))
    apply(((c(-A%*%gamma)*A)%*%Cbeta(beta,index)),2,sum) + 
      apply((Al)%*%Cbeta(beta,index)/c(Al%*%gamma),2,sum) - 
      drop(as.matrix(K(sigma,K1,K2,K3) %*% beta)) + cbetat
  }
  

  
  
  normal = function(index,pila){
    return(exp(pila[[index]] - mean(pila[[index]]))/(calc_Z(xis_c2(index,maximos),pila[[index]] - mean(pila[[index]]))))
  }
  
  normal = function(index,pila){
    return(exp(pila[[index]] - pila[[index]][pila[[index]] == max(pila[[index]])])/(calc_Z(xis_c2(index,maximos),pila[[index]] - pila[[index]][pila[[index]] == max(pila[[index]])])))
  }
  

  
  calc_x0 = function(alpha,x0,tol=1e-12) {
    x = x0
    cont = 0
    while(1) {
      g1 = calc_grad_ff(x,alpha)
      cont = cont +1
      
      H = -calc_neg_hess_ff4(x,alpha)
      
      x = x - solve(H)%*%(g1)
      if (mean((x-x0)^2 < tol)) {
        break;
      }
      if(cont == 50){
        
        break;
      } else {
        x0 = x
      }
    }
    return(x)
  }
  
  library(Matrix)
  
  

  calc_lpost32 <- function(sigma, xin) {
    cache_env <- new.env()
    
    
    if (!is.null(cache_env$x0)) {
      xin <- cache_env$x0
    }
    
  
    x0 = calc_x02(sigma,xin,A,Al,indext,K1,K2,K3,X)
    
    
    cache_env$x0 <- x0
    
    
    H <- calc_neg_hess_ff4(x0, sigma, A, Al, K1, index, K2,K3,X)
    
    diag_chol_h <- diag(chol(H))  
    
    
    ljoint <- calc_ljoint(A, x0, sigma, Al, K1, K2, K3)
    log_det <- sum(log(diag_chol_h))
    
    lpost <- ljoint + (n / 2) * log(2 * pi) - log_det
    return(lpost)
  }
  
  
  calc_lpost = function(sigma) {
    x0 = calc_x0(sigma,x0 = xin)
    H = calc_neg_hess_ff4(x0,sigma)
    chol_h = chol(H)
    f = calc_ljoint(A,x0, sigma, Al,K1,K2,K3) + (n/2)*log(2*pi)- sum(log(diag(chol_h)))
    return(list(mode = x0, Hessian = H, lpost = f))
  }
  
  calc_lpost3 = function(sigma,xin) {
    x0 = calc_x0(sigma,x0 = xin)
    chol_h = calc_neg_hess_ff4(x0,sigma)
    calc_ljoint(A,x0, sigma, Al,K1,K2,K3) + (n/2)*log(2*pi) - 0.5*(log(det(chol_h)))
  }
  
  
  
  
  xin = calc_x02(c(0,0), rep(2,length(indext)),A,Al,index,K1,K2,K3,X)

  
  #cl <- makeCluster(spec=detectCores(), type = "FORK",outfile="")
  
  
  #system.time({ d = optimParallel(par = c(1,1), fn = calc_lpost32,xin = xin, method = "L-BFGS-B",
   #                               control=list(fnscale=-1,maxit = 5000,trace = TRUE),hessian = FALSE, parallel=list(loginfo=TRUE,cl = cl))})
  
  
  #stopCluster(cl)
  
  
  
  d = optim(par = c(1,1), fn = calc_lpost32,xin = xin, method = "L-BFGS-B",
                                  control=list(fnscale=-1,maxit = 5000,trace = TRUE),hessian = FALSE)
  
  
  
  sigmamod = d$par

  mu_ini = calc_x02(sigmamod,xin,A,Al,index,K1,K2,K3,X)
  H = chol2inv(chol(calc_neg_hess_ff4(mu_ini,sigmamod, A, Al, K1, index, K2,K3,X)))
  sigma_ini = sqrt(diag(H))
  
  xis_c = function(index){
    return(seq(mu_ini[index]-7.1*sqrt(sigma_ini[index]),mu_ini[index]+7.1*sqrt(sigma_ini[index]),length = 71))
  }
  
  
  
  xis_c2 = function(index,max){
    if(index %in% which(indext == FALSE)){
      return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
    }else{
      return(seq(max[index]-2.5*sqrt(sigma_ini[index]),max[index]+2.5*sqrt(sigma_ini[index]),length = 31))
    }
  }
  
  
  calc_Z = function(alpha_vec, lpost_vec) {
    nn = length(alpha_vec)
    hh = alpha_vec[2] - alpha_vec[1]
    ww = c(1, rep(c(4,2), (nn-2)/2), c(4,1))
    return(sum(ww * exp(lpost_vec)) * hh / 3)
  }
  
  teste = function(index,mu,sigma_beta,H,sigmanovo){
    xis = xis_c(index)
    pila = vector()
    for(j in 1:length(xis)){
      mul = mu[-index] + H[-index,index]*(1/H[index,index])*(xis[j] -mu[index])
      pila[j] = (((calc_ljoint(A,append(mul,xis[j],index-1),(sigmanovo),Ad,K1,K2,K3))))
    }
    return(list(pila = pila,xis = xis))
  }
  
  #cl <- makeCluster(spec=detectCores(), type="FORK", outfile="")
  #system.time({posti =  parLapply(cl,seq(1:(length(index))),teste,mu = mu_ini,sigma_beta = sigma_ini,H=H,sigmanovo = c(sigmamod))})
  posti =  lapply(seq(1:(length(index))),teste,mu = mu_ini,sigma_beta = sigma_ini,H=H,sigmanovo = c(sigmamod))
  
  maximos_ind = sapply(lapply(posti,sapply,which.max),'[[',1)
  maximos = vector()
  for(j in 1:length(mu_ini)){
    maximos[j] = sapply(posti,'[[',2,simplify = FALSE)[[j]][maximos_ind[j]]
  }
  
  #stopCluster(cl)
  
  
  n.fine <- 501
  nx = 6
  ny =6
  delta = 0.7
  #
  
  xy.coarse <- cbind( rep( seq(sigmamod[1] -delta,sigmamod[1] + delta,length = nx), each=ny), rep(seq(sigmamod[2] -delta,sigmamod[2] + delta,length = nx),ny ) )
  
  points_t = xy.coarse
  lista =  as.list(as.data.frame(t(xy.coarse)))
  
  post_x2 =  lapply(lista, function(alpha_) {
    
    mode_ = calc_x02(alpha_,xin,A,Al,indext,K1,K2,K3,X)
    H1 <- calc_neg_hess_ff4(mode_, alpha_, A, Al,
                              K1,index,K2,K3,X)
    chol_H_ = chol(H1)
    H = chol2inv(chol_H_)
    ta = lapply(seq(1:ncol(A)), function (ii){
      Append(matrix(kronecker.prod(xis_c2(ii,maximos) -mode_[ii],H[-ii,ii]*(1/H[ii,ii])),ncol = 31,nrow = ncol(A) - 1) + mode_[-ii]
             ,xis_c2(ii,maximos),after = ii - 1,rows = TRUE)})
    pila = apply_to_list_parallel(ta,A,alpha_,Al,K1,K2,K3,X,indext)
    posti = mclapply(seq(1:ncol(A)),normal,pila = pila,mc.cores = 1)
    return(list(
      sigma = points_t,
      maximos = maximos,
      H1 = H1,
      mu_ini = mu_ini,
      sigma_ini = sigma_ini,
      alpha = alpha_, 
      x0 = mode_,
      posti_cond = posti,
      diag_sigma = drop(diag(chol2inv(chol_H_))),
      unn_log_post = calc_ljoint(A,mode_, alpha_, Al,K1,K2,K3)- sum(log(diag(chol_H_)))
    ))
  })
  
  
  return(post_x2)
  
  
}






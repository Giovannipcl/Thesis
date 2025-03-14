
############################################################3
dens_est = function(x,j,mu,sigmas_beta,post_alpha,weights){
  return(sum(dnorm(x,mu[j,],sqrt(sigmas_beta[j,]))*post_alpha*weights))}



MH = function(index,mu_post, mu,sigmas_beta,post_alpha,weights){
  x = rep(0,6000)
  x[1] = mu_post[index]     #initialize; I've set arbitrarily set this to 3
  for(i in 2:6000){
    current_x = x[i-1]
    proposed_x = current_x + rnorm(1,mean=0,sd=1)
    #proposed_x = runif(1,current_x -0.2,current_x+0.2) 
    A = min(1,dens_est(proposed_x,index,mu,sigmas_beta,post_alpha,weights)/dens_est(current_x,index,mu,sigmas_beta,post_alpha,weights))
    #A = interp1(x = xis,y = norm, xi = proposed_x,method = "spline")/interp1(x = xis,y = norm, xi = current_x,method = "spline")
    if(runif(1)<A){
      x[i] = proposed_x       # accept move with probabily min(1,A)
    } else {
      x[i] = current_x        # otherwise "reject" move, and stay where we are
    }
  }
  
  return(x[seq(2000,6000,3)])
}

generate_uniform_grid <- function(dimension, num_points, min_values, max_values) {
  # Verifica se os vetores de min e max têm o tamanho correto
  if (length(min_values) != dimension || length(max_values) != dimension) {
    stop("Os vetores min_values e max_values devem ter o mesmo comprimento que a dimensão.")
  }
  
  # Gera a grade usando expand.grid
  ranges <- lapply(1:dimension, function(i) seq(from = min_values[i], to = max_values[i], length.out = num_points))
  grid <- do.call(expand.grid, ranges)
  
  # Calcula os pesos para cada célula da grade
  cell_volume <- prod((max_values - min_values) / (num_points - 1))  # Volume de cada célula da grade
  weights <- rep(cell_volume, nrow(grid))
  
  return(list(points = grid, weights = weights,ranges = ranges))
}


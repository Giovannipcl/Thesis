## Dependencies
library("tramME")
library("ggplot2")
library("lme4")
library("survival")

## Functions
mycolors <- function(nr, type = "line") {
  cols <- list()
  cols[[1]] <- c(red = 0, green = 84, blue = 150)
  cols[[2]] <- c(red = 202, green = 108, blue = 24)
  out <- as.list(cols[[nr]])
  out$alpha <- switch(type, line = 255L, fill = 140L)
  out$maxColorValue <- 255
  do.call("rgb", out)
}

## ----sleepstudy-quantiles1, echo=FALSE, fig.width=7, fig.height=5, out.width="\\linewidth"----
plot_cquant <- function(data, newdat, pred, quantiles) {
  pred <- c(pred)
  i <- 1:nrow(newdat)
  nd <- expand.grid(quantile = factor(quantiles), i = i)
  nd <- cbind(nd, newdat[nd$i,,drop = FALSE])
  nd$predict <- pred
  newdat <- nd
  newdat$Subject <- factor(paste("Subject =", newdat$Subject))
  data$Subject <- factor(paste("Subject =", data$Subject))
  ggplot(newdat, aes(x = Days, y = predict, color = quantile)) +
    geom_line() +
    geom_point(aes(x = Days, y = Reaction), data = data, color = "black", alpha = 0.4,
               size = 0.95) +
    geom_line(aes(x = Days, y = Reaction), data = data, color = "black", alpha = 0.4) +
    facet_wrap(~ Subject, nrow = 3) +
    scale_x_continuous("Days of sleep deprivation",
                       breaks = c(0, 2, 4, 6, 8)) +
    ylab("Average reaction time (ms)") +
    labs(colour = "Estimated conditional quantiles") +
    colorspace::scale_color_discrete_diverging("Blue-Red2") +
    theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.3,
                                          linetype = "dotted"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_rect(colour = "black", fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = rel(1.2)),
          axis.text.x = element_text(size = 16),  # Tamanho dos valores do eixo x
          axis.text.y = element_text(size = 16),  # Tamanho dos valores do eixo y
          legend.text = element_text(size = 16),
          legend.title = element_text(size = rel(1.5)),
          legend.position = "bottom",
          legend.key = element_rect(fill = "transparent", colour = "transparent")) +
    guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))
}


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


calculate_cell_weights <- function(Ypontos) {
  # Número de dimensões
  n_dims <- ncol(Ypontos)
  
  # Inicializar uma lista para armazenar os volumes de cada dimensão
  volumes <- list()
  
  # Calcular os volumes para cada dimensão (usando as diferenças entre os pontos)
  for (i in 1:n_dims) {
    diffs <- diff(Ypontos[, i])  # Diferenças entre pontos em cada dimensão
    volumes[[i]] <- c(diffs, diffs[length(diffs)])  # Adiciona a última distância para manter o tamanho da grade
  }
  
  # Criar uma lista de produtos de volumes para todas as combinações de dimensões
  cell_volumes <- volumes[[1]]
  
  for (i in 2:n_dims) {
    cell_volumes <- outer(cell_volumes, volumes[[i]], FUN = "*")
  }
  
  # Retornar os pesos wY (vetor com o volume de cada célula)
  wY <- as.vector(cell_volumes)
  
  return(wY)
}

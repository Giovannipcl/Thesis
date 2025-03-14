#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]

// Assuming gamma_f, K, and Cbeta are defined in Rcpp and are available for use here



double calc_lprior(const arma::vec& sigma) {
  double sum_log_density = 0.0;
  int n = sigma.size();
  
  for (int i = 0; i < n; ++i) {
    double x = std::exp(sigma[i]);
    double log_density = -std::log(R::gammafn(1.0)) + 1.0 * std::log(0.001) - (1.0 + 1.0) * std::log(x) - 0.001 / x;
    sum_log_density += log_density;
  }
  
  return sum_log_density;
}

// [[Rcpp::export]]
arma::mat K(const arma::vec& sigma, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X) {
    arma::mat part1 = K1 / std::exp(sigma(0));
    arma::mat part2 = K2 / std::exp(sigma(1));
    arma::mat part3 = K2 / std::exp(sigma(2));
    arma::mat part4 = K3 / std::exp(sigma(3));
    
    arma::mat diag_mat = arma::eye<arma::mat>(X.n_cols, X.n_cols) * 0.00001;
    
    int total_cols = part1.n_cols + part2.n_cols + part3.n_cols + part4.n_cols + diag_mat.n_cols;
    arma::mat result = arma::zeros<arma::mat>(total_cols, total_cols);
    
    int current_col = 0;
    
    result.submat(current_col, current_col, current_col + part1.n_rows - 1, current_col + part1.n_cols - 1) = part1;
    current_col += part1.n_cols;
    
    result.submat(current_col, current_col, current_col + part2.n_rows - 1, current_col + part2.n_cols - 1) = part2;
    current_col += part2.n_cols;
    
    result.submat(current_col, current_col, current_col + part3.n_rows - 1, current_col + part3.n_cols - 1) = part3;
    current_col += part3.n_cols;
    
     result.submat(current_col, current_col, current_col + part4.n_rows - 1, current_col + part4.n_cols - 1) = part4;
    current_col += part4.n_cols;
    
    result.submat(current_col, current_col, current_col + diag_mat.n_rows - 1, current_col + diag_mat.n_cols - 1) = diag_mat;
    
    return result;
}


// [[Rcpp::export]]
arma::vec gamma_f(const arma::vec& beta, const arma::uvec& index) {
  int n = index.size();
  arma::vec gammat = beta;
  
  for (int i = 0; i < n; ++i) {
    if (index[i]) {
      gammat[i] = std::exp(beta[i]);
    }
  }
  return gammat;
}

// [[Rcpp::export]]
arma::mat Cbeta(const arma::vec& beta, const arma::uvec& index) {
  // Initialize cbetat vector with ones
  arma::vec cbetat(beta.n_elem, arma::fill::ones);
  
  // Set elements specified by index to exp(beta[index])
  for (size_t i = 0; i < index.n_elem; ++i) {
    if (index[i]) {
      cbetat[i] = std::exp(beta[i]);
    }
  }
  
  // Create a diagonal matrix from cbetat
  arma::mat diagMatrix = arma::diagmat(cbetat);
  
  return diagMatrix;
}


// [[Rcpp::export]]
arma::mat Cbetal(const arma::vec& beta, const arma::uvec& index) {
  int n = index.size();
  arma::vec cbetat(n, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    if (index[i]) {
      cbetat[i] = std::exp(beta[i]);
    }
  }
  
  arma::mat diagMatrix = arma::diagmat(cbetat);
  return diagMatrix;
}



// [[Rcpp::export]]
double calc_ljoint(const arma::mat& A, const arma::vec& beta, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index) {
  // Calculate gamma
  arma::vec gamma = gamma_f(beta, index);
  
  // Compute Cholesky decomposition of K(sigma, K1, K2)
  arma::mat chol_Q = arma::chol(K(sigma, K1,K2,K3,X));
  
  // Calculate log determinant of chol_Q
  double logdet_Q_half = arma::sum(arma::log(chol_Q.diag()));
  
  // Compute the quadratic form
  double quad_form = arma::dot(chol_Q * beta, chol_Q * beta);
  
  // Calculate the result
  double res = 0.0;
  
  // Sum of log densities of standard normal
  arma::vec A_gamma = A * gamma;
  for (size_t i = 0; i < A_gamma.size(); ++i) {
    res += R::dnorm(A_gamma[i], 0.0, 1.0, true);
  }
  // Sum of logs
  arma::vec Al_gamma = Al * gamma;
  for (size_t i = 0; i < Al_gamma.size(); ++i) {
    res += std::log(Al_gamma[i]);
  }
  
  res += logdet_Q_half;
  res += -0.5 * quad_form;
  res += calc_lprior(sigma);
  
  for (size_t i = 0; i < beta.size(); ++i) {
    if (index[i]) {
      res += beta[i];
    }
  }
  
  for (size_t i = 0; i < sigma.size(); ++i) {
    res += sigma[i];
  }
  
  return res;
}

// [[Rcpp::export]]
arma::vec calc_grad_ff(const arma::vec& beta, const arma::vec& sigma,
                       const arma::mat& A, const arma::mat& Al,
                       const arma::uvec& index, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X) {
  // Calculate gamma
  arma::vec gamma = gamma_f(beta, index);
  
  // Initialize cbetat vector
  arma::vec cbetat(index.n_elem, arma::fill::zeros);
  
  // Set elements specified by index to 1
  for (size_t i = 0; i < index.n_elem; ++i) {
    if (index[i]) {
      cbetat[i] = 1.0;
    }
  }
  
  // Compute each term
  arma::mat Cbeta_mat = Cbeta(beta, index);
  arma::vec term1 = (A * Cbeta_mat).t() * (-A * gamma);
  arma::vec term2 = (Al * Cbeta_mat).t() * (1.0 / (Al * gamma));
  arma::vec term3 = K(sigma, K1,K2,K3,X) * beta;
  
  // Sum the terms and add cbetat
  arma::vec grad = term1 + term2 - term3 + cbetat;
  
  return grad;
}



// [[Rcpp::export]]
arma::mat calc_neg_hess_ff4(const arma::vec& beta, const arma::vec& sigma, const arma::mat& A, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index) {
  // Assuming the functions Cbeta, Cbetal, gamma_f, and K are defined elsewhere and are accessible
  arma::vec gamma = gamma_f(beta, index);
  arma::mat Cbeta_beta = Cbeta(beta, index);
  arma::mat Cbetal_beta = Cbetal(beta, index);
  
  arma::mat Hess = (A * Cbeta_beta).t() * (-A * Cbeta_beta) +
    diagmat((A * Cbetal_beta).t() * (-A * gamma)) +
    (diagmat(1 / (Al * gamma)) * Al * Cbeta_beta).t() * 
    (diagmat(1 / (Al * gamma)) * (-Al * Cbeta_beta)) +
    diagmat((Al * Cbetal_beta).t() * (1 / (Al * gamma)));
  
  return -(Hess - K(sigma, K1,K2,K3,X));
}

// [[Rcpp::export]]
arma::vec calc_x0(const arma::vec& alpha, arma::vec x0, const arma::mat& A, const arma::mat& Al,
                  const arma::uvec& index, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3,  const arma::mat& X, double tol = 1e-12) {
  arma::vec x = x0;
  int cont = 0;
  const int max_iterations = 50;
  
  while (true) {
    arma::vec g1 = calc_grad_ff(x, alpha,A,Al,index,K1,K2,K3,X);
    cont += 1;
    arma::mat H = -calc_neg_hess_ff4(x, alpha,A,Al,K1,K2,K3,X, index);
    
    // Update x
    x = x - arma::solve(H, g1);
    
    // Check for convergence
    if (arma::mean(arma::square(x - x0)) < tol) {
      break;
    }
    
    // Check for maximum iterations
    if (cont == max_iterations) {
      break;
    } else {
      x0 = x;
    }
  }
  
  return x;
}

// [[Rcpp::export]]
double calc_lpost_cpp(const arma::vec& sigma, const arma::vec& xin,
                      const arma::mat& A, const arma::mat& Al,
                      const arma::uvec& index, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const size_t n) {
  // Calculate x0 using the provided function
  arma::vec x0 = calc_x0(sigma, xin, A, Al, index, K1,K2,K3,X);

  // Calculate H using the provided function (H is a sparse matrix)
  arma::mat H = calc_neg_hess_ff4(xin, sigma, A, Al, K1,K2,K3, X,index);

  
  // Perform Cholesky decomposition of H
  arma::mat chol_h;
  bool success = arma::chol(chol_h, H);
  if (!success) {
    throw std::runtime_error("Cholesky decomposition failed.");
  }

  // Calculate the log determinant of chol_h
  double logdet_h = 0.0;
  for (size_t i = 0; i < chol_h.n_rows; ++i) {
    logdet_h += std::log(chol_h(i, i));
  }

  // Calculate the log joint part
  double log_joint = calc_ljoint(A, x0, sigma, Al, K1,K2,K3, X,index);

  // Compute the posterior log likelihood
  double lpost = log_joint + (n / 2.0) * std::log(2 * M_PI) - logdet_h;

  return lpost;
}

// [[Rcpp::export]]
double calc_lpost3(const arma::vec& sigma, const arma::vec& xin, const arma::mat& A, 
                   const arma::mat& Al, const arma::uvec& index, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, size_t n) {
  // Calculate x0 using the provided function
  arma::vec x0 = calc_x0(sigma, xin, A, Al, index, K1,K2,K3,X);

  // Calculate the Hessian matrix using the provided function
  arma::mat H = calc_neg_hess_ff4(xin, sigma, A, Al, K1,K2, K3,X,index);

  // Calculate the determinant of H
  double det_H = arma::det(H);

  // Calculate the log joint part
  double log_joint = calc_ljoint(A, x0, sigma, Al, K1,K2, K3,X,index);

  // Compute the posterior log likelihood
  double lpost = log_joint + (n / 2.0) * std::log(2 * M_PI) - 0.5 * std::log(det_H);

  return lpost;
}

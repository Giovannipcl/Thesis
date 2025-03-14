#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]

// Assuming gamma_f, K, and Cbeta are defined in Rcpp and are available for use here

// Function to compute the logarithm of the multivariate gamma function
double log_multivariate_gamma(int p, double a) {
  double result = p * (p - 1) / 4.0 * log(M_PI);
  for (int i = 0; i < p; ++i) {
    result += lgamma(a + (1.0 - (i + 1)) / 2.0);
  }
  return result;
}

// Function to compute the log density of the Wishart distribution
double dwishart(const arma::mat& W, double df, const arma::mat& scale, bool log_density = false) {
  int p = W.n_rows;
  double log_det_scale, sign;
  arma::log_det(log_det_scale, sign, scale);
  
  double log_det_W, sign_W;
  arma::log_det(log_det_W, sign_W, W);
  
  double trace_term = arma::trace(arma::inv(scale) * W);
  
  double log_density_val = 0.5 * (df - p - 1) * log_det_W - 0.5 * trace_term - 0.5 * df * p * log(2) - 0.5 * df * log_det_scale - log_multivariate_gamma(p, df / 2.0);
  
  if (log_density) {
    return log_density_val;
  } else {
    return exp(log_density_val);
  }
}

// Function to compute the log density of the inverse gamma distribution
//double dinvgamma(double x) {
//  double log_density = -std::log(R::gammafn(1.0)) + 1.0 * std::log(0.001) - (1.0 + 1.0) * std::log(x) - 0.001 / x;
  
//  return log_density;
//}


//double dinvgamma(double x) {
//  double log_density = -std::log(R::gammafn(0.01)) + 0.01 * std::log(0.001) + (0.01 - 1.0) * std::log(x) - 0.001 * x;
  
//  return log_density;
//}

// Function to compute the log density of the inverse gamma distribution
double dinvgamma(double x) {
  double log_density = -std::log(R::gammafn(1.0)) + 1.0 * std::log(0.001) - (1.0 + 1.0) * std::log(x) - 0.001 / x;
  
  return log_density;
}



// [[Rcpp::export]]
double calc_lprior(const arma::vec& sigma) {
  double rho = 2 * (exp(sigma(3)) / (1 + exp(sigma(3)))) - 1;
  arma::mat W(2, 2);
  W(0, 0) = exp(sigma(1));
  W(0, 1) = rho * (sqrt(exp(sigma(1)) * exp(sigma(2))));
  W(1, 0) = W(0, 1);
  W(1, 1) = exp(sigma(2));
  
  arma::mat scale = arma::eye(2, 2);
  
  arma::mat W_inv = arma::inv(W);
  arma::mat scale_inv = arma::inv(scale);
  
  double lprior = dwishart(W_inv, 4, scale_inv, true);
  
  
  lprior += dinvgamma(exp(sigma(0))) + dinvgamma(exp(sigma(1))) + dinvgamma(1/exp(sigma(2))) + dinvgamma(1/exp(sigma(3)));
  
  
  return lprior;
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

// Function to create a block diagonal matrix from a list of matrices
arma::mat bdiag(const std::vector<arma::mat>& matrices) {
  int total_rows = 0;
  int total_cols = 0;
  
  for (const auto& mat : matrices) {
    total_rows += mat.n_rows;
    total_cols += mat.n_cols;
  }
  
  arma::mat result(total_rows, total_cols, arma::fill::zeros);
  
  int current_row = 0;
  int current_col = 0;
  
  for (const auto& mat : matrices) {
    result.submat(current_row, current_col, current_row + mat.n_rows - 1, current_col + mat.n_cols - 1) = mat;
    current_row += mat.n_rows;
    current_col += mat.n_cols;
  }
  
  return result;
}

// [[Rcpp::export]]
arma::mat K(const arma::vec& sigma, const arma::mat& K1, int nf) {
  double rho = 2 * (exp(sigma(3)) / (1 + exp(sigma(3)))) - 1;
  arma::mat W(2, 2);
  W(0, 0) = exp(sigma(1));
  W(0, 1) = rho * (sqrt(exp(sigma(1)) * exp(sigma(2))));
  W(1, 0) = W(0, 1);
  W(1, 1) = exp(sigma(2));
  
  arma::mat W_inv = arma::inv(W);
  W_inv +=  arma::eye(2, 2);
  
  arma::mat K_part1 = (1/exp(sigma(0))) * K1;
  arma::mat K_part2 = arma::diagmat(arma::vec(2).fill(1e-6));
  arma::mat K_part3 = arma::kron(W_inv, arma::eye(nf, nf));
  
  std::vector<arma::mat> matrices = {K_part1, K_part2, K_part3};
  arma::mat K = bdiag(matrices);
  
  return K;
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
double calc_ljoint(const arma::mat& A, const arma::vec& beta, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::uvec& index, double nf) {
  // Calculate gamma
  arma::vec gamma = gamma_f(beta, index);
  
  // Compute Cholesky decomposition of K(sigma, K1, K2)
  arma::mat chol_Q = arma::chol(K(sigma, K1,nf));
  
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
                       const arma::uvec& index, const arma::mat& K1, double nf) {
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
  arma::vec term3 = K(sigma, K1,nf) * beta;
  
  // Sum the terms and add cbetat
  arma::vec grad = term1 + term2 - term3 + cbetat;
  
  return grad;
}



// [[Rcpp::export]]
arma::mat calc_neg_hess_ff4(const arma::vec& beta, const arma::vec& sigma, const arma::mat& A, const arma::mat& Al, const arma::mat& K1,const arma::uvec& index, double nf) {
  // Assuming the functions Cbeta, Cbetal, gamma_f, and K are defined elsewhere and are accessible
  arma::vec gamma = gamma_f(beta, index);
  arma::mat Cbeta_beta = Cbeta(beta, index);
  arma::mat Cbetal_beta = Cbetal(beta, index);
  
  arma::mat Hess = (A * Cbeta_beta).t() * (-A * Cbeta_beta) +
    diagmat((A * Cbetal_beta).t() * (-A * gamma)) +
    (diagmat(1 / (Al * gamma)) * Al * Cbeta_beta).t() * 
    (diagmat(1 / (Al * gamma)) * (-Al * Cbeta_beta)) +
    diagmat((Al * Cbetal_beta).t() * (1 / (Al * gamma)));
  
  return -(Hess - K(sigma, K1,nf));
}

// [[Rcpp::export]]
arma::vec calc_x0(const arma::vec& alpha, arma::vec x0, const arma::mat& A, const arma::mat& Al,
                  const arma::uvec& index, const arma::mat& K1, double nf, double tol = 1e-12) {
  arma::vec x = x0;
  int cont = 0;
  const int max_iterations = 80;
  
  while (true) {
    arma::vec g1 = calc_grad_ff(x, alpha,A,Al,index,K1,nf);
    cont += 1;
    arma::mat H = -calc_neg_hess_ff4(x, alpha,A,Al,K1, index,nf);
    
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
                      const arma::uvec& index, const arma::mat& K1, const size_t n, double nf) {
  // Calculate x0 using the provided function
  arma::vec x0 = calc_x0(sigma, xin, A, Al, index, K1,nf);
  
  // Calculate H using the provided function (H is a sparse matrix)
  arma::mat H = calc_neg_hess_ff4(xin, sigma, A, Al, K1, index,nf);
  
  
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
  double log_joint = calc_ljoint(A, x0, sigma, Al, K1, index,nf);
  
  // Compute the posterior log likelihood
  double lpost = log_joint + (n / 2.0) * std::log(2 * M_PI) - logdet_h;
  
  return lpost;
}

// [[Rcpp::export]]
double calc_lpost3(const arma::vec& sigma, const arma::vec& xin, const arma::mat& A, 
                   const arma::mat& Al, const arma::uvec& index, const arma::mat& K1, size_t n, double nf) {
  // Calculate x0 using the provided function
  arma::vec x0 = calc_x0(sigma, xin, A, Al, index, K1,nf);
  
  // Calculate the Hessian matrix using the provided function
  arma::mat H = calc_neg_hess_ff4(xin, sigma, A, Al, K1, index,nf);
  
  // Calculate the determinant of H
  double det_H = arma::det(H);
  
  // Calculate the log joint part
  double log_joint = calc_ljoint(A, x0, sigma, Al, K1, index,nf);
  
  // Compute the posterior log likelihood
  double lpost = log_joint + (n / 2.0) * std::log(2 * M_PI) - 0.5 * std::log(det_H);
  
  return lpost;
}

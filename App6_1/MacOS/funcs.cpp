#include <RcppArmadillo.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// [[Rcpp::depends(RcppArmadillo,RcppParallel)]]

// Assuming gamma_f, K, and Cbeta are defined in Rcpp and are available for use here


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


arma::mat K(const arma::vec& sigma, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X) {
    // Compute exp(sigma[0]) * (K1 + K2)
    arma::mat part1 = exp(sigma[0]) * (K1 + K2);
    
    // Compute exp(sigma[1]) * K3
    arma::mat part2 = exp(sigma[1]) * K3;
    
    // Create a diagonal matrix with 10^-6 on the diagonal, with the same number of columns as X
    arma::mat diagPart = arma::diagmat(arma::vec(X.n_cols, arma::fill::ones) * 1e-6);
    
    // Combine the matrices into a block diagonal matrix
    arma::mat result = arma::join_cols(
        arma::join_rows(part1, arma::zeros(part1.n_rows, part2.n_cols + diagPart.n_cols)),
        arma::join_rows(arma::zeros(part2.n_rows, part1.n_cols), part2, arma::zeros(part2.n_rows, diagPart.n_cols)),
        arma::join_rows(arma::zeros(diagPart.n_rows, part1.n_cols + part2.n_cols), diagPart)
    );
    
    return result;
}


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
  

  
  return res;
}


// [[Rcpp::export]]
double calc_ff_cpp(const arma::mat& A, const arma::vec& beta, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index) {
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

  res += -0.5 * quad_form;
  
  for (size_t i = 0; i < beta.size(); ++i) {
    if (index[i]) {
      res += beta[i];
    }
  }
  
  return res;
}


// [[Rcpp::export]]
arma::vec calc_grad_ff_cpp(const arma::vec& beta, const arma::vec& sigma,
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
arma::mat calc_neg_hess_ff4(const arma::vec& beta, const arma::vec& sigma, const arma::mat& A, const arma::mat& Al, const arma::mat& K1,const arma::uvec& index, const arma::mat& K2,  const arma::mat& K3,const arma::mat& X) {
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
arma::vec calc_x02(const arma::vec& alpha, arma::vec x0, const arma::mat& A, const arma::mat& Al,
                  const arma::uvec& index, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, double tol = 1e-12) {
  arma::vec x = x0;
  int cont = 0;
  const int max_iterations = 50;
  
  while (true) {
    arma::vec g1 = calc_grad_ff_cpp(x, alpha,A,Al,index,K1,K2,K3,X);
    cont += 1;
    arma::mat H = -calc_neg_hess_ff4(x, alpha,A,Al,K1, index,K2,K3,X);
    
    // Imprimir o contador de iterações
    //Rcpp::Rcout << "Iteration: " << cont << std::endl;

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


double calc_lpost_cpp(const arma::vec& sigma, const arma::vec& xin,
                   const arma::mat& A, const arma::mat& Al,
                   const arma::mat& K1, const arma::uvec& index,const arma::mat& K2, const arma::mat& K3,const arma::mat& X) {
  arma::vec x0 = calc_x02(sigma, xin, A,Al,index,K1,K2,K3,X);
  arma::mat chol_h = calc_neg_hess_ff4(x0, sigma,A,Al,K1,index,K2,K3,X);
  double log_det_chol_h;
  double sign;
  arma::log_det(log_det_chol_h, sign, chol_h);
  
  double result = calc_ljoint(A, x0, sigma, Al, K1,K2,K3,X, index)  - 0.5 * log_det_chol_h;
  
  return result;
}


// [[Rcpp::export]]
arma::vec apply_calc_ljoint(const arma::mat& Teste, const arma::mat& A, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index) {
    int n_cols = Teste.n_cols;
    arma::vec results(n_cols);
    
    for (int i = 0; i < n_cols; ++i) {
        arma::vec col = Teste.col(i);
        results(i) = calc_ljoint(A, col, sigma, Al, K1, K2, K3,X,index);
    }
    
    return results;
}


// [[Rcpp::export]]
List apply_to_list(const List& ta, const arma::mat& A, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index) {
    int n = ta.size();
    List results(n);
    
    for (int i = 0; i < n; ++i) {
        arma::mat Teste = as<arma::mat>(ta[i]);
        results[i] = apply_calc_ljoint(Teste, A, sigma, Al, K1, K2, K3,X,index);
    }
    
    return results;
}


// Worker for parallel computation
struct CalcLJointWorker : public Worker {
    const std::vector<arma::mat>& ta;
    const arma::mat& A;
    const arma::vec& sigma;
    const arma::mat& Al;
    const arma::mat& K1;
    const arma::mat& K2;
    const arma::mat& K3;
    const arma::mat& X;
    const arma::uvec& index;
    std::vector<arma::vec>& results;

    CalcLJointWorker(const std::vector<arma::mat>& ta, const arma::mat& A, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3, const arma::mat& X, const arma::uvec& index, std::vector<arma::vec>& results)
        : ta(ta), A(A), sigma(sigma), Al(Al), K1(K1), K2(K2), K3(K3),X(X),index(index), results(results) {}

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            const arma::mat& Teste = ta[i];
            int n_cols = Teste.n_cols;
            arma::vec result(n_cols);
            for (int j = 0; j < n_cols; ++j) {
                arma::vec col = Teste.col(j);
                result(j) = calc_ljoint(A, col, sigma, Al, K1, K2, K3,X,index);
            }
            results[i] = result;
        }
    }
};

// [[Rcpp::export]]
List apply_to_list_parallel(const List& ta, const arma::mat& A, const arma::vec& sigma, const arma::mat& Al, const arma::mat& K1, const arma::mat& K2, const arma::mat& K3,const arma::mat& X, const arma::uvec& index) {
    int n = ta.size();
    std::vector<arma::mat> ta_vec(n);
    for (int i = 0; i < n; ++i) {
        ta_vec[i] = as<arma::mat>(ta[i]);
    }
    std::vector<arma::vec> results(n);

    CalcLJointWorker worker(ta_vec, A, sigma, Al, K1, K2, K3, X,index, results);
    parallelFor(0, n, worker);

    List res_list(n);
    for (int i = 0; i < n; ++i) {
        res_list[i] = results[i];
    }

    return res_list;
}


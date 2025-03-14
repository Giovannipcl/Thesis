#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double der_cov_log_normal2(int i, int j, const mat& Sigma, const vec& mu, const uvec& index) {
  if (!index[i - 1] && !index[j - 1]) {
    return 0;
  } else if (!index[i - 1] && index[j - 1]) {
    return((exp( + 0.5*Sigma(i-1,i-1)*(pow((Sigma(j-1,i-1))/(Sigma(i-1,i-1)),2))))*(exp(mu[j-1]  + 0.5*(Sigma(j-1,j-1) - Sigma(j-1,i-1)*(1/(Sigma(i-1,i-1)))*(Sigma(i-1,j-1)))))  - (exp(mu[j-1] + Sigma(j-1,j-1)/2)));
  } else if (index[i - 1] && index[j - 1]) {
    return(exp(mu[i-1] + mu[j-1] + 0.5*(Sigma(i-1,i-1) + Sigma(j-1,j-1)))*(exp(Sigma(i-1,j-1)) - 1));
  } else {
   return((mu[j-1] + ((Sigma(i-1,j-1))/(Sigma(j-1,j-1)))*Sigma(j-1,j-1))*(exp(+ 0.5*Sigma(j-1,j-1)*(pow((Sigma(i-1,j-1))/(Sigma(j-1,j-1)),2))))*(exp(mu[i-1] + 0.5*(Sigma(i-1,i-1) - Sigma(i-1,j-1)*(1/(Sigma(j-1,j-1)))*(Sigma(j-1,i-1)))))  - mu[j-1]*(exp(mu[i-1] + Sigma(i-1,i-1)/2)));
  }
}

// [[Rcpp::export]]
double cov_log_normal2(int i, int j, const mat& Sigma, const vec& mu, const uvec& index) {
  if (!index[i - 1] && !index[j - 1]) {
    return Sigma(i - 1, j - 1);
  } else if (!index[i - 1] && index[j - 1]) {
    return (mu[i - 1] + (Sigma(j - 1, i - 1) / Sigma(i - 1, i - 1)) * Sigma(i - 1, i - 1)) *
           exp(0.5 * Sigma(i - 1, i - 1) * pow((Sigma(j - 1, i - 1) / Sigma(i - 1, i - 1)), 2)) *
           exp(mu[j - 1] + 0.5 * (Sigma(j - 1, j - 1) - (Sigma(j - 1, i - 1) / Sigma(i - 1, i - 1)) * Sigma(i - 1, j - 1))) -
           mu[i - 1] * exp(mu[j - 1] + Sigma(j - 1, j - 1) / 2);
  } else if (index[i - 1] && index[j - 1]) {
    return exp(mu[i - 1] + mu[j - 1] + 0.5 * (Sigma(i - 1, i - 1) + Sigma(j - 1, j - 1))) *
           (exp(Sigma(i - 1, j - 1)) - 1);
  } else {
    return (mu[j - 1] + (Sigma(i - 1, j - 1) / Sigma(j - 1, j - 1)) * Sigma(j - 1, j - 1)) *
           exp(0.5 * Sigma(j - 1, j - 1) * pow((Sigma(i - 1, j - 1) / Sigma(j - 1, j - 1)), 2)) *
           exp(mu[i - 1] + 0.5 * (Sigma(i - 1, i - 1) - (Sigma(i - 1, j - 1) / Sigma(j - 1, j - 1)) * Sigma(j - 1, i - 1))) -
           mu[j - 1] * exp(mu[i - 1] + Sigma(i - 1, i - 1) / 2);
  }
}

// [[Rcpp::export]]
vec diagonalProduct(const mat& A, const mat& B) {
  int n = A.n_rows;
  int m = B.n_cols;
  
  vec result(n);
  
  for (int i = 0; i < n; ++i) {
    double sum = 0;
    for (int j = 0; j < m; ++j) {
      sum += A(i, j) * B(j, i); // Multiplying corresponding elements of A and B
    }
    result[i] = sum;
  }
  
  return result;
}

// [[Rcpp::export]]
mat diagonalMatrixMultiply(const vec& diag, const mat& dense) {
  int n = diag.n_elem;
  int m = dense.n_cols;
  mat result(n, m);
  
  for (int i = 0; i < n; ++i) {
    double diagonal_element = diag[i];
    for (int j = 0; j < m; ++j) {
      result(i, j) = diagonal_element * dense(i, j); // Multiply diagonal element with corresponding column of dense matrix
    }
  }
  
  return result;
}

// [[Rcpp::export]]
vec diagonalQuadraticForm(const mat& X, const mat& A) {
  int n = X.n_rows;
  int m = X.n_cols;
  
  vec result(m);
  
  for (int i = 0; i < m; ++i) {
    double sum = 0;
    for (int j = 0; j < n; ++j) {
      sum += X(j, i) * A(j, j) * X(j, i); // Compute diagonal element of X^TAX
      for (int k = j + 1; k < n; ++k) {
        sum += 2 * X(j, i) * A(j, k) * X(k, i); // Compute off-diagonal element of X^TAX
      }
    }
    result[i] = sum;
  }
  
  return result;
}

// [[Rcpp::export]]
vec computeExpression4(const mat& Al, const mat& Sigma_log_ti, const vec& EY) {
  int p = EY.n_elem; // Dimensionality
  int N = Al.n_rows; // Number of rows in Al
  
  // Step 1: Compute Al * EY
  vec AlEY = Al * EY;
  
  // Step 2: Compute the sum of Al * EY
  double sum_AlEY = accu(AlEY);
  
  // Step 3: Compute the element-wise reciprocal squared of the sum from step 2
  double reciprocal_squared_sum_AlEY = 1.0 / (sum_AlEY * sum_AlEY);
  
  // Step 4: Compute the transpose of Al
  mat tAl = Al.t();
  
  // Step 5: Compute the transpose of Sigma_log_ti
  mat tSigma_log_ti = Sigma_log_ti.t();
  
  // Step 6: Multiply tSigma_log_ti with tAl
  mat tSigma_log_ti_tAl = tSigma_log_ti * tAl;
  
  // Step 7: Compute the element-wise multiplication of tAl with the result from step 6
  mat element_wise_mult = tAl % tSigma_log_ti_tAl;
  
  // Step 8: Compute the row sums of the resulting matrix from step 7
  vec row_sums = sum(element_wise_mult, 1);
  
  // Step 9: Multiply the result from step 8 by 0.5 and the result from step 3
  vec result = 0.5 * row_sums * reciprocal_squared_sum_AlEY;
  
  return result;
}

// [[Rcpp::export]]
mat transposeMatrix(const mat& mat) {
  return mat.t();
}

// [[Rcpp::export]]
mat matrixMultiply(const mat& A, const mat& B) {
  return A * B;
}

// [[Rcpp::export]]
vec computeExpression(const vec& EYl, const mat& A, const vec& EY) {
  int p = EYl.n_elem; // Dimensionality

  // Compute t(A) * A
  mat AtA = A.t() * A;

  // Compute diag(EYl) * AtA
  mat diag_EYl_AtA = diagmat(EYl) * AtA;

  // Compute 2 * diag(EYl) * AtA * EY
  vec result =  diag_EYl_AtA * EY;

  return result;
}

// [[Rcpp::export]]
vec computeExpression2(const vec& EYl, const mat& Al) {
  int p = EYl.n_elem; // Dimensionality
  int N = Al.n_rows;  // Number of rows in Al

  // Step 1: Compute Al * EYl
  vec AlEYl = Al * EYl;

  // Step 2: Compute 1 / c(Al * EYl)
  vec reciprocal_AlEYl = 1.0 / AlEYl;

  // Step 3: Compute t(Al)
  mat tAl = Al.t();

  // Step 4: Compute diag(c(EYl))
  mat diag_EYl = diagmat(EYl);

  // Step 5: Perform the multiplications diag(c(EYl)) * t(Al)
  mat diag_EYl_tAl = diag_EYl * tAl;

  // Step 6: Compute the final result: (diag(c(EYl)) * t(Al)) * (1 / c(Al * EYl))
  vec result = diag_EYl_tAl * reciprocal_AlEYl;

  return result;
}

/// [[Rcpp::export]]
vec computeExpression3(const vec& EY, const vec& EYl, const mat& Al, const mat& Sigma_log) {
  int p = EY.n_elem; // Dimensionality
  int N = Al.n_rows;  // Number of rows in Al

  // Step 1: Compute Al * EY
  vec AlEY = Al * EY;

  // Step 2: Compute 1 / (c(Al * EY)^3)
  vec reciprocal_AlEY_cubed = 1.0 / pow(AlEY, 3);

  // Step 3: Compute t(Al)
  mat tAl = Al.t();

  // Step 4: Compute diag(EYl) * t(Al)
  mat diag_EYl = diagmat(EYl);
  mat diag_EYl_tAl = diag_EYl * tAl;
  


  // Step 5: Compute diagonal quadratic form diagonalQuadraticForm(t(Al), Sigma_log)
  vec diag_quadratic_form = diagonalQuadraticForm(tAl, Sigma_log);

  // Step 6: Multiply diagonal quadratic form by reciprocal_AlEY_cubed
  vec scaled_diag_quadratic_form = diag_quadratic_form % reciprocal_AlEY_cubed;

  // Step 7: Compute diagonalMatrixMultiply
  mat diag_matrix_mult = diagonalMatrixMultiply(scaled_diag_quadratic_form,diag_EYl_tAl.t());

  // Step 8: Compute column sums of the resulting matrix
  vec col_sums = sum(diag_matrix_mult, 0).t();

  return col_sums;
}

// [[Rcpp::export]]
mat diagmat_product(const vec& diag_elements, const mat& B) {
  // Create a diagonal matrix from the given diagonal elements
  mat D = diagmat(diag_elements);

  // Multiply the diagonal matrix D with the matrix B
  mat result = D * B;

  return result;
}

// [[Rcpp::export]]
mat diagmat_tcrossprod(const vec& diag_elements, const mat& A) {
  // Create a diagonal matrix from the given diagonal elements
  mat D = diagmat(diag_elements);

  // Compute the cross-product D %*% t(D)
  mat result = D * A.t();

  return result;
}

// [[Rcpp::export]]
mat tcrossprod_arma(const mat& A) {
  // Compute the cross-product A %*% t(A)
  mat result = A * A.t();

  return result;
}

// [[Rcpp::export(rng = false)]]
arma::vec gradiente_f(const arma::vec& delta, const arma::vec& mu_ini, const arma::uvec& index, 
                      const arma::mat& H, const arma::mat& A, const arma::mat& Al, const arma::mat& K_p) {
  
  // Compute mu_delta
  arma::vec mu_delta = mu_ini + delta;
  
  // Compute EY, EYl
  arma::vec EY = mu_delta;
  arma::vec EYl = arma::ones(mu_ini.size());
  for (int i = 0; i < index.size(); ++i) {
    if (index[i]) {
      EY[i] = exp(mu_delta[i] + H(i, i) / 2);
      EYl[i] = EY[i];
    }
  }
  
  // Compute Sigma_log
  int n = index.size();
  arma::mat Sigma_log(n, n, fill::zeros);
  arma::mat Sigma_log_ti(n, n, fill::zeros);
  arma::mat Sigma_log_ti2(n, n, fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      Sigma_log(i, j) = cov_log_normal2(i + 1, j + 1, H, mu_delta, index);
      Sigma_log_ti(j, i) = der_cov_log_normal2(i + 1, j + 1, H, mu_delta, index);
    }
  }
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        Sigma_log_ti2(i, j) = Sigma_log_ti(i, j)*2.0;
      }
    }
  }
  
  
  // Compute gradients
  arma::vec gr(mu_ini.size(), fill::zeros);
  
  // Compute terms individually and populate gr
  arma::mat Aq = A.t() * A;
  arma::vec AtA_Sigma_log_ti = diagonalProduct(Aq, Sigma_log_ti);
  

  gr = AtA_Sigma_log_ti + computeExpression(EYl, A, EY) - computeExpression2(EYl, Al) - 
    computeExpression3(EY, EYl, Al, Sigma_log) + computeExpression4(Al, Sigma_log_ti2, EY) + (K_p * mu_delta);

  return gr;
}


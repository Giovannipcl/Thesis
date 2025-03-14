#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// Compute matrix multiplication (A' * B)
// [[Rcpp::export(rng = false)]]
mat matrixMultiply(const mat& A, const mat& B) {
  return A.t() * B; // Compute A' * B
}

// Compute diagonal elements of a matrix
vec getDiagonal(const mat& M) {
  return diagvec(M); // Extract diagonal elements
}

// Compute the sum of elements in a vector
double sumVector(const vec& v) {
  return accu(v); // Compute sum using Armadillo
}

// Compute the sum of elements in a matrix
double sumMatrix(const mat& M) {
  return accu(M); // Compute sum using Armadillo
}

// Compute transpose of a matrix
mat transposeMatrix(const mat& M) {
  return M.t(); // Compute transpose
}


vec diagonalQuadraticForm(const mat& X, const mat& A, const uvec& index) {
  uword m = X.n_cols;
  vec result(m);

  for (uword i = 0; i < m; ++i) {
    double sum = 0.0;
    for (uword j = 0; j < X.n_rows; ++j) {
      if (index[j]) {
        sum += X(j, i) * A(j, j) * X(j, i);
        for (uword k = j + 1; k < X.n_rows; ++k) {
          if (index[k]) {
            sum += 2 * X(j, i) * A(j, k) * X(k, i);
          }
        }
      }
    }
    result[i] = sum;
  }

  return result;
}

vec calculateExpression(const mat& Al, const vec& EY, const uvec& index) {
  uvec selectedCols = find(index);
  vec result(Al.n_rows);

  for (uword i = 0; i < Al.n_rows; ++i) {
    double sum = 0.0;
    for (uword j = 0; j < selectedCols.n_elem; ++j) {
      sum += Al(i, selectedCols[j]) * EY[selectedCols[j]];
    }
    result[i] = 1 / (2 * pow(sum, 2));
  }

  return result;
}


vec calculateExpression2(const mat& Al, const vec& EY, const uvec& index) {
  uword n = Al.n_rows;
  uword selectedCols = sum(index); // Count the number of TRUE values in 'index'

  vec result(n, fill::zeros); // Initialize result vector with zeros

  vec selectedEY(selectedCols); // Vector to store selected EY values
  
  uword col_idx = 0;
  for (uword j = 0; j < Al.n_cols; ++j) {
    if (index[j]) {
      selectedEY[col_idx] = EY[j];
      for (uword i = 0; i < n; ++i) {
        result[i] += Al(i, j) * EY[j]; // Element-wise multiplication and summation
      }
      col_idx++;
    }
  }

  return result;
}


vec applyTransformation(const mat& Al, const vec& EY, const uvec& index) {
  vec multipliedResult = calculateExpression2(Al, EY, index);
  return 1 / (2 * pow(multipliedResult, 2));
}

// [[Rcpp::export(rng = false)]]
double computeExpression(const vec& mu_ini, const mat& Aq) {
  double result = 0.0;

  for (uword i = 0; i < mu_ini.n_elem; ++i) {
    vec temp = mu_ini % (Aq.col(i) * mu_ini[i]);
    result += accu(temp);
  }

  return result;
}

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

// [[Rcpp::export(rng = false)]]
double arg_min22(const vec& delta, const vec& mu_ini, const mat& H, const uvec& index,
                 const mat& A, const mat& Al, const mat& K_p) {

  vec mu_delta = mu_ini + delta;
  vec EY = mu_delta;
  
  for (uword i = 0; i < index.n_elem; ++i) {
    if (index[i]) {
      EY[i] = exp(mu_delta[i] + H(i, i) / 2);
    }
  }

  mat Sigma_log(H.n_rows, H.n_cols);
  for (uword i = 0; i < H.n_rows; ++i) {
    for (uword j = 0; j < H.n_cols; ++j) {
      Sigma_log(i, j) = cov_log_normal2(i + 1, j + 1, H, mu_delta, index);
    }
  }
  
  // Function to compute Al_EY_index
  
  uword nrows = Al.n_rows;
  vec Al_EY_index(nrows); // Create a vector to store the result

  for (uword i = 0; i < nrows; ++i) {
    double sum = 0.0; // Initialize the sum for each row
    for (uword j = 0; j < Al.n_cols; ++j) {
      if (index[j]) { // Check if the jth column should be included based on 'index'
        sum += Al(i, j) * EY[j]; // Perform element-wise multiplication and sum
      }
    }
    Al_EY_index[i] = sum; // Store the sum in the result vector
  }




  double Eyg = sum(log(Al_EY_index)) -sum(diagonalQuadraticForm(transposeMatrix(Al),Sigma_log, index) % applyTransformation(Al, EY, index));
//

  mat Aq = matrixMultiply(A, A);
  double ll_soma = -(-0.5 * (sumVector(getDiagonal(matrixMultiply(Aq, Sigma_log))) + computeExpression(EY, Aq)) + Eyg);
//
 double KLD = computeExpression(mu_delta, K_p);
//
  return ((ll_soma + 0.5 * KLD));
//    return(Eyg);
}


#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]

// Função para calcular gamma_f
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


// Função para gerar amostras de uma normal multivariada
arma::mat mvrnormArma(const arma::vec& mu, const arma::mat& sigma, int n) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu.t(), n, 1) + Y * arma::chol(sigma);
}

// Função para calcular densidade_bctm
// [[Rcpp::export]]
double densidade_bctm(const arma::uword ind,
                      const arma::mat& A,
                      const arma::mat& beta_cpo,
                      const arma::mat& Ad) {
  arma::rowvec a_row = A.row(ind);  // Linha i de A (1x52)
  arma::rowvec ad_row = Ad.row(ind);  // Linha i de Ad (1x52)

  
  // Calcular a multiplicação de a_row e beta_cpo, e ad_row e beta_cpo
  arma::mat a_row_beta_cpo = a_row * beta_cpo.t();  // Resultado: 1x1000
  arma::mat ad_row_beta_cpo = ad_row * beta_cpo.t();  // Resultado: 1x1000
  
  // Calcular a densidade usando a normal pdf
  arma::mat norm_pdf_a_row_beta_cpo = arma::normpdf(a_row_beta_cpo);  // Resultado: 1x1000
  
  // Realizar a multiplicação element-wise e calcular a densidade
  arma::vec densities = (1.0 / (norm_pdf_a_row_beta_cpo % ad_row_beta_cpo)).t();  // Resultado: 1x1000
  
  // Agora você pode usar densities para continuar seu cálculo
  return arma::mean(densities);  // Média de todas as densidades
  
  
}

// Função principal para calcular CPO e PIT
// [[Rcpp::export]]
List calculate_cpo_pit(const arma::mat& A,
                       const arma::mat& Ad,
                       const arma::mat& Aq,
                       const arma::mat& Ayq,
                       const arma::vec& yseq,
                       const arma::vec& reaction,
                       const arma::mat& mu,
                       const arma::mat& pontos2,
                       const arma::vec& post_alpha,
                       const arma::vec& weights,
                       const arma::uvec& indext,
                       const arma::mat& K1,
                       double nf) {
  
  int n_rows = 4;
  int k_values = post_alpha.n_elem;
  arma::vec CPO_values_all(n_rows, arma::fill::zeros);
  arma::vec PIT_values_all(n_rows, arma::fill::zeros);
  
  for (int i = 0; i < n_rows; ++i) {
    arma::vec CPO_values(k_values, arma::fill::zeros);
    arma::vec PIT_values(k_values, arma::fill::zeros);
    
    for (int k = 0; k < k_values; ++k) {
      // Calcular a matriz Hessiana inversa
      arma::rowvec pontos2_row = pontos2.row(k); // Selecionar a coluna k
      arma::mat sigma = arma::inv(calc_neg_hess_ff4(mu.col(k), pontos2_row.t(), A, Ad, K1, indext, nf));
      
//     // Gerar amostras de beta_cpo
     arma::mat beta_cpo = mvrnormArma(mu.col(k), sigma, 1000);
//     
//     // Aplicar gamma_f a cada coluna de beta_cpo
     for (int j = 0; j < beta_cpo.n_cols; ++j) {
       beta_cpo.col(j) = gamma_f(beta_cpo.col(j), indext);
     }
//     
//     // Selecionar valores onde yseq < reaction[i]
     arma::uvec indic_true = find(yseq < reaction[i]);
     arma::vec newyseq = yseq.elem(indic_true);
     arma::mat Aqn = Aq.rows(indic_true);
     arma::mat Aqnd = Ayq.rows(indic_true);
     arma::vec diff_newyseq = arma::diff(newyseq);
     
     // Pegue o primeiro valor da diferença, já que todos são iguais
     double diff_value = diff_newyseq[0];  // ou simplesmente diff_newyseq(0)
//     
//     // Calcular h1 e h1l
     arma::vec h1 = arma::mean(Aqn * beta_cpo.rows(0, Aqn.n_cols - 1), 1);
     arma::vec h1l = arma::mean(Aqnd * beta_cpo.rows(0, Aqn.n_cols - 1), 1);
//     
//     // Calcular h2
     arma::rowvec a_row = A.row(i);
     arma::vec h2 = arma::mean(a_row.cols(Aqn.n_cols, A.n_cols - 1) *
       beta_cpo.rows(Aqn.n_cols, A.n_cols - 1), 1);
//     
//     // Calcular fy
    arma::vec fy = arma::normpdf(h1 + arma::ones(h1.n_elem) * h2) % h1l;
//     
//     // Calcular densidade_bctm
     double fyi = densidade_bctm(i, A, beta_cpo, Ad);
//     
//     // Calcular PIT e CPO para k
     PIT_values[k] = arma::sum(diff_value * fy) * fyi;
     CPO_values[k] = post_alpha[k] * fyi * weights[0];
   }
//   
//   // Armazenar resultados para i
   CPO_values_all[i] = 1.0 / arma::sum(CPO_values);
   PIT_values_all[i] = arma::sum(PIT_values % post_alpha * weights[0]) * CPO_values_all[i];
 }
// 
  return List::create(Named("CPO_values_all") = CPO_values_all,
                      Named("PIT_values_all") = PIT_values_all);
}

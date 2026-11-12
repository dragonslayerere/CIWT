#include <RcppArmadillo.h>
#include "CIWTStatistic.h"
#include "CIWTIWT.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Categorical Interval-Wise Test (CIWT)
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame ciwt_cpp_backend(const arma::cube& R_in, const arma::cube& E_in, int B, const arma::vec& dt) {
    // 1. Run statistical permutation engine
    CIWTStatistic stat_engine(R_in, E_in, B);
    stat_engine.generate_rademacher();
    stat_engine.compute_score_statistics();

    arma::mat S_obs = stat_engine.get_S_obs();
    arma::mat S_perm = stat_engine.get_S_perm();

    // 2. Run Interval-Wise Testing adjustments
    CIWTIWT iwt_engine(S_obs, S_perm, B, dt);
    iwt_engine.compute_interval_pvalues();
    iwt_engine.compute_adjusted_pvalues();

    arma::vec p_pointwise = iwt_engine.get_pointwise_pvalues().t();
    arma::vec p_adj = iwt_engine.get_results().t();
    
    return Rcpp::DataFrame::create(
        Rcpp::Named("P_Value_Pointwise") = p_pointwise,
        Rcpp::Named("P_Value_IWT") = p_adj
    );
}
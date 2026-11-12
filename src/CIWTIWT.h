//
// Created by albor on 31/08/2026.
//

#ifndef CIWT_IWT_H
#define CIWT_IWT_H

#include "CIWTStatistic.h"
// [[Rcpp::depends(RcppArmadillo)]]

class CIWTIWT {
private:
    int B;
    int n_times;
    arma::mat S_obs;
    arma::mat S_perm;
    arma::vec dt;
    arma::mat p_unadj;
    arma::rowvec p_adj;

public:
    CIWTIWT(const arma::mat& S_obs_in, const arma::mat& S_perm_in, int B_in, const arma::vec& dt_in);

    void compute_interval_pvalues();
    void compute_adjusted_pvalues();

    arma::rowvec get_results() const;
    arma::rowvec get_pointwise_pvalues() const;
};

#endif //CIWT_IWT_H
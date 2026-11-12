//
// Created by albor on 31/08/2026.
//

#include "CIWTIWT.h"

CIWTIWT::CIWTIWT(const arma::mat& S_obs_in, const arma::mat& S_perm_in, int B_in, const arma::vec& dt_in) {
    S_obs = S_obs_in;
    S_perm = S_perm_in;
    B = B_in;
    dt = dt_in;
    n_times = S_obs.n_cols;

    p_unadj.zeros(n_times, n_times);
    p_adj.zeros(1, n_times);
}

void CIWTIWT::compute_interval_pvalues() {
    for (int start = 0; start < n_times; ++start) {
        double T_obs = 0.0;
        arma::vec T_perm(B, arma::fill::zeros);

        for (int end = start; end < n_times; ++end) {
            T_obs += S_obs(0, end) * dt(end);
            T_perm += S_perm.col(end) * dt(end);

            arma::uvec exceeds = arma::find(T_perm >= T_obs);
            p_unadj(start, end) = static_cast<double>(exceeds.n_elem) / B;
        }
    }
}

void CIWTIWT::compute_adjusted_pvalues() {
    for (int t = 0; t < n_times; ++t) {
        double max_p = 0.0;

        for (int start = 0; start <= t; ++start) {
            for (int end = t; end < n_times; ++end) {
                if (p_unadj(start, end) > max_p) {
                    max_p = p_unadj(start, end);
                }
            }
        }
        p_adj(0, t) = max_p;
    }
}

arma::rowvec CIWTIWT::get_results() const {
    return p_adj;
}

arma::rowvec CIWTIWT::get_pointwise_pvalues() const {
    return p_unadj.diag().t(); 
}
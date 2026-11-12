//
// Created by albor on 31/08/2026.
//

#include "CIWTStatistic.h"

// 1) Constructor
CIWTStatistic::CIWTStatistic(const arma::cube& R_in, const arma::cube& E_in, int B_in) {
    R = R_in;
    E = E_in;
    B = B_in;

    n_subjects = R.n_rows;
    q_states = R.n_cols;
    n_times = R.n_slices;

    S_obs.set_size(1, n_times);
    S_perm.set_size(B, n_times);
    multipliers.set_size(B, n_subjects);
}

// 2) Generate Rademacher variables
void CIWTStatistic::generate_rademacher() {
    multipliers = arma::randi<arma::mat>(B, n_subjects, arma::distr_param(0, 1)) * 2 - 1;
}

// 3) Compute empirical covariance V(t) and its Moore-Penrose pseudoinverse
arma::cube CIWTStatistic::compute_V_inverse() {
    arma::cube V_inv_cube(q_states, q_states, n_times);

    for(int t = 0; t < n_times; ++t) {
        arma::mat R_t = R.slice(t);
        arma::mat E_t = E.slice(t);

        arma::mat M_t = R_t % E_t;
        arma::mat V_t = M_t.t() * M_t;

        V_inv_cube.slice(t) = arma::pinv(V_t);
    }
    return V_inv_cube;
}

// 4 & 5) Compute Score Statistics U(t) and S(t) simultaneously
void CIWTStatistic::compute_score_statistics() {
    arma::cube V_inv_cube = compute_V_inverse();

    for(int t = 0; t < n_times; ++t) {
        arma::mat R_t = R.slice(t);
        arma::mat E_t = E.slice(t);

        arma::mat M_t = R_t % E_t;
        arma::mat V_inv = V_inv_cube.slice(t);

        // --- OBSERVED STATISTIC ---
        arma::mat U_obs = arma::sum(M_t, 0).t();
        arma::mat S_obs_t = U_obs.t() * V_inv * U_obs;
        S_obs(0, t) = S_obs_t(0, 0);

        // --- PERMUTED STATISTICS ---
        arma::mat U_all_b = multipliers * M_t;
        arma::mat perm_scores = (U_all_b * V_inv) % U_all_b;
        S_perm.col(t) = arma::sum(perm_scores, 1);
    }
}

// Getters
arma::mat CIWTStatistic::get_S_obs() const { return S_obs; }
arma::mat CIWTStatistic::get_S_perm() const { return S_perm; }
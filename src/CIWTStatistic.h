//
// Created by albor on 31/08/2026.
//

#ifndef CIWT_STATISTIC_H
#define CIWT_STATISTIC_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class CIWTStatistic {
private:
    int B;               // Number of permutations
    int n_subjects;      // Number of panelists/sequences
    int q_states;        // Number of active categories
    int n_times;         // Length of the temporal grid

    arma::cube R;        // Dimensions: n_subjects x q_states x n_times
    arma::cube E;        // Dimensions: n_subjects x q_states x n_times

    arma::mat S_obs;     // 1 x n_times to store the observed statistics
    arma::mat S_perm;    // B x n_times to store permuted statistics
    arma::mat multipliers; // B x n_subjects Rademacher variables

public:
    // 1) Constructor
    CIWTStatistic(const arma::cube& R_in, const arma::cube& E_in, int B_in);

    // 2) Generate Rademacher variables matrix (B x n)
    void generate_rademacher();

    // 3) Compute V(t) and its Moore-Penrose pseudoinverse [V(t)]^+
    arma::cube compute_V_inverse();

    // 4 & 5) Compute U(t) and S(t) simultaneously
    void compute_score_statistics();

    // Getters
    arma::mat get_S_obs() const;
    arma::mat get_S_perm() const;
};

#endif //CIWT_STATISTIC_H
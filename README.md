# Categorical Interval-Wise Test (CIWT)

[![R-CMD-check](https://img.shields.io/badge/R-CMD-check-pass.svg)](https://github.com/dragonslayerere/CIWT)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**CIWT** (`dragonslayerere/CIWT`) is an open-source R package coupled with a high-performance **C++ backend** (powered by `RcppArmadillo`) designed for non-parametric statistical inference on **categorical functional data**. 

The framework evaluates the time-varying impact of target covariates on continuous-time categorical trajectories while strictly controlling for multidimensional nuisance covariates and handling complex cross-category dependencies[cite: 3].

---

## Key Methodological Pillars

1. **Regularized Baseline Null Estimation:** Integrates point-wise Firth's penalization via the Jeffreys invariant prior to guarantee stable, finite parameter estimates under complete or quasi-complete data separation[cite: 3].
2. **Nuisance Space Orthogonalization:** Implements a functional adaptation of the **Frisch-Waugh-Lovell (FWL) theorem** via Weighted Least Squares (WLS) projections to isolate the target covariate from confounding nuisance parameters[cite: 3].
3. **Multivariate Score Statistics & Pseudoinverse Mapping:** Constructs robust multivariate score statistics that automatically adapt to rank-deficient covariance structures (e.g., strict mutual exclusivity in **TDS** paradigms) or full-rank structures (e.g., unconstrained multi-label **TCATA** paradigms) using the Moore-Penrose pseudoinverse[cite: 3, 52, 57].
4. **Rademacher Multiplier Wild Bootstrap:** Employs a time-invariant subject-level multiplier bootstrap to compute valid empirical p-values and achieve **strong asymptotic control of the Family-Wise Error Rate (FWER)** across continuous temporal domains without iterative model refitting[cite: 3, 16, 18].
5. **Analytical Continuity:** Recovers classical multivariate benchmarks—such as Pearson's $\chi^2$ and Hotelling's $T^2$ tests—as exact algebraic reductions when nuisance covariates are omitted[cite: 3, 23].

---

## Installation

You can install the development version of `CIWT` directly from GitHub:
devtools::install_github("dragonslayerere/CIWT")

## Computational Architecture  
To eliminate the prohibitive computational bottlenecks of traditional resampling over dense functional domains, CIWT decouples null parameter estimation from the bootstrap loop:   Single-Pass Estimation: Firth-penalized GLM null fitting and WLS projections are executed strictly once across the evaluation grid.   Matrix Vectorization: Precomputed Moore-Penrose pseudoinverses and Level-3 BLAS operations evaluate all $B$ bootstrap score trajectories simultaneously via vectorized matrix multiplications ($\mathcal{U}(t) = \mathcal{SM}(t)$).   Parallel Decomposition: Built-in dual-level parallelism (grid-level and batch-level CPU multi-threading) delivers linear scaling on multi-core systems.

## Quick start and Basic Usage

library(CIWT)

### Example syntax structure (refer to package documentation for full arguments)
fit <- ciwt_test(
  formula_target  = ~ X_covariate,
  formula_nuisance = ~ Z1 + Z2,
  data            = my_categorical_fd_dataframe,
  id              = "subject_id",
  time            = "time_grid",
  paradigm        = "TDS", # or "TCATA"
  B               = 10000  # Number of Rademacher multiplier bootstrap replicates
)

### Summarize and plot interval-wise adjusted p-values
summary(fit)
plot(fit, alpha = 0.05)

## Supported Experimental Paradigms

- TDS (Temporal Dominance of Sensations): Mutually exclusive categorical states ($\sum_{j} Y_j(t) = 1$) with rank-deficient covariance tracking ($\chi_{q_{obs}-1}^2$).   
- TCATA (Temporal Check-All-That-Apply): Unconstrained multi-label activation ($\sum_j Y_j(t) \in \{0, \dots, q\}$) with full-rank convergence ($\chi_{q_{obs}}^2$).   
- Paired / Longitudinal Designs: Extension for multi-session or repeated-measures gustometer panels utilizing intra-subject difference score vectors. 

## License
Distributed under the MIT License. See LICENSE for more information.


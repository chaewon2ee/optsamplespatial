# Optimal Sampling in Spatial Data

This repository contains code for a response sampling algorithm that sequentially selects which spatial locations to measure the response variable at, under measurement cost or sample size constraints. We show theoretically that the final estimator obtained from a pilot-only-based sequential design is asymptotically as accurate, in terms of the target criterion, as an oracle-optimal design that assumes the covariance parameters are known, and we verify this through simulation.

## Table of Contents
- [Overview](#overview)
- [Method](#method)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Simulation Settings](#simulation-settings)
- [Output](#output)
- [Results Interpretation](#results-interpretation)
- [Acknowledgement](#acknowledgement)

## Overview

We consider a setting with a candidate-site index set $I_N = \{1, \dots, N\}$, where each site's coordinate $s_i$, covariate vector $x_i$, and measurement cost $c_i$ are known in advance, but the response $Y_i$ is observed only after the site is actually selected into the sample. Under a fixed-size or cost-constrained sampling budget, we sequentially construct a sample $S$ that minimizes a target-weighted GLS variance criterion $L_N(S;X_N,\theta)$ for the target functional $A_N\beta_0$.

The core idea is as follows:

1. First select a spatially well-spread pilot sample and use it to estimate the covariance parameters $\theta=(\phi,\tau^2)$.
2. At each step, greedily add the candidate site that produces the largest one-step decrease in the target criterion ($\Delta_i$, or $\Delta_i/c_i$ under a cost constraint) (Algorithm 1).
3. Compute the final GLS estimate $\hat\beta$ on the resulting sample.

We show theoretically that this sequential design, built using only the pilot estimate, achieves the same target covariance on the $m_N$-scale as the oracle-optimal design that assumes the true parameter $\theta_0$ is known (Theorem 1).

## Method

| Step | Description | Related Function/Concept |
|---|---|---|
| Pilot selection | Select a pilot sample $S_{pil}$ using a response-free spatial coverage rule such as maximin | `select_pilot_maximin` |
| Pilot estimation | Measure the response on $S_{pil}$ and estimate $\hat\theta_{pil}$ | GLS-based covariance parameter estimation |
| Sequential greedy | Starting from $S_{pil}$, repeatedly add the candidate in $A_N(S)$ that maximizes $\Psi_i(S;X_N,\hat\theta_{pil})$ | `Ψ_i`, `Δ_i`, `v_i`, `r_i` (Sherman–Morrison update) |
| Final estimation | Re-estimate $\hat\theta_{fin}$ and $\hat\beta$ on the final sample $\hat S^{seq}_N$ | GLS |

The theoretical properties (consistency of the covariance parameter estimates, the oracle gap, and the asymptotic equivalence of the final GLS covariance to the oracle) are established in Sections 3–4 of the paper (Algorithm 1, Propositions 1–2, Lemma 1, Theorem 1).

## Repository Structure
.
├── pilot_only_vs_proposed_incremental.R # Simulation: Pilot-only design vs Proposed (pilot + sequential selection) design
├── optsamplespatial_abridged.pdf # Methodology paper (abridged): notation, theory, Algorithm 1
└── README.md


## Requirements
- R >= 4.0
- Runs with base R only (no external package dependencies; uses base functions such as `uniroot` and standard matrix operations)

## Installation

```bash
git clone https://github.com/<username>/<repo-name>.git
cd <repo-name>
```

## Usage

```bash
Rscript pilot_only_vs_proposed_incremental.R
```

The script runs repeated Monte Carlo simulations comparing two designs on the same candidate population (coordinates, covariates):

- **A. Pilot-only**: Draw a pilot of roughly the final sample size (≈ `n_N`) from the start, and estimate $\hat\beta$, $\hat\tau^2$, $\hat\phi$ from that pilot alone, with no sequential selection.
- **B. Proposed**: Draw a small pilot (≈ `n_pil_B_target` + 5) to estimate $\hat\theta_{pil}$, then use Algorithm 1's sequential greedy procedure to expand the sample up to `n_N` and estimate $\hat\beta$.

| Variable | Value | Meaning |
|---|---|---|
| `n_iter` | 100 | Number of Monte Carlo replicates |
| `N0` | 10,000 | Number of candidate sites |
| `n_N` | 500 | Final sample size |
| `n_pil_A_target` | 495 | Pilot target size for the Pilot-only design (actual pilot ≈ 500) |
| `n_pil_B_target` | 100 | Pilot target size for the Proposed design (actual pilot ≈ 205) |

## Simulation Settings

| Parameter | Value |
|---|---|
| Number of candidate sites (N0) | 10,000 |
| Number of covariates (p) | 3 |
| Final sample size (n_N) | 500 |
| Pilot-only actual pilot size | ≈ 500 |
| Proposed pilot size | ≈ 205 |
| β (true) | (2, 0.5, 1) |
| Spatial correlation function | Matérn ($\nu = 0.5$) |
| Correlation at reference distance (ρ*) | 0.3 |
| Nugget ratio (τ²/σ²) | 0.01 |
| Number of Monte Carlo replicates | 100 |

The spatial covariance uses a Matérn correlation function, with the range parameter $\phi$ solved via `find_phi` so that the correlation equals 0.05 at a specified distance. A small number of closely spaced "clone" point pairs are also added to include cases with strong local spatial dependence.

## Output

| File | Description |
|---|---|
| `pilot_only_100_beta_estimates.csv` | Per-replicate β estimates from the Pilot-only design |
| `proposed_20_80_beta_estimates.csv` | Per-replicate β estimates from the Proposed design |
| `pilot100_vs_proposed_theta_estimates.csv` | Per-replicate τ², φ estimates from both designs |
| `pilot100_vs_proposed_beta_summary.csv` | Bias/Variance/MSE/RMSE comparison for each component of β |
| `pilot100_vs_proposed_tau2_summary.csv` | Bias/Variance/MSE/RMSE comparison for τ² |
| `pilot100_vs_proposed_phi_summary.csv` | Bias/Variance/MSE/RMSE comparison for φ |
| `pilot100_vs_proposed_complete_summary.csv` | Combined comparison table of the three summaries above |
| `proposed_beta_MSE_improvement.csv` | MSE reduction and relative MSE reduction of Proposed vs Pilot-only, by β component |
| `pilot_only_500_vs_proposed_200_300_results.rds` | Full result object containing settings and all estimates/summaries |

## Results Interpretation

- In the `*_summary.csv` files, smaller Bias, MSE, and RMSE indicate more accurate estimation.
- A positive `relative_MSE_reduction` in `proposed_beta_MSE_improvement.csv` means that, at the same final sample size (`n_N`), the Proposed design — a small pilot expanded via sequential selection — achieves lower MSE than the Pilot-only design, which uses one large pilot and stops there. This serves as finite-sample evidence supporting the paper's theoretical result that the sequential design asymptotically approaches the oracle-optimal design.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/

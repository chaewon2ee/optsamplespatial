# optsamplespatial : Criterion-Guided Sequential Response Sampling for Spatial Data under Measurement Constraints

공간 데이터에서 측정 비용/표본 크기 제약이 있을 때, 어느 지점에서 반응변수를 측정할지 순차적으로 골라주는 표본설계(response sampling) 알고리즘 코드입니다. Pilot 표본만으로 얻는 최종 추정량이, 공분산 파라미터를 안다고 가정한 oracle-optimal 설계와 asymptotic하게 동일한 target-criterion 정확도를 가짐을 이론적으로 보이고, 이를 시뮬레이션으로 검증합니다.

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
- [Research Status](#research-status)
- [Acknowledgement](#acknowledgement)

## Overview

후보 지점(candidate site) 집합 $I_N = \{1, \dots, N\}$ 이 주어지고, 각 지점의 좌표 $s_i$, 공변량 $x_i$, 측정 비용 $c_i$ 는 사전에 알려져 있지만 반응변수 $Y_i$ 는 실제로 표본에 포함되어 측정하기 전까지는 관측되지 않는 상황을 다룹니다. 표본 크기(fixed-size) 또는 총 비용(cost-constrained) 제약 하에서, 목표 함수 $A_N\beta_0$ 에 대한 GLS 추정의 target-weighted 분산 기준 $L_N(S;X_N,\theta)$ 을 최소화하는 표본 $S$ 를 순차적으로(sequentially) 구성합니다.

핵심 아이디어는 다음과 같습니다.

1. 공간적으로 고르게 퍼진 pilot 표본을 먼저 뽑아 공분산 파라미터 $\theta=(\phi,\tau^2)$ 를 추정합니다.
2. 매 단계마다, 현재 표본에 한 지점을 추가했을 때 target criterion을 가장 많이 줄이는(one-step decrease $\Delta_i$, 비용 제약 시 $\Delta_i/c_i$) 지점을 greedy하게 추가합니다 (Algorithm 1).
3. 최종 표본으로 GLS 추정 $\hat\beta$ 를 계산합니다.

이론적으로, pilot 표본만으로 구성된 이 순차적 설계는 참 파라미터 $\theta_0$ 를 안다고 가정한 oracle-optimal 설계와 $m_N$-scale에서 동일한 target covariance를 갖는다는 것(Theorem 1)을 보입니다.

## Method

| 단계 | 설명 | 관련 함수/개념 |
|---|---|---|
| Pilot 선택 | Maximin 등 반응변수와 무관한 공간 커버리지 규칙으로 pilot 표본 $S_{pil}$ 선택 | `select_pilot_maximin` |
| Pilot 추정 | $S_{pil}$ 의 반응변수를 측정하여 $\hat\theta_{pil}$ 추정 | GLS 기반 공분산 파라미터 추정 |
| Sequential greedy | $S_{pil}$ 에서 시작해 $A_N(S)$ 의 후보 중 $\Psi_i(S;X_N,\hat\theta_{pil})$ 이 최대인 지점을 하나씩 추가 | `Ψ_i`, `Δ_i`, `v_i`, `r_i` (Sherman–Morrison 갱신) |
| 최종 추정 | 최종 표본 $\hat S^{seq}_N$ 으로 $\hat\theta_{fin}$, $\hat\beta$ 재추정 | GLS |

이론적 성질(공분산 파라미터 정합성, oracle gap, 최종 GLS 공분산의 oracle 대비 asymptotic equivalence)은 논문 3~4절(Algorithm 1, Proposition 1–2, Lemma 1, Theorem 1)에 정리되어 있습니다.

## Repository Structure

```
.
├── pilot_only_vs_proposed_incremental.R   # 시뮬레이션: Pilot-only 설계 vs Proposed(pilot+순차선택) 설계 비교
├── optsamplespatial_축약본.pdf             # 방법론 논문(축약본): 표기, 이론, Algorithm 1
└── README.md
```

## Requirements
- R >= 4.0
- Base R만으로 실행 가능 (별도 패키지 의존성 없음; `uniroot`, 행렬 연산 등 base 함수 사용)

## Installation

```bash
git clone https://github.com/<username>/<repo-name>.git
cd <repo-name>
```

## Usage

```bash
Rscript pilot_only_vs_proposed_incremental.R
```

스크립트는 같은 후보 모집단(좌표, 공변량) 위에서 두 가지 설계를 반복 시뮬레이션(Monte Carlo)하여 비교합니다.

- **A. Pilot-only** : 처음부터 최종 표본 크기만큼(≈ `n_N`)을 pilot으로 뽑고, sequential selection 없이 그 표본만으로 $\hat\beta$, $\hat\tau^2$, $\hat\phi$ 를 추정
- **B. Proposed** : 작은 pilot(≈ `n_pil_B_target` + 5)만 뽑아 $\hat\theta_{pil}$ 을 추정한 뒤, Algorithm 1의 sequential greedy로 `n_N` 까지 확장하여 $\hat\beta$ 를 추정


| 변수 | 값 | 의미 |
|---|---|---|
| `n_iter` | 100 | Monte Carlo 반복 횟수 |
| `N0` | 10,000 | 후보 지점(candidate site) 개수 |
| `n_N` | 500 | 최종 표본 크기 |
| `n_pil_A_target` | 495 | Pilot-only 설계의 pilot 목표 크기 (실제 pilot ≈ 500) |
| `n_pil_B_target` | 100 | Proposed 설계의 pilot 목표 크기 (실제 pilot ≈ 205) |

## Simulation Settings

| Parameter | Value |
|---|---|
| 후보 지점 수 (N0) | 10,000 |
| 공변량 개수 (p) | 3 |
| 최종 표본 크기 (n_N) | 500 |
| Pilot-only 실제 pilot 크기 | ≈ 500 |
| Proposed pilot 크기 | ≈ 205 |
| β (true) | (2, 0.5, 1) |
| 공간 상관함수 | Matérn ($\nu = 0.5$) |
| 유효거리 기준 상관계수 (ρ*) | 0.3 |
| 너깃 비율 (τ²/σ²) | 0.01 |
| Monte Carlo 반복 수 | 100 |

공간 공분산은 Matérn 상관함수를 사용하며, `find_phi` 로 지정한 거리에서 상관계수가 0.05가 되도록 range parameter $\phi$ 를 역산합니다. 가까운 위치에 클론(clone) 지점 쌍을 소수 추가해 국소적으로 강한 공간 의존성을 갖는 상황도 포함합니다.

## Output

| File | Description |
|---|---|
| `pilot_only_100_beta_estimates.csv` | Pilot-only 설계의 반복별 β 추정치 |
| `proposed_20_80_beta_estimates.csv` | Proposed 설계의 반복별 β 추정치 |
| `pilot100_vs_proposed_theta_estimates.csv` | 두 설계의 반복별 τ², φ 추정치 |
| `pilot100_vs_proposed_beta_summary.csv` | β 각 성분에 대한 Bias/Variance/MSE/RMSE 비교 |
| `pilot100_vs_proposed_tau2_summary.csv` | τ² 에 대한 Bias/Variance/MSE/RMSE 비교 |
| `pilot100_vs_proposed_phi_summary.csv` | φ 에 대한 Bias/Variance/MSE/RMSE 비교 |
| `pilot100_vs_proposed_complete_summary.csv` | 위 세 요약을 합친 전체 비교표 |
| `proposed_beta_MSE_improvement.csv` | Pilot-only 대비 Proposed의 β별 MSE 감소량/상대 감소율 |
| `pilot_only_500_vs_proposed_200_300_results.rds` | 설정값과 모든 추정치·요약을 포함한 전체 결과 객체 |

## Results Interpretation

- `*_summary.csv` 에서 Bias, MSE, RMSE가 작을수록 더 정확한 추정입니다.
- `proposed_beta_MSE_improvement.csv` 의 `relative_MSE_reduction` 이 양수이면, 큰 pilot 하나로 끝내는 Pilot-only 설계보다 소규모 pilot + 순차 선택으로 확장하는 Proposed 설계가 같은 최종 표본 크기(`n_N`)에서 더 낮은 MSE를 달성했다는 뜻입니다. 이는 논문의 이론적 결과(순차 설계가 oracle-optimal 설계에 asymptotic하게 근접)를 유한표본에서 뒷받침하는 근거로 사용됩니다.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/

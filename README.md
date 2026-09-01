# optsamplespatial

# Criterion-Guided Sequential Response Sampling for Spatial Data under Measurement Constraints

공간 데이터에서 측정 비용/표본 크기 제약이 있을 때, 어느 지점에서 반응변수를 측정할지 순차적으로 골라주는 표본설계 알고리즘 코드입니다.

## 방법

1. Pilot 표본을 먼저 뽑아 공간 공분산 파라미터를 추정
2. 매 단계마다 target criterion(추정 정확도)을 가장 많이 개선하는 지점을 하나씩 추가 (Sequential greedy)
3. 최종 표본으로 GLS 추정

Pilot만으로 얻는 추정량이, 알려지지 않은 파라미터를 안다고 가정한 oracle-optimal 설계와 asymptotic하게 동일한 정확도를 가짐을 이론적으로 보였습니다.

## 시뮬레이션

`pilot_only_vs_proposed.R`: pilot만 쓰는 설계 vs pilot+순차선택 설계를 같은 후보 모집단에서 비교하여 β, τ², φ의 MSE를 확인합니다.

```bash
Rscript pilot_only_vs_proposed.R
```

---
논문 *"Criterion-Guided Sequential Response Sampling for Spatial Data under Measurement Constraints"* (Algorithm 1) 검증용 코드입니다.

## Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 - http://dsplus.uos.ac.kr/

# gFreqPCA

This repo contains the source code for PCA in the graph frequency domain.

For the detail, please refer to our paper: [https://arxiv.org/pdf/2410.08422](https://arxiv.org/pdf/2410.08422).

## Description

- Code
  - `utils.R` is a code for functions used for gFreqPCA.
  - `simulations.R` is a code for simulation study.
  - `seoulmetro.R` and `G20.R` are codes for real data analyses.


- Data
  - `seoulmetro` contains data of daily number of people getting on and off the Seoul Metropolitan Subway in South Korea.
  - `stationary` contains data of hourly temperature measurements recorded in Fahrenheit across the United States on August 1, 2010.
  - `BACI_HS92_V202501` (available at [http://www.cepii.fr/anglaisgraph/bdd/baci.htm](http://www.cepii.fr/anglaisgraph/bdd/baci.htm)) and `economic` contain information on world trade and economic indicators.

## Code overview
We present a PCA method in the graph frequency domain.

## References
Kim, K. and Oh, H.-S. (2024). Principal Component Analysis in the Graph Frequency Domain. arXiv preprint arXiv:2410.08422.
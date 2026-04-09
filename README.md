# Exact Computation for the Normalized Power Prior
This repository includes my master's thesis work (temporarily) entitled **"Exact Computation for the Normalized Power Prior"** at FGV EMAp (School of Applied Mathematics, Fundação Getulio Vargas, Brazil).

This work is led by [Eduardo Adame](https://eadame.ovh) and supervised by Professors [Luiz Max Carvalho (FGV EMAp)](https://github.com/maxbiostat) and [Flávio B. Gonçalves (UFMG)](https://est.ufmg.br/~fbgoncalves/).

Special thanks to our contributors Ezequiel Braga and Professor Joseph G. Ibrahim.

## Running the R Scripts

All scripts must be run from the `software/experiments/` directory, as they use relative paths to load each other and the Stan model (`../npp_mvn.stan`).

### Prerequisites

Install the required R packages:

```r
install.packages(c("MASS", "mvtnorm", "foreach", "doParallel",
                   "tidyverse", "posterior", "cmdstanr"))
```

`cmdstanr` also requires a working CmdStan installation:

```r
cmdstanr::install_cmdstan()
```

### Fixed-iteration experiment (`run_all.R`)

Runs all estimators for a fixed 5000 iterations × 4 chains and saves posterior plots and KS summary tables.

```r
setwd("software/experiments")
source("run_all.R")
```

Outputs (written to `software/experiments/`):
- `posterior_{small,medium,large}.png` — posterior density comparisons
- `summary_{small,medium,large}.csv` — per-estimator KS statistics and quantiles

To generate the LaTeX KS table after the above:

```r
source("make_ks_table.R")   # writes ks_table.tex
```

To run a quick sanity check before the full experiment, open `run_all.R` and set `RUN_TEST <- TRUE`.

### Convergence-based experiment (`run_convergence.R`)

Runs each approximate estimator in iterative batches until its KS distance to the Stan baseline drops below 0.05, reporting chains and computation time required.

```r
setwd("software/experiments")
source("run_convergence.R")
```

Outputs (written to `software/experiments/results/{small,medium,large}/`):
- `{estimator}/checkpoint_NN/posterior.png` — density plot at each 5000-iteration checkpoint, with chains/time/KS in the legend
- `{estimator}/checkpoint_NN/convergence.csv` — convergence log up to that batch
- `{estimator}/checkpoint_NN/traces/chain_NN_eta.png` — per-chain trace plots
- `comparison_summary.csv` — convergence summary across all estimators
- `final_comparison.png` — final posterior density comparison with all methods

---

## Overview

The Bayesian framework offers a natural solution to the issue of incorporating expert knowledge into statistical analyses. A major tool for the construction of informative priors is the normalized power prior (Ibrahim, 2015). This work addresses two critical aspects:

1. We examine the power prior approach (Ibrahim, 2015), which controls historical data influence through a parameter $\eta \in [0,1]$, where $\eta=0$ discards and $\eta=1$ fully incorporates historical information.
    
2. Current methods rely on approximate inference (Carvalho, 2021) for these doubly intractable posteriors, lacking theoretical guarantees.

We aim to develop an exact MCMC algorithm (Bhattacharya, 2021; Vats, 2021; Gonçalves, 2017; Andrieu, 2009; Murray, 2012) that, thus, provides **guaranteed convergence**, enables efficient sampling from complex power prior posteriors, and offers practical advantages for clinical trial design.

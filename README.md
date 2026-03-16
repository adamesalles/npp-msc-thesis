# Exact MCMC Methods for the Normalized Power Prior
This repository includes my master's thesis work (temporarily) entitled **"Exact MCMC Methods for the Normalized Power Prior"** at FGV EMAp (School of Applied Mathematics, Fundação Getulio Vargas, Brazil).

This work is led by [Eduardo Adame](https://eadame.ovh) and supervised by Professors [Luiz Max Carvalho (FGV EMAp)](https://github.com/maxbiostat) and [Flávio B. Gonçalves (UFMG)](https://est.ufmg.br/~fbgoncalves/).

Special thanks to our contributors Ezequiel Braga and Professor Joseph G. Ibrahim.

## Overview

The Bayesian framework offers a natural solution to the issue of incorporating expert knowledge into statistical analyses. A major tool for the construction of informative priors is the normalized power prior (Ibrahim, 2015). This work addresses two critical aspects:

1. We examine the power prior approach (Ibrahim, 2015), which controls historical data influence through a parameter $\eta \in [0,1]$, where $\eta=0$ discards and $\eta=1$ fully incorporates historical information.
    
2. Current methods rely on approximate inference (Carvalho, 2021) for these doubly intractable posteriors, lacking theoretical guarantees.

We aim to develop an exact MCMC algorithm (Bhattacharya, 2021; Vats, 2021; Gonçalves, 2017; Andrieu, 2009; Murray, 2012) that, thus, provides **guaranteed convergence**, enables efficient sampling from complex power prior posteriors, and offers practical advantages for clinical trial design.

# sDCM-CaI: Spectral Dynamic Causal Modeling for Calcium Imaging

Code accompanying:

> **Spectral dynamic causal modeling for effective connectivity from spontaneous population calcium imaging**
> Hae-Jeong Park, Dongmyeong Lee, Euisun Kim, Jinseok Eo, Karl Friston

sDCM-CaI extends spectral Dynamic Causal Modeling to resting-state calcium imaging, combining a reduced neuronal state-space model, calcium-state dynamics, a fluorescence observation model, and Bayesian inversion of cross-spectral density features to estimate directed effective connectivity from population calcium signals.

## Requirements

- MATLAB with [SPM25](https://www.fil.ion.ucl.ac.uk/spm/software/spm25/) (or SPM12) on the path
- Python 3 with `numpy`, `scipy`, `matplotlib` for figure generation

## Repository structure

```
core/            Generative model and inversion machinery
  spm_fx_calcium.m                State equation (neuronal + calcium dynamics)
  spm_gx_calcium*.m                Observation models (Hill, linear variants)
  spm_dcm_calcium_priors.m         Prior specification
  spm_dcm_calcium_csd.m            Main DCM inversion wrapper
  spm_csd_calcium_mtf.m,
  spm_dcm_mtf_calcium.m            Model transfer function (spectral prediction)
  spm_dcm_csd_Q.m                  CSD precision/noise model
  spm_dcm_calcium_*.m              Result utilities
  misspecified_model_sim.m         Observation-model misspecification simulation
  run_hidden_driver_sim.m          Hidden common-driver simulation

simulations/
  run_mc_simulation_4node.m        4-node Monte Carlo recovery simulation (main-text Fig. 2)
  demo_simulation_and_inversion.m  Self-contained demo: simulates 4-node calcium data and
                                    inverts it under both the Hill and linear observation
                                    models (set DCM.options.custom_g to switch), reporting
                                    recovered connectivity and free energy for each

empirical_12n/   12-node empirical zebrafish analysis pipeline
  extract_*_signals_12n_from_raw.m, teo_kmeans_split.m
                                    Node signal extraction and TeO superficial/deep clustering
  run_sc_*_12n.m, run_sc_*_reglin_mean12n.m, run_sc_*_reghill*.m
                                    Subject-level DCM inversion, one script per compared
                                    model variant (RegLin/RegHill x Mean/PC1 x
                                    uniform/binary-SC/graded-SC prior architecture).
                                    run_sc_flat_reglin_mean12n.m is the primary,
                                    BMS-winning specification.
  run_peb_*.m, run_bms_*.m          Group-level PEB and Bayesian model selection/reduction
  run_splithalf_*.m                 Temporal and cross-subset split-half reliability analyses
  run_sensitivity_s16_*.m           Prior-covariance and prior-mean sensitivity analysis (Fig. 6)
  run_mc_reglin_mean_12n.m,
  run_mc_simulation_12node.m        12-node-scale Monte Carlo recovery simulation

figures/         Figure-generation scripts (Python)
```

## Data

The empirical analysis uses whole-brain larval zebrafish calcium-imaging recordings from Chen et al. (2018), publicly available at [github.com/xiuyechen/fishexplorer](https://github.com/xiuyechen/fishexplorer). Raw data are not redistributed here; scripts expect the extracted per-subject node time series as produced by `empirical_12n/extract_mean_signals_12n_from_raw.m`.

## Primary empirical pipeline

The primary (BMS-selected) specification used for all reported group-level results is:

1. `empirical_12n/run_sc_flat_reglin_mean12n.m` — per-subject inversion (RegLin observation model, mean population signal, uniform connectivity prior), run once per subject
2. `empirical_12n/run_peb_sc_flat_reglin_mean12n.m` — group-level PEB and BMR/BMA model reduction

## License

TBD.

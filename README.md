# Log-Transformed RTH for Single-Diode PV Parameter Extraction

MATLAB code, numerical outputs, and result figures supporting the article:

> Salma El Aamery and Rachid El Bachtiri, *Log-Transformed Red-Tailed Hawk Optimization for Single-Diode Photovoltaic Parameter Extraction*, **Next Energy**, 2026. https://doi.org/10.1016/j.nxener.2026.100879

## What this repository contains

```text
code/                 MATLAB scripts and optimisation implementations
data/                 Benchmark-data provenance and data-use notes
results/              Run-level outputs, tables, and saved MATLAB workspaces
figures/eps/          Publication-quality EPS result figures
figures/pdf-preview/  PDF previews of the EPS figures
supplementary/        Supplementary-material source and PDF
```

## Main scripts

Run these scripts from the repository root in MATLAB:

| Script | Purpose |
| --- | --- |
| `code/Main_PV_RTH.m` | Shell SM55 RTH extraction experiment |
| `code/Main_RTC_France.m` | RTC France RTH extraction experiment |
| `code/Main_Comparison_PV.m` | Initial comparison with benchmark algorithms |
| `code/Benchmark_EqualFE.m` | Matched function-evaluation experiment for Shell SM55, RTC France, and Photowatt-PWP201 |
| `code/Postprocess_EqualFE.m` | Exports equal-FE summaries and figures |
| `code/Ablation_LogI0_RTH.m` | Linear versus `log10(I0)` search-variable ablation |
| `code/Generate_ParamScatter.m` | Run-level parameter-distribution figures |

## Method summary

The decision vector is `[Iph, z, Rs, Rsh, n]`, where `z = log10(I0)` and the physical saturation current is recovered as `I0 = 10^z`. The single-diode current is evaluated implicitly at each measured voltage with Newton-Raphson iteration. RTH minimises the current RMSE and is compared with PSO, GWO, WOA, HHO, AVOA, and MGO.

The matched benchmark uses 30 independent runs and a common budget of 149,050 objective-function evaluations per run. The raw run-level RMSE values are provided in `results/EqualFE_run_RMSE.csv`.

## Data

The benchmark I-V values used by the MATLAB experiments are embedded in the experiment scripts and reported in the accompanying supplementary material. Please cite the original benchmark sources given in the article when reusing them.

## Reproducing the numerical outputs

1. Open MATLAB and set the current folder to this repository.
2. Run `code/Benchmark_EqualFE.m` to reproduce the main matched-budget experiment.
3. Run `code/Postprocess_EqualFE.m` to export summary tables and figures.
4. Run `code/Ablation_LogI0_RTH.m` for the `log10(I0)` ablation.

Stochastic trajectories can vary when executed with a different MATLAB version or random-seed sequence. The archived numerical outputs in `results/` are the results reported in the article.

## Citation

Please cite both the article above and this repository. Citation metadata are available in `CITATION.cff`.

## License

This repository is released under the [MIT License](LICENSE).

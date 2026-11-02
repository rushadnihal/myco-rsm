# myco-rsm

**Response Surface Methodology (RSM) and multi-objective media optimization for submerged mycelial fermentation.**

An R Shiny decision-support app for upstream bioprocess scientists who adapt filamentous-fungus strains to variable agro-industrial sidestreams (dairy whey, brewer's spent grain, and similar). Instead of one-variable-at-a-time (OFAT) testing, it uses sequential Central Composite (CCD) and Box-Behnken (BBD) designs to capture non-linear interactions, such as carbon-to-nitrogen (C:N) overflow or shear-induced pellet lysis, with far fewer bioreactor runs.

> **Project goal:** reduce wet-lab experimentation by an estimated 50-70% relative to OFAT. This is the design target, not a result measured by this repository.

## Features

| Area | What it does |
|---|---|
| **Design generator** | Randomized, blocked, rotatable CCD and Box-Behnken run-sheets in natural units, with CSV download |
| **Blocking** | Batch-to-batch drift is absorbed by a fixed `Block` effect, so it does not bias the surface coefficients |
| **Model fit** | Second-order polynomial (first-order, two-way interaction and pure quadratic terms) via `rsm`, with a formal lack-of-fit F-test against pure error from replicated centers |
| **Canonical analysis** | Stationary point, eigenvalues of the quadratic matrix, surface classification (maximum, minimum, saddle, ridge), and a check that the optimum lies inside the tested region |
| **Steepest ascent** | Path of steepest ascent from a first-order fit when the optimum lies outside the tested region |
| **Multi-objective optimizer** | Derringer-Suich desirability (maximize, minimize, target) with weights and optional media cost; returns a setpoint with 95% prediction intervals |
| **Surfaces** | Interactive contour and 3D response-surface plots (plotly) over any two factors |
| **CSV upload** | Validates a user's design (levels, replicated centers, model size) before fitting |

## Quick start

```r
install.packages(c("shiny", "bslib", "rsm", "plotly", "DT"))
shiny::runApp()          # from the repository root
```

The app opens with three datasets to choose from:

1. **ChemReact benchmark**: the 14-run, 2-factor blocked design from `rsm::ChemReact`.
2. **Simulated 4-factor demo**: a 30-run rotatable CCD with *synthetic* biomass and protein responses.
3. **Upload CSV**: your own data.

### Upload format

- One row per run; numeric columns for each factor (natural units) and each response.
- Each factor needs **3** (BBD) or **5** (CCD) distinct, symmetric levels. The center and step size are inferred from them.
- At least **2 replicated center points** (ideally 2 or more per block), so lack of fit can be tested.
- An optional block column (any labels). It is treated as a categorical factor.
- More runs than model parameters.

`data/example_ccd4_runsheet.csv` shows a generated run-sheet (design only; add your response columns after the runs are done).

## Methods

Predictors are analysed in coded units, `x_i = (X_i - X_i0) / ΔX_i`, and mapped back to natural units for reporting.

Fitted model (with blocks):

```
Y = β0 + Σ βi·xi + Σ βii·xi² + Σ Σ βij·xi·xj + Block + ε
```

- **Stationary point:** `x_s = -½ B⁻¹ b`, where `b` is the vector of first-order coefficients and `B` the symmetric matrix of quadratic coefficients.
- **Classification:** all eigenvalues of `B` negative means a maximum, all positive means a minimum, mixed signs means a saddle, and a near-zero eigenvalue means a ridge.
- **Lack of fit:** F-test of the second-order model against a one-mean-per-design-point model. Pure error comes from replicated centers.
- **Desirability:** `D = (d1^w1 · d2^w2 ··· dm^wm)^(1/Σw)`. A single response with `d = 0` forces `D = 0`.
- **Optimizer:** multi-start Nelder-Mead, penalized to stay inside the design region (radius of the farthest design run, in coded units).
- **Prediction intervals:** 95% intervals for one new run in the reference (first) block.

## Verification

The numerical results below were reproduced in R on the built-in data and are checked by the app itself (Optimum tab, ChemReact dataset).

**ChemReact benchmark**

| Check | Result |
|---|---|
| Block 1 alone, FO + TWI: lack of fit | F(1, 2) = 190.02, p = 0.0052, so curvature is significant and Block 2 is needed |
| Full blocked model: block shift (B2 − B1) | −4.46 |
| Full blocked model: lack of fit | F(3, 4) = 0.53, p = 0.685 (second-order model adequate) |
| Stationary point (coded) | (0.372, 0.334) |
| Stationary point (natural) | 86.86 min, 176.67 °C |
| Eigenvalues | −0.92, −1.32 (true maximum) |
| Radius from center | 0.50 (coded), inside the region (limit 1.414) |
| Predicted yield at optimum | 84.37 % (Block 1 basis) |

**Design properties** (4 factors)

- Rotatable CCD: 30 runs in 2 blocks (19 + 11), axial distance α = 2. Scaled prediction variance is identical in every direction at a given radius.
- Box-Behnken: 33 runs in 3 blocks. It is also rotatable for 4 factors but has higher prediction variance at equal radius, and its region only reaches radius √2.

**Recovery test** (simulated data from a known surface, 200 repeated experiments on the 30-run CCD)

- Mean estimated optimum (0.400, 0.301, −0.200, 0.103) against a true (0.4, 0.3, −0.2, 0.1); run-to-run SD 0.016-0.025 coded.
- Lack-of-fit p > 0.05 in 95.5% of runs (expected ≈ 95%).

## Repository layout

```
myco-rsm/
├── app.R                              # single-file Shiny app (functions + UI + server)
├── analysis/
│   └── benchmark_walkthrough.R        # step-by-step analysis that produced the numbers above
├── data/
│   └── example_ccd4_runsheet.csv      # generated 30-run CCD run-sheet (no responses)
├── README.md
└── .gitignore
```

`analysis/benchmark_walkthrough.R` is meant to be run top to bottom from the repository root; each section builds on the previous one.

## Deploying to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
rsconnect::deployApp(appDir = ".", appFiles = "app.R")
```

## Limitations

- The 4-factor demo dataset is **simulated**. It exercises the pipeline and is not experimental data. The published 4-factor fermentation dataset (Mendeley Data, DOI 10.17632/bdhc5xcb6h.3) is not bundled here.
- Models are second-order only, with a single fixed block effect. Responses with stronger non-linearity need a larger design or a different model.
- The steepest-ascent path comes from a first-order fit and extrapolates linearly. It overshoots when the surface curves, so use only the first steps before running new experiments.
- Desirability limits and weights are user-chosen and drive the optimum; D values from different weightings are not comparable with each other.
- The optional media cost assumes the first two factors are concentrations in g/L and that prices are linear.
- The upload validator expects symmetric 3- or 5-level designs.

## References

- Lenth, R. V. (2009). Response-surface methods in R, using `rsm`. *Journal of Statistical Software*, 32(7), 1-17.
- Derringer, G., & Suich, R. (1980). Simultaneous optimization of several response variables. *Journal of Quality Technology*, 12(4), 214-219.
- Myers, R. H., Montgomery, D. C., & Anderson-Cook, C. M. (2016). *Response Surface Methodology* (4th ed.). Wiley.

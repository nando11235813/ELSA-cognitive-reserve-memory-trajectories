# Bayesian smooth-transition differential equation mixed-effects model

This repository contains a **self-contained R example using simulated data only**. No ELSA data or other restricted records are included.

## Files

- `DEM_example_rjags_snow.R`: simulates longitudinal data for 50 subjects, writes the JAGS model to `dem.bug`, and fits it using four independent parallel MCMC chains.
- `dem.bug`: generated when the script is run; it need not be committed separately.

## Dependencies

R packages `rjags`, `snow`, and `coda`, as well as an installed JAGS system library.

## Expected input: `data`

The script creates a long-format `data.frame` named `data` with **one row per assessment**, having columns:

| Column | Meaning |
| --- | --- |
| `id` | Integer subject identifier, consecutive from 1 to the number of subjects. Repeated across visits. |
| `age` | Age or time at assessment, numeric. The example uses ages 50–90. |
| `y` | Longitudinal continuous outcome, numeric. |
| `x1` | Numeric subject-level predictor, constant for all visits of the subject. |
| `x2` | Binary subject-level predictor (0/1), constant for all visits of the subject. |

The simulation uses 50 subjects, each with **4–8 measurements**. Both `x1` and `x2` independently enter the subject-level distributions of *Level*, *Onset*, and *Rate*, with six independent regression coefficients. The observation model assumes conditional Gaussian residuals.

To fit other data, **replace only Section 1 (simulation)** with an import or construction of a long-format `data.frame` with these five columns. Set `n <- length(unique(data$id))` and ensure `id` is coded as 1, ..., `n` (or recode it). Subjects need not have the same number or ages of measurements. Covariates must be time-invariant; categorical variables with more than two categories require dummy coding and additional coefficient terms in the JAGS model.

## Model and assumptions

- Memory/outcome trajectories are flat until the individual onset parameter, undergo a smooth transition of duration `dt`, then follow exponential change with individual relative rate `b2_i`.
- `b0_i`, `CP_i`, and `b2_i` contain independent normal random effects. Separate fixed-effect coefficients are modeled for the numerical (`x1`) and binary (`x2`) predictors.
- The priors and scale parameters are the **illustrative values from the supplied analysis code**, not universally non-informative priors. The population Level prior assumes a 0–20 scale, and the onset prior assumes ages 60–80; users of other outcomes should adapt them.
- The Gaussian likelihood does not constrain outcomes to a fixed range. The pre-onset level is not truncated to be positive at the subject level in the original specification.
- Numerical stability may deteriorate if transition duration `dt` approaches zero; inspect sampling diagnostics.

## Execution

Run `Rscript DEM_example_rjags_snow.R` in a directory where you can write `dem.bug`. The program runs four independent MCMC chains in parallel using `snow` and prints summaries, Gelman–Rubin diagnostics, and effective sample sizes. MCMC execution requires a working local R/JAGS installation and may be time-consuming.

**Note:** The code has not been executed here because R/JAGS was not available in the authoring environment. Verify it locally before committing it to a public repository.

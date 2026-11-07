# Experimental two-way CMG backend for fixest

**Status: experimental proposal.** This is not upstream fixest. Native MAP remains the default unless a backend is selected per call or through session defaults set with `setFixest_estimation()`, and fixest uses cmgr only when CMG is selected.

This source adds an optional backend for demeaning with two fixed effects. It builds on the `demeaner` constructor and dispatcher interface proposed in [fixest PR #647](https://github.com/lrberge/fixest/pull/647). The numerical work is done by [cmgr](https://github.com/johannes-schmieder/CMG/tree/main/cmgr), a separate R package kept in the CMG repository.

This file appears on two branches of [johannes-schmieder/fixest-cmg](https://github.com/johannes-schmieder/fixest-cmg). Both report fixest 0.14.2.9000 and suggest cmgr 0.0.3 or later.

- `codex/cmg-two-way` is the full working branch, and the release's fixest source package comes from it.
- `codex/cmg-pr647` is the minimal review patch on top of #647.

The branches place `demeaner` at different positions in the argument list, so always pass it by name.

## Install

Follow the recipe in the [`cmgr-v0.0.3-preview.1` release](https://github.com/johannes-schmieder/CMG/releases/tag/cmgr-v0.0.3-preview.1). It installs cmgr and this fixest source into a separate library, leaving your regular fixest alone. It needs a Rust toolchain.

## Usage

```r
.libPaths(c(path.expand("~/R/cmg-preview"), .libPaths()))
library(fixest)

# Small low-mobility panel. To use your own data, replace d and the names
# y, x, worker and firm.
set.seed(1)
n_workers <- 2000L; n_years <- 6L; n_firms <- 100L
worker <- rep(seq_len(n_workers), each = n_years)
home   <- sample(rep_len(seq_len(n_firms), n_workers))
firm   <- home[worker]
moved  <- runif(length(worker)) < 0.05
firm[moved] <- sample.int(n_firms, sum(moved), replace = TRUE)
x <- rnorm(length(worker))
y <- 0.5 * x + rnorm(n_workers)[worker] + rnorm(n_firms)[firm] + rnorm(length(worker))
d <- data.frame(y, x, worker, firm)

fit <- feols(y ~ x | worker + firm, d, demeaner = CmgDemeaner(), nthreads = 2)
```

fixest still handles the formula, sample selection, singleton removal, collinearity, estimation, covariance and prediction. The backend only replaces the step that residualizes columns on the two fixed effects.

**Select the backend per call in this preview**, as above.

## Exact worker elimination

A worker with few employers doesn't need to be part of the graph that CMG solves. Given the firm effects, that worker's effect is just the worker's weighted mean of the outcome net of those firm effects. cmgr substitutes this into the system, which removes the worker exactly and connects the worker's firms directly. After solving the smaller system, it recovers the worker's effect.

This happens by default, and only for levels of the larger fixed-effect block (normally workers) with **at most five distinct neighbors**. Higher-degree workers stay in the solved graph. Nothing is sampled or approximated, no observations are dropped, and CMG's own numerics are unchanged. In linear algebra this is a partial exact Schur reduction, done as preprocessing inside cmgr.

`schur = FALSE` is an optional comparison and troubleshooting switch; the default is `TRUE`. It disables the algebraic reduction, so CMG solves the full worker-firm graph. The model and the observations stay the same, and CMG stays in use.

```r
full_graph <- feols(y ~ x | worker + firm, d, nthreads = 2,
                    demeaner = CmgDemeaner(schur = FALSE))
```

## Constructors

Each constructor returns a configuration object for the `demeaner` argument of `feols()`, `feols.fit()` and `demean()`.

```r
MapDemeaner(fixef_tol = 1e-06, fixef_maxiter = 10000, extraProj = 0, iter_warmup = 15,
            iter_projAfterAcc = 40, iter_grandAcc = 4)

LsmrDemeaner(fixef_atol = 1e-06, fixef_btol = 1e-06, fixef_maxiter = 10000,
             preconditioner = "auto", local_size = NULL)

CmgDemeaner(fixef_tol = 1e-8, fixef_maxiter = 10000L, schur = TRUE,
            tol = NULL, maxiter = NULL)
```

| Constructor | Solver | Needs | Fixed effects |
|---|---|---|---|
| `MapDemeaner()` | fixest's native alternating projections | nothing extra | as native fixest |
| `LsmrDemeaner()` | LSMR from [withinr](https://github.com/py-econometrics/within) | withinr with `LsmrOptions()` | two or more; one-way uses native demeaning; no varying slopes |
| `CmgDemeaner()` | CMG-PCG from cmgr | cmgr | exactly two categorical intercepts |

- **Controls go in the constructor.** Passing `fixef.tol`, `fixef.iter` or `fixef.algo` together with a non-NULL `demeaner` is an error.
- **Missing packages.** `LsmrDemeaner()` and `CmgDemeaner()` report a missing package as soon as they are called.
- **CMG aliases.** `tol` and `maxiter` are compatibility aliases for `fixef_tol` and `fixef_maxiter`. Supplying both names for one control is an error.
- **withinr tolerance.** withinr accepts one tolerance, so fixest passes the larger of `fixef_atol` and `fixef_btol`. `preconditioner = "auto"` selects withinr's additive Schwarz preconditioner, which applies its own approximate Schur step.
- **Threads.** `nthreads` controls native fixest and cmgr. withinr uses a process-wide Rayon pool with no per-call setting, so set `RAYON_NUM_THREADS` before starting R.
- **Tolerances aren't interchangeable.** MAP, LSMR and CMG stop on different criteria, so equal numbers don't mean equal accuracy.
- **Nonconvergence.** CMG and withinr raise an error, including in multiple-model fits, and never silently switch solvers. MAP keeps its existing warning and `convStatus` behavior.
- **Interrupts.** Cancelling a CMG solve joins all native work and then raises an ordinary R error. `try()` and error handlers can catch it, so a loop that wraps fits in `try()` may continue after an interrupt.

## Global defaults, prepared environments and saved fits

**Global selection.** `setFixest_estimation(demeaner = ...)` accepts these constructors, but it isn't recommended for CMG in this preview. A global CMG selector applies to every `feols()` call, so ordinary calls with no fixed effects, one, or more than two fail by design, with no silent fallback. It can't be combined with `fixef.tol` or `fixef.iter` in the same call, and `setFixest_estimation(reset = TRUE)` clears it. It doesn't change the inner loop of `feglm()` or `fepois()`.

**Prepared environments.** Set the backend and its controls when you create the environment, and don't supply the backend again when you reuse it:

```r
prepared <- feols(y ~ x | worker + firm, d, demeaner = CmgDemeaner(), only.env = TRUE)
fit <- feols(env = prepared)
```

Supplying native controls such as `fixef.tol` when reusing an optional-backend environment is an error. When reusing a single-model native environment, explicit native tolerance and iteration overrides are honored. For a `split` environment, set these controls when you create the environment, because overrides supplied on reuse don't reach the split fits.

**Saved fits.** `demean(fit)` follows these rules:

| Call | Result |
|---|---|
| `demean(fit)` | Returns the saved columns if the fit used `demeaned = TRUE`; otherwise recomputes with the fit's own backend |
| `demean(fit, demeaner = ctl)` | Recomputes with `ctl` |
| `demean(fit, demeaner = NULL)` | Recomputes with native MAP |
| `demean(fit, tol = , iter = , fixef.algo = )` | Recomputes with native MAP; on a fit from an optional backend, controls you omit take native defaults |
| `demean(fit, demeaner = ctl, tol = )` | Error |

An optional backend's tolerance is never copied into MAP. A saved CMG configuration without a `schur` field is treated as `schur = TRUE`; an explicit `FALSE` stays `FALSE`.

**Other `demean()` input.** Bare vectors and `data.frame` input with `as.matrix = TRUE` are normalized for optional backends. Native `demean()` without an explicit backend keeps its existing positive-tolerance and trusted-input validation.

## What stays native

**On a CMG fit, some results use native demeaning internally:**

- IV diagnostics
- named fixed-effect extraction

Their residualization doesn't necessarily reach the same accuracy as the CMG projection, and no numerical equivalence is promised. An internal already-demeaned, no-FE IV regression also stays native, so it doesn't inherit a global optional backend.

**Compared with pristine upstream fixest**, native MAP remains the default, but native results aren't all guaranteed identical. The shared interface patch forwards intended controls, so explicit native tolerances now reach split and multiple-model calls, which can change those results. Fits also gain saved backend metadata.

## Supported and unsupported

**Supported with `CmgDemeaner()`:** two categorical intercept effects in weighted or unweighted OLS, including FE-only fits, IV, multiple outcomes, stepwise regressors, split samples, offsets, fixest's covariance choices, `feols.fit()` and `demean()`.

**Rejected with an explicit error:**

- no fixed effects, one effect, or three or more effects
- varying slopes
- stepwise fixed-effect specifications
- an optional backend passed to GLM/PPML functions

**Not available in this preview:** a memory cap, a cache across calls, and IRLS support.

**Platform:** tested only on Apple Silicon (M2 Ultra) with Homebrew R 4.6.1, Apple clang 21 and Rust 1.85.1. CRAN builds of R, Linux and Windows aren't qualified. See the [cmgr README](https://github.com/johannes-schmieder/CMG/tree/main/cmgr) for limits.

## How to assess it

Compare methods on your own data, with the same sample and thread count:

```r
on  <- feols(y ~ x | worker + firm, d, demeaner = CmgDemeaner(fixef_tol = 1e-10))
off <- feols(y ~ x | worker + firm, d, demeaner = CmgDemeaner(fixef_tol = 1e-10, schur = FALSE))
map <- feols(y ~ x | worker + firm, d, fixef.tol = 1e-10, fixef.iter = 20000)

max(abs(coef(on) - coef(off)))
max(abs(resid(on) - resid(off)))
max(abs(coef(on) - coef(map)))
max(abs(resid(on) - resid(map)))
```

On a low-mobility graph, native residuals can differ noticeably from CMG's even when coefficients agree closely, so look at both. Use `system.time()` around complete calls to include setup. Speed depends on your graph: the published benchmark covers only six short synthetic panels, and no result here supports a general claim for long panels or many high-degree workers.

From a clone of this repository, with the packages installed, run `tests/demeaner_tests.R`, `tests/demeaner_boundary_tests.R`, `tests/cmg_tests.R` and `tests/cmg_sample_tests.R`. They use independent weighted-QR references and cover:

- the shared interface and full CMG integration
- columns collinear with worker, firm and additive effects, at scales 1 and 1e5, in both CMG modes
- actual elimination after NA, zero-weight and singleton cleanup
- outcomes with different missingness, including fixed-effect levels that disappear
- one and two threads

withinr checks run if it is installed.

## Links

- [Integration proposal](https://github.com/johannes-schmieder/CMG/blob/cmgr-v0.0.3-preview.1/cmgr/PROPOSAL.md)
- [Benchmark results](https://github.com/johannes-schmieder/CMG/blob/cmgr-v0.0.3-preview.1/cmgr/benchmarks/results/2026-11-07/README.md)
- [Minimal patch against #647](https://github.com/johannes-schmieder/fixest-cmg/compare/0752916d36e7376439817712ddf835b662946171...codex/cmg-pr647)
- [CMG-only part of that patch](https://github.com/johannes-schmieder/fixest-cmg/compare/codex/pr647-interface-fixes...codex/cmg-pr647)

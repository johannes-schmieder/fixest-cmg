# Experimental two-way CMG backend

This branch adds an optional `cmgr` backend at the existing fixed-effect
residualization boundary. It is based on upstream fixest `cf946133` (0.14.2).
Its constructors and explicit dispatcher follow fixest PR #647 (head 0752916).
The proposal remains a draft; this branch adds CMG as a third backend case.
It also forwards resolved MAP controls from the prepared environment to the
native solver, so constructor settings apply to actual OLS solves.

```r
library(fixest)
fit <- feols(y ~ x1 + x2 | worker + firm, data,
             demeaner = CmgDemeaner(fixef_tol = 1e-8), nthreads = 4)
```


## Optional exact worker elimination

Elimination is opt-in, so the existing full-graph path remains available:

```r
cmg_only <- feols(y ~ x1 + x2 | worker + firm, data,
                  demeaner = CmgDemeaner(schur = FALSE), nthreads = 4)
cmg_schur <- feols(y ~ x1 + x2 | worker + firm, data,
                   demeaner = CmgDemeaner(schur = TRUE), nthreads = 4)
```

`schur = TRUE` selects the FE block with more levels (normally workers; ties
select the first block in the prepared ID order). It eliminates levels with
at most five distinct neighbors after combining repeated worker--firm
observations. Higher-degree workers remain as vertices of the system actually
solved. Stayers create no firm edges, but their observations and within-worker
variation remain in the regression. No sampling is used, and CMG's source,
options and PCG implementation are unchanged. This option is available through
both `CmgDemeaner()` and `cmgr::demean_fixest()`.

Use the same tolerance, threads, sample and formula for paired on/off timings.
Both modes check the original weighted FE moments. The low-level function's
`reduction` diagnostics report the selected block and eliminated/retained level
counts; named `timing` entries report the actual graph size and native setup and
projection times. Include setup in complete-call comparisons. Turning
elimination off isolates its incremental contribution conditional on CMG; it
does not create a Schur-without-CMG solver. Native fixest remains a third
comparison with its own, numerically different convergence criterion.

Install the companion `cmgr` source package first. The local sources are
`/Users/johannes/Git/CMG/cmgr` and `/Users/johannes/Git/fixest-cmg`. Build and use
an isolated R library. The current local installation is at
`~/research/fereg/fixest-cmg-current/library`. Reinstall after changing cmgr's
source or bundled CMG revision; updating a Git checkout does not update R's
installed package.

```sh
mkdir -p "$HOME/research/fereg/fixest-cmg-current/library"
R CMD INSTALL --library="$HOME/research/fereg/fixest-cmg-current/library" /Users/johannes/Git/CMG/cmgr
R_LIBS_USER="$HOME/research/fereg/fixest-cmg-current/library" R CMD INSTALL \
  --library="$HOME/research/fereg/fixest-cmg-current/library" /Users/johannes/Git/fixest-cmg
```

In a fresh R session:

```r
.libPaths(c(path.expand("~/research/fereg/fixest-cmg-current/library"), .libPaths()))
library(fixest)
readLines(system.file("CMG-SOURCE", package = "cmgr"))
```

The installed `CMG-SOURCE` record identifies the exact source commit used by a
build. The November 7 monorepo build records `1f577d5`; its numerical CMG code
is unchanged from `e5613e8` and retains the default 32-level cap. CMGr is a
separately installable R package within the CMG repository.
Archived benchmark libraries are historical evidence, not the current installation.

Exactly two categorical intercept fixed effects are supported, including
weighted OLS, FE-only, IV, multiple outcomes, stepwise regressors, split samples,
offsets, `feols.fit()` and `demean()`. Unsupported effects and GLM/PPML requests
fail explicitly. Named-effect recovery remains native. Configure tolerance and
iteration limits through `CmgDemeaner`; the native controls cannot also be supplied.
Set the backend when creating an `only.env` environment, then reuse that environment
without resupplying the backend argument.

Calls without `demeaner` take the original native path and do not require cmgr.
No existing positional argument has moved. No fixest C++ source is changed.

The backend returns original-unit weighted residuals, shares a hierarchy within
each demeaning call, supports interruption, and fails on numerical errors rather
than silently switching solvers. It has no cache across calls and does not yet
offer a memory cap. CMG and native MAP tolerance values have different meanings;
equal values are not an accuracy-matched benchmark.

Run `tests/cmg_tests.R` with both packages installed for the focused integration
checks. The companion package contains independent dense-QR projection checks.
The experiment directory retains timing scripts, numerical comparisons, an
interrupt test, and default-path comparisons against pristine upstream.
The current platform evidence is Apple Silicon macOS, not release qualification.

## Shared backend interface and withinr

```r
native <- feols(y ~ x | worker + firm, data, demeaner = MapDemeaner())
lsmr <- feols(y ~ x | worker + firm, data,
              demeaner = LsmrDemeaner(fixef_atol = 1e-8, fixef_btol = 1e-8,
                                      preconditioner = "additive"))
cmg <- feols(y ~ x | worker + firm, data,
             demeaner = CmgDemeaner(fixef_tol = 1e-8, schur = TRUE))
```

`fixef_maxiter` sets the iteration cap for all three constructors. CMG's former
`tol` and `maxiter` names remain compatibility aliases; an alias and its canonical
name cannot be supplied together. Older serialized CMG configurations are
normalized when used. Native tolerance/acceleration arguments cannot accompany
any non-NULL `demeaner`. `setFixest_estimation(demeaner = ...)` accepts the same
objects. Constructors for optional backends report missing packages immediately.

PR #647 supplies an explicit MAP/LSMR dispatcher, not a generic plugin registry.
Installing CMGr alone therefore does not add CMG support to that PR: a thin
constructor, configuration handling and a CMG dispatcher case are still needed.
The withinr adapter accepts both the older `solve_batch()` entry point and
`within_solve_batch()` from within PR #84. Nonconvergence is an error, including
in multiple-model fits; it is not silently accepted as a completed projection.

Withinr uses a process-wide Rayon pool and has no per-call thread argument in
these bindings. Set `RAYON_NUM_THREADS` before starting R for reproducible
comparisons. `nthreads` controls native fixest and CMGr work but does not resize
withinr's pool. Withinr's one-way case uses native demeaning, as in #647.

Run `tests/demeaner_tests.R` for the common interface and independent weighted-QR
checks, including withinr when installed; run `tests/cmg_tests.R` for the broader
CMG on/off integration checks. The native C++ demeaning implementation is unchanged.

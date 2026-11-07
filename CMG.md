# Experimental two-way CMG backend

This branch adds an optional `cmgr` backend at the existing fixed-effect
residualization boundary. It is based on upstream fixest `cf946133` (0.14.2).
It is an independent local implementation; fixest PR #647 remains the place
to coordinate a general backend interface.

```r
library(fixest)
fit <- feols(y ~ x1 + x2 | worker + firm, data,
             demeaner = CmgDemeaner(tol = 1e-8), nthreads = 4)
```

Install the companion `cmgr` source package first. The local sources are
`/Users/johannes/Git/cmgr` and `/Users/johannes/Git/fixest-cmg`. Build and use
an isolated R library. The current local installation is at
`~/research/fereg/fixest-cmg-current/library`. Reinstall after changing cmgr's
source or bundled CMG revision; updating a Git checkout does not update R's
installed package.

```sh
mkdir -p "$HOME/research/fereg/fixest-cmg-current/library"
R CMD INSTALL --library="$HOME/research/fereg/fixest-cmg-current/library" /Users/johannes/Git/cmgr
R_LIBS_USER="$HOME/research/fereg/fixest-cmg-current/library" R CMD INSTALL \
  --library="$HOME/research/fereg/fixest-cmg-current/library" /Users/johannes/Git/fixest-cmg
```

In a fresh R session:

```r
.libPaths(c(path.expand("~/research/fereg/fixest-cmg-current/library"), .libPaths()))
library(fixest)
readLines(system.file("CMG-SOURCE", package = "cmgr"))
```

The current build uses CMG `e5613e8` and its default 32-level cap. The dated
`fixest-cmg-20261106/implementation` and `main-review-20261107` libraries preserve
previous experiments; they are not the current installation.

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

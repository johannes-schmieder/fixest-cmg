\name{CmgDemeaner}
\alias{CmgDemeaner}
\title{Configure the optional two-way CMG backend}
\usage{CmgDemeaner(fixef_tol = 1e-8, fixef_maxiter = 10000L, schur = FALSE,
  tol = NULL, maxiter = NULL)}
\arguments{
\item{schur}{Enable exact low-degree Schur elimination. Defaults to FALSE, preserving the full CMG graph.}
\item{fixef_tol}{Relative weighted fixed-effect moment tolerance.}
\item{fixef_maxiter}{Maximum PCG iterations per residualized column.}
\item{tol,maxiter}{Compatibility aliases. Do not supply an alias and its corresponding canonical argument together.}
}
\description{Configure explicit CMG residualization for two categorical intercept
fixed effects in OLS or \code{demean()}. Requires the optional experimental cmgr
package. The native backend remains the default.}
\details{Supply the result as \code{demeaner = CmgDemeaner()} to \code{feols},
\code{feols.fit} or \code{demean}. Put convergence controls in this constructor;
do not also supply the native tolerance, iteration or acceleration controls.
The tolerance differs from the native MAP iterate-change tolerance. Unsupported
FE dimensions, varying slopes, stepwise FE specifications and GLM requests fail
explicitly. Native fixed-effect extraction and prediction remain available.

With \code{schur = TRUE}, eliminate levels of the fixed-effect block
with more levels (normally workers; the first block breaks ties) if they have
at most five distinct neighbors. Repeated worker--firm observations are combined
before counting degree. Higher-degree levels remain in the actual solved graph.
This is exact partial elimination, without sampling. CMG and its PCG settings
are unchanged. No observations are dropped.}
\value{A configuration of class \code{fixest_demeaner_cmg}.}

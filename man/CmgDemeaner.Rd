\name{CmgDemeaner}
\alias{CmgDemeaner}
\title{Configure the optional two-way CMG backend}
\usage{CmgDemeaner(tol = 1e-8, maxiter = 10000L)}
\arguments{
\item{tol}{Relative weighted fixed-effect moment tolerance.}
\item{maxiter}{Maximum PCG iterations per residualized column.}
}
\description{Configure explicit CMG residualization for two categorical intercept
fixed effects in OLS or \code{demean()}. Requires the optional experimental cmgr
package. The native backend remains the default.}
\details{Supply the result as \code{demeaner = CmgDemeaner()} to \code{feols},
\code{feols.fit} or \code{demean}. Put convergence controls in this constructor;
do not also supply the native tolerance, iteration or acceleration controls.
The tolerance differs from the native MAP iterate-change tolerance. Unsupported
FE dimensions, varying slopes, stepwise FE specifications and GLM requests fail
explicitly. Native fixed-effect extraction and prediction remain available.}
\value{A configuration of class \code{fixest_demeaner_cmg}.}

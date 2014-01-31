\name{jackstraw}
\alias{jackstraw}
\title{Non-Parametric Jackstraw for Principal Component Analysis (PCA)}
\usage{
  jackstraw(dat, PC = NULL, r = NULL, s = NULL, B = NULL,
    covariate = NULL, verbose = TRUE, seed = NULL)
}
\arguments{
  \item{dat}{a data matrix with \code{m} rows as variables
  and \code{n} columns as observations.}

  \item{PC}{a numeric vector of principal components of
  interest. Choose a subset of \code{r} significant PCs to
  be used.}

  \item{r}{a number (a positive integer) of significant
  principal components. See \link{permutationPA} and other
  methods.}

  \item{s}{a number (a positive integer) of ``synthetic''
  null variables. Out of \code{m} variables, \code{s}
  variables are independently permuted.}

  \item{B}{a number (a positive integer) of resampling
  iterations. There will be a total of \code{s*B} null
  statistics.}

  \item{covariate}{a data matrix of covariates with
  corresponding \code{n} observations.}

  \item{verbose}{a logical indicator as to whether to print
  the progress.}

  \item{seed}{a seed for the random number generator.}
}
\value{
  \code{jackstraw} returns a list consisting of
  \item{p.value}{the \code{m} p-values of association tests
  between variables and their principal components}
  \item{obs.stat}{the observed F-test statistics}
  \item{null.stat}{the \code{s*B} null F-test statistics}
}
\description{
  Estimates statistical significance of association between
  variables and their principal components (PCs).
}
\details{
  This function computes \code{m} p-values of linear
  association between \code{m} variables and their PCs. Its
  resampling strategy accounts for the over-fitting
  characteristics due to direct computation of PCs from the
  observed data and protects against an anti-conservative
  bias.

  Provide the data matrix, with \code{m} variables as rows
  and \code{n} observations as columns. Given that there
  are \code{r} significant PCs, this function tests for
  linear association between \code{m} varibles and their
  \code{r} PCs.

  You could specify a subset of significant PCs that you
  are interested in (\code{PC}). If \code{PC} is given,
  then this function computes statistical significance of
  association between \code{m} variables and \code{PC},
  while adjusting for other PCs (i.e., significant PCs that
  are not your interest). For example, if you want to
  identify variables associated with 1st and 2nd PCs, when
  your data contains three significant PCs, set \code{r=3}
  and \code{PC=c(1,2)}.

  Please take a careful look at your data and use
  appropriate graphical and statistical criteria to
  determine a number of significant PCs, \code{r}. The
  number of significant PCs depends on the data structure
  and the context. In a case when you fail to specify
  \code{r}, it will be estimated from a permutation test
  (Buja and Eyuboglu, 1992) using a function
  \link{permutationPA}.

  If \code{s} is not supplied, \code{s} is set to about
  10\% of \code{m} variables. If \code{B} is not supplied,
  \code{B} is set to \code{m*10/s}.
}
\examples{
set.seed(1234)
## simulate data from a latent variable model: Y = BL + E
B = c(rep(1,50),rep(-1,50), rep(0,900))
L = rnorm(20)
E = matrix(rnorm(1000*20), nrow=1000)
dat = B \%*\% t(L) + E
dat = t(scale(t(dat), center=TRUE, scale=TRUE))

## apply the jackstraw
out = jackstraw(dat, PC=1, r=1)

## Use optional arguments
## For example, set s and B for a balance between speed of the algorithm and accuracy of p-values
## out = jackstraw(dat, PC=1, r=1, s=10, B=1000, seed=5678)
}
\author{
  Neo Christopher Chung \email{nc@princeton.edu}
}
\references{
  Chung and Storey (2013) Statistical Significance of
  Variables Driving Systematic Variation in
  High-Dimensional Data. arXiv:1308.6013 [stat.ME]
  \url{http://arxiv.org/abs/1308.6013}

  More information available at
  \url{http://www.jackstraw.info}
}
\seealso{
  \link{permutationPA}
}


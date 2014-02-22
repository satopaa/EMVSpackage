\name{EMVS}
\alias{EMVS}
\title{Create an EMVS-object}
\description{
  This function creates an EMVS-object.
}
\usage{
EMVS(	
		Y,
		X,
		v0s,
		v1 = NULL,
		type,
		beta_init,
		sigma_init,
		epsilon,
		temperature = NULL,
		Z,
		mu,
		Sigma,
		p,
		a,
		b,
		a_v1,
		b_v1,
		v1_g
)
}
\arguments{
  \item{Y}{Response vector}
  \item{X}{Covariates}
}
\value{
  An EMVS-object with the following elements: 
    \item{}
}
\details{}
%\references{
%}
\seealso{
  \code{\link{EMVSpackage}}
}
\examples{}
\author{Ville Satopaa \email{satopaa@wharton.upenn.edu}}
\keyword{variable selection}

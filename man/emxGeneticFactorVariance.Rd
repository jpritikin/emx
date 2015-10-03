\name{emxGeneticFactorVariance}
\alias{emxGeneticFactorVariance}

\title{Creates a variance matrix accoring to the Genetic Factor Model}

\description{
This function creates a variance matrix according to the genetic factor model
}

\usage{
emxGeneticFactorVariance(x, name, values=.8, free=TRUE, lbound=NA, ubound=NA)
}


\arguments{
   \item{x}{character vector. The names of the variables used for the matrix.}
   \item{name}{character. The name of the variance matrix created.}
   \item{values}{numeric vector. The starting values for the lower triangular matrix.}
   \item{free}{logical vector. Whether the lower triangular elements are free.}
   \item{lbound}{numeric vector. Lower bounds on free parameters.}
   \item{ubound}{numeric vector. Upper bounds on free parameters.}
}

%\details{
%TODO add more information on Genetic factor model
%
%}

\value{
    A list with two components. The first component is the factor loadings matrix. The second component is an MxAlgebra, the result of which is the variance matrix implied by the factor loadings.
}

\seealso{
\link{emxCholeskyVariance}
}

%\references{
%
%}

\examples{
   
# Create a loadings matrix
require(emx)
amat <- emxCholeskyVariance(x, 'A')
ahmat <- emxRelatednessMatrix(2, c(1, .5, 1), name='AH')
ab <- emxKroneckerVariance(ahmat, amat[[2]], 'AB')

}


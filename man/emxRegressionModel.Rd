

\name{emxRegressionModel}
\alias{emxRegressionModel}

\title{Create a regression model}

\description{
This function creates a regression model as an MxModel object.
}

\usage{
emxRegressionModel(model, data, run, ...)
}


\arguments{
   \item{model}{formula.  See Details.}
   \item{data}{data used for the model}
   \item{run}{logical. Whether to run the model before returning.}
   \item{...}{Further named arguments to be passed to \code{lm} for the formula}
}

\details{
The \code{model} argument is a formula identical to what is used in \code{link{lm}}.


}

\value{
    An MxModel.
}

\seealso{
\link{lm}
}

\references{

}

\examples{
   
# Example
require(emx)
data(myRegDataRaw)
myrdr <- myRegDataRaw
myrdr[1, 4] <- NA

run <- emxRegressionModel(y~1+x*z, data=myrdr, run=TRUE)
summary(run)

summary(lm(y~1+x*z, data=myrdr))



}

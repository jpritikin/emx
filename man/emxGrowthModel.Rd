

\name{emxGrowthModel}
\alias{emxGrowthModel}

\title{Create a latent growth curve model}

\description{
This function creates a latent growth curve model as an MxModel object.
}

\usage{
emxGrowthModel(model, data, name, run=FALSE, identification, use, ordinal)
}


\arguments{
   \item{model}{character or numeric. See Details.}
   \item{data}{data used for the model}
   \item{name}{character.  Optional name of the model created.}
   \item{run}{logical. Whether to run the model before returning.}
   \item{identification}{Not yet implemented.  How the model is identified.  Currently idngored.}
   \item{use}{character vector.  The names of the variables to use.}
   \item{ordinal}{character vector.  The names of the ordinal variables.}
}

\details{
The \code{model} argument can be either a character or a number that tells the kind of growth curve to make.  If it is a character it currently must be one of "Intercept", "Linear", "Quadratic", "Cubic", "Quartic", or "Quintic", and it produces a polynomial growth curve of the corresponding type.  If it is a number, the function produces a polynomial growth curve of the corresponding order.  Zero is an intercept only, one is linear, two is quadratic; and so on.
}

\value{
    An MxModel
}

\seealso{
\link{emxFactorModel,emxGrowthModel}
}

\references{

}

\examples{
   
# Example
require(emx)
data(myLongitudinalData)
mod <- emxGrowthModel('Linear', data=myLongitudinalData, use=names(myLongitudinalData), run=TRUE)

}

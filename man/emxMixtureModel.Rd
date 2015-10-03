

\name{emxMixtureModel}
\alias{emxMixtureModel}

\title{Create a mixture model}

\description{
This function creates a mxiture model as an MxModel object.
}

\usage{
emxMixtureModel(model, data, run=FALSE, p=NA, ...)
}


\arguments{
   \item{model}{list. The MxModel objects that compose the mixture.}
   \item{data}{data used for the model}
   \item{run}{logical. Whether to run the model before returning.}
   \item{p}{character.  Name of the mixing proportions matrix.}
   \item{...}{Further Mx Objects passed into the mixture model.}
}

\details{
The \code{model} argument is list of MxModel objects.  These are the classes over which the mixture model operates.

}

\value{
    An MxModel.
}

\seealso{
\link{emxLoadings}
}

%\references{
%
%}

\examples{
   
# Example
require(emx)
data(myFADataRaw)
xmap1 <- list(F1=paste0('x', 1:6), F2=paste0('y', 1:3), F3=paste0('z', 1:3))
mod1 <- emxFactorModel(xmap1, data=myFADataRaw, name='m1')

xmap2 <- list(F1=c(paste0('x', 1:6), paste0('y', 1:3), paste0('z', 1:3)))
mod2 <- emxFactorModel(xmap2, data=myFADataRaw, name='m2')

mod <- emxMixtureModel(list(mod1, mod2), data=myFADataRaw, run=TRUE)
summary(mod)


}

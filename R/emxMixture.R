

emxMixtureModel <- function(models, p=NA, data){
	modelNames <- sapply(models, slot, 'name')
	# Check that all models have ML fit functions
	for(i in 1:length(models)){
		if(!(class(models[[i]]$fitfunction) %in% 'MxFitFunctionML')){
			stop('Cannot proceed because model ', OpenMx::omxQuotes(modelNames[i]), ' does not have an ML fit function.')
		}
		models[[i]]$fitfunction <- OpenMx::mxFitFunctionML(vector=TRUE, rowDiagnostics=models[[i]]$fitfunction$rowDiagnostics)
	}
	model <- OpenMx::mxModel(model='MixtureModel', models, data)
	#if(single.na(p)){
		p1 <- OpenMx::mxMatrix('Full', nrow=length(models), ncol=1, lbound=0, ubound=1, values=1/length(models), free=TRUE, name='theSmallPMatrix', labels=paste0('unscaledProportionForModel', modelNames))
		theSmallPMatrix <- NULL
		p2 <- OpenMx::mxAlgebra((1/sum(theSmallPMatrix)) %x% theSmallPMatrix, name='theScaledPMatrix')
		p <- 'theScaledPMatrix'
		model <- OpenMx::mxModel(model, p1, p2)
	#}
	theAlg1 <- OpenMx::mxAlgebraFromString(paste('cbind(', paste(modelNames, '.fitfunction', sep='', collapse=', '), ') %*% ', p, sep=''), name='theMixtureFitVector')
	theMixtureFitVector <- NULL
	theAlg2 <- OpenMx::mxAlgebra(-2*sum(log(theMixtureFitVector)), name='theMixtureFit')
	theFit <- OpenMx::mxFitFunctionAlgebra('theMixtureFit')
	model <- OpenMx::mxModel(model, theAlg1, theAlg2, theFit)
}

#v2 <- OpenMx::mxMixtureModel(list(class1, class2), data=dataRaw)
#v2r <- OpenMx::mxRun(v2)
#summary(v2r)


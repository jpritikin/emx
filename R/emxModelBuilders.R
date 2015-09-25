




#------------------------------------------------------------------------------
emxFactorModel <- function(model, data, name, run=FALSE, identification, use, ordinal){
	if(missing(name)){name <- 'Model'}
	if(missing(use)){use <- sort(unique(unlist(model)))}
	latents <- names(model)
	manifests <- use
	if( nrow(data) == ncol(data) && all(data == t(data)) ){
		data <- data[use, use]
		bdata <- mxData(data, 'cov')
	} else {
		data <- data[,use]
		if(missing(ordinal)){
			if (is.data.frame(data)) {
				ordinalCols <- sapply(data, is.ordered)
				binaryCols <- sapply(data, is.binary)
			} else {
				ordinalCols <- rep(FALSE, numVar)
				binaryCols <- rep(FALSE, numVar)
			}
		} else {
			if(!all(ordinal %in% use)){stop('Some of the ordinal variables are not among those being used.  Either include the ordinal variables in use, or exclude the ordinal variables not being used.')}
			ordinalCols <- rep(FALSE, ncol(data))
			ordinalCols[match(ordinal, use)] <- TRUE
			binaryCols <- sapply(data, is.binary)
		}
		if(!any(is.na(data)) && !any(ordinalCols)){
			bdata <- mxData(cov(data), 'cov', means=colMeans(data), numObs=nrow(data))
		} else {
			bdata <- mxData(data, 'raw')
		}
	}
	mmat <- emxMeans(x=use, free=!ordinalCols)
	lmat <- emxLoadings(x=model)
	rmat <- emxResiduals(x=use, free=!ordinalCols, values=1)
	ka <- emxMeans(x=latents, free=FALSE, name='LatentMeans')
	ph <- emxCovariances(x=latents, type='corr', name='LatentVariances')
	bmodel <- mxModel(name=name,
		lmat, rmat, mmat, ka, ph,
		bdata,
		mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name')),
		mxFitFunctionML()
	)
	if(any(ordinalCols)){
		threshList <- emxThresholds(data, ordinalCols)
		bmodel <- mxModel(bmodel, threshList,
			mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name'),
			thresholds=slot(threshList[[3]], 'name'))
		)
	}
	if(run){return(mxRun(bmodel))}
	return(bmodel)
}








#------------------------------------------------------------------------------
# Growth model

GrowthBasisMatrix <- function(numTimes, deltaT=1, order=2, steps=NA){
	if(all(is.na(steps))){
		m <- numTimes
		steps <- (0:(m-1))*deltaT
	}
	L <- outer(steps, 0:order, "^")
	return(L)
}


resolveGrowthModel <- function(x){
	nums <- c("Intercept", "Linear", "Quadratic", "Cubic", "Quartic", "Quintic")
	if(length(x) != 1){
		stop('wtf am I supposed to do with that?')
	}
	if(is.character(x)){
		if( x %in% nums ){
			return( match(x, nums) - 1 )
		} else {stop(
			paste('I do not understand what model is meant by', omxQuotes(x), '. ', 
			'The model argument must be a number or one of the following: ', 
			omxQuotes(nums), sep=''))
		}
	}
	else if( round(x) == x && x >= 0){
		return(x)
	}
	else{stop('I do not understand your model arument.')}
}

emxGrowthModel <- function(model, data, name, run=FALSE, identification, use, ordinal){
	if(missing(name)){name <- 'Model'}
	if(missing(use)){stop('You must specify use for growth models')}
	order <- resolveGrowthModel(model)
	latents <- paste('F', 0:order, sep='')
	manifests <- use
	if( nrow(data) == ncol(data) && all(data == t(data)) ){
		bdata <- mxData(data, 'cov')
	} else {
		if(!any(is.na(data))){
			bdata <- mxData(cov(data), 'cov', means=colMeans(data), numObs=nrow(data))
		} else {
		bdata <- mxData(data, 'raw')
		}
	}
	mmat <- emxMeans(x=use, values=0, free=FALSE, type='equal')
	lval <- GrowthBasisMatrix(length(use), order=order)
	lmat <- mxMatrix('Full', length(use), length(latents), FALSE, values=lval, name='Loadings', dimnames=list(manifests, latents))
	rmat <- emxResiduals(x=use, type='identical')
	ka <- emxMeans(x=latents, free=TRUE, name='LatentMeans')
	ph <- emxCovariances(x=latents, type='full', name='LatentVariances')
	bmodel <- mxModel(name=name,
		lmat, rmat, mmat, ka, ph,
		bdata,
		mxExpectationLISREL(
			LX=slot(lmat, 'name'),
			TD=slot(rmat, 'name'),
			TX=slot(mmat, 'name'),
			KA=slot(ka, 'name'),
			PH=slot(ph, 'name')),
		mxFitFunctionML()
	)
	if(run){return(mxRun(bmodel))}
	return(bmodel)
}




#------------------------------------------------------------------------------
# Regression with FIML

emxRegressionModel <- function(model, data, run, ...){
	theFormu <- model #a formula as passed to lm
	theFrame <- lm(theFormu, data=data, na.action=na.pass, method='model.frame', ...)
	theTerms <- attr(theFrame, 'terms')
	theRespo <- model.response(theFrame, 'numeric')
	theMatri <- model.matrix(theTerms, theFrame)
	mnam <- gsub(":", "_x_", colnames(theMatri))
	colnames(theMatri) <- mnam
	whichInt <- match('(Intercept)', colnames(theMatri))
	if(!is.na(whichInt)){
		theMatri <- theMatri[, -whichInt]
	}
	
	namRespo <- names(theFrame)[1]
	namTerms <- colnames(theMatri)
	namAll <- c(namRespo, namTerms)
	cdat <- cbind(theRespo, theMatri)
	colnames(cdat) <- namAll
	regPaths <- mxPath(from=namTerms, to=namRespo, free=TRUE, values=.8, labels=paste(namRespo, '_on_', namTerms, sep=''))
	predVars <- mxPath(from=namTerms, arrows=2, free=TRUE, values=1, labels=paste('var', namTerms, sep='_'))
	residVar <- mxPath(from=namRespo, arrows=2, free=TRUE, values=1, labels=paste('residVar', namRespo, sep='_'))
	theInter <- mxPath(from='one', to=namRespo, labels='Intercept', values=0, free=!is.na(whichInt))
	theMeans <- mxPath(from='one', to=namTerms, labels=paste('mean', namTerms, sep='_'))
	theData <- mxData(cdat, 'raw')

	model <- mxModel('model', type='RAM', manifestVars=namAll, theInter, regPaths, residVar, theMeans, predVars, theData)
	
	if(run==TRUE){
		model <- mxRun(model)
	}
	return(model)
}





`cvpl` <-
function(x, surv.time, surv.event, strata.cox=NULL, nfold=1, setseed=NULL, na.rm=FALSE, verbose=FALSE) {

	require(survival)
	
	data <- as.matrix(data);
	
	#remove NA values
	cc.ix <- complete.cases(data, surv.time, surv.event, strata.cox);
	surv.time <- surv.time[cc.ix];
	surv.event <- surv.event[cc.ix];
	data <- as.matrix(data[cc.ix, ,drop=FALSE]);
	strata.cox <- strata.cox[cc.ix];
	nr <- sum(cc.ix);
	if (!all(cc.ix) && !na.rm) { stop("NA values are present!"); }
	if(verbose) { cat(sprintf("%i cases (%i cases are removed due to NA values)\n", nr, sum(!cc.ix))); }
	
	print("function not implemented yet")
	return(-1)
}


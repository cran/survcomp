`D.index` <-
function(x, surv.time, surv.event, weights, strat, alpha=0.05, na.rm=FALSE, ...) {
	require(survival)
	require(SuppDists)
	if(!missing(weights)) {
		if(length(weights) != length(x)) { stop("bad length for parameter weights!") }
	} else { weights <- rep(1,  length(x)) }
	if(!missing(strat)) {
		if(length(strat) != length(x)) { stop("bad length for parameter strat!") }
	} else { strat <- rep(1,  length(x)) }
	cc.ix <- complete.cases(x, surv.time, surv.event, weights, strat)
	if(all(!cc.ix)) {
		data <- list("x"=x, "z"=rep(NA, length(x)), "surv.time"=surv.time, "surv.event"=surv.event, "weights"=weights, "strat"=strat)
		return(list("d.index"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=0, "coxm"=NA, "data"=data))	
	}
	if(any(!cc.ix) & !na.rm) { stop("NA values are present!") }
	sx <- x[cc.ix]
	oo <- order(sx, decreasing=FALSE)
	sx <- sx[oo]
	stime <- surv.time[cc.ix][oo]
	sevent <- surv.event[cc.ix][oo]
	sweights <- weights[cc.ix][oo]
	sstrat <- strat[cc.ix][oo]	
	kap <- sqrt(8/pi)
	z <- kap^-1 * normOrder(N=length(sx))
	#ties?
	dup <- duplicated(sx)
	if(any(dup)) {
		udup <- unique(sx[dup])
		for(i in 1:length(udup)) { z[sx == udup[i]] <- mean(z[sx == udup[i]]) }
	}
	z2 <- x
	z2[!cc.ix] <- NA
	z2[cc.ix] <- z[match(1:length(sx), oo)]
	data <- list("x"=x, "z"=z2, "surv.time"=surv.time, "surv.event"=surv.event, "weights"=weights, "strat"=strat)
	#fit the cox model		
	options(warn=2)
	rr <- try(coxph(Surv(stime, sevent) ~ strata(sstrat) + z, weights=sweights, ...))
	options(warn=0)
	if(class(rr) == "try-error") {
		res <- list("d.index"=NA, "coef"=NA, "se"=NA, "lower"=NA, "upper"=NA, "p.value"=NA, "n"=sum(cc.ix), "coxm"=NA, "data"=data)
	} else {
		dicoef <- rr$coefficients
		dise <- sqrt(drop(rr$var))
		names(dicoef) <- names(dise) <- NULL
		res <- list("d.index"=exp(dicoef), "coef"=dicoef, "se"=dise, "lower"=exp(dicoef - qnorm(alpha / 2, lower.tail=FALSE) * dise), "upper"=exp(dicoef + qnorm(alpha / 2, lower.tail=FALSE) * dise), "p.value"=pchisq((dicoef / dise)^2, df=1, lower.tail=FALSE), "n"=rr$n, "coxm"=rr, "data"=data)
	}
	
	return(res)
}
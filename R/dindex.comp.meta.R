`dindex.comp.meta` <-
function(list.dindex1, list.dindex2, hetero=FALSE) {

	if(length(list.dindex1) != length(list.dindex2)) { stop("the concordance indices are computed from different number of samples!") }

	n <- 0
	x1 <- x1.se <- x2 <- x2.se <- corz <- corz.se <- NULL
	for(i in 1:length(list.dindex1)) {
		nn <- list.dindex1[[i]]$n
		if(nn != list.dindex2[[i]]$n) { stop("the number of samples to compute the concordance indices is not the same!") }
		n <- n + nn
		x1 <- c(x1, list.dindex1[[i]]$coef)
		x1.se <- c(x1.se, list.dindex1[[i]]$se)
		x2 <- c(x2, list.dindex2[[i]]$coef)
		x2.se <- c(x2.se, list.dindex2[[i]]$se)
		cort <- cor(list.dindex1[[i]]$data$z, list.dindex2[[i]]$data$z, method="pearson", use="complete.obs")
		corz <- c(corz, fisherz(cort, inv=FALSE))
		if(nn > 3) { corz.se <- c(corz.se, 1 / sqrt(nn - 3)) } else { corz.se <- c(corz.se, NA) }
	}
	x1.meta <- combine.est(x=x1, x.se=x1.se, hetero=hetero, na.rm=TRUE)
	x2.meta <- combine.est(x=x2, x.se=x2.se, hetero=hetero, na.rm=TRUE) 
	r <- fisherz(combine.est(x=corz, x.se=corz.se, na.rm=TRUE, hetero=hetero)$estimate, inv=TRUE)

	if(abs(r) < 1) {
		t.stat <- (x1.meta$estimate - x2.meta$estimate) / sqrt(x1.meta$se^2 + x2.meta$se^2 - 2 * r * x1.meta$se * x2.meta$se)
		diff.ci.p <- pt(q=t.stat, df=n - 1, lower.tail=FALSE)
	} else { diff.ci.p <- 1 }
	return(list("p.value"=diff.ci.p, "dindex1"=exp(x1.meta$estimate), "dindex2"=exp(x2.meta$estimate)))
}

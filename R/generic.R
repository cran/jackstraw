fstat.adjust = function(dat, LV, ALV, decomposition, ext.covariate, n) {
	# Calculate F-statistics and parametric significance when there exists ALV
	
	if(is.null(ext.covariate)) {
		model.alt = model.matrix(seq(n) ~ decomposition$v[,LV] + decomposition$v[,ALV])
		model.null = model.matrix(seq(n) ~ decomposition$v[,ALV])
	}
	if(!is.null(ext.covariate)) {
		model.alt = model.matrix(seq(n) ~ decomposition$v[,LV] + decomposition$v[,ALV] + ext.covariate[,-1])
		model.null = model.matrix(seq(n) ~ decomposition$v[,ALV] + ext.covariate[,-1])
	}
	
	RSS.alt = RSS(dat, model.alt)
	RSS.null = RSS(dat, model.null)
	
	fstat = (RSS.null - RSS.alt)/(ncol(model.alt)-ncol(model.null)) / (RSS.alt/(n-ncol(model.alt)))
	fstat.pval = 1-pf(fstat, ncol(model.alt)-ncol(model.null), n-ncol(model.alt))
	
	return(list(fstat=fstat, p.value=fstat.pval))
}

fstat = function(dat, LV, decomposition, ext.covariate, n) {
	# Calculate F-statistics and parametric significance

	if(is.null(ext.covariate)) {
		model.alt = model.matrix(seq(n) ~ decomposition$v[,LV])
		model.null = model.matrix(seq(n) ~ 1)
	}
	if(!is.null(ext.covariate)) {
		model.alt = model.matrix(seq(n) ~ decomposition$v[,LV] + ext.covariate[,-1])
		model.null = model.matrix(seq(n) ~ ext.covariate[,-1])
	}
	
	RSS.alt = RSS(dat, model.alt)
	RSS.null = RSS(dat, model.null)
	
	fstat = (RSS.null - RSS.alt)/(ncol(model.alt)-ncol(model.null)) / (RSS.alt/(n-ncol(model.alt)))
	fstat.pval = 1-pf(fstat, ncol(model.alt)-ncol(model.null), n-ncol(model.alt))
	
	return(list(fstat=fstat, p.value=fstat.pval))
}

RSS = function(dat,mod){
	# Calculate residual sum of squares, comparing to an alternative model

	if(is.vector(dat)) {
		m = 1
		n = length(dat)
	} else {
		m = dim(dat)[1]	
		n = dim(dat)[2]
	}
	
	Id = diag(n)
	res = dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))	#Residuals of the model
	rss = res^2 %*% rep(1,n)
	rsst = (dat-mean(dat))^2 %*% rep(1,n)
	r2.total = 1 - (rss/rsst)	#R^2 explained by the model
	residual.r2 = rss/rsst		#R^2 not explained by the model
	
	return(residual.r2)
}

getp = function(lr,lr0) {
	# Get resampled p-values, pulling across variables (e.g., genes)
	# lr: observed statistics
	# lr0: null statistics (i.e. from resampled residuals)
	
	m = length(lr)
	v = c(rep(TRUE,m),rep(FALSE,length(lr0)))
	v = v[rev(order(c(lr,lr0)))]
	u = 1:length(v)
	w = 1:m
	p = ((u[v==TRUE]-w)+1)/(length(lr0)+2)
	p = p[rank(-lr)]
	return(p)
}
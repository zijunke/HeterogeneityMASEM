# vector to matrix
v2m <- function(vec,p,corr= T){
	M = matrix(0,p,p)
	M[lower.tri(M)] = vec
	M = M + t(M)
	if(corr==TRUE){
		diag(M) = 1
	}else{
		diag(M) = diag(M)/2
	}
	return(M)
}

# impute missing values in covariance / correlation matrices of each study
# to obtain a rough estimate of the covariance matrix of covariance / correlation matrix
# weighted average correlation
Mimpute <- function(R,N,missing){
	if(is.null(missing)){
		return(R)
	}else{
		na.pos = which(is.na(R),arr.ind = TRUE)
		mu.N = mean(N)
		Rbar = apply(R,2,mean,na.rm = TRUE)# Becker's mean r
		
		for(coli in unique(na.pos[,2])){
			id = na.pos[(na.pos[,2] == coli),1]
			R[id,coli] = Rbar[coli]
		}
		return(R)
	}
}

# change the coordinating system of a vectorized matrix to the coordinating system of 
# the original matrix
# e.g., from vS to S, the former uses one coordinate (vil), whereas the latter uses two (j,k).
Get.vi2jk <- function(p,diag.incl=FALSE,byrow=FALSE){
	A = matrix(1,p,p)
	if(diag.incl ==FALSE){
		pp = p*(p-1)/2
		vi2jk <- matrix(NA,pp,3)
		vi2jk[,3] <- 1:pp
		if(byrow == FALSE){
			vi2jk[,1:2] <- which(lower.tri(A)==1,arr.ind = TRUE)
		}else{
			vi2jk[,1:2] <- which(upper.tri(A)==1,arr.ind = TRUE)
		}
		colnames(vi2jk) = c('j','k','vi')
	}else{
		pp = p*(p+1)/2
		vi2jk <- matrix(NA,pp,3)
		vi2jk[,3] <- 1:pp
		if(byrow == FALSE){
			vi2jk[,1:2] <- which(lower.tri(A,diag = TRUE)==1,arr.ind = TRUE)
		}else{
			vi2jk[,1:2] <- which(upper.tri(A,diag = TRUE)==1,arr.ind = TRUE)
		}
		colnames(vi2jk) = c('j','k','vi')		
	}
	return(vi2jk)
}

# change the coordinating system of a matrix to the coordinating system of 
# the corresponding vectorized matrix
# e.g., from S to vS, the former uses two coordinates (j,k), whereas the latter uses only one (vil).
Get.jk2vi <- function(vi2jk,p,diag.incl=FALSE){
	jk2vi = matrix(0,p,p)
	jk2vi[vi2jk[,1:2]] = vi2jk[,3]
	if(diag.incl){
		jk2vi = jk2vi + t(jk2vi)
		diag(jk2vi) = diag(jk2vi)/2
	}else{
		pp = p*(p-1)/2
		jk2vi = jk2vi + t(jk2vi) + diag(rep(pp+1,p))
	}
	return(jk2vi)
}

jkvil <- function(p){
	vi2jk	= Get.vi2jk(p)
	j	= vi2jk[,1]
	k	= vi2jk[,2]
	vil	= Get.jk2vi(vi2jk,p)
	return(list(j=j,k=k,vil=vil))
}

# compute the covariance matrix of correlation matrix
# based on Steiger (1980)
Corr.Cov <- function(vR,N,index.list){
	nvR	= length(vR)
	vR	= c(vR,1)
	NvR.cov = matrix(NA,nvR,nvR)
	j = index.list$j
	k = index.list$k
	vil = index.list$vil

	for(vi in 1:nvR){
		NvR.cov[vi,vi] = (1-(vR[vi])^2)^2
	}
	for(vi in 1:(nvR-1)){
	for(vj in (vi+1):nvR){
		NvR.cov[vi,vj] = ((vR[vil[j[vi],j[vj]]]-vR[vi]*vR[vil[k[vi],j[vj]]])*(vR[vil[k[vi],k[vj]]]-vR[vil[k[vi],j[vj]]]*vR[vj])
		 +(vR[vil[j[vi],k[vj]]]-vR[vil[j[vi],j[vj]]]*vR[vj])*(vR[vil[k[vi],j[vj]]]-vR[vi]*vR[vil[j[vi],j[vj]]])
		 +(vR[vil[j[vi],j[vj]]]-vR[vil[j[vi],k[vj]]]*vR[vj])*(vR[vil[k[vi],k[vj]]]-vR[vi]*vR[vil[j[vi],k[vj]]])
		 +(vR[vil[j[vi],k[vj]]]-vR[vi]*vR[vil[k[vi],k[vj]]])*(vR[vil[j[vj],k[vi]]]-vR[vil[k[vi],k[vj]]]*vR[vj]))/2
		NvR.cov[vj,vi] <- NvR.cov[vi,vj]
	}
	}

	vR.cov = NvR.cov/(N)
	vR.cov = as.matrix(nearPD(vR.cov,posd.tol = 1e-5)$mat)
	return(vR.cov)
}

# Use average correlation vector to compute V_psi
Vj <- function(vR.bar,N,pp,Nstudy,index.list){

	mu.N = mean(N)
	S.vR.bar = Corr.Cov(vR.bar,mu.N,index.list)
	inv.S.vR.bar = solve(S.vR.bar)
	tau.vR = array(NA,dim = c(Nstudy,pp,pp))
	S.vR = array(NA,dim = c(Nstudy,pp,pp))
	for(i in 1:Nstudy){
		S.vR[i,,]<- S.vR.bar/N[i]*mu.N
		tau.vR[i,,] <- inv.S.vR.bar/mu.N*N[i]
	}	
	return(list(S.vR = S.vR,tau.vR = tau.vR))
}

wbugs <-function(data,initsl,prm,mfn,prm.del.conv,
	nchains=1,niter=60000,nburnin=30000,nthin=1,wd,
	diagm){
# data: a named list of the data in the likelihood model for OpenBUGS
# initsl: a list with nchains elements; each element is a list of starting values
# prm: vector of names of the parameters to save
# mfn: the file name of the likelihood model for OpenBUGS
# diagm: name of the convergence diagnostic method; either 'Geweke' or 'Gelman'
# The function checks convergence every niter-nburnin iterations
	
	fit = bugs(data,initsl,prm,mfn,
		n.chains=nchains,n.iter=niter,n.burnin=nburnin,n.thin=1,
		debug=F,saveExec=T,working.directory = wd)

	for(tryi in 2:20){
		print(paste0('Iteration: ',tryi*(niter-nburnin)))
		fit.coda = read.openbugs(stem="",thin = nthin,quiet=T)
		del.id = na.omit(match('ppp',varnames(fit.coda)))
		if(is.null(prm.del.conv)==0){
			for(i in 1:length(prm.del.conv)){
			   del.id = c(del.id,grep(prm.del.conv[i],varnames(fit.coda)))
			}
		}
		#print(summary(fit.coda),3)
		if(diagm=='Geweke'){
			if(length(del.id)>0){
				tmp.conv = geweke.diag(fit.coda[,-del.id])[[1]]$z
			}else{ tmp.conv = geweke.diag(fit.coda)[[1]]$z }
			crit = (sum((abs(tmp.conv)>1.96),na.rm = T)==0)
		}else if(diagm=='Gelman'){
			if(length(del.id)>0){
				tmp.conv = gelman.diag(fit.coda)$psrf[-del.id,2]
			}else{ tmp.conv = gelman.diag(fit.coda)$psrf[,2] }
			crit = (sum((tmp.conv>1.1),na.rm = T)==0)
		}
		if(crit){
			print(tmp.conv)
			print(summary(fit.coda),3)	
			break
		}else{	
			fit = bugs(data,initsl,prm,mfn,
			n.chains=nchains,n.iter=niter-nburnin+1,n.burnin=1,n.thin=1,
			restart=T,saveExec=T,working.directory = wd)
		}
	}
	ppp.id = match('ppp',prm)
	sel = NA
	if(is.na(ppp.id)){
		nprm = length(prm)
		for(i in 1:nprm){ 
			sel = c(sel,grep(prm[i],rownames(summary(fit.coda)$quantiles)))
		}
	}else{
		prm = prm[-ppp.id]
		nprm = length(prm)
		for(i in 1:nprm){ 
			sel = c(sel,grep(prm[i],rownames(summary(fit.coda)$quantiles))) 
		}
	}
	sel = sel[-1]
	sel = unique(sel)
	
	if(is.na(ppp.id)){
		est = round(summary(fit.coda)$quantiles[sel,'50%'],3)
	}else{
		ppp = summary(fit.coda)$statistics['ppp','Mean']
		namesl = c(rownames(summary(fit.coda)$quantiles)[sel],'ppp')
		est = round(c(summary(fit.coda)$quantiles[sel,'50%'],ppp),3)
		names(est) = namesl
	}
	psd = round(summary(fit.coda)$statistics[sel,'SD'],3)
	if(diagm=='Geweke'){
		CIl = round(HPDinterval(fit.coda,prob = .95)[[1]][sel,1],3)
		CIu = round(HPDinterval(fit.coda,prob = .95)[[1]][sel,2],3)
	}else if(diagm=='Gelman'){
		fit.coda.l = do.call(rbind,fit.coda)
		HPDCI = HPDinterval(mcmc(fit.coda.l),prob = .95)
		CIl = HPDCI[sel,1]
		CIu = HPDCI[sel,2]
	}
	conv = round(c(tryi,tmp.conv),3)
	return(list(est=est,psd=psd,CIl=CIl,CIu=CIu,conv=conv,
		DIC=fit$DIC,fit.coda=fit.coda))
}

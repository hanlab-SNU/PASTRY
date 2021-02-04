#!/usr/bin/env Rscript

## PASTRY and LS simulation code
## This code only consider overlapping controls between the studies(all) of the meta-analysis

## Now load or install & load all
if(!require(optparse)){
    install.packages("optparse")
    library(optparse)
}

## Option argparse
option_list = list(
    make_option(c("-n","--N"), type="integer", default=5,
              help=":# of studies combined (default : 5)", metavar="integer"), # meta-analysis studies	
    make_option(c("-r", "--RRs"), type="numeric", default=1.3, # Risk Ratio
              help=":Risk Ratio (default : 1.3)", metavar="float"),
    make_option(c("-k", "--nrep"), type="integer", default=5000, # simulation N
              help=":# of Simulation (default : 5000)", metavar="integer"),                
    make_option(c("-m", "--maf"), type="numeric", default=0.3, # MAF
              help=":MAF ", metavar="float"),         
    make_option(c("-s", "--shared"), type="integer", default=1000,
              help=":# of shared control (default : 1000)", metavar="integer"),  # shared control N
    make_option(c("-l","--cont"), type="integer", default=1000,
              help=":# of control (default : 1000)", metavar="integer"),  # control 
    make_option(c("-e","--case"), type="integer", default=1000,
              help=":# of case (default : 1000)", metavar="integer"),            # case N
    make_option(c("-t","--thrs"), type="numeric", default=5E-8,
            help=":Threshold (default : 5e-08)", metavar="float")          # threshold
); 

 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sim_PASTRY <- function(RRs=opt$RRs, nrep=opt$nrep, ncase=opt$case, ncont=opt$cont, ncontshared=opt$shared, nstudy=opt$N, t=opt$thrs, maf=opt$maf) {

	## Study-specific controls: ncontspe (Type: nstudy sized vector)
	## Shared controls: ncontshared (Type: nstudy sized vector
        ## Number of cases: ncase (Type: nstudy sized vector)
        ## Number of controls: ncontspe + ncontshared
	## Minor allele frequency: MAF (default=0.3)
        ## Prevalence: preval (default:0.00)

 	m=nstudy
 	maf=maf
	preval=0.000
	casenshared=0
	casensplit=casenshared/m
	contnsplit=ncontshared/m
    RRs=rep(RRs, nstudy)
    ncase=rep(ncase,nstudy)
    ncont=rep(ncont,nstudy)
    ncontspe=ncont-rep(ncontshared,nstudy)


	cat(" * Check the parameter : \n")
	cat("  - Studies # for meta-analysis : ",nstudy,"\n")
	cat("  - Risk Ratio : ",RRs,"\n")
	cat("  - MAF : ",maf,"\n")
	cat("  - Simulation # : ",nrep,"\n")
	cat("  - Case # : ",ncase,"\n")
	cat("  - Non-shared cont # : ",ncontspe,"\n")
	cat("  - Shared cont # : ",ncontshared,"\n")
	cat("  - Threshold : ",t,"\n\n")

	## Case/control frequencies
	casemaf=(RRs*maf) / ((RRs-1)*maf + 1)
	contmaf=(maf-preval*casemaf)/(1-preval)
    
	## Lin_correlation matrix
	R = diag(m)

	for (i in 1:m) {
		for (j in 1:m) {
			if (i == j) {
				next
			} else {
				first=ncontshared*sqrt(((ncase[i]*ncase[j]))/((ncontspe[i]+ncontshared)*(ncontspe[j]+ncontshared)))
				sec=casenshared*sqrt(((ncontspe[i]+ncontshared)*(ncontspe[j]+ncontshared))/((ncase[i])*(ncase[j])))
				rij=first/sqrt((ncase[i]+ncontspe[i]+ncontshared)*(ncase[j]+ncontspe[j]+ncontshared))
				R[i,j] = R[j,i] = rij
			}
		}
	}

	invR = solve(R)
    #cat(paste("\n\n Lin correlation matrix : \n ",R))

	## FE subroutine
	FE <- function(beta,stders){
		inv.var=stders**(-2)
		inv.var.beta=(inv.var %*% beta)/sum(inv.var)
		inv.var.beta.var = 1/sum(inv.var)
		new.z = inv.var.beta / sqrt(inv.var.beta.var)
		p.FE = 2*pnorm(-abs(new.z),lower.tail=sign)
		list(z=new.z,p=p.FE,beta=inv.var.beta,se=sqrt(inv.var.beta.var))
	}

	##### Overlapping Subjects #####

	pvalue.LS.FE = rep(0, nrep)
	pvalue.PASTRY.FE = rep(0, nrep)
	pvalue.split.FE = rep(0, nrep)

	# m = 2 two cohort study
	i = 1; j = 2
	
	## Now generate studies
	record=2
	for (eachrep in 1:nrep) {
		betas.nonde = rep(0, m)
		stder.nonde = rep(0, m) 
		betas.split = rep(0, m)
		stder.split = rep(0, m) 
		casemac=rep(0,m)
		contmac=rep(0,m)
		contsplitmac=rep(0,m)

		for (i in 1:m) {
			#case-specific
			casemac[i]=rbinom(1, ncase[i]*2, casemaf[i])
			#cont-specific
			contmac[i]=rbinom(1, ncontspe[i]*2, contmaf[i])
			#cont-shared
			contsplitmac[i]=rbinom(1, 2*contnsplit, contmaf[i])
		}

		contsharedmac=sum(contsplitmac)
		beta.nonde = NULL
		stderr.nonde = NULL
		beta.split = NULL
		stderr.split = NULL

		for (i in 1:m) {
			## non-splitting
			n11=casemac[i]
			n10=2*ncase[i] - casemac[i]
			n01=contmac[i]+contsharedmac
			n00=2*ncontspe[i] + 2*ncontshared - contmac[i] - contsharedmac
			betas.nonde[i] = log(n11*n00/n01/n10)
  			stder.nonde[i] = sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)
			beta.nonde = c(beta.nonde, betas.nonde[i])
			stderr.nonde = c(stderr.nonde, stder.nonde[i])

			## splitting
			n11=casemac[i]
			n10=2*ncase[i] - casemac[i]
			n01=contmac[i]+contsplitmac[i]
			n00=2*ncontspe[i] + 2*contnsplit - contmac[i] - contsplitmac[i]
			betas.split[i] = log(n11*n00/n01/n10)
			stder.split[i] = sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)
			beta.split = c(beta.split, betas.split[i])
			stderr.split = c(stderr.split, stder.split[i])
		}

	    ## Lin
		invstdmat <- solve(diag(stder.nonde))
		invsigma <- invstdmat %*% invR %*% invstdmat
		unit <- t(rep(1,nstudy))
		invvar <- unit %*% invsigma %*% t(unit)
		beta_num <- unit %*% invsigma %*% (betas.nonde)
		ls_beta <- beta_num/invvar
		ls_var <- 1/invvar
		pvalue.LS.FE[eachrep] = 2*pnorm(-abs(ls_beta/sqrt(ls_var)),lower.tail=sign)
		
        if (eachrep == 1) {cat(" * PASTRY statistics are generating... for each SNP. \n")
        start_time <- Sys.time()}
        
		######################
		i = 1; j = 2
		ntotali = ncontspe[i] + ncase[i] + ncontshared
        ntotalj = ncontspe[j] + ncase[j] + ncontshared

		p0cont=1-contmaf[i]
        p1cont=contmaf[i]
        p0case=1-casemaf[i]
        p1case=casemaf[i]

		x_0 = matrix(c(1,0,0,0),nrow=2,ncol=2)
        x_1 = matrix(c(1,1,1,1),nrow=2,ncol=2)

		info_i = matrix(c(0,0,0,0),nrow=2,ncol=2)
        info_j = matrix(c(0,0,0,0),nrow=2,ncol=2)

		up_covmat = matrix(c(0,0,0,0),nrow=2,ncol=2)
		######################

		#x <- t(betas.nonde)
        exp_alpha_i = as.numeric(exp(log(ncase[i]/(ncontspe[i]+ncontshared)) - (ls_beta*((ncase[i]*casemaf[i]+(ncontspe[i]+ncontshared)*contmaf[i])/ntotali))))
		exp_alpha_j = as.numeric(exp(log(ncase[j]/(ncontspe[j]+ncontshared)) - (ls_beta*((ncase[j]*casemaf[j]+(ncontspe[j]+ncontshared)*contmaf[j])/ntotalj))))

		x0i = exp_alpha_i ; x1i = as.numeric(exp_alpha_i*exp(ls_beta))
		x0j = exp_alpha_j ; x1j = as.numeric(exp_alpha_j*exp(ls_beta))

        ### Information Matrix for PASTRY
        info_i = 2*ncase[i]*(p0case*(x0i/(1+x0i)**2)*x_0 + p1case*(x1i/(1+x1i)**2)*x_1) + 2*(ncontspe[i]+ncontshared)*(p0cont*(x0i/(1+x0i)**2)*x_0 + p1cont*(x1i/(1+x1i)**2)*x_1)
        info_j = 2*ncase[j]*(p0case*(x0j/(1+x0j)**2)*x_0 + p1case*(x1j/(1+x1j)**2)*x_1) + 2*(ncontspe[j]+ncontshared)*(p0cont*(x0j/(1+x0j)**2)*x_0 + p1cont*(x1j/(1+x1j)**2)*x_1)

        ### Covarinace Matrix for PASTRY
        cov00=p0cont*(x0i*x0j/((1+x0i)*(1+x0j)))*x_0
        cov11=p1cont*(x1i*x1j/((1+x1i)*(1+x1j)))*x_1

        covmatrix=2*ncontshared*(cov00+cov11)
        invi = solve(info_i) ; invj = solve(info_j);
		covfinal=invi %*% covmatrix %*% invj

		#PASTRY correlation
		cor=(covfinal[2,2]/sqrt(invi[2,2]*invj[2,2]))

		#depends on study #
		PASTRY.cor=matrix(rep(cor,nstudy*nstudy),nrow=nstudy,ncol=nstudy)
		diag(PASTRY.cor)=1

 		#################################################

		## new decoupling
		invsigmaPASTRY <- invstdmat %*% solve(PASTRY.cor) %*% invstdmat
		invvarPASTRY <- unit %*% invsigmaPASTRY %*% t(unit)
	    pastry_var <- 1/invvarPASTRY
		beta_num <- unit %*% invsigmaPASTRY %*% (betas.nonde)	
        pastry_beta <- beta_num/invvarPASTRY
		
		pvalue.split.FE[eachrep] = FE(beta.split,stder.split)$p ## pvalue calculation by FE
		pvalue.PASTRY.FE[eachrep] = 2*pnorm(-abs(pastry_beta/sqrt(pastry_var)),lower.tail=sign)

		if(eachrep%%50000==0){cat(eachrep)}
		
}
        end_time <- Sys.time()
		cat(" * # Simulation : ",nrep," \n * Simulation Results : \n  LS : ",sum(pvalue.LS.FE<t)/nrep,", PASTRY",sum(pvalue.PASTRY.FE<t)/nrep,", split_FE",sum(pvalue.split.FE<t)/nrep,"\n")
        cat("\n")
        cat("PASTRY running time in R : ", end_time - start_time, ". \n")
}

cat("================ PASTRY v1.0 -- Simulation ================ \n")
sim_PASTRY()
cat("================ PASTRY Simulation END ================")
cat("\n\n")
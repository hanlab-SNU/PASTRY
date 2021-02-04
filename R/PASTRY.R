#!/usr/bin/env Rscript

## Now load or install & load all
packages = c("data.table", "optparse", "pracma")

package.load <- lapply(packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      suppressMessages(library(x,character.only=TRUE, quietly = T))
      }
  }
)

## Option argparse
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
              help=":input files prefix name (.out)", metavar="character"), # input files prefix
    make_option(c("-m", "--maf"), type="character", default="./",
              help=":maf files directory (casemaf.out, contmaf.out)", metavar="character"), # input maf file directory path          
    make_option(c("-o", "--out"), type="character", default=NULL, 
              help=":output file prefix name (.out)", metavar="character"), # output files prefix
    make_option(c("-n","--N"), type="integer", default=NULL,
              help=":# of studies combined", metavar="integer"), # meta-analysis studies
    make_option(c("-s", "--shared"), type="integer", default=NULL,
              help=":# of shared control", metavar="integer"),  # shared control N
    make_option(c("-l","--cont"), type="integer", default=NULL,
              help=":# of control", metavar="integer"),          # control N
    make_option(c("-e","--case"), type="integer", default=NULL,
              help=":# of case", metavar="integer")            # case N
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if ( length(opt$input) == 0 | length(opt$N) == 0 | length(opt$shared) == 0 | length(opt$case) == 0 | length(opt$cont) == 0) {
    stop("Please provide PASTRY parameters (check the required parameters using 'Rscript PASTRY.R -h')")
}

## Function define
pastry <- function(files=opt$input,maf=opt$maf,outdir=opt$out,N=opt$N, ncase=opt$case, ncont=opt$cont, ncontshared=opt$shared){

        files <- paste(rep(files,N),seq(1,N),sep="")
        
        casemaf <- as.matrix(read.table(paste(maf,"casemaf.out",sep=""),header=T))
        contmaf <- as.matrix(read.table(paste(maf,"contmaf.out",sep=""),header=T))

        #maf=as.matrix(fread(maf,header=T)) # for (.frq.cc files)
        #casemaf=as.numeric(maf[,"MAF_A"]) # for (.frq.cc files)
        #contmaf=as.numeric(maf[,"MAF_U"]) # for (.frq.cc files)

        ncase <- rep(ncase,N)
        ncont <- rep(ncont,N)
        ncontspe <- ncont-ncontshared
        casenshared <- 0  # only consider overlapping controls

        # Load File names and snp #
        cat(paste("You have ",N," Input files. \n", sep=""))

        study <- as.matrix(fread(paste(files[1],".out",sep="")))
        nsnp <- length(study[,1])

       cat(paste("You have ",nsnp," SNPs in each studies. \n \n", sep="")) 

        # Matrix Initialize
        betas <- matrix(0,N,nsnp)
        ses <- matrix(0,N,nsnp)
        ls_stat <- rep(0,nsnp)
        pastry_stat <- rep(0,nsnp)
        unit <- t(rep(1,N)) ## e^t = unit

        betas[1,] <- as.numeric(study[,1])
        ses[1,] <- as.numeric(study[,2])

        ls_stat <- rep(0,nsnp)
        pastry_stat <- rep(0,nsnp)

        # Beta, SE Files loading
        for (l in 2:N){
        study <- as.matrix(fread(paste(files[l],".out",sep="")))
        betas[l,] <- as.numeric(study[,1])
        ses[l,] <- as.numeric(study[,2])
        }

        ## NA handling
        betas[is.na(betas)] <- 0
        ses[is.na(ses)] <- 0

        ## Lin_correlation matrix
        R <- diag(N)
        for (i in 1:N) {
            for (j in 1:N) {
                if (i == j) {next}
                else {
                    first <- ncontshared*sqrt(((ncase[i]*ncase[j]))/((ncontspe[i]+ncontshared)*(ncontspe[j]+ncontshared)))
                    sec <- casenshared*sqrt(((ncontspe[i]+ncontshared)*(ncontspe[j]+ncontshared))/((ncase[i])*(ncase[j])))
                    rij <- first/sqrt((ncase[i]+ncontspe[i]+ncontshared)*(ncase[j]+ncontspe[j]+ncontshared))
                    R[i,j] = R[j,i] = rij
                }
            }
        }

        invR <- solve(R)
        cat(" * Lin Correlation Matrix generated. \n")

        ## Generating Statistics for LS & PASTRY
        for (k in 1:nsnp){
            for (i in 1:N){
                for (j in 1:N){
                    if (i != j){

                        ntotali <- ncont[i]+ncase[i]
                        ntotalj <- ncont[j]+ncase[j]

                        p0cont <- 1-contmaf[k]
                        p1cont <- contmaf[k]

                        p0case <- 1-casemaf[k]
                        p1case <- casemaf[k]

                        x_0 <- matrix(c(1,0,0,0),nrow=2,ncol=2)
                        x_1 <- matrix(c(1,1,1,1),nrow=2,ncol=2)

                        info_i <- matrix(c(0,0,0,0),nrow=2,ncol=2)
                        info_j <- matrix(c(0,0,0,0),nrow=2,ncol=2)

                        up_covmat <- matrix(c(0,0,0,0),nrow=2,ncol=2)

                        ## LS method
                        if (k == 1 && i == 1 && j == 2){cat(" * LS statistics are generating... for each SNP. \n")}
                        invstdmat <- pinv(diag(ses[,k]))
                        invsigma <- invstdmat %*% invR %*% invstdmat
                        invvar <- unit %*% invsigma %*% t(unit)

                        beta_num <- unit %*% invsigma %*% (betas[,k])

                        ls_beta <- beta_num/invvar
                        ls_var <- 1/invvar

                        ls_stat[k] <- ls_beta/sqrt(ls_var)

                        exp_alpha_i <- as.numeric(exp(log(ncase[i]/(ncontspe[i]+ncontshared)) - (ls_beta*((ncase[i]*casemaf[k]+(ncontspe[i]+ncontshared)*contmaf[k])/ntotali))))
                        exp_alpha_j <- as.numeric(exp(log(ncase[j]/(ncontspe[j]+ncontshared)) - (ls_beta*((ncase[j]*casemaf[k]+(ncontspe[j]+ncontshared)*contmaf[k])/ntotalj))))

                        x0i <- exp_alpha_i
                        x1i <- as.numeric(exp_alpha_i*exp(ls_beta))

                        x0j <- exp_alpha_j
                        x1j <- as.numeric(exp_alpha_j*exp(ls_beta))

                        if (k == 1 && i == 1 && j == 2) {cat(" * PASTRY statistics are generating... for each SNP. \n")
                        start_time <- Sys.time()}

                        ### Information Matrix for PASTRY
                        info_i <- 2*ncase[i]*(p0case*(x0i/(1+x0i)**2)*x_0 + p1case*(x1i/(1+x1i)**2)*x_1) + 2*(ncontspe[i]+ncontshared)*(p0cont*(x0i/(1+x0i)**2)*x_0 + p1cont*(x1i/(1+x1i)**2)*x_1)
                        info_j <- 2*ncase[j]*(p0case*(x0j/(1+x0j)**2)*x_0 + p1case*(x1j/(1+x1j)**2)*x_1) + 2*(ncontspe[j]+ncontshared)*(p0cont*(x0j/(1+x0j)**2)*x_0 + p1cont*(x1j/(1+x1j)**2)*x_1)

                        ### Covarinace Matrix for PASTRY
                        cov00 <- p0cont*(x0i*x0j/((1+x0i)*(1+x0j)))*x_0
                        cov11 <- p1cont*(x1i*x1j/((1+x1i)*(1+x1j)))*x_1

                        covmatrix <- 2*ncontshared*(cov00+cov11)

                        invi <- pinv(info_i)
                        invj <- pinv(info_j)
                        covfinal <- invi %*% covmatrix %*% invj
                    }
                }
            }
            
        ### Correlation matrix for PASTRY
        cor <- (covfinal[2,2]/sqrt(invi[2,2]*invj[2,2]))
        PASTRY.cor <- matrix(rep(cor,N*N),nrow=N,ncol=N)
        diag(PASTRY.cor) <- 1

        if (k == 1 && i == 1 && j == 2){print(PASTRY.cor)}
        invsigmaPASTRY <- invstdmat %*% solve(PASTRY.cor) %*% invstdmat
        invvarpastry <- unit %*% invsigmaPASTRY %*% t(unit)

        pastry_var <- 1/invvarpastry
        beta_num <- unit %*% invsigmaPASTRY %*% (betas[,k])

        pastry_beta <- beta_num/invvarpastry
        pastry_stat[k] <- pastry_beta/sqrt(pastry_var)
    } 
    end_time <- Sys.time()

out="PASTRY_"
cat("\n\nAll done. Writing your results... \n")
cat("Your PASTRY output files are generated. \n")
fwrite(format(as.data.frame(ls_stat),digits=3),paste(outdir,"_LS_zstat.r.out",sep=""),sep="\n",row.names=FALSE)
fwrite(format(as.data.frame(pastry_stat),digits=3),paste(outdir,"_PASTRY_zstat.r.out",sep=""),sep="\n",row.names=FALSE)
#fwrite(as.data.frame(casemaf),paste(outdir,"_",out,"casemaf.out",sep=""),sep="\n",row.names=FALSE) # for (.frq.cc files)
#fwrite(as.data.frame(contmaf),paste(outdir,"_",out,"contmaf.out",sep=""),sep="\n",row.names=FALSE) # for (.frq.cc files)
cat("\n")
cat(paste("PASTRY running time in R : ", end_time - start_time, ". \n", sep=""))
}


## Main

cat("\n")
cat("================ PASTRY v 1.0 ================ \n")
cat("==== v.1.0 Seoul National University Hanlab ==== \n")
cat("\n")

pastry()

cat("================ PASTRY END ================")
cat("\n\n")
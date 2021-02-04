#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, time
import numpy as np
import scipy.stats as stats
from numpy import savetxt
from numpy import random as np_random
import pandas as pd


def pastry_sim(RRs, maf, nrep, k, ncase, ncont, ncontshared, t):

    nsnp=k
    betas=np.zeros((k,nsnp))

    unit=np.ones(k)
    unit_t=np.transpose(unit)
    preval=0.000

    RRs=np.repeat(RRs, k)
    nstudy=np.repeat(ncont+ncase,k)
    ncont=np.repeat(ncont,k)
    ncase=np.repeat(ncase,k)
    ncontspe=ncont-np.repeat(ncontshared,k)
    contnsplit=ncontshared/k
    
    print(" * Check the parameter : ")
    print("  - Studies # for meta-analysis : ", k)
    print("  - Risk Ratio : ",RRs)
    print("  - MAF : ",maf)
    print("  - Simulation # : ",nrep)
    print("  - Case # : ",ncase)
    print("  - Non-shared cont # : ",ncontspe)
    print("  - Shared cont # : ",ncontshared)
    print("  - Threshold : ",t,"\n")

    ## Case/control frequencies
    casemaf=(RRs*maf) / ((RRs-1)*maf + 1)
    contmaf=(maf-preval*casemaf)/(1-preval)

    ## Generating Correlation Matrix
    ls_R=np.diagflat(np.repeat(1,k)).astype(float)
    
    for i in range(k):
        for j in range(k):
            if i!=j:
                num=ncontshared*np.sqrt((ncase[i]*ncase[j])/(ncont[i]*ncont[j]))
                den=np.sqrt(nstudy[i]*nstudy[i])
                ls_R[i,j]=num/den
	
    inv_ls_R=np.linalg.inv(ls_R)

    ## FE subroutine
    def FE(beta, stders):
        inv_var=stders**(-2)
        inv_var_beta=(np.dot(inv_var, beta)/sum(inv_var))
        inv_var_beta_var=1/sum(inv_var)
        z=inv_var_beta/np.sqrt(inv_var_beta_var)
        p_FE=stats.norm.pdf(abs(z))*2
        return(pd.DataFrame(data={'z':[z],'p':[p_FE],'beta':[inv_var_beta],'se':[np.sqrt(inv_var_beta_var)]}))

    pval_LS_FE=np.zeros(nrep)
    pval_PASTRY_FE=np.zeros(nrep)
    pval_split_FE=np.zeros(nrep)

    # k = 2 two cohort study
    i = 1; j = 2
	
    ## Now generate studies
    record=2
    beta_nonde=[]
    stderr_nonde=[]
    beta_split=[]
    stderr_split=[]

    for eachrep in range(nrep):
        betas_nonde=np.zeros(k)
        stder_nonde=np.zeros(k)
        betas_split=np.zeros(k)
        stder_split=np.zeros(k)
        casemac=np.zeros(k)
        contmac=np.zeros(k)
        contsplitmac=np.zeros(k)

        for i in range(k):
            #case-specific
            casemac[i]=np.random.binomial(n=ncase[i]*2, p=casemaf[i], size=1)
            #cont-specific
            contmac[i]=np.random.binomial(n=ncontspe[i]*2, p=contmaf[i], size=1)
            #cont-shared
            contsplitmac[i]=np.random.binomial(n=2*contnsplit, p=contmaf[i], size=1)

        contsharedmac=np.sum(contsplitmac)

        for i in range(k):
            ## non-splitting
            n11=casemac[i]
            n10=2*ncase[i]-casemac[i]
            n01=contmac[i]+contsharedmac
            n00=2*ncontspe[i]+2*ncontshared-contmac[i]-contsharedmac
            betas_nonde[i]=np.log(n11*n00/n01/n10)
            stder_nonde[i]=np.sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)
            beta_nonde.append(betas_nonde[i])
            stderr_nonde.append(stder_nonde[i])

            ## splitting
            n11=casemac[i]
            n10=2*ncase[i]-casemac[i]
            n01=contmac[i]+contsplitmac[i]
            n00=2*ncontspe[i]+2*contnsplit-contmac[i]-contsplitmac[i]
            betas_split[i]=np.log(n11*n00/n01/n10)
            stder_split[i]=np.sqrt(1/n11 + 1/n00 + 1/n01 + 1/n10)
            beta_split.append(betas_split[i])
            stderr_split.append(stder_split[i])

        ## Lin
        i = 1; j = 2

        p0cont=1-contmaf[i]
        p1cont=contmaf[i]
        p0case=1-casemaf[i]
        p1case=casemaf[i]

        invstdmat=np.linalg.pinv(np.diagflat(stder_nonde))
        ls_invsigma=np.dot(np.dot(invstdmat,inv_ls_R),invstdmat)
        ls_invvar=np.dot(np.dot(unit,ls_invsigma),unit_t)
        ls_beta_num=np.dot(np.dot(unit,ls_invsigma),betas_nonde)
        ls_beta=ls_beta_num/ls_invvar
        ls_var=1/ls_invvar
        ls_stat=ls_beta/np.sqrt(ls_var)
        pval_LS_FE[eachrep]=stats.norm.pdf(abs(ls_stat))*2

        if eachrep == 1:
           print(" * PASTRY statistics are generating... for each SNP. \n")
           start=time.time()

        ntotali=ncont[i]+ncase[i]
        ntotalj=ncont[j]+ncase[j]
        
        exp_alpha_i = np.exp(np.log(ncase[i]/(ncont[i]))-(ls_beta*(((ncase[i]*casemaf[i])+(ncont[i]*contmaf[i]))/ntotali)))
        exp_alpha_j = np.exp(np.log(ncase[j]/(ncont[j]))-(ls_beta*(((ncase[j]*casemaf[j])+(ncont[j]*contmaf[j]))/ntotalj)))

        x0_i = exp_alpha_i
        x1_i = exp_alpha_i*np.exp(ls_beta)

        x0_j = exp_alpha_j
        x1_j = exp_alpha_j*np.exp(ls_beta)

        ## Information Matrix
        x0=np.array([1,0,0,0]).reshape(2,2)
        x1=np.array([1,1,1,1]).reshape(2,2)

        up_covmat=np.array([0,0,0,0]).reshape(2,2)

        info_i=2*ncase[i]*(p0case*(x0_i/(1+x0_i)**2)*x0 + p1case*(x1_i/(1+x1_i)**2)*x1) + 2*(ncont[i])*(p0cont*(x0_i/(1+x0_i)**2)*x0 + p1cont*(x1_i/(1+x1_i)**2)*x1)
        info_j=2*ncase[j]*(p0case*(x0_j/(1+x0_j)**2)*x0 + p1case*(x1_j/(1+x1_j)**2)*x1) + 2*(ncont[j])*(p0cont*(x0_j/(1+x0_j)**2)*x0 + p1cont*(x1_j/(1+x1_j)**2)*x1)

        ### Covarinace Matrix for PASTRY
        cov00=p0cont*(x0_i*x0_j/((1+x0_i)*(1+x0_j)))*x0
        cov11=p1cont*(x1_i*x1_j/((1+x1_i)*(1+x1_j)))*x1

        covmatrix=2*ncontshared*(cov00+cov11)

        inv_i = np.linalg.pinv(info_i)
        inv_j = np.linalg.pinv(info_j)

        covfinal=np.dot(np.dot(inv_i,covmatrix),inv_j)

        ### Correlation matrix for PASTRY
        cor=covfinal[1,1]/np.sqrt(np.dot(inv_i[1,1],inv_j[1,1]))
        pastry_R=np.tile(cor,(k,k))
        np.fill_diagonal(pastry_R, 1)

        ## new decoupling
        pastry_invsigma=np.dot(np.dot(invstdmat,np.linalg.inv(pastry_R)),invstdmat)
        pastry_invvar=np.dot(np.dot(unit,pastry_invsigma),unit_t)
        pastry_var=1/pastry_invvar
        beta_num=np.dot(np.dot(unit,pastry_invsigma),betas_nonde)	
        pastry_beta=beta_num/pastry_invvar
        pastry_stat=pastry_beta/np.sqrt(pastry_var)

        ## pvalue calculation by FE
        pval_split_FE[eachrep]=FE(betas_split, stder_split)['p']
        pval_PASTRY_FE[eachrep]=stats.norm.pdf(abs(pastry_stat))*2

    end=time.time()
    print(" * # Simulation : ",nrep," \n * Simulation Results : \n  LS : ",sum(pval_LS_FE<t)/nrep,", PASTRY",sum(pval_PASTRY_FE<t)/nrep,", split_FE",sum(pval_split_FE<t)/nrep,"\n")
    print("Done.\n")
    print(" * PASTRY Simulation running time in Python : "+str(end - start))
    return()

# PASTRY module run as the main program.
if __name__ == "__main__":

    print("\n")
    print("================ PASTRY_sim v 1.0 ================ ")
    print("===== v.1.0 Seoul National University Hanlab ===== \n")

    parser = argparse.ArgumentParser()
    parser.set_defaults(method = pastry_sim)

    parser.add_argument("-n","--N", type=int, default=5, help=":# of studies combined", metavar="integer")
    parser.add_argument("-r", "--RRs", type=float, default=1.3, help=":Risk Ratio (default : 1.3)")
    parser.add_argument("-m", "--maf", type=float, default=0.3, help=":maf files directory (casemaf.out, contmaf.out) ")      
    parser.add_argument("-p", "--nrep", type=int, default=5000, help=":# of Simulation (default : 5000)")
    parser.add_argument("-s", "--shared", type=int, default=1000, help=":# of shared control", metavar="integer")
    parser.add_argument("-l", "--cont", type=int, default=1000, help=":# of control (default : 1000)")
    parser.add_argument("-e","--case", type=int, default=1000, help=":# of case (default : 1000)", metavar="integer")
    parser.add_argument("-t","--thrs", type=float, default=5E-8, help=":# of control", metavar="integer")
    options = parser.parse_args()

    pastry_sim(RRs=options.RRs, maf=options.maf, nrep=options.nrep, k=options.N, ncase=options.case, ncont=options.cont, ncontshared=options.shared, t=options.thrs)
    print("PASTRY Simulation Done. ========================= ")
    print("================================================= \n ")
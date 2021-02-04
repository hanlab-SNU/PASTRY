#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, time
import numpy as np
import datatable as dt
import scipy.stats as stats
from numpy import savetxt

def pastry(indir, maf, outdir, k, ncase, ncont, ncontshared):

	# K : Combining Studies (# of files)
	# X : Genotype of K study
	# Y : Phenotype of k study
	# nsnp : # of SNP for one study

	# Importing input files and setting parameters
	print("Your input files count for meta-analysis : ",k)

	casemaf=dt.fread(maf+"casemaf.out").to_numpy()
	contmaf=dt.fread(maf+"contmaf.out").to_numpy()
	nsnp=len(casemaf[:,0])
	#maf=dt.fread(mafpath).to_pandas() # for (.frq.cc files)
	#casemaf=np.array(maf["MAF_A"]) # for (.frq.cc files)
	#contmaf=np.array(maf["MAF_U"]) # for (.frq.cc files)
	unit=np.repeat(1,k).reshape(1,k)
	unit_t=np.transpose(unit)

	nstudy=np.repeat(ncont+ncase,k)
	ncont=np.repeat(ncont,k)
	ncase=np.repeat(ncase,k)

	ls_stat=np.zeros(nsnp)
	pastry_stat=np.zeros(nsnp)

	betas=np.zeros((k,nsnp))
	ses=np.zeros((k,nsnp))

	for i in range(1,k+1):
		study=dt.fread(indir+str(i)+".out").to_numpy()
		betas[i-1]=study[:,0]
		ses[i-1]=study[:,1]

	# Na handling
	betas=np.nan_to_num(betas)
	ses=np.nan_to_num(ses)

	# LS
	## Generating Correlation Matrix
	ls_R=np.diagflat(np.repeat(1,k)).astype(float)
	pastry_R=np.diagflat(np.repeat(1,k)).astype(float)
	
	for i in range(k):
		for j in range(k):
			if i!=j:
				num=ncontshared*np.sqrt((ncase[i]*ncase[j])/(ncont[i]*ncont[j])) #numerator
				den=np.sqrt(nstudy[i]*nstudy[i]) #denominator
				ls_R[i,j]=num/den
	
	inv_ls_R=np.linalg.inv(ls_R)


	invstdmat={}

	ls_beta=np.zeros(nsnp)
	ls_stat=np.zeros(nsnp)

	pastry_beta=np.zeros(nsnp)
	pastry_stat=np.zeros(nsnp)

	## Generate Statistics
	for n in range(nsnp):
		p0cont=1-contmaf[n]
		p1cont=contmaf[n]

		p0case=1-casemaf[n]
		p1case=casemaf[n]

		invstdmat[n]=np.linalg.pinv(np.diagflat(ses[:,n]))
		ls_invsigma=np.dot(np.dot(invstdmat[n],inv_ls_R),invstdmat[n])
		ls_invvar=np.dot(np.dot(unit,ls_invsigma),unit_t)
		ls_beta_num=np.dot(np.dot(unit,ls_invsigma),betas[:,n])

		ls_beta[n]=ls_beta_num/ls_invvar
		ls_var=1/ls_invvar

		ls_stat[n]=ls_beta[n]/np.sqrt(ls_var)

	# PASTRY
	start = time.time()
	for n in range(nsnp):
		for i in range(k):
			for j in range(k):
				if i!=j:
					p0cont=1-contmaf[n]
					p1cont=contmaf[n]

					p0case=1-casemaf[n]
					p1case=casemaf[n]

					ntotali=ncont[i]+ncase[i]
					ntotalj=ncont[j]+ncase[j]

					exp_alpha_i = np.exp(np.log(ncase[i]/(ncont[i]))-(ls_beta[n]*(((ncase[i]*casemaf[n])+(ncont[i]*contmaf[n]))/ntotali)))
					exp_alpha_j = np.exp(np.log(ncase[j]/(ncont[j]))-(ls_beta[n]*(((ncase[j]*casemaf[n])+(ncont[j]*contmaf[n]))/ntotalj)))

					x0_i = exp_alpha_i
					x1_i = exp_alpha_i*np.exp(ls_beta[n])

					x0_j = exp_alpha_j
					x1_j = exp_alpha_j*np.exp(ls_beta[n])

					## Information Matrix
					x0=np.array([1,0,0,0]).reshape(2,2)
					x1=np.array([1,1,1,1]).reshape(2,2)

					up_covmat=np.array([0,0,0,0]).reshape(2,2)

					info_i=2*ncase[i]*(p0case*(x0_i/(1+x0_i)**2)*x0 + p1case*(x1_i/(1+x1_i)**2)*x1) + 2*(ncont[i])*(p0cont*(x0_i/(1+x0_i)**2)*x0 + p1cont*(x1_i/(1+x1_i)**2)*x1)
					info_j=2*ncase[j]*(p0case*(x0_j/(1+x0_j)**2)*x0 + p1case*(x1_j/(1+x1_j)**2)*x1) + 2*(ncont[j])*(p0cont*(x0_j/(1+x0_j)**2)*x0 + p1cont*(x1_j/(1+x1_j)**2)*x1)

					## Covarinace Matrix
					#cov00=p0cont*(x0_i*x0_j/(((1+x0_i)*(1+x0_j))**2))*x0
					#cov11=p1cont*(x1_i*x1_j/(((1+x1_i)*(1+x1_j))**2))*x1

					cov00=p0cont*(x0_i*x0_j/((1+x0_i)*(1+x0_j)))*x0
					cov11=p1cont*(x1_i*x1_j/((1+x1_i)*(1+x1_j)))*x1

					covmatrix=2*ncontshared*(cov00+cov11)

					inv_i = np.linalg.pinv(info_i)
					inv_j = np.linalg.pinv(info_j)

					covfinal=np.dot(np.dot(inv_i,covmatrix),inv_j)

					### Correlation matrix for PASTRY
					cor=covfinal[1,1]/np.sqrt(np.dot(inv_i[1,1],inv_j[1,1]))
					pastry_R[i,j]=cor
					pastry_R[j,i]=cor

		pastry_invsigma=np.dot(np.dot(invstdmat[n],np.linalg.inv(pastry_R)),invstdmat[n])
		pastry_invvar=np.dot(np.dot(unit,pastry_invsigma),unit_t)

		pastry_var=1/pastry_invvar
		beta_num=np.dot(np.dot(unit,pastry_invsigma),betas[:,n])

		pastry_beta[n]=beta_num/pastry_invvar
		pastry_stat[n]=pastry_beta[n]/np.sqrt(pastry_var)

	end = time.time()
	print("Your output is writing...\n")
	savetxt("{}_PASTRY_zstat.py.out".format(outdir),pastry_stat,fmt='%.3f',header="PASTRY_Z",delimiter="\n")
	savetxt("{}_LS_zstat.py.out".format(outdir),ls_stat,fmt='%.3f',header="LS_Z",delimiter="\n")
	#savetxt("{}_contmaf.py.out".format(outdir),contmaf,delimiter="\n") # for (.frq.cc files)
	#savetxt("{}_casemaf.py.out".format(outdir),casemaf,delimiter="\n") # for (.frq.cc files)
	print("Done.\n")
	print(" * PASTRY processed time : "+str(end - start))

	return()

# PASTRY module run as the main program.
if __name__ == "__main__":


	print("\n")
	print("================ PASTRY v 1.0 ================ \n")
	print("==== v.1.0 Seoul National University Hanlab ==== \n")
	print("\n")

	parser = argparse.ArgumentParser()
	parser.set_defaults(method = pastry)
	parser.add_argument("-i", "--input", type=str,help=":input files prefix name (.out)") # input files prefix
	parser.add_argument("-m", "--maf", type=str, default="./",help=":maf files directory (casemaf.out, contmaf.out) ") # input maf file directory path          
	parser.add_argument("-o", "--out", type=str, help=":output file prefix name (.out)") # output files prefix
	parser.add_argument("-n","--N", type=int,help=":# of studies combined", metavar="integer") # meta-analysis studies
	parser.add_argument("-s", "--shared", type=int,help=":# of shared control", metavar="integer")  # shared control N
	parser.add_argument("-e","--case", type=int,help=":# of case", metavar="integer")            # case N
	parser.add_argument("-l","--cont", type=int,help=":# of control", metavar="integer")          # control N
	options = parser.parse_args()

	if options.input is None:
		raise ValueError('input file is required. -h for help message.')
	if options.maf is None:
		raise ValueError('frequency file is required. -h for help message.')
	if options.out is None:
		raise ValueError('Result file prefix is required. -h for help message.')
	if options.N is None:
		raise ValueError('# N is required. -h for help message.')
	if options.shared is None:
		raise ValueError('# of shared is required. -h for help message.')
	if options.case is None:
		raise ValueError('# of control is required. -h for help message.')
	if options.cont is None:
		raise ValueError('# of control is required. -h for help message.')

	pastry(indir=options.input, maf=options.maf, outdir=options.out, k=options.N, ncase=options.case, ncont=options.cont, ncontshared=options.shared)
	print("PASTRY Analysis & Writing Results Done. ========= ")
	print("================================================= \n ")
PASTRY `v1.0`
===================

Introduction
------------
`PASTRY`  is a software package for avoiding **P**ower **AS**ymme**TRY** in meta-analysis.
To uncover the associated variants with small effects in Genome-wide association studies (GWAS), we need a large sample size to obtain statistical power. If multiple studies participating in a meta-analysis utilize the same public dataset as controls, the summary statistics from these studies become correlated and lead to increased false positives. Fortunately, [Lin and Sullivan (LS) method](https://www.ncbi.nlm.nih.gov/pubmed/20004761) proposed the correlation estimator based on the shared and unshared sample sizes, along with an optimal test statistic to account for the correlations. However, we identified that can lead to unbalanced power for detecting protective alleles (OR<1) and risk alleles (OR>1) and this power asymmetry problem occurred because the standard correlation estimator did not exactly predict the true correlation. We developed a method, called `PASTRY`, estimating the accurate correlation to overcome the power asymmetry problem.

Below we briefly describe short instructions for using the software.

#<!Diagram Figure>#

Instructions
-------------
### Downloading the package
`PASTRY` supports Fortran and R.

#### Using Fortran,
In order to download `PASTRY`, you can clone this repository via the commands. Also, when you use R, you can download it from Bioconductor.

```
$ git clone https://github.com/hanlab-SNU/PASTRY.git
$ cd PASTRY
```

#### Using R,
To install the bioconductor `PASTRY` package, type the following in an R command window:

```
## Install Bioconductor core packages
> if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
> BiocManager::install(c("PASTRY"))
> library(PASTRY)
```

### Installing required dependencies
If you are using Fortran, make sure [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries) is installed.
Also, you need [lapack](http://www.netlib.org/lapack/) package for Fortran 90.

```
# You can use 'yum, brew, dnf' instead of 'apt-get' according to each linux version / Mac OS
$ sudo apt-get install gcc gfortran libblas-dev liblapack-dev
$ gfortran --version
```

### Input data format
You need control, case frequency files(`casemaf.out`/`contmaf.out`) and easily processed n study files(`studyN.out`) for `PASTRY` input.

#### Making `.out` files
From the result from plink association studies(`.assos.logistic`), you can extract beta and standard error (7, 8 columns) from it.
```
# If you have n studies, make 'study1.out, study2.out,...,studyn.out'
$  awk '{print $7,$8}' study1name.assoc.logistic > study1.out
$  awk '{print $7,$8}' study2name.assoc.logistic > study2.out
....
$  awk '{print $7,$8}' studyNname.assoc.logistic > studyN.out
```

#### Making `.frq.cc` files
```
# You need plink frequency files for making maf files.
$ awk '{print $5}' studyname.frq.cc > casemaf.out
$ awk '{print $6}' studyname.frq.cc > contmaf.out
```

### Executing the codes

#### Using Fortran,
```
# You can give the path for n studies input (if it is not the same folder) and 
$ chmod +x run_pastry.sh
$ ./run_pastry.sh -n N 
  ( + optional arguments : -i "input path" / -o "output path and prefix" / -m "maf prefix" )
```

#### Using R,
```
# If all input files exists and you prefer to create the PASTRY output in R current working directory, 
# you don't need to give the arguments except  for the study number n.

> PASTRY(n=N, input="input_path", maf="", output="output_path&prefix")

or

> PASTRY(n=N) # All files in current working directory

```

### Output
You can get the z statistics for the LS and PASTRY method (`Zstat_LS_xxx.out` and `Zstat_PASTRY_xxx.out`).

License
---------
This project is licensed under the terms of the XXXX license.

Citation
----------
If you use the software `PASTRY`, please cite [Kim et al. Achieving balanced power for detecting risk and protective alleles in meta-analysis of association studies with overlapping subjects. (under review) (2019)](www.)

Reference
------------
1. [Lin and Sullivan (LS) method](https://www.ncbi.nlm.nih.gov/pubmed/20004761) | Lin,D.Y. and Sullivan,P.F., Meta-analysis of genome-wide associationstudies with overlapping subjects. Am. J. Hum. Genet., 85, 862–872. (2009) doi:10.1016/j.ajhg.2009.11.001.

Support
----------
This software was implemented by Emma E. Kim, Chloe Soohyun Jang and Buhm Han. Please contact [hanlab.snu@gmail.com](mailto:hanlab.snu@gmail.com)

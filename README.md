
This is the code computing the likelihoods of the SDSS III DR12 data using the redshift dependent Alcock-Paczynski method.

The code post-process the Planck MCMC chains and create the CMB+BAO+SNIa+H0+AP joint likelihood chains:
  /src/AP_MCMC
The created chains are used to make the likelihood contours

Some hints:

  1. 2-point correlation and covariance matrice files necessary for the likelihood analysis are provided
      2pCFs/
      covmats/
  
  2. You may want to modify the values of following variables
      filedir,  covmatdir,  chisqdir (AP_tools.f90)
      mcmcdir, MCMCfilestr, inputMCMCfile, outputMCMCfile (AP_MCMC.f90)

If any questions, feel free to contact Dr. Xiao-Dong Li: xiaodongli@kias.re.kr

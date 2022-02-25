# AH-BCS
R code for simulation experiments in the additive hazards model for bivariate current status data

Manuscript: "Sieve estimation of the additive hazards model with bivariate current status data" (Under Review)

Authors: Ce Zhang, Riyadh Rustam Al-Mosawi, Dipankar Bandyopadhyay, Haiwu Huang and Xuewen Lu

Packages used in these codes
-	numDeriv()     
This package is used to compute the numerical gradient the hessian matrix
-	Alabama()    
This package is used to implement constrOptim.nl() function for constraint optimization.
-	Copula()      
This function is used to compute and generate different copula models
User-defined function:
-	log.like() is the function of negative of  log-likelihood 
-	bern() is the function of computing the coefficients of Bernstein polynomial
-	hin.b1b2() is the monotone constraint function on b1 and b2
-	hin.phi() is the constraint function on phi
-	phi.val() is the function to evaluate the phi function
-	Gen.Data() the function for generating data of bivariate copula
-	Init.Val() function used to compute the initial values of b1,b2 and phi 
-	Prof.Inf()  function used  to compute the standard error using observed likelihood approach
-	Copula.fun(): This function to define copula name, copula function, copula generation and the association parameter. 
-	Types of copula that are used are Clayton, FGM, Frank, Gumbel, Gaussian with tau=0.25 & 0.5; for all of the types except tau=0.22 for FGM


Running the code: 

Step 1: Define the copula type, and the value of Kendall's measure "tau". 
For example, see lines 279-280 in code Biv-Sieve.R. Here, 

Copula.Name="clayton"   
Copula.tau=0.25        

implies, utilizing the Clayton copula, with tau=0.25. 

Step 2: Then, simply copy & paste the R code: Biv-Sieve.R

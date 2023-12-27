nR_S1 <- c(100,50,20,10,5,1)
nR_S2 <- c(3,7,8,12,27,89)
# SETUP-----------
nRatings <- length(nR_S1)/2 # number of ratings on confidence scale
# find actual type 2 FAR and HR (data to be fit)
  # I_nR and C_nR are rating trial counts for incorrect and correct trials
  # element i corresponds to # (in)correct w/ rating i
I_nR_rS2 = nR_S1[(nRatings+1):length(nR_S1)];  # each side of confidence scale
I_nR_rS1 = nR_S2[nRatings:1];
C_nR_rS2 = nR_S2[(nRatings+1):length(nR_S2)];
C_nR_rS1 = nR_S1[nRatings:1];

t2FAR_rS2 <- vector(length=length(I_nR_rS2)-1)
t2HR_rS2  <- vector(length=length(C_nR_rS2)-1)
t2FAR_rS1 <- vector(length=length(I_nR_rS1)-1)
t2HR_rS1  <- vector(length=length(C_nR_rS1)-1)

for(i in 2:nRatings){
  t2FAR_rS2[i-1] <- sum(I_nR_rS2[i:length(I_nR_rS2)])/sum(I_nR_rS2) #false alarm rate = for each stimulus: (incorrect count above each level of confidence) / n(incorrect ratings)
  t2HR_rS2[i-1]  <- sum(C_nR_rS2[i:length(C_nR_rS2)])/sum(C_nR_rS2)
  
  t2FAR_rS1[i-1] <- sum(I_nR_rS1[i:length(I_nR_rS1)])/sum(I_nR_rS1)
  t2HR_rS1[i-1]  <- sum(C_nR_rS1[i:length(C_nR_rS1)])/sum(C_nR_rS1)    
}

t2FAR_rS1 <- t2FAR_rS1[length(t2FAR_rS1):1]; # reverse rs1 FAR and HR so that...
t2HR_rS1 <- t2HR_rS1[length(t2HR_rS1):1];

# set up constraints
  #these are for feeding into the fmincon function https://search.r-project.org/CRAN/refmans/pracma/html/fmincon.html
  #which is an R implementation based on the MATLAB function of the same name: https://uk.mathworks.com/help/optim/ug/fmincon.html
  #fmincon() finds the minimum of a constrained nonlinear multivariable function
  #video tutorial: https://www.youtube.com/watch?v=_Il7GQdL3Sk ; written version: https://apmonitor.com/che263/index.php/Main/MatlabOptimization
nCriteria <- 2*nRatings-1 # number of confidence criteria > c
nCritPerResp <- (nCriteria-1)/2; # number of confidence criteria for each response

# A and b are 'linear ineqality constraints of the form A x <= b.'
  #	i.e. when finding a solution to the problem (fitting a distribution to the data), it will accept the solution that is closest to that constraint.
A <- matrix(0L, nrow=nCritPerResp-1, ncol=nCritPerResp+1) #	A is a matrix, run each row at a time, where each col is compared to the same index of b. 
b <- -.001 # B is -.001 followed by nCritPerResp-1 -1.0s if nCritPerResp>2
# This means that the distribution must start somewhere near 0 and the cumulative probability must increase for each level.
if(nCritPerResp>2){ b[2] <- -1.0 }
offset <- 1;
for(crit in 1:(nCritPerResp-1)){
  A[crit,c(offset+crit,offset+crit+1)] <- c(1,-1); #
  if(nCritPerResp>2){ b[crit+1] <- -1.0 } #
}

LB_rS1 <- -10 # meta-d' rS1
LB_rS1[2:(nCritPerResp+1)] <- rep(-20,nCritPerResp) # criteria lower than t1c #check shouldn't be 1s
UB_rS1 <- 10
UB_rS1[2:(nCritPerResp+1)] <- rep(0,nCritPerResp)
LB_rS2 <- -10
LB_rS2[2:(nCritPerResp+1)] <- rep(0,nCritPerResp)
UB_rS2 <- 10
UB_rS2[2:(nCritPerResp+1)] <- rep(20,nCritPerResp)

# set up guess
ratingHR  <- vector(length=(nRatings*2)-1)
ratingFAR <- vector(length=(nRatings*2)-1)

# (reversed) cumulative probabilities: Proportion of responses at that confidence rating or higher (i.e. high confidence responses) for each stim presentation
for(c in 2:(nRatings*2)){ # starting at the second one, the first would be = 1, will be d' later, calculated below
  ratingHR[c-1] <- sum(nR_S2[c:length(nR_S2)]) / sum(nR_S2) #ratingHR is the nR_S2 responses above each confidence rating
  ratingFAR[c-1] <- sum(nR_S1[c:length(nR_S1)]) / sum(nR_S1) #NOTE: NR_S1 and S2 are ordered in the same direction i.e. [(r=1,c=3),(r=1,c=2)...]
}

t1_index <- nRatings # as ratingHR uses nR_S2[2:length]: for each stim presentation, this is the first S1 response
t2_index <- setdiff(1:(2*nRatings-1), t1_index) # indexes for each response type above the first one

##CONFUSED ABOUT CALCULATION OF D' HERE
# qnorm gets the inverse Cumulative Distribution Function of the Normal Distribution (iCDFN)
  # ratingXX[t1_index] gives proportion of responses for nR_S1 and 2 where r=S1
  # consdier the ordering and length of the raw data: this is also the centre of the distribution
  # i.e. the proportion of ratings that are correct.
s <-1 # ratio of S1 and S2 SDs, recommended to set to 1
d1 <- (1/s)*qnorm(ratingHR[t1_index])-qnorm(ratingFAR[t1_index]);
meta_d1_rS1 <- d1;
meta_d1_rS2 <- d1;

c1 <- (-1/(1+s))*(qnorm(ratingHR)+qnorm(ratingFAR));
t1c1 <- c1[t1_index];
t2c1 <- c1[t2_index];

constant_criterion_rS1 <- 'meta_d1_rS1 * (t1c1 / d1)'
constant_criterion_rS2 <- 'meta_d1_rS2 * (t1c1 / d1)'
guess_rS1 <- c(meta_d1_rS1, t2c1[1:nCritPerResp]-eval(parse(text=constant_criterion_rS1)));
guess_rS2 <- c(meta_d1_rS2, t2c1[(nCritPerResp+1):length(t2c1)]- eval(parse(text=constant_criterion_rS1)));

# fit for S1 responses ------------------
#install.packages("pracma") install.packages("NlcOptim")
library(NlcOptim)
library(pracma) #https://search.r-project.org/CRAN/refmans/pracma/html/fmincon.html

fitM_rS1_logL <- function(input){
  # set up parameters
  meta_d1_rS1 <- input[1]; 
  t2c1        <- input[2:length(input)];
  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS1/2; S1sd = 1
  S2mu <- meta_d1_rS1/2; S2sd = S1sd/s
  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd)
  f <- 1-pnorm(0,S1mu,S1sd)
 
  shift_c1 <- (1/(1+s)) * (qnorm(h)+qnorm(f)) #this is the value of c1 midway b/t S1 and S2
  S1mu <- S1mu + shift_c1 #shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1
  # adjust so that t1c1 = 0
  S1mu = S1mu - eval(parse(text=constant_criterion_rS1))
  S2mu = S2mu - eval(parse(text=constant_criterion_rS1));
  t1c1 = 0;
  ## set up MLE analysis
  nC_rS1 <- vector(length=nRatings)
  nI_rS1 <- vector(length=nRatings)
  for(ic in 1:nRatings){ # get type 2 response counts
    nC_rS1[ic] <- nR_S1[ic] # S1 responses
    nI_rS1[ic] <- nR_S2[ic]
  }
  # get type 2 probabilities
  C_area_rS1 <- pnorm(t1c1,S1mu,S1sd);
  I_area_rS1 <- pnorm(t1c1,S2mu,S2sd);
  t2c1x <- c(-Inf,t2c1,t1c1)
  
  prC_rS1 <- vector(length=nRatings)
  prI_rS1 <- vector(length=nRatings)
  for(ic in 1:nRatings){
    prC_rS1[ic] <- (pnorm(t2c1x[ic+1],S1mu,S1sd) - pnorm(t2c1x[ic],S1mu,S1sd)) /C_area_rS1;
    prI_rS1[ic] <- (pnorm(t2c1x[ic+1],S2mu,S2sd) - pnorm(t2c1x[ic],S2mu,S2sd)) /I_area_rS1;   
  }
  
  logL = 0;
  for(ic in 1:nRatings){
    logL <- logL + nC_rS1[ic]*log(prC_rS1[ic]) + nI_rS1[ic]*log(prI_rS1[ic])
  }
  if(is.nan(logL)){logL=-Inf}
  logL <- -logL
  return(logL)
}

minconNMF <- fmincon(x0=guess_rS1,fn=fitM_rS1_logL,A=A,b=b,lb=LB_rS1,ub=UB_rS1,maxfeval=100000)

x <- minconNMF$par
f <- minconNMF$value

meta_d1_rS1 <- x[1];
meta_c1_rS1 <- eval(parse(text=constant_criterion_rS1));
t2c1_rS1  <- x[2:length(x)] + eval(parse(text=constant_criterion_rS1));
logL_rS1  <- -f;

# find model-estimated t2FAR and t2HR for S1 responses---------------------
# find the estimated type 2 FAR and HR
S1mu <- -meta_d1_rS1/2; S1sd = 1
S2mu <-  meta_d1_rS1/2; S2sd = S1sd/s 
# adjust so that everything is centered on t1c1=0
h <- 1-pnorm(0,S2mu,S2sd);
f <- 1-pnorm(0,S1mu,S1sd);
shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)) # this is the value of c1 midway b/t S1 and S2
S1mu <- S1mu + shift_c1 # shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
S2mu <- S2mu + shift_c1
C_area_rS1 <- pnorm(meta_c1_rS1,S1mu,S1sd);
I_area_rS1 <- pnorm(meta_c1_rS1,S2mu,S2sd);
est_t2FAR_rS1 <- vector(length=length(t2c1_rS1))
est_t2HR_rS1 <-  vector(length=length(t2c1_rS1))
for(i in 1:length(t2c1_rS1)){
  t2c1_lower <- t2c1_rS1[i]
  I_FAR_area_rS1 <- pnorm(t2c1_lower,S2mu,S2sd)
  C_HR_area_rS1  <- pnorm(t2c1_lower,S1mu,S1sd)
  est_t2FAR_rS1[i] <- I_FAR_area_rS1/I_area_rS1
  est_t2HR_rS1[i]  <- C_HR_area_rS1/C_area_rS1
}

# fit for S2 responses------------------
fitM_rS2_logL <- function(input){
  # set up parameters
  meta_d1_rS2 <- input[1]
  t2c1        <- input[2:length(input)]
  # find the estimated type 2 FAR and HR
  S1mu <- -meta_d1_rS2/2; S1sd = 1
  S2mu <- meta_d1_rS2/2; S2sd = S1sd/s
  # adjust so that everything is centered on t1c1=0
  h <- 1-pnorm(0,S2mu,S2sd)
  f <- 1-pnorm(0,S1mu,S1sd)

  shift_c1 <- (1/(1+s)) * (qnorm(h)+qnorm(f)) #this is the value of c1 midway b/t S1 and S2
  S1mu <- S1mu + shift_c1 #shift S1 and S2mu so that they lie on an axis for 0 --> c1=0
  S2mu <- S2mu + shift_c1

  # adjust so that t1c1 = 0
  S1mu = S1mu - eval(parse(text=constant_criterion_rS2))
  S2mu = S2mu - eval(parse(text=constant_criterion_rS2));
  t1c1 = 0;
  ## set up MLE analysis
  nC_rS2 <- vector(length=nRatings)
  nI_rS2 <- vector(length=nRatings)
  for(ic in 1:nRatings){ # get type 2 response counts
    nC_rS2[ic] <- nR_S2[ic+nRatings] # S1 responses
    nI_rS2[ic] <- nR_S1[ic+nRatings]
  }
  # get type 2 probabilities
  C_area_rS2 <- 1-pnorm(t1c1,S2mu,S2sd);
  I_area_rS2 <- 1-pnorm(t1c1,S1mu,S1sd);
  t2c1x <- c(t1c1,t2c1,Inf)
  
  prC_rS2 <- vector(length=nRatings)
  prI_rS2 <- vector(length=nRatings)
  for(ic in 1:nRatings){
    prC_rS2[ic] = ((1-pnorm(t2c1x[ic],S2mu,S2sd)) - (1-pnorm(t2c1x[ic+1],S2mu,S2sd))) /C_area_rS2;
    prI_rS2[ic] = ((1-pnorm(t2c1x[ic],S1mu,S1sd)) - (1-pnorm(t2c1x[ic+1],S1mu,S1sd))) /I_area_rS2;   
  }
  
  logL = 0;
  for(ic in 1:nRatings){
    logL <- logL + nC_rS2[ic]*log(prC_rS2[ic]) + nI_rS2[ic]*log(prI_rS2[ic])
  }
  if(is.nan(logL)){logL=-Inf}
  logL <- -logL
  return(logL)
}

minconNMF <- fmincon(x0=guess_rS2,fn=fitM_rS2_logL,A=A,b=b,lb=LB_rS2,ub=UB_rS2,maxfeval=100000)

x <- minconNMF$par
f <- minconNMF$value

meta_d1_rS2 <- x[1];
meta_c1_rS2 <- eval(parse(text=constant_criterion_rS2));
t2c1_rS2  <- x[2:length(x)] + eval(parse(text=constant_criterion_rS2));
logL_rS2  <- -f;



# find estimated t2FAR and t2HR ------------
S1mu <- -meta_d1_rS2/2; S1sd <- 1;
S2mu <-  meta_d1_rS2/2; S2sd <- S1sd/s;
#adjust so that everything is centered on t1c1<-0
h <- 1-pnorm(0,S2mu,S2sd);
f <- 1-pnorm(0,S1mu,S1sd);
shift_c1 <- (-1/(1+s)) * (qnorm(h)+qnorm(f)); # this is the value of c1 midway b/t S1 and S2
S1mu <- S1mu + shift_c1; # shift S1 and S2mu so that they lie on an axis for 0 --> c1<-0
S2mu <- S2mu + shift_c1;
C_area_rS2 <- 1-pnorm(meta_c1_rS2,S2mu,S2sd);
I_area_rS2 <- 1-pnorm(meta_c1_rS2,S1mu,S1sd);
est_t2FAR_rS2 <- vector(length=length(t2c1_rS2))
est_t2HR_rS2 <-  vector(length=length(t2c1_rS2))
for(i in 1:length(t2c1_rS2)){
  t2c1_upper <- t2c1_rS2[i];
  I_FAR_area_rS2 <- 1-pnorm(t2c1_upper,S1mu,S1sd);
  C_HR_area_rS2  <- 1-pnorm(t2c1_upper,S2mu,S2sd);
  est_t2FAR_rS2[i] <- I_FAR_area_rS2 / I_area_rS2;
  est_t2HR_rS2[i] <- C_HR_area_rS2 / C_area_rS2;
}

# Output------
SDT_s_convert <- function(in1,s,cont1,cont2){
  # out = SDT_s_convert(in1,s,cont1,cont2)
  # in1 - input var to be converted
  # s - sd(S1)/sd(S2)
  # cont1 - identity of in1
  # cont2 - identity of out
  # cont1 and cont2 can be these tokens
  # 'da','d1','d2','ca','c1',c2'
  
  # convert d'
  # s = d2 / d1
  # da = sqrt(2./(1+s.^2)) .* d2
  if(cont1=='da'){
    da = in1;
    d2 = da / sqrt(2/(1+s^2));
    d1 = d2 / s;
  } else if(cont1=='d1'){
    d1 = in1;
    d2 = d1*s;
    da = sqrt(2/(1+s^2))*d2;
  } else if(cont1=='d2'){
    d2 = in1;
    d1 = d2 / s;
    da = sqrt(2/(1+s^2)) * d2;
  } else if(cont1=='ca'){
    ca = in1;
    c1 = ( sqrt(1+s^2) / sqrt(2)*s ) * ca;
    c2 = c1 * s;
  } else if(cont1=='c1'){
    c1 = in1;
    c2 = c1 * s;
    ca = ( sqrt(2)*s / sqrt(1+s^2) ) * c1;
  } else if(cont1=='c2'){
    c2 = in1;
    c1 = c2 / s;
    ca = ( sqrt(2)*s / sqrt(1+s^2) ) * c1;
  }
    # convert c
    # s = c2 / c1
    # ca = ( sqrt(2).*s ./ sqrt(1+s.^2) ) .* c1;

  out <- eval(parse(text=cont2))
  return(out)
}

# type 1 params
fit.da         = SDT_s_convert(d1,   s,'d1','da');
fit.t1ca       = SDT_s_convert(t1c1, s,'c1','ca');
fit.s          = s;

# type 2 fits for rS1
fit.meta_da_rS1   = SDT_s_convert(meta_d1_rS1, s,'d1','da');
fit.t1ca_rS1      = SDT_s_convert(meta_c1_rS1, s,'c1','ca');
fit.t2ca_rS1      = SDT_s_convert(t2c1_rS1,    s,'c1','ca');
fit.M_ratio_rS1   = fit.meta_da_rS1 / fit.da;
fit.M_diff_rS1    = fit.meta_da_rS1 - fit.da;
fit.logL_rS1      = logL_rS1;
fit.obs_HR2_rS1   = t2HR_rS1;
fit.est_HR2_rS1   = est_t2HR_rS1;
fit.obs_FAR2_rS1  = t2FAR_rS1;
fit.est_FAR2_rS1  = est_t2FAR_rS1;

# type 2 fits for rS2
fit.meta_da_rS2   = SDT_s_convert(meta_d1_rS2, s,'d1','da');
fit.t1ca_rS2      = SDT_s_convert(meta_c1_rS2, s,'c1','ca');
fit.t2ca_rS2      = SDT_s_convert(t2c1_rS2,    s,'c1','ca');
fit.M_ratio_rS2   = fit.meta_da_rS2 / fit.da;
fit.M_diff_rS2    = fit.meta_da_rS2 - fit.da;
fit.logL_rS2      = logL_rS2;
fit.obs_HR2_rS2   = t2HR_rS2;
fit.est_HR2_rS2   = est_t2HR_rS2;
fit.obs_FAR2_rS2  = t2FAR_rS2;
fit.est_FAR2_rS2  = est_t2FAR_rS2;

# S1 units
fit.S1units.d1   = d1;
fit.S1units.t1c1 = t1c1;
fit.S1units.meta_d1_rS1 = meta_d1_rS1;
fit.S1units.t1c1_rS1    = meta_c1_rS1;
fit.S1units.t2c1_rS1    = t2c1_rS1;
fit.S1units.meta_d1_rS2 = meta_d1_rS2;
fit.S1units.t1c1_rS2    = meta_c1_rS2;
fit.S1units.t2c1_rS2    = t2c1_rS2;

# Breath data version ------
nR_B9 <- c(100,50) #FAs and Hits


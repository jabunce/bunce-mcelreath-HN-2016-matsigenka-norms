#Fits all 19 IRT models and constructs appendix Table B1 with
#parameter posterior estimates and model WAIC weights.
#This can take a few hours to run.

#Before beginning, download the data file Manu_interviews_31oct16.csv into
#a local directory on your machine. I'd recommend creating a new folder and
#putting the data file in there. Then open R from within that local directory
#and run this script. The figures will be created as pdfs in that same
#directory.

#On some machines, copying and pasting the entire script to run
#in the terminal results in errors. If you get errors, try running
#the individual models one at a time.

rm (list = ls(all=TRUE))
library(rstan)
library(rethinking)
library(stargazer)


rstan_options(auto_write = TRUE) #to let stan make a copy of the model and use multiple cores
options(mc.cores = parallel::detectCores())

#Read the data from the csv data file into R:
Interviews.raw <- read.csv(
	file=
	"./Manu_interviews_31oct16.csv", 
	header=TRUE)

#Check the variable names and dimensions in the data frame Interviews.raw
names(Interviews.raw)
dim(Interviews.raw)


#working dataframe
d <- Interviews.raw
dim(d)
names(d)


for ( i in 1:length(d$Question) ) { #flip coding of some questions to polarize latent axis as in text
  if ( d$Question[i] == 6 || 
       d$Question[i] == 17 ||
       d$Question[i] == 19 ||
       d$Question[i] == 22 ||
       d$Question[i] == 29 ||
       d$Question[i] == 30 ) {

      if (d$Response[i]==1) {
        d$Response[i] <- 0
      } else if (d$Response[i]==0) {
        d$Response[i] <- 1
      } #else

  } #if
} #for i


#replace Person ID with a new consecutive ID number
d$newID <- rep(0, nrow(d)) #initialize newID column of d
for ( i in 1:length(d$ID) ) {
	for ( m in 1:length(unique(d$ID)) ) {
		if ( d$ID[i] == unique(d$ID)[m] ) {
			d$newID[i] <- m
		} #if	
	} #m loop
} #i loop

#replace Question with a new consecutive QID number
d$QID <- rep(0, nrow(d)) #initialize QID column of d
for ( j in 1:length(d$Question) ) {
  for ( k in 1:length(unique(d$Question)) ) {
    if ( d$Question[j] == unique(d$Question)[k] ) {
      d$QID[j] <- k
    } #if 
  } #k loop
} #j loop


d <- d[order(d$newID),] #reorder by consecutive newID
dim(d)
names(d)


J <- length(unique(d$newID))	#number of people
K <- length(unique(d$QID))		#number of questions
N <- nrow(d)			#total number of responses
jj <- d$newID			#vector of person IDs
kk <- d$QID			#vector of question IDs
y <- d$Response			#vector of responses


quest_names <- c("1.wife hunts/works", "2.daughter babysits", "3.not wear dead hat",
               "4.hit students", "5.no questions", "6.post flu",
               "7.pot to needy", "8.good nonbaptized no heaven", "9.bad baptized heaven",
               "10.stop work to visit", "11.cheap mean store", "12.xcousin marriage",
               "13.arranged marriage", "14.laborer party")


#WAIC function
waic <- function(stanfit){
  log_lik <- extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
    p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
    pointwise=pointwise, total=total, se=se))
}

#colVars function needed by waic function
colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
                     sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}


#######base model: m1, rand effect for person ############################

model_code_1 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J];	// location of people (differing from the mean), i.e., random effect for person

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1);
  b1 ~ normal(0,1); //identifying prior for location and scale

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] );	//random effect for person

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  } //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] );

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
				            gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_1 <- list(
	J = length(unique(d$newID)),	#number of people
	K = length(unique(d$Question)),		#number of questions
	N = nrow(d),			#total number of responses
	jj = d$newID,			#vector of person IDs
	kk = d$QID,			#vector of question IDs
	y = d$Response			#vector of responses
)

start_1 <- list(
  b0=0, b1 = as.array(rep(0, times=J)), 
  beta = as.array(rep(0, times=K)),
  gamma= as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)

set.seed(1)
m1 <- stan( model_code=model_code_1, data=data_list_1,
            init=list(start_1, start_1, start_1, start_1), 
            	          iter=4000 , chains=4, 
                        control=list(adapt_delta=0.99) )

print(m1, pars=c("b0","b1","beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post1 <- extract.samples( m1 )
str(post1)


#look at all traces
pdf(file="./traces_m1.pdf", 
  height=3, width=8)
par(mfrow=c(2,1))

traceplot(m1, pars="b0", inc_warmup=T)
for ( z1 in 1:J ){
  print(traceplot(m1, pars=paste("b1[", z1, "]", sep=""), inc_warmup=T ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m1, pars=paste("beta[", z2, "]", sep=""), inc_warmup=T ))
  }
for ( z3 in 1:K ){
  print(traceplot(m1, pars=paste("gamma[", z3, "]", sep=""), inc_warmup=T ))
}
traceplot(m1, pars="sigma_beta", inc_warmup=T)
traceplot(m1, pars="sigma_gamma", inc_warmup=T)
graphics.off()


#######base model: m2, rand effect for person, Ethnicity predictor ############################

model_code_2 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained close to 0
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  } //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_2 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi   #vector of ethnicity codes
)

start_2 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m2 <- stan( model_code=model_code_2, data=data_list_2,
            init=list(start_2, start_2, start_2, start_2), 
                        iter=4000 , chains=4, 
                        control=list(adapt_delta=0.99) )

print(m2, pars=c("b0","b1","bMachi","beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post2 <- extract.samples( m2 )
str(post2)


#look at all traces
pdf(file="./traces_m2.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m2, pars="b0")
traceplot(m2, pars="bMachi")
for ( z1 in 1:J ){
  print(traceplot(m2, pars=paste("b1[", z1, "]", sep="") ))
  } 
for ( z2 in 1:(K) ){
  print(traceplot(m2, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m2, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m2, pars="sigma_beta")
traceplot(m2, pars="sigma_gamma")
graphics.off()



#######base model: m3, Ethnicity, Sex, Eth x Sex ############################

model_code_3 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> Male[N]; // predictor for sex
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}

model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained close to 0
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + bMale*Male[n] +
                                bMaleMach*Machi[n]*Male[n] );

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n]+ bMale*Male[n] +
                                bMaleMach*Machi[n]*Male[n] );

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_3 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  Male = d$Sex.1male    #vector of sex codes
)

start_3 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, bMale=0, bMaleMach=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m3 <- stan( model_code=model_code_3, data=data_list_3,
            init=list(start_3, start_3, start_3, start_3), 
                        iter=4000 , chains=4, 
                        control=list(adapt_delta=0.99) )

print(m3, pars=c("b0","b1","bMachi","bMale","bMaleMach", 
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post3 <- extract.samples( m3 )
str(post3)


#look at all traces
pdf(file="./traces_m3.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m3, pars="b0")
traceplot(m3, pars="bMachi")
traceplot(m3, pars="bMale")
traceplot(m3, pars="bMaleMach")
for ( z1 in 1:J ){
  print(traceplot(m3, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m3, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K){
  print(traceplot(m3, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m3, pars="sigma_beta")
traceplot(m3, pars="sigma_gamma")
graphics.off()

 

#######base model: m4, Ethnicity, Age ############################

model_code_4 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained positive
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_4 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50
)

start_4 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m4 <- stan( model_code=model_code_4, data=data_list_4,
            init=list(start_4, start_4, start_4, start_4), 
                        iter=5000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m4, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post4 <- extract.samples( m4 )
str(post4)


#look at all traces
pdf(file="./traces_m4.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m4, pars="b0")
traceplot(m4, pars="bMachi")
traceplot(m4, pars="badol")
traceplot(m4, pars="bmat")
traceplot(m4, pars="bold")
for ( z1 in 1:J ){
  print(traceplot(m4, pars=paste("b1[", z1, "]", sep="") ))
  } 
for ( z2 in 1:(K) ){
  print(traceplot(m4, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m4, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m4, pars="sigma_beta")
traceplot(m4, pars="sigma_gamma")
graphics.off()



#######base model: m5, Ethnicity, Age, Sex, Sex x Eth ############################

model_code_5 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1);
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_5 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male    #vector of sex codes
)

start_5 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m5 <- stan( model_code=model_code_5, data=data_list_5,
            init=list(start_5, start_5, start_5, start_5), 
                        iter=4000 , chains=4,
                        control=list(adapt_delta=0.99) )

print(m5, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post5 <- extract.samples( m5 )
str(post5)


#look at all traces
pdf(file="./traces_m5.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m5, pars="b0")
traceplot(m5, pars="bMachi")
traceplot(m5, pars="badol")
traceplot(m5, pars="bmat")
traceplot(m5, pars="bold")
traceplot(m5, pars="bMale")
traceplot(m5, pars="bMaleMach")
for ( z1 in 1:J ){
  print(traceplot(m5, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m5, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m5, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m5, pars="sigma_beta")
traceplot(m5, pars="sigma_gamma")
graphics.off()




#######base model: m6, Ethnicity, Age, Sex, Sex x Eth, School.mest x Machi ############################

model_code_6 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> School[N]; // school experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bSchool;
  real bMachiSchool;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained positive
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_6 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  School = d$School.mest
)

start_6 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bSchool=0, bMachiSchool=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m6 <- stan( model_code=model_code_6, data=data_list_6,
            init=list(start_6, start_6, start_6, start_6), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m6, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bSchool", "bMachiSchool",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post6 <- extract.samples( m6 )
str(post6)


#look at all traces
pdf(file="./traces_m6.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m6, pars="b0")
traceplot(m6, pars="bMachi")
traceplot(m6, pars="badol")
traceplot(m6, pars="bmat")
traceplot(m6, pars="bold")
traceplot(m6, pars="bMale")
traceplot(m6, pars="bMaleMach")
traceplot(m6, pars="bSchool")
traceplot(m6, pars="bMachiSchool")
for ( z1 in 1:J ){
  print(traceplot(m6, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m6, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m6, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m6, pars="sigma_beta")
traceplot(m6, pars="sigma_gamma")
graphics.off()



#######base model: m7, Ethnicity, Age, Sex, Sex x Eth, ExEmpMest x Machi ############################

model_code_7 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> EmpMest[N]; // work experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bEmpMest;
  real bMachiEmpMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_7 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  EmpMest = d$ExEmpMest
)

start_7 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bEmpMest=0, bMachiEmpMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m7 <- stan( model_code=model_code_7, data=data_list_7,
            init=list(start_7, start_7, start_7, start_7), 
                        iter=5000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m7, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bEmpMest", "bMachiEmpMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post7 <- extract.samples( m7 )
str(post7)


#look at all traces
pdf(file="./traces_m7.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m7, pars="b0")
traceplot(m7, pars="bMachi")
traceplot(m7, pars="badol")
traceplot(m7, pars="bmat")
traceplot(m7, pars="bold")
traceplot(m7, pars="bMale")
traceplot(m7, pars="bMaleMach")
traceplot(m7, pars="bEmpMest")
traceplot(m7, pars="bMachiEmpMest")
for ( z1 in 1:J ){
  print(traceplot(m7, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m7, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m7, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m7, pars="sigma_beta")
traceplot(m7, pars="sigma_gamma")
graphics.off()


#######base model: m8, Ethnicity, Age, Sex, Sex x Eth, ExComMest x Machi ############################

model_code_8 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_8 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  ComMest = d$ExComMest
)

start_8 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m8 <- stan( model_code=model_code_8, data=data_list_8,
            init=list(start_8, start_8, start_8, start_8), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m8, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bComMest", "bMachiComMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post8 <- extract.samples( m8 )
str(post8)


#look at all traces
pdf(file="./traces_m8.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m8, pars="b0")
traceplot(m8, pars="bMachi")
traceplot(m8, pars="badol")
traceplot(m8, pars="bmat")
traceplot(m8, pars="bold")
traceplot(m8, pars="bMale")
traceplot(m8, pars="bMaleMach")
traceplot(m8, pars="bComMest")
traceplot(m8, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m8, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m8, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m8, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m8, pars="sigma_beta")
traceplot(m8, pars="sigma_gamma")
graphics.off()


#######base model: m9, Eth, Age, Sex, Sex x Eth, School.mest x Machi, ExEmpMest x Machi ############################

model_code_9 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> School[N]; // school experience with mestizos
  int<lower=0,upper=1> EmpMest[N]; // work experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bSchool;
  real bMachiSchool;
  real bEmpMest;
  real bMachiEmpMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained positive
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_9 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  School = d$School.mest,
  EmpMest = d$ExEmpMest
)

start_9 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bSchool=0, bMachiSchool=0,
  bEmpMest=0, bMachiEmpMest=0,
  beta = as.array(rep(0, times=K)),
  gamma= as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m9 <- stan( model_code=model_code_9, data=data_list_9,
            init=list(start_9, start_9, start_9, start_9), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m9, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bSchool", "bMachiSchool",
                "bEmpMest", "bMachiEmpMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post9 <- extract.samples( m9 )
str(post9)


#look at all traces
pdf(file="./traces_m9.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m9, pars="b0")
traceplot(m9, pars="bMachi")
traceplot(m9, pars="badol")
traceplot(m9, pars="bmat")
traceplot(m9, pars="bold")
traceplot(m9, pars="bMale")
traceplot(m9, pars="bMaleMach")
traceplot(m9, pars="bSchool")
traceplot(m9, pars="bMachiSchool")
traceplot(m9, pars="bEmpMest")
traceplot(m9, pars="bMachiEmpMest")
for ( z1 in 1:J ){
  print(traceplot(m9, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m9, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m9, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m9, pars="sigma_beta")
traceplot(m9, pars="sigma_gamma")
graphics.off()



#######base model: m10, Eth, Age, Sex, Sex x Eth, School.mest x Machi, ExComMest x Machi ############################

model_code_10 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> School[N]; // school experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bSchool;
  real bMachiSchool;
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_10 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  School = d$School.mest,
  ComMest = d$ExComMest
)

start_10 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bSchool=0, bMachiSchool=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m10 <- stan( model_code=model_code_10, data=data_list_10,
            init=list(start_10, start_10, start_10, start_10), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m10, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bSchool", "bMachiSchool",
                "bComMest", "bMachiComMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post10 <- extract.samples( m10 )
str(post10)


#look at all traces
pdf(file="./traces_m10.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m10, pars="b0")
traceplot(m10, pars="bMachi")
traceplot(m10, pars="badol")
traceplot(m10, pars="bmat")
traceplot(m10, pars="bold")
traceplot(m10, pars="bMale")
traceplot(m10, pars="bMaleMach")
traceplot(m10, pars="bSchool")
traceplot(m10, pars="bMachiSchool")
traceplot(m10, pars="bComMest")
traceplot(m10, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m10, pars=paste("b1[", z1, "]", sep="") ))
  } 
for ( z2 in 1:(K) ){
  print(traceplot(m10, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m10, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m10, pars="sigma_beta")
traceplot(m10, pars="sigma_gamma")
graphics.off()



#######base model: m11, Eth, Age, Sex, Sex x Eth, ExEmpMest x Machi, ExComMest x Machi ############################

model_code_11 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bEmpMest;
  real bMachiEmpMest;
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained positive
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_11 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  EmpMest = d$ExEmpMest,
  ComMest = d$ExComMest
)

start_11 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bEmpMest=0, bMachiEmpMest=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m11 <- stan( model_code=model_code_11, data=data_list_11,
            init=list(start_11, start_11, start_11, start_11), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m11, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bEmpMest", "bMachiEmpMest",
                "bComMest", "bMachiComMest",
                "beta","gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post11 <- extract.samples( m11 )
str(post11)


#look at all traces
pdf(file="./traces_m11.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m11, pars="b0")
traceplot(m11, pars="bMachi")
traceplot(m11, pars="badol")
traceplot(m11, pars="bmat")
traceplot(m11, pars="bold")
traceplot(m11, pars="bMale")
traceplot(m11, pars="bMaleMach")
traceplot(m11, pars="bEmpMest")
traceplot(m11, pars="bMachiEmpMest")
traceplot(m11, pars="bComMest")
traceplot(m11, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m11, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m11, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m11, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m11, pars="sigma_beta")
traceplot(m11, pars="sigma_gamma")
graphics.off()



#######base model: m12, Eth, Age, Sex, Sex x Eth, ExEmpMest x Machi, ExComMest x Machi,
############################ School.mest x Machi


model_code_12 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> adol[N]; // adolescent
  int<lower=0,upper=1> mat[N]; // mature
  int<lower=0,upper=1> old[N]; // old
  int<lower=0,upper=1> Male[N]; // predictor for sex
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
  int<lower=0,upper=1> School[N]; // school experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real badol;
  real bmat;
  real bold;
  real bMale; // effect of maleness on person location in latent space
  real bMaleMach; //interaction effect
  real bEmpMest;
  real bMachiEmpMest;
  real bComMest;
  real bMachiComMest;
  real bSchool;
  real bMachiSchool;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}

model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); //constrained positive
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  badol ~ normal(0,1);
  bmat ~ normal(0,1);
  bold ~ normal(0,1);
  bMale ~ normal(0,1);
  bMaleMach ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] );

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] );

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_12 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  adol = d$adol.less20,
  mat = d$mat.20to50,
  old = d$old.over50,
  Male = d$Sex.1male,    #vector of sex codes
  EmpMest = d$ExEmpMest,
  ComMest = d$ExComMest,
  School = d$School.mest
)

start_12 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0, badol=0, bmat=0, bold=0,
  bMale=0, bMaleMach=0,
  bEmpMest=0, bMachiEmpMest=0,
  bComMest=0, bMachiComMest=0,
  bSchool=0, bMachiSchool=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m12 <- stan( model_code=model_code_12, data=data_list_12,
            init=list(start_12, start_12, start_12, start_12), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m12, pars=c("b0","b1","bMachi",
                "badol", "bmat", "bold",
                "bMale", "bMaleMach",
                "bEmpMest", "bMachiEmpMest",
                "bComMest", "bMachiComMest",
                "bSchool", "bMachiSchool",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post12 <- extract.samples( m12 )
str(post12)


#look at all traces
pdf(file="./traces_m12.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m12, pars="b0")
traceplot(m12, pars="bMachi")
traceplot(m12, pars="badol")
traceplot(m12, pars="bmat")
traceplot(m12, pars="bold")
traceplot(m12, pars="bMale")
traceplot(m12, pars="bMaleMach")
traceplot(m12, pars="bEmpMest")
traceplot(m12, pars="bMachiEmpMest")
traceplot(m12, pars="bComMest")
traceplot(m12, pars="bMachiComMest")
traceplot(m12, pars="bSchool")
traceplot(m12, pars="bMachiSchool")
for ( z1 in 1:J ){
  print(traceplot(m12, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m12, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z3 in 1:K){
  print(traceplot(m12, pars=paste("gamma[", z3, "]", sep="") ))
}
traceplot(m12, pars="sigma_beta")
traceplot(m12, pars="sigma_gamma")
graphics.off()



#######base model: m13, Eth, School.mest x Machi ############################

model_code_13 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> School[N]; // school experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bSchool;
  real bMachiSchool;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_13 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  School = d$School.mest
)

start_13 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bSchool=0, bMachiSchool=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m13 <- stan( model_code=model_code_13, data=data_list_13,
            init=list(start_13, start_13, start_13, start_13), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m13, pars=c("b0","b1","bMachi",
                "bSchool", "bMachiSchool",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post13 <- extract.samples( m13 )
str(post13)


#look at all traces
pdf(file="./traces_m13.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m13, pars="b0")
traceplot(m13, pars="bMachi")
traceplot(m13, pars="bSchool")
traceplot(m13, pars="bMachiSchool")
for ( z1 in 1:J ){
  print(traceplot(m13, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m13, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m13, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m13, pars="sigma_beta")
traceplot(m13, pars="sigma_gamma")
graphics.off()



#######base model: m14, Eth, ExEmpMest x Machi ############################

model_code_14 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bEmpMest;
  real bMachiEmpMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1);
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_14 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  EmpMest = d$ExEmpMest
)

start_14 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bEmpMest=0, bMachiEmpMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m14 <- stan( model_code=model_code_14, data=data_list_14,
            init=list(start_14, start_14, start_14, start_14), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m14, pars=c("b0","b1","bMachi",
                "bEmpMest", "bMachiEmpMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post14 <- extract.samples( m14 )
str(post14)


#look at all traces
pdf(file="./traces_m14.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m14, pars="b0")
traceplot(m14, pars="bMachi")
traceplot(m14, pars="bEmpMest")
traceplot(m14, pars="bMachiEmpMest")
for ( z1 in 1:J ){
  print(traceplot(m14, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m14, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m14, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m14, pars="sigma_beta")
traceplot(m14, pars="sigma_gamma")
graphics.off()


#######base model: m15, Eth, ExComMest x Machi ############################

model_code_15 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_15 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  ComMest = d$ExComMest
)

start_15 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m15 <- stan( model_code=model_code_15, data=data_list_15,
            init=list(start_15, start_15, start_15, start_15), 
                        iter=5000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m15, pars=c("b0","b1","bMachi",
                "bComMest", "bMachiComMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post15 <- extract.samples( m15 )
str(post15)


#look at all traces
pdf(file="./traces_m15.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m15, pars="b0")
traceplot(m15, pars="bMachi")
traceplot(m15, pars="bComMest")
traceplot(m15, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m15, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m15, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m15, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m15, pars="sigma_beta")
traceplot(m15, pars="sigma_gamma")
graphics.off()



#######base model: m16, Eth, School.mest x Machi, ExEmpMest x Machi ############################

model_code_16 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> School[N]; // school experience with mestizos
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of people (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bSchool;
  real bMachiSchool;
  real bEmpMest;
  real bMachiEmpMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_16 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  School = d$School.mest,
  EmpMest = d$ExEmpMest
)

start_16 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bSchool=0, bMachiSchool=0,
  bEmpMest=0, bMachiEmpMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m16 <- stan( model_code=model_code_16, data=data_list_16,
            init=list(start_16, start_16, start_16, start_16), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m16, pars=c("b0","b1","bMachi",
                "bSchool", "bMachiSchool",
                "bEmpMest", "bMachiEmpMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post16 <- extract.samples( m16 )
str(post16)


#look at all traces
pdf(file="./traces_m16.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m16, pars="b0")
traceplot(m16, pars="bMachi")
traceplot(m16, pars="bSchool")
traceplot(m16, pars="bMachiSchool")
traceplot(m16, pars="bEmpMest")
traceplot(m16, pars="bMachiEmpMest")
for ( z1 in 1:J ){
  print(traceplot(m16, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m16, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m16, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m16, pars="sigma_beta")
traceplot(m16, pars="sigma_gamma")
graphics.off()



#######base model: m17, Eth, School.mest x Machi, ExComMest x Machi ############################

model_code_17 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> School[N]; // school experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bSchool;
  real bMachiSchool;
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1);
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_17 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  School = d$School.mest,
  ComMest = d$ExComMest
)

start_17 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bSchool=0, bMachiSchool=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m17 <- stan( model_code=model_code_17, data=data_list_17,
            init=list(start_17, start_17, start_17, start_17), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m17, pars=c("b0","b1","bMachi",
                "bSchool", "bMachiSchool",
                "bComMest", "bMachiComMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post17 <- extract.samples( m17 )
str(post17)


#look at all traces
pdf(file="./traces_m17.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m17, pars="b0")
traceplot(m17, pars="bMachi")
traceplot(m17, pars="bSchool")
traceplot(m17, pars="bMachiSchool")
traceplot(m17, pars="bComMest")
traceplot(m17, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m17, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m17, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m17, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m17, pars="sigma_beta")
traceplot(m17, pars="sigma_gamma")
graphics.off()



#######base model: m18, Eth, ExEmpMest x Machi, ExComMest x Machi ############################

model_code_18 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bEmpMest;
  real bMachiEmpMest;
  real bComMest;
  real bMachiComMest;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1); 
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_18 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  EmpMest = d$ExEmpMest,
  ComMest = d$ExComMest
)

start_18 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bEmpMest=0, bMachiEmpMest=0,
  bComMest=0, bMachiComMest=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m18 <- stan( model_code=model_code_18, data=data_list_18,
            init=list(start_18, start_18, start_18, start_18), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m18, pars=c("b0","b1","bMachi",
                "bEmpMest", "bMachiEmpMest",
                "bComMest", "bMachiComMest",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post18 <- extract.samples( m18 )
str(post18)


#look at all traces
pdf(file="./traces_m18.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m18, pars="b0")
traceplot(m18, pars="bMachi")
traceplot(m18, pars="bEmpMest")
traceplot(m18, pars="bMachiEmpMest")
traceplot(m18, pars="bComMest")
traceplot(m18, pars="bMachiComMest")
for ( z1 in 1:J ){
  print(traceplot(m18, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m18, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z2 in 1:K ){
  print(traceplot(m18, pars=paste("gamma[", z2, "]", sep="") ))
}
traceplot(m18, pars="sigma_beta")
traceplot(m18, pars="sigma_gamma")
graphics.off()


#######base model: m19, Eth, ExEmpMest x Machi, ExComMest x Machi, School.mest x Machi

model_code_19 <- '

data {
  int<lower=1> J; // number of interviewees
  int<lower=1> K; // number of questions
  int<lower=1> N; // number of answers to questions (observations)
  int<lower=1,upper=J> jj[N]; // interviewee ID for observation n
  int<lower=1,upper=K> kk[N]; // question for observation n
  int<lower=0,upper=1> y[N]; // prob of positive response for obs n
  int<lower=0,upper=1> Machi[N]; // predictor for ethnicity
  int<lower=0,upper=1> EmpMest[N]; // employment experience with mestizos
  int<lower=0,upper=1> ComMest[N]; // commerce experience with mestizos
  int<lower=0,upper=1> School[N]; // school experience with mestizos
}

parameters {

  real beta[K];
  real gamma[K];

  real b0;     // mean interviewee location in latent dimension (mean ability intercept)
  real b1[J]; // location of everyone else (differing from the mean), i.e., random effect for person

  real bMachi; // effect of machiness on person location in latent space
  real bEmpMest;
  real bMachiEmpMest;
  real bComMest;
  real bMachiComMest;
  real bSchool;
  real bMachiSchool;

  real<lower=0> sigma_beta; // scale of question
  real<lower=0> sigma_gamma; // scale of discrimination
}


model {
  vector[N] params;
  real alpha;
  
  b0 ~ normal(0,1);
  b1 ~ normal(0,1); //identifying prior for location and scale

  bMachi ~ normal(0,1);
  bEmpMest ~ normal(0,1);
  bMachiEmpMest ~ normal(0,1);
  bComMest ~ normal(0,1);
  bMachiComMest ~ normal(0,1);
  bSchool ~ normal(0,1);
  bMachiSchool ~ normal(0,1);

  beta ~ normal(0,sigma_beta);
  gamma ~ normal(0,sigma_gamma);
  
  sigma_beta ~ exponential(1); //exponential(beta), where here beta = lambda = 1/mean
  sigma_gamma ~ exponential(1); //or use uniform(0,5), both prevent ceiling effect


  for (n in 1:N) {

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] <- ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] <- bernoulli_logit_log( y[n],    //needed for waic function, see Stan manual pg340
                    gamma[kk[n]]*(alpha[jj[n]] - beta[kk[n]]) );
  } //for
}
'

data_list_19 <- list(
  J = length(unique(d$newID)),  #number of people
  K = length(unique(d$Question)),   #number of questions
  N = nrow(d),      #total number of responses
  jj = d$newID,     #vector of person IDs
  kk = d$QID,     #vector of question IDs
  y = d$Response,      #vector of responses
  Machi = d$Machi,   #vector of ethnicity codes
  EmpMest = d$ExEmpMest,
  ComMest = d$ExComMest,
  School = d$School.mest
)

start_19 <- list(
  b0=0, b1 = as.array(rep(0, times=J)),
  bMachi=0,
  bEmpMest=0, bMachiEmpMest=0,
  bComMest=0, bMachiComMest=0,
  bSchool=0, bMachiSchool=0,
  beta = as.array(rep(0, times=K)),
  gamma = as.array(rep(1, times=K)),
  sigma_beta=1, sigma_gamma=1
)


m19 <- stan( model_code=model_code_19, data=data_list_19,
            init=list(start_19, start_19, start_19, start_19), 
                        iter=4000 , chains=4,  
                        control=list(adapt_delta=0.99) )

print(m19, pars=c("b0","b1","bMachi",
                "bEmpMest", "bMachiEmpMest",
                "bComMest", "bMachiComMest",
                "bSchool", "bMachiSchool",
                "beta", "gamma",
                "sigma_beta", "sigma_gamma"),  
      probs = c(0.025,0.975), digits_summary=2)

#for dotplots below 
post19 <- extract.samples( m19 )
str(post19)


#look at all traces
pdf(file="./traces_m19.pdf") #, 
  #height=10, width=8)
par(mfrow=c(2,1))

traceplot(m19, pars="b0")
traceplot(m19, pars="bMachi")
traceplot(m19, pars="bEmpMest")
traceplot(m19, pars="bMachiEmpMest")
traceplot(m19, pars="bComMest")
traceplot(m19, pars="bMachiComMest")
traceplot(m19, pars="bSchool")
traceplot(m19, pars="bMachiSchool")
for ( z1 in 1:J ){
  print(traceplot(m19, pars=paste("b1[", z1, "]", sep="") ))
  }
for ( z2 in 1:(K) ){
  print(traceplot(m19, pars=paste("beta[", z2, "]", sep="") ))
  }
for ( z3 in 1:K ){
  print(traceplot(m19, pars=paste("gamma[", z3, "]", sep="") ))
}
traceplot(m19, pars="sigma_beta")
traceplot(m19, pars="sigma_gamma")
graphics.off()




####################### WAIC model comparison
model_waics <- c(waic(m1)$waic,
                 waic(m2)$waic,
                 waic(m3)$waic,
                 waic(m4)$waic,
                 waic(m5)$waic,
                 waic(m6)$waic,
                 waic(m7)$waic,
                 waic(m8)$waic,
                 waic(m9)$waic,
                 waic(m10)$waic,
                 waic(m11)$waic,
                 waic(m12)$waic,
                 waic(m13)$waic,
                 waic(m14)$waic,
                 waic(m15)$waic,
                 waic(m16)$waic,
                 waic(m17)$waic,
                 waic(m18)$waic,
                 waic(m19)$waic)
model_pwaics <- c(waic(m1)$p_waic, #effective number of parameters
                 waic(m2)$p_waic,
                 waic(m3)$p_waic,
                 waic(m4)$p_waic,
                 waic(m5)$p_waic,
                 waic(m6)$p_waic,
                 waic(m7)$p_waic,
                 waic(m8)$p_waic,
                 waic(m9)$p_waic,
                 waic(m10)$p_waic,
                 waic(m11)$p_waic,
                 waic(m12)$p_waic,
                 waic(m13)$p_waic,
                 waic(m14)$p_waic,
                 waic(m15)$p_waic,
                 waic(m16)$p_waic,
                 waic(m17)$p_waic,
                 waic(m18)$p_waic,
                 waic(m19)$p_waic)



compare(m19, m17, m16, m12) #last column gives standard error of the difference in waic b/t models
x <- compare(m19,m17,m16,m12)
str(x)
x@dSE #gives matrix of differences in waic between all models
x #last column gives standard error of differences from the best model
#waic difference b/t m17 and m19 is a lot smaller than standard error of waic difference
#so m19 is probably slightly overfit, and effect size of labor is likely small
#In other words: predictive power of labor is small if you already know edu and com.
#but use m19 anyway because it has all the predictors of interest.


waic_diffs <- rep(0, length(model_waics))
for ( i in 1:length(model_waics) ) {
  waic_diffs[i] <- model_waics[i] - min(model_waics) 
} #for i

model_weights <- rep(0, length(model_waics))
for ( j in 1:length(model_waics) ) {
  model_weights[j] <- exp(-0.5*waic_diffs[j])/sum(exp(-0.5*waic_diffs)) #McElreath rethinking pg199
} #for j

WAIC_sum <- cbind(model_waics, model_pwaics, waic_diffs, model_weights)
colnames(WAIC_sum) <- c("WAIC", "pWAIC", "dWAIC", "weight")
rownames(WAIC_sum) <- c("m1", "m2", "m3", "m4", "m5",
                        "m6", "m7", "m8", "m9", "m10",
                        "m11", "m12", "m13", "m14", "m15",
                        "m16", "m17", "m18", "m19")
WAIC_sum <- WAIC_sum[order(WAIC_sum[,1]),]
print(WAIC_sum, digits=2)


#compare posterior parameter estimates for models 19 and 12

#labor
HPDI(post19$bEmpMest - post12$bEmpMest, prob=0.95)
HPDI(post19$bMachiEmpMest - post12$bMachiEmpMest, prob=0.95)




#### table of model outputs

cnames <- c("Model", "Ethn", "Sex", "Ethn X Sex", "Adol", "Adult", "Elder",
                    "Educ", "Ethn X Educ", "Labor", "Ethn X Labor",
                    "Comm", "Ethn X Comm", "WAIC weight")
rnames <- c(
                      "m1", " ",
                      "m2", " ",
                      "m3", " ",
                      "m4", " ",
                      "m5", " ",
                      "m6", " ",
                      "m7", " ",
                      "m8", " ",
                      "m9", " ",
                      "m10", " ",
                      "m11", " ",
                      "m12", " ",
                      "m13", " ",
                      "m14", " ",
                      "m15", " ",
                      "m16", " ",
                      "m17", " ",
                      "m18", " ",
                      "m19", " "
                  )

digs <- 2 #number of rounding digits
digs2 <- 2 #rounding digits for waic weights
cred <- 0.95 #HPDI
options(digits=10)

eth_col <- c(

  #Ethn
  " ",                              #m1
  " ",
  round(mean(post2$bMachi), digs),  #m2
  paste( round(as.vector(HPDI(post2$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post2$bMachi, prob=cred))[2], digs) ),
  round(mean(post3$bMachi), digs),  #m3
  paste( round(as.vector(HPDI(post3$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post3$bMachi, prob=cred))[2], digs) ),
  round(mean(post4$bMachi), digs),  #m4
  paste( round(as.vector(HPDI(post4$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post4$bMachi, prob=cred))[2], digs) ),
  round(mean(post5$bMachi), digs),  #m5
  paste( round(as.vector(HPDI(post5$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$bMachi, prob=cred))[2], digs) ),
  round(mean(post6$bMachi), digs),  #m6
  paste( round(as.vector(HPDI(post6$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bMachi, prob=cred))[2], digs) ),
  round(mean(post7$bMachi), digs),  #m7
  paste( round(as.vector(HPDI(post7$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bMachi, prob=cred))[2], digs) ),
  round(mean(post8$bMachi), digs),  #m8
  paste( round(as.vector(HPDI(post8$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bMachi, prob=cred))[2], digs) ),
  round(mean(post9$bMachi), digs),  #m9
  paste( round(as.vector(HPDI(post9$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bMachi, prob=cred))[2], digs) ),
  round(mean(post10$bMachi), digs),  #m10
  paste( round(as.vector(HPDI(post10$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bMachi, prob=cred))[2], digs) ),
  round(mean(post11$bMachi), digs),  #m11
  paste( round(as.vector(HPDI(post11$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bMachi, prob=cred))[2], digs) ),
  round(mean(post12$bMachi), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMachi, prob=cred))[2], digs) ),
  round(mean(post13$bMachi), digs),  #m13
  paste( round(as.vector(HPDI(post13$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post13$bMachi, prob=cred))[2], digs) ),
  round(mean(post14$bMachi), digs),  #m14
  paste( round(as.vector(HPDI(post14$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post14$bMachi, prob=cred))[2], digs) ),
  round(mean(post15$bMachi), digs),  #m15
  paste( round(as.vector(HPDI(post15$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post15$bMachi, prob=cred))[2], digs) ),
  round(mean(post16$bMachi), digs),  #m16
  paste( round(as.vector(HPDI(post16$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post16$bMachi, prob=cred))[2], digs) ),
  round(mean(post17$bMachi), digs),  #m17
  paste( round(as.vector(HPDI(post17$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post17$bMachi, prob=cred))[2], digs) ),
  round(mean(post18$bMachi), digs),  #m18
  paste( round(as.vector(HPDI(post18$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post18$bMachi, prob=cred))[2], digs) ),
  round(mean(post19$bMachi), digs),  #m19
  paste( round(as.vector(HPDI(post19$bMachi, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bMachi, prob=cred))[2], digs) )
)

sex_col <- c(

  #Sex
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  round(mean(post3$bMale), digs),  #m3
  paste( round(as.vector(HPDI(post3$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post3$bMale, prob=cred))[2], digs) ),
  " ",                              #m4
  " ",
  round(mean(post5$bMale), digs),  #m5
  paste( round(as.vector(HPDI(post5$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$bMale, prob=cred))[2], digs) ),
  round(mean(post6$bMale), digs),  #m6
  paste( round(as.vector(HPDI(post6$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bMale, prob=cred))[2], digs) ),
  round(mean(post7$bMale), digs),  #m7
  paste( round(as.vector(HPDI(post7$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bMale, prob=cred))[2], digs) ),
  round(mean(post8$bMale), digs),  #m8
  paste( round(as.vector(HPDI(post8$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bMale, prob=cred))[2], digs) ),
  round(mean(post9$bMale), digs),  #m9
  paste( round(as.vector(HPDI(post9$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bMale, prob=cred))[2], digs) ),
  round(mean(post10$bMale), digs),  #m10
  paste( round(as.vector(HPDI(post10$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bMale, prob=cred))[2], digs) ),
  round(mean(post11$bMale), digs),  #m11
  paste( round(as.vector(HPDI(post11$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bMale, prob=cred))[2], digs) ),
  round(mean(post12$bMale), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMale, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMale, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  " ",                              #m16
  " ",
  " ",                              #m17
  " ",
  " ",                              #m18
  " ",
  " ",                              #m19
  " "
)

ethsex_col <- c(

#Eth X Sex
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  round(mean(post3$bMaleMach), digs),  #m3
  paste( round(as.vector(HPDI(post3$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post3$bMaleMach, prob=cred))[2], digs) ),
  " ",                              #m4
  " ",
  round(mean(post5$bMaleMach), digs),  #m5
  paste( round(as.vector(HPDI(post5$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post6$bMaleMach), digs),  #m6
  paste( round(as.vector(HPDI(post6$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post7$bMaleMach), digs),  #m7
  paste( round(as.vector(HPDI(post7$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post8$bMaleMach), digs),  #m8
  paste( round(as.vector(HPDI(post8$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post9$bMaleMach), digs),  #m9
  paste( round(as.vector(HPDI(post9$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post10$bMaleMach), digs),  #m10
  paste( round(as.vector(HPDI(post10$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post11$bMaleMach), digs),  #m11
  paste( round(as.vector(HPDI(post11$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bMaleMach, prob=cred))[2], digs) ),
  round(mean(post12$bMaleMach), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMaleMach, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMaleMach, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  " ",                              #m16
  " ",
  " ",                              #m17
  " ",
  " ",                              #m18
  " ",
  " ",                              #m19
  " "
)

adol_col <- c(

#Adol
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  round(mean(post4$badol), digs),  #m4
  paste( round(as.vector(HPDI(post4$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post4$badol, prob=cred))[2], digs) ),
  round(mean(post5$badol), digs),  #m5
  paste( round(as.vector(HPDI(post5$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$badol, prob=cred))[2], digs) ),
  round(mean(post6$badol), digs),  #m6
  paste( round(as.vector(HPDI(post6$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$badol, prob=cred))[2], digs) ),
  round(mean(post7$badol), digs),  #m7
  paste( round(as.vector(HPDI(post7$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$badol, prob=cred))[2], digs) ),
  round(mean(post8$badol), digs),  #m8
  paste( round(as.vector(HPDI(post8$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$badol, prob=cred))[2], digs) ),
  round(mean(post9$badol), digs),  #m9
  paste( round(as.vector(HPDI(post9$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$badol, prob=cred))[2], digs) ),
  round(mean(post10$badol), digs),  #m10
  paste( round(as.vector(HPDI(post10$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$badol, prob=cred))[2], digs) ),
  round(mean(post11$badol), digs),  #m11
  paste( round(as.vector(HPDI(post11$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$badol, prob=cred))[2], digs) ),
  round(mean(post12$badol), digs),  #m12
  paste( round(as.vector(HPDI(post12$badol, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$badol, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  " ",                              #m16
  " ",
  " ",                              #m17
  " ",
  " ",                              #m18
  " ",
  " ",                              #m19
  " "
)


adult_col <- c(

#Adult
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  round(mean(post4$bmat), digs),  #m4
  paste( round(as.vector(HPDI(post4$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post4$bmat, prob=cred))[2], digs) ),
  round(mean(post5$bmat), digs),  #m5
  paste( round(as.vector(HPDI(post5$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$bmat, prob=cred))[2], digs) ),
  round(mean(post6$bmat), digs),  #m6
  paste( round(as.vector(HPDI(post6$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bmat, prob=cred))[2], digs) ),
  round(mean(post7$bmat), digs),  #m7
  paste( round(as.vector(HPDI(post7$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bmat, prob=cred))[2], digs) ),
  round(mean(post8$bmat), digs),  #m8
  paste( round(as.vector(HPDI(post8$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bmat, prob=cred))[2], digs) ),
  round(mean(post9$bmat), digs),  #m9
  paste( round(as.vector(HPDI(post9$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bmat, prob=cred))[2], digs) ),
  round(mean(post10$bmat), digs),  #m10
  paste( round(as.vector(HPDI(post10$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bmat, prob=cred))[2], digs) ),
  round(mean(post11$bmat), digs),  #m11
  paste( round(as.vector(HPDI(post11$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bmat, prob=cred))[2], digs) ),
  round(mean(post12$bmat), digs),  #m12
  paste( round(as.vector(HPDI(post12$bmat, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bmat, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  " ",                              #m16
  " ",
  " ",                              #m17
  " ",
  " ",                              #m18
  " ",
  " ",                              #m19
  " "
)


eld_col <- c(

#Elder
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  round(mean(post4$bold), digs),  #m4
  paste( round(as.vector(HPDI(post4$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post4$bold, prob=cred))[2], digs) ),
  round(mean(post5$bold), digs),  #m5
  paste( round(as.vector(HPDI(post5$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post5$bold, prob=cred))[2], digs) ),
  round(mean(post6$bold), digs),  #m6
  paste( round(as.vector(HPDI(post6$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bold, prob=cred))[2], digs) ),
  round(mean(post7$bold), digs),  #m7
  paste( round(as.vector(HPDI(post7$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bold, prob=cred))[2], digs) ),
  round(mean(post8$bold), digs),  #m8
  paste( round(as.vector(HPDI(post8$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bold, prob=cred))[2], digs) ),
  round(mean(post9$bold), digs),  #m9
  paste( round(as.vector(HPDI(post9$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bold, prob=cred))[2], digs) ),
  round(mean(post10$bold), digs),  #m10
  paste( round(as.vector(HPDI(post10$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bold, prob=cred))[2], digs) ),
  round(mean(post11$bold), digs),  #m11
  paste( round(as.vector(HPDI(post11$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bold, prob=cred))[2], digs) ),
  round(mean(post12$bold), digs),  #m12
  paste( round(as.vector(HPDI(post12$bold, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bold, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  " ",                              #m16
  " ",
  " ",                              #m17
  " ",
  " ",                              #m18
  " ",
  " ",                              #m19
  " "
)


edu_col <- c(

  #Educ
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  round(mean(post6$bSchool), digs),  #m6
  paste( round(as.vector(HPDI(post6$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bSchool, prob=cred))[2], digs) ),
  " ",                              #m7
  " ",
  " ",                              #m8
  " ",
  round(mean(post9$bSchool), digs),  #m9
  paste( round(as.vector(HPDI(post9$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bSchool, prob=cred))[2], digs) ),
  round(mean(post10$bSchool), digs),  #m10
  paste( round(as.vector(HPDI(post10$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bSchool, prob=cred))[2], digs) ),
  " ",                              #m11
  " ",
  round(mean(post12$bSchool), digs),  #m12
  paste( round(as.vector(HPDI(post12$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bSchool, prob=cred))[2], digs) ),
  round(mean(post13$bSchool), digs),  #m13
  paste( round(as.vector(HPDI(post13$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post13$bSchool, prob=cred))[2], digs) ),
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  round(mean(post16$bSchool), digs),  #m16
  paste( round(as.vector(HPDI(post16$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post16$bSchool, prob=cred))[2], digs) ),
  round(mean(post17$bSchool), digs),  #m17
  paste( round(as.vector(HPDI(post17$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post17$bSchool, prob=cred))[2], digs) ),
  " ",                              #m18
  " ",
  round(mean(post19$bSchool), digs),  #m19
  paste( round(as.vector(HPDI(post19$bSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bSchool, prob=cred))[2], digs) )
)


ethedu_col <- c(

  #Ethn X Educ
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  round(mean(post6$bMachiSchool), digs),  #m6
  paste( round(as.vector(HPDI(post6$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post6$bMachiSchool, prob=cred))[2], digs) ),
  " ",                              #m7
  " ",
  " ",                              #m8
  " ",
  round(mean(post9$bMachiSchool), digs),  #m9
  paste( round(as.vector(HPDI(post9$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bMachiSchool, prob=cred))[2], digs) ),
  round(mean(post10$bMachiSchool), digs),  #m10
  paste( round(as.vector(HPDI(post10$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bMachiSchool, prob=cred))[2], digs) ),
  " ",                              #m11
  " ",
  round(mean(post12$bMachiSchool), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMachiSchool, prob=cred))[2], digs) ),
  round(mean(post13$bMachiSchool), digs),  #m13
  paste( round(as.vector(HPDI(post13$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post13$bMachiSchool, prob=cred))[2], digs) ),
  " ",                              #m14
  " ",
  " ",                              #m15
  " ",
  round(mean(post16$bMachiSchool), digs),  #m16
  paste( round(as.vector(HPDI(post16$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post16$bMachiSchool, prob=cred))[2], digs) ),
  round(mean(post17$bMachiSchool), digs),  #m17
  paste( round(as.vector(HPDI(post17$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post17$bMachiSchool, prob=cred))[2], digs) ),
  " ",                              #m18
  " ",
  round(mean(post19$bMachiSchool), digs),  #m19
  paste( round(as.vector(HPDI(post19$bMachiSchool, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bMachiSchool, prob=cred))[2], digs) )
)


lab_col <- c(

  #Labor
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  " ",                              #m6
  " ",
  round(mean(post7$bEmpMest), digs),  #m7
  paste( round(as.vector(HPDI(post7$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bEmpMest, prob=cred))[2], digs) ),
  " ",                              #m8
  " ",
  round(mean(post9$bEmpMest), digs),  #m9
  paste( round(as.vector(HPDI(post9$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bEmpMest, prob=cred))[2], digs) ),
  " ",                              #m10
  " ",
  round(mean(post11$bEmpMest), digs),  #m11
  paste( round(as.vector(HPDI(post11$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bEmpMest, prob=cred))[2], digs) ),
  round(mean(post12$bEmpMest), digs),  #m12
  paste( round(as.vector(HPDI(post12$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bEmpMest, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  round(mean(post14$bEmpMest), digs),  #m14
  paste( round(as.vector(HPDI(post14$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post14$bEmpMest, prob=cred))[2], digs) ),
  " ",                              #m15
  " ",
  round(mean(post16$bEmpMest), digs),  #m16
  paste( round(as.vector(HPDI(post16$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post16$bEmpMest, prob=cred))[2], digs) ),
  " ",                              #m17
  " ",
  round(mean(post18$bEmpMest), digs),  #m18
  paste( round(as.vector(HPDI(post18$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post18$bEmpMest, prob=cred))[2], digs) ),
  round(mean(post19$bEmpMest), digs),  #m19
  paste( round(as.vector(HPDI(post19$bEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bEmpMest, prob=cred))[2], digs) )
)


ethlab_col <- c(

  #Ethn X Labor
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  " ",                              #m6
  " ",
  round(mean(post7$bMachiEmpMest), digs),  #m7
  paste( round(as.vector(HPDI(post7$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post7$bMachiEmpMest, prob=cred))[2], digs) ),
  " ",                              #m8
  " ",
  round(mean(post9$bMachiEmpMest), digs),  #m9
  paste( round(as.vector(HPDI(post9$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post9$bMachiEmpMest, prob=cred))[2], digs) ),
  " ",                              #m10
  " ",
  round(mean(post11$bMachiEmpMest), digs),  #m11
  paste( round(as.vector(HPDI(post11$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bMachiEmpMest, prob=cred))[2], digs) ),
  round(mean(post12$bMachiEmpMest), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMachiEmpMest, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  round(mean(post14$bMachiEmpMest), digs),  #m14
  paste( round(as.vector(HPDI(post14$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post14$bMachiEmpMest, prob=cred))[2], digs) ),
  " ",                              #m15
  " ",
  round(mean(post16$bMachiEmpMest), digs),  #m16
  paste( round(as.vector(HPDI(post16$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post16$bMachiEmpMest, prob=cred))[2], digs) ),
  " ",                              #m17
  " ",
  round(mean(post18$bMachiEmpMest), digs),  #m18
  paste( round(as.vector(HPDI(post18$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post18$bMachiEmpMest, prob=cred))[2], digs) ),
  round(mean(post19$bMachiEmpMest), digs),  #m19
  paste( round(as.vector(HPDI(post19$bMachiEmpMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bMachiEmpMest, prob=cred))[2], digs) )
)


com_col <- c(

  #Comm
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  " ",                              #m6
  " ",
  " ",                              #m7
  " ",
  round(mean(post8$bComMest), digs),  #m8
  paste( round(as.vector(HPDI(post8$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bComMest, prob=cred))[2], digs) ),
  " ",                              #m9
  " ",
  round(mean(post10$bComMest), digs),  #m10
  paste( round(as.vector(HPDI(post10$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bComMest, prob=cred))[2], digs) ),
  round(mean(post11$bComMest), digs),  #m11
  paste( round(as.vector(HPDI(post11$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bComMest, prob=cred))[2], digs) ),
  round(mean(post12$bComMest), digs),  #m12
  paste( round(as.vector(HPDI(post12$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bComMest, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  round(mean(post15$bComMest), digs),  #m15
  paste( round(as.vector(HPDI(post15$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post15$bComMest, prob=cred))[2], digs) ),
  " ",                              #m16
  " ",
  round(mean(post17$bComMest), digs),  #m17
  paste( round(as.vector(HPDI(post17$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post17$bComMest, prob=cred))[2], digs) ),
  round(mean(post18$bComMest), digs),  #m18
  paste( round(as.vector(HPDI(post18$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post18$bComMest, prob=cred))[2], digs) ),
  round(mean(post19$bComMest), digs),  #m19
  paste( round(as.vector(HPDI(post19$bComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bComMest, prob=cred))[2], digs) )
)


ethcom_col <- c(

  #Ethn X Comm
  " ",                              #m1
  " ",
  " ",                              #m2
  " ",
  " ",                              #m3
  " ",
  " ",                              #m4
  " ",
  " ",                              #m5
  " ",
  " ",                              #m6
  " ",
  " ",                              #m7
  " ",
  round(mean(post8$bMachiComMest), digs),  #m8
  paste( round(as.vector(HPDI(post8$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post8$bMachiComMest, prob=cred))[2], digs) ),
  " ",                              #m9
  " ",
  round(mean(post10$bMachiComMest), digs),  #m10
  paste( round(as.vector(HPDI(post10$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post10$bMachiComMest, prob=cred))[2], digs) ),
  round(mean(post11$bMachiComMest), digs),  #m11
  paste( round(as.vector(HPDI(post11$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post11$bMachiComMest, prob=cred))[2], digs) ),
  round(mean(post12$bMachiComMest), digs),  #m12
  paste( round(as.vector(HPDI(post12$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post12$bMachiComMest, prob=cred))[2], digs) ),
  " ",                              #m13
  " ",
  " ",                              #m14
  " ",
  round(mean(post15$bMachiComMest), digs),  #m15
  paste( round(as.vector(HPDI(post15$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post15$bMachiComMest, prob=cred))[2], digs) ),
  " ",                              #m16
  " ",
  round(mean(post17$bMachiComMest), digs),  #m17
  paste( round(as.vector(HPDI(post17$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post17$bMachiComMest, prob=cred))[2], digs) ),
  round(mean(post18$bMachiComMest), digs),  #m18
  paste( round(as.vector(HPDI(post18$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post18$bMachiComMest, prob=cred))[2], digs) ),
  round(mean(post19$bMachiComMest), digs),  #m19
  paste( round(as.vector(HPDI(post19$bMachiComMest, prob=cred))[1], digs), ", ",
         round(as.vector(HPDI(post19$bMachiComMest, prob=cred))[2], digs) )
)


waic_col <- c(

  #WAIC weight
  signif(model_weights[1], digs2),     #m1
  " ",
  signif(model_weights[2], digs2),     #m2
  " ",
  signif(model_weights[3], digs2),     #m3
  " ",
  signif(model_weights[4], digs2),     #m4
  " ",
  signif(model_weights[5], digs2),     #m5
  " ",
  signif(model_weights[6], digs2),     #m6
  " ",
  signif(model_weights[7], digs2),     #m7
  " ",
  signif(model_weights[8], digs2),     #m8
  " ",
  signif(model_weights[9], digs2),     #m9
  " ",
  signif(model_weights[10], digs2),     #m10
  " ",
  signif(model_weights[11], digs2),     #m11
  " ",
  signif(model_weights[12], digs2),     #m12
  " ",
  signif(model_weights[13], digs2),     #m13
  " ",
  signif(model_weights[14], digs2),     #m14
  " ",
  signif(model_weights[15], digs2),     #m15
  " ",
  signif(model_weights[16], digs2),     #m16
  " ",
  signif(model_weights[17], digs2),     #m17
  " ",
  signif(model_weights[18], digs2),     #m18
  " ",
  signif(model_weights[19], digs2),     #m19
  " "
)


tab_mat <- cbind(rnames, eth_col, sex_col, ethsex_col, adol_col, 
                                adult_col, eld_col, edu_col, ethedu_col,
                                lab_col, ethlab_col, com_col, ethcom_col,
                                waic_col)

dim(tab_mat)

colnames(tab_mat) <- cnames

tab_mat[1:20,1:7]

stargazer(tab_mat, summary=FALSE, rownames=FALSE, type="text",
          out="./table.txt")

#Full script. Fits all models and creates all figures. It can take a few hours to run.


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
library(lattice)
library(stargazer)


rstan_options(auto_write = TRUE) #to let stan make a copy of the model and use multiple cores
options(mc.cores = parallel::detectCores())

#Read the data from the csv data file into R:
Interviews.raw <- read.csv(
	file=
  #"./Desktop/IRT_analysis/Manu_interviews_31oct16.csv",
	"./Matsi_IRT/Manu_interviews_31oct16.csv", 
	header=TRUE)

#Check the variable names and dimensions in the data frame Interviews.raw
names(Interviews.raw)
dim(Interviews.raw)


#working dataframe
d <- Interviews.raw
dim(d)
names(d)


#flip coding of some questions to polarize latent axis as in text

Flip <- ifelse(d$Question == 6 | #index of questions to flip
               d$Question == 17 |
               d$Question == 19 |
               d$Question == 22 |
               d$Question == 29 |
               d$Question == 30 , 1 , 0 )

d$ResponseFlipped <- ifelse( Flip == 1, ifelse( d$Response==1, 0, 1), d$Response )


#check flipping
d[1:20, c("Question", "Response", "ResponseFlipped")]

pdf(file="./Flip_check.pdf", 
height=4, width=4)
par(mfrow=c(1,1))

plot( jitter(d$ResponseFlipped) ~ jitter(d$Response), 
  col=ifelse(Flip==1,"red","black") )

graphics.off()


#replace Person ID with a new consecutive ID number
d$newID <- as.numeric( factor(d$ID, levels=unique(d$ID)) )
unique(d[,c("ID","newID")])


#replace Question with a new consecutive QID number
d$QID <- as.numeric( factor(d$Question, levels=unique(d$Question)) )
unique(d[,c("Question","QID")])


d <- d[order(d$newID),] #reorder by consecutive newID
dim(d)
names(d)


#number of people answering each question (Questions 9 and 12 not asked in Boca Manu)
table( d$QID, d$Machi )



J <- length(unique(d$newID))	#number of people
K <- length(unique(d$QID))		#number of questions
N <- nrow(d)			#total number of responses
jj <- d$newID			#vector of person IDs
kk <- d$QID			#vector of question IDs
y <- d$ResponseFlipped			#vector of responses


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


###########################raw proportions, Figure 2

props1 <- data.frame(
  Quest=c(1:K) ,
  Matsi=table( d$QID[which(d$Machi==1)], d$ResponseFlipped[which(d$Machi==1)] )[,2] /
    rowSums( table( d$QID[which(d$Machi==1)], d$ResponseFlipped[which(d$Machi==1)]) ),

  Mest=table( d$QID[which(d$Machi==0)], d$ResponseFlipped[which(d$Machi==0)] )[,2] /
    rowSums( table( d$QID[which(d$Machi==0)], d$ResponseFlipped[which(d$Machi==0)]) )
)

colnames(props1)[c(1:3)] <- c("Quest", "mach_yes", "mest_yes")


props1 <- cbind(props, quest_names)
print(props1, digits=2)

props <- props1[order(props1$mach_yes),] #reorder by machi percentage


#point altitude caluclations
x0_5 <- props[props$Quest==5, "mest_yes"]
y0_5 <- props[props$Quest==5, "mach_yes"]
b_5 <- y0_5 - x0_5
b_5_round <- round(b_5,2)

x0_11 <- props[props$Quest==11, "mest_yes"]
y0_11 <- props[props$Quest==11, "mach_yes"]
b_11 <- y0_11 - x0_11
b_11_round <- round(b_11,2)

pdf(file="./Fig2_raw_props.pdf", 
height=4, width=4)
par(mfrow=c(1,1))

xyplot( props$mach_yes ~ props$mest_yes,
  main="",
  xlim=c(-0.05, 1.05),
  ylim=c(-0.05, 1.05),
  xlab=list("Mestizo Proportion Inter-dependence", cex=0.8),
  ylab=list("Matsigenka Proportion Inter-dependence", cex=0.8),
  scales=list( tck=c(1,1), cex=0.75, alternating=1 ),
  panel = function (x, y) {
          panel.segments( x0=x0_5,
                          y0=y0_5,
                          x1=x0_5, 
                          y1=x0_5, col="black", lwd=0.5 )
          panel.segments( x0=x0_11,
                          y0=y0_11,
                          x1=x0_11, 
                          y1=x0_11, col="black", lwd=0.5 )
          panel.xyplot( x=props$mest_yes, y=props$mach_yes,
                         pch = 16, cex=1.8, col = "white" )
          panel.xyplot( x=props$mest_yes, y=props$mach_yes,
                         pch = 1, cex=1.8, col = "black" )
          panel.segments( x0=-0.05, y0=-0.05, #diag line y=x
                     x1=1.05, y1=1.05, col="black", lwd=0.5 )
          ltext(x=props$mest_yes,
                 y=props$mach_yes,
                 labels=props$Quest, pos=1, offset=-0.35, cex=0.7, col="black")
          ltext(x=x0_5 - 0.04,
                 y=0.35,
                 labels=b_5_round, pos=1, offset=0, cex=0.65, col="black")
          ltext(x=x0_11 - 0.05,
                 y=0.53,
                 labels=b_11_round, pos=1, offset=0, cex=0.65, col="black")
  } #panel

) #xyplot
graphics.off()



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

    alpha[jj[n]] = ( b0 + b1[jj[n]] );

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
	y = d$ResponseFlipped			#vector of responses
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



################## Appendix Fig B1
#plot logistic response functions for a given question and latent trait estimates of indivs

questions <- quest_names #from above

#look at all plots
pdf(file="./FigB1_logistic_m1.pdf", 
  height=10, width=12)
par(mfrow=c(4,4), oma=c(5,7,5,5), mar=c(3,2,2,2))

gamma <- post1$gamma

for (Quest in 1:K) {
  plot( t( d[which(d$QID == Quest), "ResponseFlipped"] ) ~
        t(  mean(post1$b0) +
            colMeans(post1$b1[,d[which(d$QID == Quest), "newID"]]) ),
      col="red", ylim=c(0,1), xlim=c(-2.5,2.5), 
      ylab="",
      xlab="",
      main=c(questions[Quest]),
      cex.main=1.2,
      cex.axis=1.2,
      xaxp=c(-2, 2, 2),
      yaxp=c(0, 1, 2)
      )

  InverseLogit <- function(x) 1/(1+exp(-1*x)) # logit^-1 = logistic, undoing the logit function in the model
  curve( InverseLogit( mean(gamma[,Quest])*
    (x - mean(post1$beta[,Quest])) ),
    from=-2.5, to=2.5, add=T)
  abline( v= ( mean(post1$beta[,Quest]) ) , #Bafumi et al. 2005, pg 174-175
          col=col.alpha("black",0.2) )
} #for
mtext(text="Prob(Response = 1)", side=2, outer=TRUE, line=2, cex=2, las=3, adj=0.5)
mtext(text="Latent Axis", side=1, outer=TRUE, line=1, cex=2, las=1, adj=0.5)
graphics.off()



############################Fig 3 and Appendix Fig C1: base model with just random effect

####create plotting matrix
#get many more samples from the stan samples
b0_samp <- sample(post1$b0, size=10000, replace=T)
b1_samp <- matrix(nrow=10000, ncol=J)
for (z6 in 1:J) {
  b1_samp[,z6] <- sample(post1$b1[,z6], size=10000, replace=T)
}

machiID <- unique(d[which(d$Machi == 1), "newID"])
length(machiID) #number of machis in dataset
mestizoID <- unique(d[which(d$Machi == 0), "newID"])
length(mestizoID) #number of mestizos in dataset
machi <- matrix(data= cbind(machiID, 1), nrow <- length(machiID), ncol <- 2)
mestizo <- matrix(data= cbind(mestizoID, 0), nrow <- length(mestizoID), ncol <- 2)
machi_dum <- rbind(machi, mestizo)
machi_dum <- machi_dum[order(machi_dum[,1]),] #order based on newID
machi <- machi_dum[,2]

School.mestID <- unique(d[which(d$School.mest == 1), "newID"])
notSchool.mestID <- unique(d[which(d$School.mest == 0), "newID"])
School.mest <- matrix(data= cbind(School.mestID, 1), nrow <- length(School.mestID), ncol <- 2)
notSchool.mest <- matrix(data= cbind(notSchool.mestID, 0), nrow <- length(notSchool.mestID), ncol <- 2)
School.mest_dum <- rbind(School.mest, notSchool.mest)
School.mest_dum <- School.mest_dum[order(School.mest_dum[,1]),]
School.mest <- School.mest_dum[,2]

WorkMestID <- unique(d[which(d$ExEmpMest == 1), "newID"])
notWorkMestID <- unique(d[which(d$ExEmpMest == 0), "newID"])
WorkMest <- matrix(data= cbind(WorkMestID, 1), nrow <- length(WorkMestID), ncol <- 2)
notWorkMest <- matrix(data= cbind(notWorkMestID, 0), nrow <- length(notWorkMestID), ncol <- 2)
WorkMest_dum <- rbind(WorkMest, notWorkMest)
WorkMest_dum <- WorkMest_dum[order(WorkMest_dum[,1]),]
WorkMest <- WorkMest_dum[,2]

ComMestID <- unique(d[which(d$ExComMest == 1), "newID"])
notComMestID <- unique(d[which(d$ExComMest == 0), "newID"])
ComMest <- matrix(data= cbind(ComMestID, 1), nrow <- length(ComMestID), ncol <- 2)
notComMest <- matrix(data= cbind(notComMestID, 0), nrow <- length(notComMestID), ncol <- 2)
ComMest_dum <- rbind(ComMest, notComMest)
ComMest_dum <- ComMest_dum[order(ComMest_dum[,1]),]
ComMest <- ComMest_dum[,2]


abils <- mean(b0_samp) +  colMeans(b1_samp)
o <- c(J:1)
abils[o]  #inverse order for plotting
labs <- c(1:J)
ID <- unique(d[,"ID"])

#data matrix for plotting
all_1 <- as.data.frame( cbind(labs, ID, abils,
                            machi,
                            School.mest, WorkMest, ComMest
                            ) ) 
#sort by abils
all <- all_1[order(all_1$abils),]
all <- cbind(all, c(1:J))
names(all)
colnames(all)[8] <- c("ord")


#######################Appendix Fig C1: individual location highlights
pdf(file="./FigC1_m1.pdf",
    height=2, width=6)
dotplot( all$abils ~ all$ord , data=all,
         horizontal=F ,
         ylim=c(-2.5, 2.5),
         xlim=c(-0.5, 162),
         main=NULL, 
         cex=0.4, xlab=NULL, ylab="Latent Axis",
         scales=list(x=list(draw=F),
                     tck=c(1,0), cex=0.6,
                     alternating=1),
         panel = function (x, y) {
           panel.segments( x0=-0.5, y0=mean(all$abils), #mean line
                     x1=162, y1=mean(all$abils), lty=1, lwd=0.5, col=grey(0.25) )
           panel.xyplot(x, y, pch = 1, cex=0.3, col = "black") #grey(0.75))
           panel.xyplot(x, y, pch = 16, cex=0.25, col = "white") 
           for (p in 1:J) {
             if (all$machi[p] == 1) {
               panel.points(all$ord[p],
                        all$abils[p],
                        pch=16, cex=0.35, col="black") #grey(0.75))
             } #if
           } #for
 
           for ( p1 in 1:J ) { #work w/o school: 15 years
             if ( any(all$ID[p1] == c(4,29)) ) {
               panel.points(all$ord[p1],
                            all$abils[p1],
                            pch=16, cex=0.45, lwd=0.5, col="magenta")
             } #if
           } #for

           for ( p1 in 1:J ) { #school & work 10-15 years
             if ( any(all$ID[p1] == c(297,304,109,19,146,155,165,101,307)) ) {
               panel.points(all$ord[p1],
                            all$abils[p1],
                            pch=16, cex=0.45, lwd=0.5, col="orange")
             } #if
           } #for

           for ( p1 in 1:J ) { #unusual indiv
             if ( any(all$ID[p1] == c(55)) ) {
               panel.points(all$ord[p1],
                            all$abils[p1],
                            pch=16, cex=0.45, lwd=0.5, col="turquoise")
             } #if
           } #for

           for ( p1 in 1:J ) { #extreme machis who want to work
             if ( any(all$ID[p1] == c(38, 61)) ) {
               panel.points(all$ord[p1],
                            all$abils[p1],
                            pch=16, cex=0.45, lwd=0.5, col="green")
             } #if
           } #for

         } #panel
) #dotplot
graphics.off()



###############################Fig 3 individual location by interaction experience
jit <- jitter( rep(0,nrow(all)), factor=0, amount=0.15 )
all_jit <- all
all_jit$machi_jit <- all_jit$machi + jit
all$machi[1:10]
all_jit$machi_jit[1:10]

pdf(file="./Fig3_m1.pdf", 
    height=3.5, width=3.5)
xyplot( all_jit$abils ~ all_jit$machi,
  xlim=c(-0.6, 4.6),
  ylim=c(-2.7, 1.7),
  cex=0.4,
  xlab="",
  ylab="Latent Axis",
  scales=list( draw=T, tck=c(1,0), cex=0.6, #alternating=1,
              x=list(draw=F)
              ),
  panel = function (x, y) {
          #mestizos
          panel.xyplot( x=all_jit$machi_jit[all_jit$machi==0],
                        y=all_jit$abils[all_jit$machi==0],
                        pch = 1, cex=0.5, col = "black",
                        jitter.x = F, jitter.y = FALSE,
                        amount = 0.1)

          #all machis
          panel.xyplot( x=all_jit$machi_jit[all_jit$machi==1],
                        y=all_jit$abils[all_jit$machi==1],
                        pch = 16, cex=0.5, col = "black",
                        jitter.x = F, jitter.y = FALSE,
                        amount = 0.1)

          panel.xyplot( x=all_jit$machi_jit[all_jit$ComMest==1 & all_jit$machi==1] + 1,
                        y=all_jit$abils[all_jit$ComMest==1 & all_jit$machi==1],
                        pch = 16, cex=0.5, col = "black",
                        jitter.x = F, jitter.y = FALSE,
                        amount = 0.1) 

          panel.xyplot( x=all_jit$machi_jit[all_jit$WorkMest==1 & all_jit$machi==1] + 2,
                        y=all_jit$abils[all_jit$WorkMest==1 & all_jit$machi==1],
                        pch = 16, cex=0.5, col = "black",
                        jitter.x = F, jitter.y = FALSE,
                        amount = 0.1)

          panel.xyplot( x=all_jit$machi_jit[all_jit$School.mest==1 & all_jit$machi==1] + 3,
                        y=all_jit$abils[all_jit$School.mest==1 & all_jit$machi==1],
                        pch = 16, cex=0.5, col = "black",
                        jitter.x = F, jitter.y = FALSE,
                        amount = 0.1)

          ltext(x=0,
                 y=-2.4,
                 labels="Mestizos", pos=1, offset=0, cex=0.7, col="black")
          ltext(x=2.5,
                 y=-2.3,
                 labels="Matsigenka", pos=1, offset=0, cex=0.7, col="black")
          ltext(x=1,
                 y=-1.9,
                 labels="All", pos=1, offset=0, cex=0.7, col="black")
          ltext(x=2,
                 y=-1.9,
                 labels="w/ Com", pos=1, offset=0, cex=0.7, col="black")
          ltext(x=3,
                 y=-1.9,
                 labels="w/ Lab", pos=1, offset=0, cex=0.7, col="black")
          ltext(x=4,
                 y=-1.9,
                 labels="w/ Edu", pos=1, offset=0, cex=0.7, col="black")
          panel.segments( x0=0.7,
                          y0=-2.2,
                          x1=4.3, 
                          y1=-2.2, col="black", lwd=1 )
          panel.segments( x0=0.7,
                          y0=-2.2,
                          x1=0.7, 
                          y1=-2.1, col="black", lwd=1 )
          panel.segments( x0=4.3,
                          y0=-2.2,
                          x1=4.3, 
                          y1=-2.1, col="black", lwd=1 )

          panel.segments( x0=-0.4,
                          y0=mean(all_jit$abils[all_jit$machi==0]),
                          x1=0.4, 
                          y1=mean(all_jit$abils[all_jit$machi==0]),
                          col="black", lwd=0.5 )
          panel.segments( x0=0.6,
                          y0=mean(all_jit$abils[all_jit$machi==1]),
                          x1=1.4, 
                          y1=mean(all_jit$abils[all_jit$machi==1]),
                          col="black", lwd=0.5 )
          panel.segments( x0=1.6,
                          y0=mean(all_jit$abils[all_jit$ComMest==1 & all_jit$machi==1]),
                          x1=2.4, 
                          y1=mean(all_jit$abils[all_jit$ComMest==1 & all_jit$machi==1]),
                          col="black", lwd=0.5 )
          panel.segments( x0=2.6,
                          y0=mean(all_jit$abils[all_jit$WorkMest==1 & all_jit$machi==1]),
                          x1=3.4, 
                          y1=mean(all_jit$abils[all_jit$WorkMest==1 & all_jit$machi==1]),
                          col="black", lwd=0.5 )
          panel.segments( x0=3.6,
                          y0=mean(all_jit$abils[all_jit$School.mest==1 & all_jit$machi==1]),
                          x1=4.4, 
                          y1=mean(all_jit$abils[all_jit$School.mest==1 & all_jit$machi==1]),
                          col="black", lwd=0.5 )

  } #panel
) #xyplot
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + bMale*Male[n] +
                                bMaleMach*Machi[n]*Male[n] );

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n]+ bMale*Male[n] +
                                bMaleMach*Machi[n]*Male[n] );

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n]);

    params[n] = gamma[kk[n]]*(alpha - beta[kk[n]]);

  }; //for

  y ~ bernoulli_logit(params);
}

generated quantities {       //for computing waic
  vector[N] log_lik;
  vector[J] alpha;

  for (n in 1:N){

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      badol*adol[n] + bmat*mat[n] + bold*old[n] +
                      bMale*Male[n] + bMaleMach*Machi[n]*Male[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] );

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] +
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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

    alpha[jj[n]] = ( b0 + b1[jj[n]] + bMachi*Machi[n] + 
                      bEmpMest*EmpMest[n] + bMachiEmpMest*Machi[n]*EmpMest[n] + 
                      bComMest*ComMest[n] + bMachiComMest*Machi[n]*ComMest[n] +
                      bSchool*School[n] + bMachiSchool*Machi[n]*School[n]);

    log_lik[n] = bernoulli_logit_lpmf( y[n] |    //needed for waic function, see Stan manual pg479
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
  y = d$ResponseFlipped,      #vector of responses
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


########################Fig B2: plot discrimination posterior estimates for m19

#Richard's density plot function
denschart1 <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), 
    color = "grey20", colorHPDI ="grey60", HPDI=0.9,
    gcolor = par("fg"), lcolor = "gray", xlim = range(unlist(x)), 
    main = NULL, xlab = NULL, ylab = NULL, height=0.7 , border=NA, adjust=1, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.list(x)) 
        stop("'x' must be a list of vectors or matrices")
    n <- length(x)
    glabels <- NULL
    if (is.list(x)) {
        if (is.null(labels)) 
            labels <- names(x)
        if (is.null(labels)) 
            labels <- as.character(1L:n)
        labels <- rep_len(labels, n)
        #if (is.null(groups)) 
        #    groups <- col(x, as.factor = TRUE)
        #glabels <- levels(groups)
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1L:n
        y <- o #c(1,2,3, 4.5,5.5,6.5, 8,9,10) #o    #space out in groups of three
        ylim <- c(0.2, max(y)+0.7) #n + 1)
    }
    else {
        # sub-groups, so need more rows
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        #mtext(labs, side = 2, line = -1 #0.4, #loffset,           #### y-labels
        #      at = y, adj = 1, 
        #      col = "black", las = 2, cex = cex, ...)
        text(labels=labs, x=-0.5, y=y, pos=2, adj=1)
    }

    abline(v=0, lty=2, lwd=0.75) #dotted vertical line at 0

    #abline(h = y, lty = "dotted", col = lcolor)
    #points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)


    # draw densities at each y offset
    for ( i in 1:n ) {
        a <- density( x[[i]] , adjust=adjust )
        a$y <- a$y/max(a$y) * height + y[i] - 0.3
        polygon( a$x , a$y , col=color , border=border )
        Cuts <- HPDI( x[[i]] , HPDI)
        XX <- a$x[which(a$x > Cuts[1] & a$x < Cuts[2])]
        YY <- a$y[which(a$x > Cuts[1] & a$x < Cuts[2])]
        polygon( c(min(XX), XX, max(XX)), c(min(a$y), YY, min(a$y)),
          col=colorHPDI, border=NA )
    }

    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                cex = pt.cex/cex, ...)
        }
    }
    axis(side=1, at=c(0,1,2))
    box()
    #title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}


disc_list <- list(
                    post19$gamma[,1],
                    post19$gamma[,2],
                    post19$gamma[,3],
                    post19$gamma[,4],
                    post19$gamma[,5],
                    post19$gamma[,6],
                    post19$gamma[,7],
                    post19$gamma[,8],
                    post19$gamma[,9],
                    post19$gamma[,10],
                    post19$gamma[,11],
                    post19$gamma[,12],
                    post19$gamma[,13],
                    post19$gamma[,14]
                  )
names(disc_list) <- quest_names

str(disc_list)
disc_list <- rev(disc_list) #reverse order for plotting
str(disc_list)

#FigB2: dotplot of discrimination estimates for m19
pdf(file="./FigB2_disc_dens.pdf", 
height=5, width=6.5)
par(mfrow=c(1,1))

denschart1( disc_list , adjust=1 , color="black",
          colorHPDI=grey(0.45),
          HPDI=0.9,
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(disc_list))-4, 2.5) 
 )

title(xlab="Discrimination (Slope)", cex.lab=0.8, adj=0.95,
      line=2.3)

graphics.off()



#######################Fig B3: plot parameter posterior estimates for m19

#Richard's density plot function
denschart2 <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), 
    color = "grey20", colorHPDI ="grey60", HPDI=0.9, 
    gcolor = par("fg"), lcolor = "gray", xlim = range(unlist(x)), 
    main = NULL, xlab = NULL, ylab = NULL, height=0.7 , border=NA, adjust=1, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.list(x)) 
        stop("'x' must be a list of vectors or matrices")
    n <- length(x)
    glabels <- NULL
    if (is.list(x)) {
        if (is.null(labels)) 
            labels <- names(x)
        if (is.null(labels)) 
            labels <- as.character(1L:n)
        labels <- rep_len(labels, n)
        #if (is.null(groups)) 
        #    groups <- col(x, as.factor = TRUE)
        #glabels <- levels(groups)
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1L:n
        y <- o #c(1,2,3, 4.5,5.5,6.5, 8,9,10) #o    #space out in groups of three
        ylim <- c(0.2, max(y)+0.7) #n + 1)
    }
    else {
        # sub-groups, so need more rows
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        #mtext(labs, side = 2, line = -1 #0.4, #loffset,           #### y-labels
        #      at = y, adj = 1, 
        #      col = "black", las = 2, cex = cex, ...)
        text(labels=labs, x=-6, y=y, pos=2, adj=1)
    }

    abline(v=0, lty=2, lwd=0.75) #dotted vertical line at 0

    #abline(h = y, lty = "dotted", col = lcolor)
    #points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)

    # draw densities at each y offset
    for ( i in 1:n ) {
        a <- density( x[[i]] , adjust=adjust )
        a$y <- a$y/max(a$y) * height + y[i] - 0.3
        polygon( a$x , a$y , col=color , border=border )
        Cuts <- HPDI( x[[i]] , HPDI)
        XX <- a$x[which(a$x > Cuts[1] & a$x < Cuts[2])]
        YY <- a$y[which(a$x > Cuts[1] & a$x < Cuts[2])]
        polygon( c(min(XX), XX, max(XX)), c(min(a$y), YY, min(a$y)),
          col=colorHPDI, border=NA )
    }

    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                cex = pt.cex/cex, ...)
        }
    }
    axis(side=1, at=c(-2,0,2,4,6))
    box()
    #title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}


par_list <- list(
                    post19$b0,
                    post19$bMachi,
                    post19$bComMest,
                    post19$bMachiComMest,
                    post19$bEmpMest,
                    post19$bMachiEmpMest,
                    post19$bSchool,
                    post19$bMachiSchool
                  )
names(par_list) <- c(
                      "Intercept", "Matsi",
                      "Commerce", "Matsi X Commerce",
                      "Labor", "Matsi X Labor",
                      "Education", "Matsi X Education"
                      )
str(par_list)
par_list <- rev(par_list) #reverse order for plotting
str(par_list)

#FigB3: dotplot of parameter estimates
pdf(file="./FigB3_params_dens.pdf", 
height=5, width=6)
par(mfrow=c(1,1))

denschart2( par_list , adjust=1 , color="black",
          colorHPDI=grey(0.45),
          HPDI=0.9,
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(par_list))-7.5, max(unlist(par_list)) ) 
 )
title(xlab="Effect Size on Latent Axis", cex.lab=0.8, adj=0.95,
      line=2.3)

graphics.off()


##########################Fig 4: Calculate contrasts for m19

post0 <- post19 #posterior estimates of complete model

best_str_machi_none <- (post0$b0) +     #mean location of Matsi with no inter-ethnic experience
                          (post0$bMachi)

best_str_machi_just_sch <- (post0$b0) + #mean location of Matsi with just school experience
                          (post0$bMachi) +
                          (post0$bSchool) + 
                          (post0$bMachiSchool)

best_mest_sch <- (post0$b0) +   #mean location of Mestizo with just school experience
                (post0$bSchool)

#Matsi without inter-ethnic experience minus Matsi with just school experience
best_str_machi_wo_sch_machi_w_sch <- best_str_machi_none - best_str_machi_just_sch

#Matsi without inter-ethnic experience minus Mestizo with just school experience
best_str_machi_wo_sch_mest <- best_str_machi_none - best_mest_sch

#Matsi with just school experience minus Mestizo with just school experience
best_str_machi_w_sch_mest <- best_str_machi_just_sch - best_mest_sch



best_str_machi_just_emp <- (post0$b0) + #mean location of Matsi with just employment experience
                          (post0$bMachi) +
                          (post0$bEmpMest) + 
                          (post0$bMachiEmpMest)

best_mest_emp <- (post0$b0) + #mean location of Mestizo with just employment experience
                (post0$bEmpMest)

#Matsi without inter-ethnic experience minus Matsi with just employment experience
best_str_machi_wo_emp_machi_w_emp <- best_str_machi_none - best_str_machi_just_emp

#Matsi without inter-ethnic experience minus Mestizo with just employment experience
best_str_machi_wo_emp_mest <- best_str_machi_none - best_mest_emp

#Matsi with just employment experience minus Mestizo with just employment experience
best_str_machi_w_emp_mest <- best_str_machi_just_emp - best_mest_emp



best_str_machi_just_com <- (post0$b0) + #mean location of Matsi with just commerce experience
                          (post0$bMachi) +
                          (post0$bComMest) + 
                          (post0$bMachiComMest)

best_mest_com <- (post0$b0) + #mean location of Mestizo with just commerce experience
                (post0$bComMest)

#Matsi without inter-ethnic experience minus Matsi with just commerce experience
best_str_machi_wo_com_machi_w_com <- best_str_machi_none - best_str_machi_just_com

#Matsi without inter-ethnic experience minus Mestizo with just commerce experience
best_str_machi_wo_com_mest <- best_str_machi_none - best_mest_com

#Matsi with just commerce experience minus Mestizo with just commerce experience
best_str_machi_w_com_mest <- best_str_machi_just_com - best_mest_com


cred <- 0.90
contrasts_sch <- NULL       #table of school contrasts and HPDI limits
contrasts_sch <- cbind( 
                      c( mean(best_str_machi_wo_sch_machi_w_sch),
                          mean(best_str_machi_wo_sch_mest),
                          mean(best_str_machi_w_sch_mest)
                          ),
                      rbind(t(as.vector(HPDI(best_str_machi_wo_sch_machi_w_sch, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_wo_sch_mest, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_w_sch_mest, prob=cred)))
                          )
                      )

contrasts_sch <- cbind(c(1:3), contrasts_sch) ######## change according to # of parameters
rownames(contrasts_sch) <- c("Matsi w/o exp - Matsi w/ edu", 
                         "Matsi w/o exp - Mest w/ edu",
                         "Matsi w/ edu - Mest w/ edu")
colnames(contrasts_sch) <- c("index", "m19",
                             "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_sch, digit=2)
contrasts_sch <- contrasts_sch[order(-contrasts_sch[,"index"]),]


contrasts_emp <- NULL       #table of employment contrasts and HPDI limits
contrasts_emp <- cbind( 
                      c( mean(best_str_machi_wo_emp_machi_w_emp),
                          mean(best_str_machi_wo_emp_mest),
                          mean(best_str_machi_w_emp_mest)
                          ),
                      rbind(t(as.vector(HPDI(best_str_machi_wo_emp_machi_w_emp, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_wo_emp_mest, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_w_emp_mest, prob=cred)))
                          )
                      )

contrasts_emp <- cbind(c(1:3), contrasts_emp) ######## change according to # of parameters
rownames(contrasts_emp) <- c("Matsi w/o exp - Matsi w/ lab", 
                         "Matsi w/o exp - Mest w/ lab",
                         "Matsi w/ lab - Mest w/ lab")
colnames(contrasts_emp) <- c("index", "m19",
                             "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_emp, digit=2)
contrasts_emp <- contrasts_emp[order(-contrasts_emp[,"index"]),]


contrasts_com <- NULL         #table of commerce contrasts and HPDI limits
contrasts_com <- cbind( 
                      c( mean(best_str_machi_wo_com_machi_w_com),
                          mean(best_str_machi_wo_com_mest),
                          mean(best_str_machi_w_com_mest)
                          ),
                      rbind(t(as.vector(HPDI(best_str_machi_wo_com_machi_w_com, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_wo_com_mest, prob=cred))),
                          t(as.vector(HPDI(best_str_machi_w_com_mest, prob=cred)))
                          )
                      )

contrasts_com <- cbind(c(1:3), contrasts_com) ######## change according to # of parameters
rownames(contrasts_com) <- c("Matsi w/o exp - Matsi w/ com", 
                         "Matsi w/o exp - Mest w/ com",
                         "Matsi w/ com - Mest w/ com")
colnames(contrasts_com) <- c("index", "m19",
                              "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_com, digit=2)
contrasts_com <- contrasts_com[order(-contrasts_com[,"index"]),]


#combine sub-tables of contrasts
contrasts_all <- rbind(contrasts_com, contrasts_emp, contrasts_sch)
contrasts_all[,"index"] <- c( length(contrasts_all[,"index"]):1 )
print(contrasts_all, digit=3)
contrasts_all <- contrasts_all[c(1,4,7,2,5,8,3,6,9),] #new order



#Richard's density plot function
denschart3 <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), 
    color = "grey20", colorHPDI ="grey60", HPDI=0.9, 
    gcolor = par("fg"), lcolor = "gray", xlim = range(unlist(x)), 
    main = NULL, xlab = NULL, ylab = NULL, height=0.7 , border=NA, adjust=1, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (!is.list(x)) 
        stop("'x' must be a list of vectors or matrices")
    n <- length(x)
    glabels <- NULL
    if (is.list(x)) {
        if (is.null(labels)) 
            labels <- names(x)
        if (is.null(labels)) 
            labels <- as.character(1L:n)
        labels <- rep_len(labels, n)
        #if (is.null(groups)) 
        #    groups <- col(x, as.factor = TRUE)
        #glabels <- levels(groups)
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1L:n
        y <- c(1,2,3, 4.5,5.5,6.5, 8,9,10) #o    #space out in groups of three
        ylim <- c(0.2, max(y)+0.7) #n + 1)
    }
    else {
        # sub-groups, so need more rows
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        #mtext(labs, side = 2, line = -1 #0.4, #loffset,           #### y-labels
        #      at = y, adj = 1, 
        #      col = "black", las = 2, cex = cex, ...)
        text(labels=labs, x=-2.2, y=y, pos=2, adj=1)
    }

    abline(v=0, lty=2, lwd=0.75) #dotted vertical line at 0

    #abline(h = y, lty = "dotted", col = lcolor)
    #points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)

    # draw densities at each y offset
    for ( i in 1:n ) {
        a <- density( x[[i]] , adjust=adjust )
        a$y <- a$y/max(a$y) * height + y[i] - 0.3
        polygon( a$x , a$y , col=color , border=border )
        Cuts <- HPDI( x[[i]] , HPDI)
        XX <- a$x[which(a$x > Cuts[1] & a$x < Cuts[2])]
        YY <- a$y[which(a$x > Cuts[1] & a$x < Cuts[2])]
        polygon( c(min(XX), XX, max(XX)), c(min(a$y), YY, min(a$y)),
          col=colorHPDI, border=NA )
    }

    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                cex = pt.cex/cex, ...)
        }
    }
    axis(side=1, at=c(-2,0,2,4,6))
    box()
    #title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}


cont_list <- list(
                    best_str_machi_wo_com_mest,
                    best_str_machi_wo_emp_mest,
                    best_str_machi_wo_sch_mest,

                    best_str_machi_wo_com_machi_w_com,
                    best_str_machi_wo_emp_machi_w_emp,
                    best_str_machi_wo_sch_machi_w_sch,
                    
                    best_str_machi_w_com_mest,
                    best_str_machi_w_emp_mest,
                    best_str_machi_w_sch_mest
                  )
names(cont_list) <- c(
                      "Matsi w/o exp - Mest w/ com",
                      "Matsi w/o exp - Mest w/ lab",
                      "Matsi w/o exp - Mest w/ edu",

                      "Matsi w/o exp - w/ com",
                      "Matsi w/o exp - w/ lab",
                      "Matsi w/o exp - w/ edu",
                      
                      "Matsi w/ com - Mest w/ com",
                      "Matsi w/ lab - Mest w/ lab",
                      "Matsi w/ edu - Mest w/ edu"
                      )
str(cont_list)
cont_list <- rev(cont_list) #reverse order for plotting
str(cont_list)

#Fig4: dotplot of contrasts
pdf(file="./Fig4_contrasts_dens.pdf", 
height=5, width=6.75)
par(mfrow=c(1,1))

denschart3( cont_list , adjust=1 , color="black",
          colorHPDI=grey(0.45),
          HPDI=0.9,
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(cont_list))-7.5, 7)#max(unlist(cont_list)) ) 
 )
#abline(v=0)
lines( list( x=c(min(unlist(cont_list))-8.15,7.7), y=c(3.75,3.75) ) )
lines( list( x=c(min(unlist(cont_list))-8.15,7.7), y=c(7.25,7.25) ) )
#box(bty="o")
title(#main="Contrast Estimates",
      xlab="Contrast on Latent Axis", cex.lab=0.8, adj=0.87,
      line=2.3)
text(x=6.5, y=10.3, cex=1.3, labels="A")
text(x=6.5, y=6.8, cex=1.3, labels="B")
text(x=6.5, y=3.3, cex=1.3, labels="C")

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

stargazer(tab_mat, summary=FALSE, rownames=FALSE, type="latex",
          out="./table.tex")


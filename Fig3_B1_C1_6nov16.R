#Fits IRT model m1 and makes Fig 3 and appendix Figs B1 and C1

#Before beginning, download the data file Manu_interviews_31oct16.csv into
#a local directory on your machine. I'd recommend creating a new folder and
#putting the data file in there. Then open R from within that local directory
#and run this script. The figures will be created as pdfs in that same
#directory.

rm (list = ls(all=TRUE))
library(rstan)
library(rethinking)
library(lattice)


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



################## Appendix Fig B1
#plot logistic response functions for a given question and latent trait estimates of indivs

questions <- quest_names #from above

#look at all plots
pdf(file="./FigB1_logistic_m1.pdf", 
  height=10, width=12)
par(mfrow=c(4,4), oma=c(5,7,5,5), mar=c(3,2,2,2))

gamma <- post1$gamma

for (Quest in 1:K) {
  plot( t( d[which(d$QID == Quest), "Response"] ) ~
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



############################Fig 3 and Appendix Fig C1

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
#labs <- c(1:J)
ID <- unique(d[,"ID"])
newID <- unique(d[,"newID"])


#data matrix for plotting
all_1 <- as.data.frame( cbind(newID, ID, abils,
                            machi,
                            School.mest, WorkMest, ComMest
                            ) ) 
#sort by abils
all <- all_1[order(all_1$abils),]
all <- cbind(all, c(1:J))
names(all)
colnames(all)[8] <- c("ord")
print(all, digits=2)


######################Appendix Fig C1: individual location highlights
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



###############################Fig 3: individual location by interaction experience
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


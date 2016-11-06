#Fits IRT model m19 and makes Fig 4 and appendix Figs B2 and B3


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



########################Fig B2: plot discrimination posterior estimates for m19

#Richard's density plot function
denschart <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), color = "gray", 
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

denschart( disc_list , adjust=1 , color= grey(0.25), #"slateblue" ,
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(disc_list))-4, 2.5) 
 )

title(xlab="Slope", cex.lab=0.8, adj=0.85,
      line=2.3)

graphics.off()


#######################Fig B3: plot parameter posterior estimates for m19

#Richard's density plot function
denschart <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), color = "gray", 
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

denschart( par_list , adjust=1 , color= grey(0.25),
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(par_list))-7.5, max(unlist(par_list)) ) 
 )
title(xlab="Effect Size on Latent Axis", cex.lab=0.8, adj=0.95,
      line=2.3)

graphics.off()


##########################Fig 4: Calculate contrasts for m19

post0 <- post19 #posterior estimates of best-fitting model

best_str_machi_none <- (post0$b0) + 
                          (post0$bMachi)

best_str_machi_just_sch <- (post0$b0) + 
                          (post0$bMachi) +
                          (post0$bSchool) + 
                          (post0$bMachiSchool)

best_mest_sch <- (post0$b0) + 
                (post0$bSchool)

best_str_machi_wo_sch_machi_w_sch <- best_str_machi_none - best_str_machi_just_sch
best_str_machi_wo_sch_mest <- best_str_machi_none - best_mest_sch
best_str_machi_w_sch_mest <- best_str_machi_just_sch - best_mest_sch


best_str_machi_none <- (post0$b0) + 
                          (post0$bMachi)

best_str_machi_just_emp <- (post0$b0) + 
                          (post0$bMachi) +
                          (post0$bEmpMest) + 
                          (post0$bMachiEmpMest)

best_mest_emp <- (post0$b0) + 
                (post0$bEmpMest)

best_str_machi_wo_emp_machi_w_emp <- best_str_machi_none - best_str_machi_just_emp
best_str_machi_wo_emp_mest <- best_str_machi_none - best_mest_emp
best_str_machi_w_emp_mest <- best_str_machi_just_emp - best_mest_emp


best_str_machi_none <- (post0$b0) + 
                          (post0$bMachi)

best_str_machi_just_com <- (post0$b0) + 
                          (post0$bMachi) +
                          (post0$bComMest) + 
                          (post0$bMachiComMest)

best_mest_com <- (post0$b0) + 
                (post0$bComMest)

best_str_machi_wo_com_machi_w_com <- best_str_machi_none - best_str_machi_just_com
best_str_machi_wo_com_mest <- best_str_machi_none - best_mest_com
best_str_machi_w_com_mest <- best_str_machi_just_com - best_mest_com


cred <- 0.90
contrasts_sch <- NULL
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
rownames(contrasts_sch) <- c("Matsi w/o edu - Matsi w/ edu", 
                         "Matsi w/o edu - Mestizo",
                         "Matsi w/ edu - Mestizo")
colnames(contrasts_sch) <- c("index", "m19",
                             "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_sch, digit=2)
contrasts_sch <- contrasts_sch[order(-contrasts_sch[,"index"]),]


contrasts_emp <- NULL
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
rownames(contrasts_emp) <- c("Matsi w/o lab - Matsi w/ lab", 
                         "Matsi w/o lab - Mestizo",
                         "Matsi w/ lab - Mestizo")
colnames(contrasts_emp) <- c("index", "m19",
                             "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_emp, digit=2)
contrasts_emp <- contrasts_emp[order(-contrasts_emp[,"index"]),]


contrasts_com <- NULL
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
rownames(contrasts_com) <- c("Matsi w/o com - Matsi w/ com", 
                         "Matsi w/o com - Mestizo",
                         "Matsi w/ com - Mestizo")
colnames(contrasts_com) <- c("index", "m19",
                              "strm19HPDI_low", "strm19HPDI_high")
print(contrasts_com, digit=2)
contrasts_com <- contrasts_com[order(-contrasts_com[,"index"]),]



contrasts_all <- rbind(contrasts_com, contrasts_emp, contrasts_sch)
contrasts_all[,"index"] <- c( length(contrasts_all[,"index"]):1 )
print(contrasts_all, digit=3)
contrasts_all <- contrasts_all[c(1,4,7,2,5,8,3,6,9),] #new order



#Richard's density plot function
denschart <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, bg = par("bg"), color = "gray", 
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
                      "Matsi w/o com - Mestizo",
                      "Matsi w/o lab - Mestizo",
                      "Matsi w/o edu - Mestizo",

                      "Matsi w/o com - w/ com",
                      "Matsi w/o lab - w/ lab",
                      "Matsi w/o edu - w/ edu",
                      
                      "Matsi w/ com - Mestizo",
                      "Matsi w/ lab - Mestizo",
                      "Matsi w/ edu - Mestizo"
                      )
str(cont_list)
cont_list <- rev(cont_list) #reverse order for plotting
str(cont_list)

#Fig4: dotplot of contrasts
pdf(file="./Fig4_contrasts_dens.pdf", 
height=5, width=6)
par(mfrow=c(1,1))

denschart( cont_list , adjust=1 , color= grey(0.25), #"slateblue" ,
          border=NA, yaxt="n",
          cex=0.8, height=0.7,
          xlim=range( min(unlist(cont_list))-7.5, 7)#max(unlist(cont_list)) ) 
 )
#abline(v=0)
lines( list( x=c(min(unlist(cont_list))-8.15,7.7), y=c(3.75,3.75) ) )
lines( list( x=c(min(unlist(cont_list))-8.15,7.7), y=c(7.25,7.25) ) )
#box(bty="o")
title(#main="Contrast Estimates",
      xlab="Contrast on Latent Axis", cex.lab=0.8, adj=0.9,
      line=2.3)
text(x=6.5, y=10.3, cex=1.3, labels="A")
text(x=6.5, y=6.8, cex=1.3, labels="B")
text(x=6.5, y=3.3, cex=1.3, labels="C")

graphics.off()



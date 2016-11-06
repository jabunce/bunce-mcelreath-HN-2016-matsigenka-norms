#Calculates response raw proportions and makes Fig 2 


#Before beginning, download the data file Manu_interviews_31oct16.csv into
#a local directory on your machine. I'd recommend creating a new folder and
#putting the data file in there. Then open R from within that local directory
#and run this script. The figures will be created as pdfs in that same
#directory.

rm (list = ls(all=TRUE))


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

J <- length(unique(d$newID))  #number of people
K <- length(unique(d$QID))    #number of questions
N <- nrow(d)      #total number of responses
jj <- d$newID     #vector of person IDs
kk <- d$QID     #vector of question IDs
y <- d$Response     #vector of responses


###########################raw proportions, Figure 2

Quest <- c(1:K) #question to plot

num_mach_yes <- 0
num_mach_no <- 0
num_mest_yes <- 0
num_mest_no <- 0

props <- as.data.frame( cbind(rep(0, K),
                              rep(0, K),
                              rep(0, K)) )
colnames(props)[c(1:3)] <- c("Quest", "mach_yes",
                                    "mest_yes")

for (l in 1:K) {
  resp <- d[which(d$QID == Quest[l]),
              c("Response", "Machi", "newID")] #responses to question
  for (i in 1:length(resp$newID)) {
    if ((resp$Machi[i] == 1) && (resp$Response[i] == 1)) {
      num_mach_yes <- num_mach_yes + 1
    } else if ((resp$Machi[i] == 0) && (resp$Response[i] == 1)) {
       num_mest_yes <- num_mest_yes + 1
    } else if ((resp$Machi[i] == 1) && (resp$Response[i] == 0)) {
       num_mach_no <- num_mach_no + 1
   } else if ((resp$Machi[i] == 0) && (resp$Response[i] == 0)) {
      num_mest_no <- num_mest_no + 1
    }#if
  }#for i

  props[Quest[l],] <- c(Quest[l],
           num_mach_yes/(num_mach_yes + num_mach_no),
           num_mest_yes/(num_mest_yes + num_mest_no))
  
  #reset temp variables
  num_mach_yes <- 0
  num_mach_no <- 0
  num_mest_yes <- 0
  num_mest_no <- 0
  resp <- NULL
} #for l


quest_names <- c("1.wife hunts/works", "2.daughter babysits", "3.not wear dead hat",
               "4.hit students", "5.no questions", "6.post flu",
               "7.pot to needy", "8.good nonbaptized no heaven", "9.bad baptized heaven",
               "10.stop work to visit", "11.cheap mean store", "12.xcousin marriage",
               "13.arranged marriage", "14.laborer party")

props <- cbind(props, quest_names)
print(props, digits=2)

props <- props[order(props$mach_yes),] #reorder by machi percentage


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
  main="Questions",
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


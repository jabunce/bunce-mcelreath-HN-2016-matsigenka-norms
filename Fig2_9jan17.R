#Calculates response raw proportions and makes Fig 2


#Before beginning, download the data file Manu_interviews_31oct16.csv into
#a local directory on your machine. I'd recommend creating a new folder and
#putting the data file in there. Then open R from within that local directory
#and run this script. The figures will be created as pdfs in that same
#directory.


rm (list = ls(all=TRUE))
library(lattice)

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



###########################raw proportions, Figure 2

props1 <- data.frame(
  Quest=c(1:K) ,
  Matsi=table( d$QID[which(d$Machi==1)], d$ResponseFlipped[which(d$Machi==1)] )[,2] /
    rowSums( table( d$QID[which(d$Machi==1)], d$ResponseFlipped[which(d$Machi==1)]) ),

  Mest=table( d$QID[which(d$Machi==0)], d$ResponseFlipped[which(d$Machi==0)] )[,2] /
    rowSums( table( d$QID[which(d$Machi==0)], d$ResponseFlipped[which(d$Machi==0)]) )
)

colnames(props1)[c(1:3)] <- c("Quest", "mach_yes", "mest_yes")


props1 <- cbind(props1, quest_names)
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

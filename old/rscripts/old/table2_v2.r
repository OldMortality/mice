# 
# Verifies the results in table2-v2.rmd
#

## get the data 
source("~/Documents/VITAMIN_D/rscripts/createlongdata3.r")
d <- longdata


##
## Mothers
##

## vit d deficiency
# row 1
table(d[which(d$season==0 & d$period=='AN1'),"ohd"])
table(d[which(d$season==0 & d$period=='AN1'),"ohd"]<50)
table(d[which(d$season==0 & d$period=='AN2'),"ohd"]<50)
table(d[which(d$season==0 & d$period=='AN3'),"ohd"]<50)
r <- table(d[which(d$season==0 & d$period=='PN1'),"ohd"]<50)
n <- length(which(d$season==0 & d$period=='PN1') )
r / n
r <- table(d[which(d$season==0 & d$period=='PN2'),"ohd"]<50)
n <- length(which(d$season==0 & d$period=='PN2') )
r / n

# row 2
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='AN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='AN2') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN3'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==1 & d$period=='AN3') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='PN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==1 & d$period=='PN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='PN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='PN2') )
r / n

# row 2
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='AN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='AN2') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='AN3'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==1 & d$period=='AN3') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='PN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==1 & d$period=='PN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==1 & d$period=='PN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==1 & d$period=='PN2') )
r / n

# row 3
r <- table(d[which(!is.na(d$ohd) & d$season==2 & d$period=='AN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==2 & d$period=='AN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==2 & d$period=='AN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==2 & d$period=='AN2') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==2 & d$period=='AN3'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==2 & d$period=='AN3') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==2 & d$period=='PN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==2 & d$period=='PN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==2 & d$period=='PN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==2 & d$period=='PN2') )
r / n

# row 4
r <- table(d[which(!is.na(d$ohd) & d$season==3 & d$period=='AN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==3 & d$period=='AN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==3 & d$period=='AN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==3 & d$period=='AN2') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==3 & d$period=='AN3'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==3 & d$period=='AN3') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==3 & d$period=='PN1'),"ohd"]<50)
n <- length(which(!is.na(d$ohd)& d$season==3 & d$period=='PN1') )
r / n
r <- table(d[which(!is.na(d$ohd) & d$season==3 & d$period=='PN2'),"ohd"]<50)
n <- length(which(!is.na(d$ohd) & d$season==3 & d$period=='PN2') )
r / n

## PTH
options(digits=4)

# row 0
# row 1
p <- d[which(d$period=='AN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which(d$period=='AN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which(d$period=='AN3'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which(d$period=='PN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which(d$period=='PN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))

# row 1
p <- d[which( d$season==0 & d$period=='AN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==0 & d$period=='AN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==0 & d$period=='AN3'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==0 & d$period=='PN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==0 & d$period=='PN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))

# row 2
p <- d[which( d$season==1 & d$period=='AN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==1 & d$period=='AN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==1 & d$period=='AN3'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==1 & d$period=='PN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==1 & d$period=='PN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))

# row 3
p <- d[which( d$season==2 & d$period=='AN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==2 & d$period=='AN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==2 & d$period=='AN3'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==2 & d$period=='PN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==2 & d$period=='PN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))

# row 4
p <- d[which( d$season==3 & d$period=='AN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==3 & d$period=='AN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==3 & d$period=='AN3'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==3 & d$period=='PN1'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- d[which( d$season==3 & d$period=='PN2'),"pth"]
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))


##
## baby
##
d <- read.csv("~/Documents/VITAMIN_D/data/mdl1.csv",na.strings=c("","NA"),stringsAsFactors=FALSE,
                   header=TRUE)

# row 0
mean(as.numeric(d$Cord.PTH),na.rm=T)
sqrt(var(as.numeric(d$Cord.PTH),na.rm=T))
mean(as.numeric(d$Mth.5.baby.PTH),na.rm=T)
sqrt(var(as.numeric(d$Mth.5.baby.PTH),na.rm=T))

# we can't use season.of.delivery, because we have defined seasons differently
# row 1-4 pth
p <- as.numeric(d[which( d$Month.of.delivery==1 | d$Month.of.delivery==2 | d$Month.of.delivery==3),"Cord.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$Month.of.delivery==4 | d$Month.of.delivery==5 | d$Month.of.delivery==6),"Cord.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$Month.of.delivery==7 | d$Month.of.delivery==8 | d$Month.of.delivery==9),"Cord.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$Month.of.delivery==10 | d$Month.of.delivery==11 | d$Month.of.delivery==12),"Cord.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))

#vit-d
# row 0
mean(as.numeric(d$Cord.Baby.Total.25OHD),na.rm=T)
sqrt(var(as.numeric(d$Cord.Baby.Total.25OHD),na.rm=T))
mean(as.numeric(d$Mth.5.baby.PTH),na.rm=T)
sqrt(var(as.numeric(d$Mth.5.baby.PTH),na.rm=T))

p <- as.numeric(d[which( d$Month.of.delivery==1 | d$Month.of.delivery==2 | d$Month.of.delivery==3),"Cord.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$Month.of.delivery==4 | d$Month.of.delivery==5 | d$Month.of.delivery==6),"Cord.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$Month.of.delivery==7 | d$Month.of.delivery==8 | d$Month.of.delivery==9),"Cord.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$Month.of.delivery==10 | d$Month.of.delivery==11 | d$Month.of.delivery==12),"Cord.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))

## final column PN2 baby
source('~/Documents/VITAMIN_D/rscripts/getseason.r')

to_be_removed <- which(d$Mth.5.Baby.Date == "Declined M5" | 
                         d$Mth.5.Baby.Date == "Declined  M5")

d <- d[-to_be_removed,]

d$x1 <- as.Date(d$Mth.5.Baby.Date,format="%d/%m/%y")
d$x2 <- as.POSIXlt(d$x1)


d$m5.baby.season <- NA
for (i in 1:dim(d)[1]) {
  d$m5.baby.season[i] <- as.character(get_season(d[i,"x2"]$mon))
  
}
table(d$m5.baby.season)
p <- as.numeric(d[which( d$m5.baby.season==0 ),"Mth.5.baby.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$m5.baby.season==1 ),"Mth.5.baby.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$m5.baby.season==2 ),"Mth.5.baby.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))
p <- as.numeric(d[which( d$m5.baby.season==3 ),"Mth.5.baby.PTH"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))


p <- as.numeric(d[,"Mth.5.Baby.Total.25OHD"])
mean(p,na.rm=T)
sqrt(var(p,na.rm=T))


p <- as.numeric(d[which( d$m5.baby.season=='0'),"Mth.5.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$m5.baby.season=='1'),"Mth.5.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$m5.baby.season=='2'),"Mth.5.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))
p <- as.numeric(d[which( d$m5.baby.season=='3'),"Mth.5.Baby.Total.25OHD"])
sum(p<50,na.rm=T)/length(which(!is.na(p)))


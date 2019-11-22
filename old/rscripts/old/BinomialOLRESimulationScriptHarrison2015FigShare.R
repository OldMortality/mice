library(lme4)
library(MuMIn)
#Required for MuMIn
options(na.action = "na.fail")   #  prevent fitting models to different datasets


# Function To Generate Data from Beta-Binomial Distribution	
	#This function generates an outcome variable from 2 covariates - number of flowers visited, and body size
	# The slope for 'flower vistis' is positive but varies in strength, the slope for previous breeding is negative and has a weak effect
	# Probabilties are derived from a linear predictor (B0 + b1*flowers + b2*previous breed)
	# Probs are then scaled to 'a' and 'b' for the beta distribution using the value for theta, and binomial counts are drawn from this distribution
data.betabinom<-function(n.pops,n.ind,theta,intercept.mu,intercept.sd,beta1slope,beta2slope,clutch){

			sample<- n.pops * n.ind	
			popid<-gl(n.pops,n.ind )
			
			#Draw Flower Visits
			#flowervisits<-sample(seq(1,5,1),sample,replace=T)
			flowervisits<-rpois(sample,1.5)
			
			#Draw Clutch Size (deprecated as all individuals now have the same 'C')
			#clutch<-sample(seq(1,4,1),sample,replace=T)
			#clutch<-rpois(sample,nclutch)


			#Body Size
			body<-rnorm(sample,20,3)
			bodyc<-body-mean(body)
			
			#Draw Random Intecepts
			alphas<-rnorm(n.pops,intercept.mu,intercept.sd)
			
			#Generate Proportion of Clutch Hatched - Overdispersed Data
			linpred.od<- alphas[popid] + beta1slope * flowervisits + beta2slope*body 
			probability.od<- plogis( linpred.od)
			a<-probability.od / theta
			b<- (1-probability.od)/theta
			n.hatched.od<- rbinom(sample,clutch,prob=rbeta(sample,a,b))
			
			dat<-data.frame(hatched=n.hatched.od,clutch=clutch,flowervisits=flowervisits,bodysize=body,popid=popid,obs=seq(sample),a=a,b=b,bodysize.c=bodyc)
			return(dat)
			}


# Function To Generate Data from  Overdispersed Binomial Distribution	
	#This function generates an outcome variable from 2 covariates - number of flowers visited, and body size
	# The slope for 'flower vistis' is positive but varies in strength, the slope for previous breeding is negative and has a weak effect
	# Probabilties are derived from a linear predictor (B0 + b1*flowers + b2*previous breed) but also add random extra-Binomial noise from a Normal distribution with mean 0 and a fixed SD

data.eps<-function(n.pops,n.ind,epsilon,intercept.mu,intercept.sd,beta1slope,beta2slope,clutch){

			sample<- n.pops * n.ind	
			popid<-gl(n.pops,n.ind )
			
			#Draw Flower Visits
			#flowervisits<-sample(seq(1,5,1),sample,replace=T)
			flowervisits<-rpois(sample,1.5)
			
			#Draw Clutch Size (deprecated as all individuals have fixed 'C')
			#clutch<-sample(seq(1,4,1),sample,replace=T)
			#clutch<-rpois(sample,nclutch)


			#Body Size
			body<-rnorm(sample,20,3)
			bodyc<-body-mean(body)
			
			#Draw Random Intecepts
			alphas<-rnorm(n.pops,intercept.mu,intercept.sd)
			
			#Generate Proportion of Clutch Hatched - Overdispersed Data
			linpred.od<- alphas[popid] + beta1slope * flowervisits + beta2slope*body + rnorm(sample,0,epsilon)
			probability.od<- plogis( linpred.od)
			n.hatched.od<- rbinom(sample,clutch,prob=probability.od)
			
			dat<-data.frame(hatched=n.hatched.od,clutch=clutch,flowervisits=flowervisits,bodysize=body,popid=popid,obs=seq(sample),bodysize.c=bodyc)
			return(dat)
			}

  
##########
#Parameters
##########

intercept.sample<--1
intercept.sd.sample<-0.5
n.pops.sample<-10
n.ind.sample<-20
nclutch<-5

#Slopes 
beta1.strong<-0.6
beta2<- -0.01

#Overdispersion - BetaBinomial Data
simtheta<-c(0.1,1,2)

#Overdispersion - Binomial Logit-Normal Data
simepsilon<-c(0.1,1.5,3)

#Sample Size - Levels of Random Effect (both kinds of data)
simpops<-c(3,5,20)

#Sample Size - Binomial Trials (both kinds of data)
simclutch<-c(2,4,10)




#####################################
# Simulations to Test Effect of i) Overdispersion Levels ii) Sample Size 
#####################################

datastore<-function(nrows,ncols){
	x<-list(slope1=matrix(NA,nrows,ncols),slope2=matrix(NA,nrows,ncols),interceptmean=matrix(NA,nrows,ncols),overdispersion=matrix(NA,nrows,ncols),olre=matrix(NA,nrows,ncols),popsd=matrix(NA,nrows,ncols))
	return(x)
}


################# Rep 1 Beta Binomial With Weak Slope and Both Effects in the Model
simrep<-1000


od.betabinom<-datastore(nrows=length(simtheta),ncols=simrep)
od.eps<-datastore(nrows=length(simtheta),ncols=simrep)

pops.betabinom<-datastore(nrows=length(simtheta),ncols=simrep)
pops.eps<-datastore(nrows=length(simtheta),ncols=simrep)

clutch.betabinom<-datastore(nrows=length(simtheta),ncols=simrep)
clutch.eps<-datastore(nrows=length(simtheta),ncols=simrep)


for(k in 1:3){
		for(j in 1:simrep){
		
		#Beta Binomial Overdispersion	
			b1<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=simtheta[k],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
			m1<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=b1)

	
			od.betabinom $slope1[k,j]<-fixef(m1)[2]
			od.betabinom $slope2[k,j]<-fixef(m1)[3]
			od.betabinom $interceptmean[k,j]<-fixef(m1)[1] - mean(b1$bodysize)*fixef(m1)[3]
			od.betabinom $olre[k,j]<-sqrt(VarCorr(m1)$obs)[1]
			od.betabinom $popsd[k,j]<-sqrt(VarCorr(m1)$popid)[1]
		
		#Epsilon Overdispersion
		e1<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=simepsilon[k],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
		m2<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=e1)

			od.eps$slope1[k,j]<-fixef(m2)[2]
			od.eps$slope2[k,j]<-fixef(m2)[3]
			od.eps $interceptmean[k,j]<-fixef(m2)[1] - mean(e1$bodysize)*fixef(m2)[3]
			od.eps$olre[k,j]<-sqrt(VarCorr(m2)$obs)[1]
			od.eps$popsd[k,j]<-sqrt(VarCorr(m2)$popid)[1]
		
		#Beta Binomial Pops
		b2<-data.betabinom(n.pops=simpops[k],n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
		m3<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=b2)

			pops.betabinom $slope1[k,j]<-fixef(m3)[2]
			pops.betabinom $slope2[k,j]<-fixef(m3)[3]
			pops.betabinom $interceptmean[k,j]<-fixef(m3)[1] - mean(b2$bodysize)*fixef(m3)[3]
			pops.betabinom $olre[k,j]<-sqrt(VarCorr(m3)$obs)[1]
			pops.betabinom $popsd[k,j]<-sqrt(VarCorr(m3)$popid)[1]

		#Epsilon Pops
		e2<-data.eps(n.pops=simpops[k],n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
		m4<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=e2)
	
			pops.eps$slope1[k,j]<-fixef(m4)[2]
			pops.eps $slope2[k,j]<-fixef(m4)[3]
			pops.eps $interceptmean[k,j]<-fixef(m4)[1] - mean(e2$bodysize)*fixef(m4)[3]
			pops.eps $olre[k,j]<-sqrt(VarCorr(m4)$obs)[1]
			pops.eps $popsd[k,j]<-sqrt(VarCorr(m4)$popid)[1]
		
		#Beta Binomial Trials
		b3<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[k])
		m5<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=b3)

			clutch.betabinom $slope1[k,j]<-fixef(m5)[2]
			clutch.betabinom $slope2[k,j]<-fixef(m5)[3]
			clutch.betabinom $interceptmean[k,j]<-fixef(m5)[1] - mean(b3$bodysize)*fixef(m5)[3]
			clutch.betabinom $olre[k,j]<-sqrt(VarCorr(m5)$obs)[1]
			clutch.betabinom $popsd[k,j]<-sqrt(VarCorr(m5)$popid)[1]

		#Epsilon Trials
		e3<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[k])
		m6<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=e3)

			clutch.eps$slope1[k,j]<-fixef(m6)[2]
			clutch.eps $slope2[k,j]<-fixef(m6)[3]
			clutch.eps $interceptmean[k,j]<-fixef(m6)[1] - mean(e3$bodysize)*fixef(m6)[3]
			clutch.eps $olre[k,j]<-sqrt(VarCorr(m6)$obs)[1]
			clutch.eps $popsd[k,j]<-sqrt(VarCorr(m6)$popid)[1]
						
					}
				}


######################
#
#	Simulation Means
#
######################

odbmean<-lapply(od.betabinom,function(x){apply(x,1,mean)})
odemean<-lapply(od.eps,function(x){apply(x,1,mean)})
odb95<-lapply(od.betabinom,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})
ode95<-lapply(od.eps,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})

pbmean<-lapply(pops.betabinom,function(x){apply(x,1,mean)})
pemean<-lapply(pops.eps,function(x){apply(x,1,mean)})
pb95<-lapply(pops.betabinom,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})
pe95<-lapply(pops.eps,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})

cbmean<-lapply(clutch.betabinom,function(x){apply(x,1,mean)})
cemean<-lapply(clutch.eps,function(x){apply(x,1,mean)})
cb95<-lapply(clutch.betabinom,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})
ce95<-lapply(clutch.eps,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})



	#Proportion of Times Wrong Slope on Beta2
	odbprop<-apply(od.betabinom$slope2,1,function(x){mean(x>0)})
	odeprop<-apply(od.eps$slope2,1,function(x){mean(x>0)})
	odbprop95<-apply(replicate(1000,apply(apply(od.betabinom$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	odeprop95<-apply(replicate(1000,apply(apply(od.eps$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	
	pbprop<-apply(pops.betabinom$slope2,1,function(x){mean(x>0)})
	peprop<-apply(pops.eps$slope2,1,function(x){mean(x>0)})
	pbprop95<-apply(replicate(1000,apply(apply(pops.betabinom$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	peprop95<-apply(replicate(1000,apply(apply(pops.eps$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	
	cbprop<-apply(clutch.betabinom$slope2,1,function(x){mean(x>0)})
	ceprop<-apply(clutch.eps$slope2,1,function(x){mean(x>0)})
	cbprop95<-apply(replicate(1000,apply(apply(clutch.betabinom$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	ceprop95<-apply(replicate(1000,apply(apply(clutch.eps$slope2,1,function(x){sample(x,ncol(od.betabinom$slope2),replace=T)}),2,function(y){mean(y>0)})),1,function(z){quantile(z,c(0.025,0.975))})
	
	pmatrix<-pmatrixl95<-pmatrixu95<-matrix(0,nrow=2,ncol=9)
	pmatrix[1,1:3]<-odbprop;pmatrix[2,1:3]<-odeprop		
	pmatrix[1,4:6]<-pbprop;pmatrix[2,4:6]<-peprop	
	pmatrix[1,7:9]<-cbprop;pmatrix[2,7:9]<-ceprop	
	dimnames(pmatrix)=list(c("betabinomial","OD binomial"),c("od1","od2","od3","re1","re2","re3","clutch1","clutch2","clutch3"))
	
	pmatrixl95[1,1:3]<-odbprop95[1,1:3];pmatrixl95[2,1:3]<-odeprop95[1,1:3]	
	pmatrixl95[1,4:6]<-pbprop95[1,1:3];pmatrixl95[2,4:6]<-peprop95[1,1:3]	
	pmatrixl95[1,7:9]<-cbprop95[1,1:3];pmatrixl95[2,7:9]<-ceprop95[1,1:3]		
	
	pmatrixu95[1,1:3]<-odbprop95[2,1:3];pmatrixu95[2,1:3]<-odeprop95[2,1:3]	
	pmatrixu95[1,4:6]<-pbprop95[2,1:3];pmatrixu95[2,4:6]<-peprop95[2,1:3]	
	pmatrixu95[1,7:9]<-cbprop95[2,1:3];pmatrixu95[2,7:9]<-ceprop95[2,1:3]		
	



######################
#
#	Separate Trial for 3 pops but holding sample size at ~ 200
#
######################

pops.betabinom2<-datastore(nrows=1,ncols=simrep)
pops.eps2<-datastore(nrows=1,ncols=simrep)


for(j in 1:simrep){
	#Beta Binomial Pops
		b2<-data.betabinom(n.pops=3,n.ind=67,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
		m3<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=b2)

			pops.betabinom2 $slope1[1,j]<-fixef(m3)[2]
			pops.betabinom2 $slope2[1,j]<-fixef(m3)[3]
			pops.betabinom2 $interceptmean[1,j]<-fixef(m3)[1] - mean(b2$bodysize)*fixef(m3)[3]
			pops.betabinom2 $olre[1,j]<-sqrt(VarCorr(m3)$obs)[1]
			pops.betabinom2 $popsd[1,j]<-sqrt(VarCorr(m3)$popid)[1]
			
		#Epsilon Pops
		e2<-data.eps(n.pops=3,n.ind=67,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
		m4<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=e2)
	
			pops.eps2$slope1[1,j]<-fixef(m4)[2]
			pops.eps2 $slope2[1,j]<-fixef(m4)[3]
			pops.eps2 $interceptmean[1,j]<-fixef(m4)[1] - mean(e2$bodysize)*fixef(m4)[3]
			pops.eps2 $olre[1,j]<-sqrt(VarCorr(m4)$obs)[1]
			pops.eps2 $popsd[1,j]<-sqrt(VarCorr(m4)$popid)[1]

			
}

pbmean2<-lapply(pops.betabinom2,function(x){apply(x,1,mean)})
pb95.2<-lapply(pops.betabinom2,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})

pemean<-lapply(pops.eps2,function(x){apply(x,1,mean)})
pe95<-lapply(pops.eps2,function(x){apply(x,1,function(y){quantile(y,c(0.025,0.975),na.rm=T)})})

######################
#
#	SPAMM fits
#
######################	

library('spaMM')
library(glmmADMB)
spamm.beta<-spamm.eps<-matrix(NA,ncol=4,nrow=1000)

for(k in 1:1000){
b2<-data.betabinom(n.pops=10,n.ind=20,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
e2<-data.eps(n.pops=10,n.ind=20,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)

s1<-HLfit(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family=binomial(),rand.family=Beta(),HLmethod="HL(0,0)",data=b2)
s2<-HLfit(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family=binomial(),rand.family=Beta(),HLmethod="HL(0,0)",data=e2)
# g1<-glmmadmb(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family="betabinomial",data=b2)
# g2<-glmmadmb(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family="betabinomial",data=e2)
# spamm.beta[k,2]<-fixef(g1)[2]; spamm.beta[k,1]<-fixef(g1)[1]
# spamm.eps[k,2]<-fixef(g2)[2]; spamm.eps[k,1]<-fixef(g2)[1]
# spamm.beta[k,4]<-sqrt(VarCorr(g1)[[1]]);spamm.beta[k,3]<-fixef(g1)[3]
# spamm.eps[k,4]<-sqrt(VarCorr(g2)[[1]]);spamm.eps[k,3]<-fixef(g2)[3]

spamm.beta[k,1:3]<-fixef(s1); spamm.beta[k,4]<-(s1$lambda)
spamm.eps[k,1:3]<-fixef(s2);  spamm.eps[k,4]<-(s2$lambda)

	}

spb.means<-apply(spamm.beta,2,mean)
spe.means<-apply(spamm.eps,2,mean)

spb.95<-apply(spamm.beta,2,function(x){quantile(x,c(0.025,0.975))})
spe.95<-apply(spamm.eps,2,function(x){quantile(x,c(0.025,0.975))})

###### glmmADMB fits
install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")),type="source")

library(glmmADMB)


g1<-glmmadmb(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family="betabinomial",data=b2)
summary(g1)

######################
#
#	AIC Comparisons
#
######################	

e2<-data.eps(n.pops=10,n.ind=20,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
g2<-glmmadmb(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family="betabinomial",data=e2)
ob1<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid)+(1|obs),family=binomial,data=e2)
AIC(g2,ob1)

######################
#
#	Can You Identify Overdispersion in Bernoulli Models? - NO
#
######################	
	
nbern<-1000
epsod<-betaod<-numeric(nbern)

for(k in 1:nbern){
			b1<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=1)
			m1<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family=binomial,data=b1)
			betaod[k]<-od.point(m1)

			e1<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=1)
			m1<-glmer(cbind(hatched,clutch-hatched) ~ flowervisits + bodysize.c  + (1|popid),family=binomial,data=e1)
			epsod[k]<-od.point(m1)			
}

mean(betaod); quantile(betaod,c(0.025,0.975))
mean(epsod); quantile(epsod,c(0.025,0.975))


######################
#
#	Bayesian Beta-Binomial Models
#
######################
library(runjags)


######################
# Beta Binomial Data : Bayesian Analogue of Top Model for Prediction
######################

bayes.beta<-"model{
	
	# Population Intercepts 
	
		for (k in 1:npop){
		alpha.pop[k]~dnorm(mu.pop,tau.pop)
		}
		
	#Intercept Priors	
		mu.pop~dnorm(0,0.001)
		sd.pop~dunif(0,10)
		tau.pop<-pow(sd.pop,-2)
					
	#Slope Priors
	b.flowervisits~dnorm(0,0.0001)
	b.body~dnorm(0,0.001)
	
	#Theta
	theta~dgamma(0.001,0.001)

		
	#Model
	for (n in 1:nsample){
		
		#Linear Predictor + Logit Trick	
		logit.p[n]<-alpha.pop[pop[n]]  + b.flowervisits*flowervisits[n] + b.body*bodysize[n]
		mean.p[n]<-exp(logit.p[n]) / (1 + exp(logit.p[n]))
		
		#Beta Parameters, sampling from Theta
		a[n]<-(mean.p[n] / theta)
		b[n]<-((1-mean.p[n]) / theta)
		
		#Response Comes From a Binomial Distribution with mean p
		hatch[n]~dbetabin(a[n],b[n],clutch[n])			
		
		
							}

		}"


############ Global Bayesian Parameters

#Model Inits
clutchinits<-function(){dump.format(list(b.body=rnorm(1,-0.5,0.5),b.flowervisits=rnorm(1,-0.5,0.5),mu.pop=rnorm(1,-0.5,0.5),sd.pop=runif(1,0.1,1)))}
clutchparams<-c("mu.pop","b.body","b.flowervisits","theta","sd.pop","dic")

clutchburn=2000
clutchsamp=1000
clutchthin=20


############ Data
odbeta1<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=simtheta[1],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
odbeta2<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=simtheta[2],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
odbeta3<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=simtheta[3],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)

pbeta1<-data.betabinom(n.pops=simpops[1],n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
pbeta2<-data.betabinom(n.pops=simpops[2],n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
pbeta3<-data.betabinom(n.pops=simpops[3],n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)

cbeta1<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[1])
cbeta2<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[2])
cbeta3<-data.betabinom(n.pops=n.pops.sample,n.ind=n.ind.sample,theta=2,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[3])

odeps1<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=simepsilon[1],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
odeps2<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=simepsilon[2],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
odeps3<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=simepsilon[3],intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)

peps1<-data.eps(n.pops=simpops[1],n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
peps2<-data.eps(n.pops=simpops[2],n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)
peps3<-data.eps(n.pops=simpops[3],n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=nclutch)

ceps1<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[1])
ceps2<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[2])
ceps3<-data.eps(n.pops=n.pops.sample,n.ind=n.ind.sample,epsilon=3,intercept.mu=intercept.sample,intercept.sd=intercept.sd.sample,beta1slope=beta1.strong,beta2slope=beta2,clutch=simclutch[3])


############ JAGS Data

#Overdispersion
odbeta1.data<-dump.format(list(pop=as.numeric(odbeta1 $popid),flowervisits= odbeta1 $flowervisits,bodysize= odbeta1 $bodysize.c,clutch= odbeta1 $clutch,hatch= odbeta1 $hatched,nsample=nrow(odbeta1),npop=max(as.numeric(odbeta1 $popid)),bodysizemean=mean(odbeta1 $bodysize)))
odbeta2.data<-dump.format(list(pop=as.numeric(odbeta2 $popid),flowervisits= odbeta2 $flowervisits,bodysize= odbeta2 $bodysize.c,clutch= odbeta2 $clutch,hatch= odbeta2 $hatched,nsample=nrow(odbeta2),npop=max(as.numeric(odbeta2 $popid)),bodysizemean=mean(odbeta2 $bodysize)))
odbeta3.data<-dump.format(list(pop=as.numeric(odbeta3 $popid),flowervisits= odbeta3 $flowervisits,bodysize= odbeta3 $bodysize.c,clutch= odbeta3 $clutch,hatch= odbeta3 $hatched,nsample=nrow(odbeta3),npop=max(as.numeric(odbeta3 $popid)),bodysizemean=mean(odbeta3 $bodysize)))

#Pops
pbeta1.data<-dump.format(list(pop=as.numeric(pbeta1 $popid),flowervisits= pbeta1 $flowervisits,bodysize= pbeta1 $bodysize.c,clutch= pbeta1 $clutch,hatch= pbeta1 $hatched,nsample=nrow(pbeta1),npop=max(as.numeric(pbeta1 $popid)),bodysizemean=mean(pbeta1 $bodysize)))
pbeta2.data<-dump.format(list(pop=as.numeric(pbeta2 $popid),flowervisits= pbeta2 $flowervisits,bodysize= pbeta2 $bodysize.c,clutch= pbeta2 $clutch,hatch= pbeta2 $hatched,nsample=nrow(pbeta2),npop=max(as.numeric(pbeta2 $popid)),bodysizemean=mean(pbeta2 $bodysize)))
pbeta3.data<-dump.format(list(pop=as.numeric(pbeta3 $popid),flowervisits= pbeta3 $flowervisits,bodysize= pbeta3 $bodysize.c,clutch= pbeta3 $clutch,hatch= pbeta3 $hatched,nsample=nrow(pbeta3),npop=max(as.numeric(pbeta3 $popid)),bodysizemean=mean(pbeta3 $bodysize)))

#Clutch
cbeta1.data<-dump.format(list(pop=as.numeric(cbeta1 $popid),flowervisits= cbeta1 $flowervisits,bodysize= cbeta1 $bodysize.c,clutch= cbeta1 $clutch,hatch= cbeta1 $hatched,nsample=nrow(cbeta1),npop=max(as.numeric(cbeta1 $popid)),bodysizemean=mean(cbeta1 $bodysize)))
cbeta2.data<-dump.format(list(pop=as.numeric(cbeta2 $popid),flowervisits= cbeta2 $flowervisits,bodysize= cbeta2 $bodysize.c,clutch= cbeta2 $clutch,hatch= cbeta2 $hatched,nsample=nrow(cbeta2),npop=max(as.numeric(cbeta2 $popid)),bodysizemean=mean(cbeta2 $bodysize)))
cbeta3.data<-dump.format(list(pop=as.numeric(cbeta3 $popid),flowervisits= cbeta3 $flowervisits,bodysize= cbeta3 $bodysize.c,clutch= cbeta3 $clutch,hatch= cbeta3 $hatched,nsample=nrow(cbeta3),npop=max(as.numeric(cbeta3 $popid)),bodysizemean=mean(cbeta3 $bodysize)))

# Epsilon Overdispersion
odeps1.data<-dump.format(list(pop=as.numeric(odeps1 $popid),flowervisits= odeps1 $flowervisits,bodysize= odeps1 $bodysize.c,clutch= odeps1 $clutch,hatch= odeps1 $hatched,nsample=nrow(odeps1),npop=max(as.numeric(odeps1 $popid)),bodysizemean=mean(odeps1 $bodysize)))
odeps2.data<-dump.format(list(pop=as.numeric(odeps2 $popid),flowervisits= odeps2 $flowervisits,bodysize= odeps2 $bodysize.c,clutch= odeps2 $clutch,hatch= odeps2 $hatched,nsample=nrow(odeps2),npop=max(as.numeric(odeps2 $popid)),bodysizemean=mean(odeps2 $bodysize)))
odeps3.data<-dump.format(list(pop=as.numeric(odeps3 $popid),flowervisits= odeps3 $flowervisits,bodysize= odeps3 $bodysize.c,clutch= odeps3 $clutch,hatch= odeps3 $hatched,nsample=nrow(odeps3),npop=max(as.numeric(odeps3 $popid)),bodysizemean=mean(odeps3 $bodysize)))

#Epsilon Pops
peps1.data<-dump.format(list(pop=as.numeric(peps1 $popid),flowervisits= peps1 $flowervisits,bodysize= peps1 $bodysize.c,clutch= peps1 $clutch,hatch= peps1 $hatched,nsample=nrow(peps1),npop=max(as.numeric(peps1 $popid)),bodysizemean=mean(peps1 $bodysize)))
peps2.data<-dump.format(list(pop=as.numeric(peps2 $popid),flowervisits= peps2 $flowervisits,bodysize= peps2 $bodysize.c,clutch= peps2 $clutch,hatch= peps2 $hatched,nsample=nrow(peps2),npop=max(as.numeric(peps2 $popid)),bodysizemean=mean(peps2 $bodysize)))
peps3.data<-dump.format(list(pop=as.numeric(peps3 $popid),flowervisits= peps3 $flowervisits,bodysize= peps3 $bodysize.c,clutch= peps3 $clutch,hatch= peps3 $hatched,nsample=nrow(peps3),npop=max(as.numeric(peps3 $popid)),bodysizemean=mean(peps3 $bodysize)))

#Epsilon Clutch
ceps1.data<-dump.format(list(pop=as.numeric(ceps1 $popid),flowervisits= ceps1 $flowervisits,bodysize= ceps1 $bodysize.c,clutch= ceps1 $clutch,hatch= ceps1 $hatched,nsample=nrow(ceps1),npop=max(as.numeric(ceps1 $popid)),bodysizemean=mean(ceps1 $bodysize)))
ceps2.data<-dump.format(list(pop=as.numeric(ceps2 $popid),flowervisits= ceps2 $flowervisits,bodysize= ceps2 $bodysize.c,clutch= ceps2 $clutch,hatch= ceps2 $hatched,nsample=nrow(ceps2),npop=max(as.numeric(ceps2 $popid)),bodysizemean=mean(ceps2 $bodysize)))
ceps3.data<-dump.format(list(pop=as.numeric(ceps3 $popid),flowervisits= ceps3 $flowervisits,bodysize= ceps3 $bodysize.c,clutch= ceps3 $clutch,hatch= ceps3 $hatched,nsample=nrow(ceps3),npop=max(as.numeric(ceps3 $popid)),bodysizemean=mean(ceps3 $bodysize)))



############# Run Models
odbeta1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odbeta1.data)
odbeta2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odbeta2.data)
odbeta3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odbeta3.data)

pbeta1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=pbeta1.data)
pbeta2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=pbeta2.data)
pbeta3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=pbeta3.data)

cbeta1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=cbeta1.data)
cbeta2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=cbeta2.data)
cbeta3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=cbeta3.data)


odeps1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odeps1.data)
odeps2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odeps2.data)
odeps3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=odeps3.data)

peps1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=peps1.data)
peps2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=peps2.data)
peps3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=peps3.data)

ceps1.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=ceps1.data)
ceps2.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=ceps2.data)
ceps3.jags<-run.jags(model=bayes.beta,monitor=clutchparams,burnin=clutchburn,sample=clutchsamp,thin=clutchthin,inits=c(clutchinits(),clutchinits()),data=ceps3.data)



								######################################
								#
								#		Plots
								#
								######################################

#Read in Data - Need to Replace With Your FIle Paths Here and Name the Output Files Above Accordingly
odbmeans<-read.csv('OD Beta Binomial Means.csv',header=T)
odbCI<-read.csv('OD Beta Binomial CIs.csv',header=T)
odemeans<-read.csv('OD Overdispersed Means.csv',header=T)
odeCI<-read.csv('OD Overdispersed CIs.csv',header=T)

pbmeans<-read.csv('Population Beta Binomial Means.csv',header=T)
pemeans<-read.csv('Population Overdispersed Means.csv',header=T)
pbCI<-read.csv('Population Beta Binomial CIs.csv',header=T)
peCI<-read.csv('Population Overdispersed CIs.csv',header=T)

cbmeans<-read.csv("Clutch Beta Binomial Means.csv",header=T)
cemeans<-read.csv("Clutch Overdispersed Means.csv",header=T)
cbCI<-read.csv("Clutch Beta Binomial CIs.csv",header=T)
ceCI<-read.csv("Clutch Overdispersed CIs.csv",header=T)

		###################
		# Effect of Overdispersion
		###################


sl1.odb.bayes<-c(odbeta1.jags$summary[[1]][3,1],odbeta2.jags$summary[[1]][3,1],odbeta3.jags$summary[[1]][3,1])
sl1.odb.bayesl95<-c(odbeta1.jags$summary[[2]][3,1],odbeta2.jags$summary[[2]][3,1],odbeta3.jags$summary[[2]][3,1])
sl1.odb.bayesu95<-c(odbeta1.jags$summary[[2]][3,5],odbeta2.jags$summary[[2]][3,5],odbeta3.jags$summary[[2]][3,5])

sl1.ode.bayes<-c(odeps1.jags$summary[[1]][3,1],odeps2.jags$summary[[1]][3,1],odeps3.jags$summary[[1]][3,1])
sl1.ode.bayesl95<-c(odeps1.jags$summary[[2]][3,1],odeps2.jags$summary[[2]][3,1],odeps3.jags$summary[[2]][3,1])
sl1.ode.bayesu95<-c(odeps1.jags$summary[[2]][3,5],odeps2.jags$summary[[2]][3,5],odeps3.jags$summary[[2]][3,5])

sl1.lCI<-as.matrix(data.frame(s1.odb.l95=as.numeric(odbCI[1,2:4]),s1.ode.l95=as.numeric(odeCI[1,2:4]),s1.odb.bayes.l95=sl1.odb.bayesl95,s1.ode.bayes.l95=sl1.ode.bayesl95))
sl1.uCI<-as.matrix(data.frame(s1.odb.u95=as.numeric(odbCI[2,2:4]),s1.ode.u95=as.numeric(odeCI[2,2:4]),s1.odb.bayes.u95=sl1.odb.bayesu95,s1.ode.bayes.u95=sl1.ode.bayesu95 ))

sl1<-as.matrix(data.frame(s1odb=odbmeans$slope1,s1ode=odemeans$slope1,s1odb.bayes=sl1.odb.bayes,s1.ode.bayes=sl1.ode.bayes))
sl1.x<-matrix(c(0.85,0.95,1.05,1.15,1.85,1.95,2.05,2.15,2.85,2.95,3.05,3.15),nrow=3,ncol=4,byrow=T)

#Mu 
int.odb.bayes<-c(odbeta1.jags$summary[[1]][1,1],odbeta2.jags$summary[[1]][1,1],odbeta3.jags$summary[[1]][1,1])
int.odb.bayesl95<-c(odbeta1.jags$summary[[2]][1,1],odbeta2.jags$summary[[2]][1,1],odbeta3.jags$summary[[2]][1,1])
int.odb.bayesu95<-c(odbeta1.jags$summary[[2]][1,5],odbeta2.jags$summary[[2]][1,5],odbeta3.jags$summary[[2]][1,5])
int.ode.bayes<-c(odeps1.jags$summary[[1]][1,1],odeps2.jags$summary[[1]][1,1],odeps3.jags$summary[[1]][1,1])
int.ode.bayesl95<-c(odeps1.jags$summary[[2]][1,1],odeps2.jags$summary[[2]][1,1],odeps3.jags$summary[[2]][1,1])
int.ode.bayesu95<-c(odeps1.jags$summary[[2]][1,5],odeps2.jags$summary[[2]][1,5],odeps3.jags$summary[[2]][1,5])
int.lCI<-as.matrix(data.frame(int.odb.l95=as.numeric(odbCI[1,8:10]),int.ode.l95=as.numeric(odeCI[1,8:10]),int.odb.bayes.l95=int.odb.bayesl95,int.ode.bayes.l95=int.ode.bayesl95))
int.uCI<-as.matrix(data.frame(int.odb.u95=as.numeric(odbCI[2,8:10]),int.ode.u95=as.numeric(odeCI[2,8:10]),int.odb.bayes.u95=int.odb.bayesu95,int.ode.bayes.u95=int.ode.bayesu95 ))
mu1<-as.matrix(data.frame(int.odb=odbmeans$interceptmean,int.ode=odemeans$interceptmean,int.odb.bayes=int.odb.bayes,s1.ode.bayes=int.ode.bayes))

#Theta Pop 
sd.odb.bayes<-c(odbeta1.jags$summary[[1]][5,1],odbeta2.jags$summary[[1]][5,1],odbeta3.jags$summary[[1]][5,1])
sd.odb.bayesl95<-c(odbeta1.jags$summary[[2]][5,1],odbeta2.jags$summary[[2]][5,1],odbeta3.jags$summary[[2]][5,1])
sd.odb.bayesu95<-c(odbeta1.jags$summary[[2]][5,5],odbeta2.jags$summary[[2]][5,5],odbeta3.jags$summary[[2]][5,5])
sd.ode.bayes<-c(odeps1.jags$summary[[1]][5,1],odeps2.jags$summary[[1]][5,1],odeps3.jags$summary[[1]][5,1])
sd.ode.bayesl95<-c(odeps1.jags$summary[[2]][5,1],odeps2.jags$summary[[2]][5,1],odeps3.jags$summary[[2]][5,1])
sd.ode.bayesu95<-c(odeps1.jags$summary[[2]][5,5],odeps2.jags$summary[[2]][5,5],odeps3.jags$summary[[2]][5,5])
sd.lCI<-as.matrix(data.frame(sd.odb.l95=as.numeric(odbCI[1,17:19]),sd.ode.l95=as.numeric(odeCI[1,17:19]),sd.odb.bayes.l95=sd.odb.bayesl95,sd.ode.bayes.l95=sd.ode.bayesl95))
sd.uCI<-as.matrix(data.frame(sd.odb.u95=as.numeric(odbCI[2,17:19]),sd.ode.u95=as.numeric(odeCI[2,17:19]),sd.odb.bayes.u95=sd.odb.bayesu95,sd.ode.bayes.u95=sd.ode.bayesu95 ))
sd1<-as.matrix(data.frame(sd.odb=odbmeans$popsd,sd.ode=odemeans$popsd,sd.odb.bayes=sd.odb.bayes,s1.ode.bayes=sd.ode.bayes))

#################Graph
par(las=1,cex.axis=2.4,lwd=1.5,mar=c(6,6,2,2),mfrow=c(1,3),oma=c(2,5,2,2))
colmat<-matrix(rep(c("cornsilk","lightblue"),each=3,times=2),ncol=4,nrow=3)
shapemat<-matrix(rep(c(21,23),each=6,times=2),nrow=3,ncol=4)

#Slope1
plot(sl1.x,sl1,xlab="",xaxt="n",type="n",ylim=c(0,3.5),ylab="")
abline(h=0.6,col="blue",lty=8)
arrows(sl1.x,sl1,sl1.x,sl1.lCI,angle=90,length=0.05)
arrows(sl1.x,sl1,sl1.x,sl1.uCI,angle=90,length=0.05)
points(sl1.x,sl1,pch=shapemat,bg=colmat,cex=2.5)
text(2,3.3,expression(beta[prey]),cex=4)
axis(1,at=seq(3),labels=c("0.1 / 0.1", "1.5 / 1" , "3 / 2"),tick=F,lty=0)
#mtext(expression(paste(sigma[epsilon], " / ", phi)),side=1,line=5,cex= 2)
mtext("Parameter Estimate",side=2,line=5,cex=2.5,las=0)


#Mu
plot(sl1.x,mu1,xlab="",xaxt="n",type="n",ylim=c(-12,5),ylab="")
abline(h=-1,col="blue",lty=8)
arrows(sl1.x,mu1,sl1.x,int.lCI,angle=90,length=0.05)
arrows(sl1.x,mu1,sl1.x,int.uCI,angle=90,length=0.05)
points(sl1.x,mu1,pch=shapemat,bg=colmat,cex=2.5)
text(2,4.1,expression(mu[pop]),cex=4)
axis(1,at=seq(3),labels=c("0.1 / 0.1", "1.5 / 1" , "3 / 2"),tick=F,lty=0)
mtext(expression(paste(sigma[epsilon], " / ", phi)),side=1,line=6.5,cex= 3.5)
legend("bottomleft",c("BetaBin. / OLRE","Overdispersed Bin. / OLRE", "BetaBin. / BetaBin.", "Overdispersed Bin. / BetaBin."),pch=c(21,21,23,23),pt.bg=colmat[1,],cex=1.8,pt.cex=3,title="Overdispersion / Model",bty="n",title.col="blue")

#SD
plot(sl1.x,sd1,xlab="",xaxt="n",type="n",ylim=c(0,3),ylab="")
abline(h=0.5,col="blue",lty=8)
arrows(sl1.x,sd1,sl1.x,sd.lCI,angle=90,length=0.05)
arrows(sl1.x,sd1,sl1.x,sd.uCI,angle=90,length=0.05)
points(sl1.x,sd1,pch=shapemat,bg=colmat,cex=2.5)
text(2,2.8,expression(sigma[pop]),cex=4)
axis(1,at=seq(3),labels=c("0.1 / 0.1", "1.5 / 1" , "3 / 2"),tick=F,lty=0)
#mtext(expression(paste(sigma[epsilon], " / ", theta)),side=1,line=5,cex= 2)


		###################
		#Effect of Levle sof Random Effect
		###################


sl1.pb.bayes<-c(pbeta1.jags$summary[[1]][3,1],pbeta2.jags$summary[[1]][3,1],pbeta3.jags$summary[[1]][3,1])
sl1.pb.bayesl95<-c(pbeta1.jags$summary[[2]][3,1],pbeta2.jags$summary[[2]][3,1],pbeta3.jags$summary[[2]][3,1])
sl1.pb.bayesu95<-c(pbeta1.jags$summary[[2]][3,5],pbeta2.jags$summary[[2]][3,5],pbeta3.jags$summary[[2]][3,5])

sl1.pe.bayes<-c(peps1.jags$summary[[1]][3,1],peps2.jags$summary[[1]][3,1],peps3.jags$summary[[1]][3,1])
sl1.pe.bayesl95<-c(peps1.jags$summary[[2]][3,1],peps2.jags$summary[[2]][3,1],peps3.jags$summary[[2]][3,1])
sl1.pe.bayesu95<-c(peps1.jags$summary[[2]][3,5],peps2.jags$summary[[2]][3,5],peps3.jags$summary[[2]][3,5])

sl1.p.lCI<-as.matrix(data.frame(s1.pb.l95=as.numeric(pbCI[1,2:4]),s1.pe.l95=as.numeric(peCI[1,2:4]),s1.pb.bayes.l95=sl1.pb.bayesl95,s1.pe.bayes.l95=sl1.pe.bayesl95))
sl1.p.uCI<-as.matrix(data.frame(s1.pb.u95=as.numeric(pbCI[2,2:4]),s1.pe.u95=as.numeric(peCI[2,2:4]),s1.pb.bayes.u95=sl1.pb.bayesu95,s1.pe.bayes.u95=sl1.pe.bayesu95 ))

sl1.p<-as.matrix(data.frame(s1pb=pbmeans$slope1,s1pe=pemeans$slope1,s1pb.bayes=sl1.pb.bayes,s1.pe.bayes=sl1.pe.bayes))
sl1.x<-matrix(c(0.85,0.95,1.05,1.15,1.85,1.95,2.05,2.15,2.85,2.95,3.05,3.15),nrow=3,ncol=4,byrow=T)

#Mu 
int.pb.bayes<-c(pbeta1.jags$summary[[1]][1,1],pbeta2.jags$summary[[1]][1,1],pbeta3.jags$summary[[1]][1,1])
int.pb.bayesl95<-c(pbeta1.jags$summary[[2]][1,1],pbeta2.jags$summary[[2]][1,1],pbeta3.jags$summary[[2]][1,1])
int.pb.bayesu95<-c(pbeta1.jags$summary[[2]][1,5],pbeta2.jags$summary[[2]][1,5],pbeta3.jags$summary[[2]][1,5])
int.pe.bayes<-c(peps1.jags$summary[[1]][1,1],peps2.jags$summary[[1]][1,1],peps3.jags$summary[[1]][1,1])
int.pe.bayesl95<-c(peps1.jags$summary[[2]][1,1],peps2.jags$summary[[2]][1,1],peps3.jags$summary[[2]][1,1])
int.pe.bayesu95<-c(peps1.jags$summary[[2]][1,5],peps2.jags$summary[[2]][1,5],peps3.jags$summary[[2]][1,5])
int.p.lCI<-as.matrix(data.frame(int.pb.l95=as.numeric(pbCI[1,8:10]),int.pe.l95=as.numeric(peCI[1,8:10]),int.pb.bayes.l95=int.pb.bayesl95,int.pe.bayes.l95=int.pe.bayesl95))
int.p.uCI<-as.matrix(data.frame(int.pb.u95=as.numeric(pbCI[2,8:10]),int.pe.u95=as.numeric(peCI[2,8:10]),int.pb.bayes.u95=int.pb.bayesu95,int.pe.bayes.u95=int.pe.bayesu95 ))
mu1.p<-as.matrix(data.frame(int.pb=pbmeans$interceptmean,int.pe=pemeans$interceptmean,int.pb.bayes=int.pb.bayes,s1.pe.bayes=int.pe.bayes))

#Theta Pop 
sd.pb.bayes<-c(pbeta1.jags$summary[[1]][5,1],pbeta2.jags$summary[[1]][5,1],pbeta3.jags$summary[[1]][5,1])
sd.pb.bayesl95<-c(pbeta1.jags$summary[[2]][5,1],pbeta2.jags$summary[[2]][5,1],pbeta3.jags$summary[[2]][5,1])
sd.pb.bayesu95<-c(pbeta1.jags$summary[[2]][5,5],pbeta2.jags$summary[[2]][5,5],pbeta3.jags$summary[[2]][5,5])
sd.pe.bayes<-c(peps1.jags$summary[[1]][5,1],peps2.jags$summary[[1]][5,1],peps3.jags$summary[[1]][5,1])
sd.pe.bayesl95<-c(peps1.jags$summary[[2]][5,1],peps2.jags$summary[[2]][5,1],peps3.jags$summary[[2]][5,1])
sd.pe.bayesu95<-c(peps1.jags$summary[[2]][5,5],peps2.jags$summary[[2]][5,5],peps3.jags$summary[[2]][5,5])
sd.p.lCI<-as.matrix(data.frame(sd.pb.l95=as.numeric(pbCI[1,17:19]),sd.pe.l95=as.numeric(peCI[1,17:19]),sd.pb.bayes.l95=sd.pb.bayesl95,sd.pe.bayes.l95=sd.pe.bayesl95))
sd.p.uCI<-as.matrix(data.frame(sd.pb.u95=as.numeric(pbCI[2,17:19]),sd.pe.u95=as.numeric(peCI[2,17:19]),sd.pb.bayes.u95=sd.pb.bayesu95,sd.pe.bayes.u95=sd.pe.bayesu95 ))
sd1.p<-as.matrix(data.frame(sd.pb=pbmeans$popsd,sd.pe=pemeans$popsd,sd.pb.bayes=sd.pb.bayes,s1.pe.bayes=sd.pe.bayes))





par(las=1,cex.axis=2.4,lwd=1.5,mar=c(6,6,2,2),mfrow=c(1,3),oma=c(2,5,2,2))
colmat<-matrix(rep(c("cornsilk","lightblue"),each=3,times=2),ncol=4,nrow=3)
shapemat<-matrix(rep(c(21,23),each=6,times=2),nrow=3,ncol=4)

#Slope1
plot(sl1.x,sl1.p,xlab="",xaxt="n",type="n",ylim=c(-0.2,5),ylab="")
abline(h=0.6,col="blue",lty=8)
arrows(sl1.x,sl1.p,sl1.x,sl1.p.lCI,angle=90,length=0.05)
arrows(sl1.x,sl1.p,sl1.x,sl1.p.uCI,angle=90,length=0.05)
points(sl1.x,sl1.p,pch=shapemat,bg=colmat,cex=2.5)
text(2,4.8,expression(beta[prey]),cex=4)
axis(1,at=seq(3),labels=c(3,5,20),tick=F,lty=0)
#mtext("n Populations for Random Intercept",side=1,line=5,cex= 1.5)
mtext("Parameter Estimate",side=2,line=5,cex=2.5,las=0)


#Mu
plot(sl1.x,mu1.p,xlab="",xaxt="n",type="n",ylim=c(-22,10),ylab="")
abline(h=-1,col="blue",lty=8)
arrows(sl1.x,mu1.p,sl1.x,int.p.lCI,angle=90,length=0.05)
arrows(sl1.x,mu1.p,sl1.x,int.p.uCI,angle=90,length=0.05)
points(sl1.x,mu1.p,pch=shapemat,bg=colmat,cex=2.5)
text(2,9.2,expression(mu[pop]),cex=4)
axis(1,at=seq(3),labels=c(3,5,20),tick=F,lty=0)
mtext("n Populations for Random Intercept",side=1,line=5,cex= 2.5)
legend("bottomright",c("BetaBin. / OLRE","Overdispersed Bin. / OLRE", "BetaBin. / BetaBin.", "Overdispersed Bin. / BetaBin."),pch=c(21,21,23,23),pt.bg=colmat[1,],cex=1.8,pt.cex=3,title="Overdispersion / Model",bty="n",title.col="blue")

#SD
plot(sl1.x,sd1.p,xlab="",xaxt="n",type="n",ylim=c(0,7),ylab="")
abline(h=0.5,col="blue",lty=8)
arrows(sl1.x,sd1.p,sl1.x,sd.p.lCI,angle=90,length=0.05)
arrows(sl1.x,sd1.p,sl1.x,sd.p.uCI,angle=90,length=0.05)
points(sl1.x,sd1.p,pch=shapemat,bg=colmat,cex=2.5)
text(2,6.8,expression(sigma[pop]),cex=4)
axis(1,at=seq(3),labels=c(3,5,20),tick=F,lty=0)
#mtext("n Populations for Random Intercept",side=1,line=5,cex= 1.5)




		###################
		#Effect of Binomial Sample Size (Clutch)
		###################


sl1.cb.bayes<-c(cbeta1.jags$summary[[1]][3,1],cbeta2.jags$summary[[1]][3,1],cbeta3.jags$summary[[1]][3,1])
sl1.cb.bayesl95<-c(cbeta1.jags$summary[[2]][3,1],cbeta2.jags$summary[[2]][3,1],cbeta3.jags$summary[[2]][3,1])
sl1.cb.bayesu95<-c(cbeta1.jags$summary[[2]][3,5],cbeta2.jags$summary[[2]][3,5],cbeta3.jags$summary[[2]][3,5])

sl1.ce.bayes<-c(ceps1.jags$summary[[1]][3,1],ceps2.jags$summary[[1]][3,1],ceps3.jags$summary[[1]][3,1])
sl1.ce.bayesl95<-c(ceps1.jags$summary[[2]][3,1],ceps2.jags$summary[[2]][3,1],ceps3.jags$summary[[2]][3,1])
sl1.ce.bayesu95<-c(ceps1.jags$summary[[2]][3,5],ceps2.jags$summary[[2]][3,5],ceps3.jags$summary[[2]][3,5])

sl1.c.lCI<-as.matrix(data.frame(s1.cb.l95=as.numeric(cbCI[1,2:4]),s1.ce.l95=as.numeric(ceCI[1,2:4]),s1.cb.bayes.l95=sl1.cb.bayesl95,s1.ce.bayes.l95=sl1.ce.bayesl95))
sl1.c.uCI<-as.matrix(data.frame(s1.cb.u95=as.numeric(cbCI[2,2:4]),s1.ce.u95=as.numeric(ceCI[2,2:4]),s1.cb.bayes.u95=sl1.cb.bayesu95,s1.ce.bayes.u95=sl1.ce.bayesu95 ))

sl1.c<-as.matrix(data.frame(s1cb=cbmeans$slope1,s1ce=cemeans$slope1,s1cb.bayes=sl1.cb.bayes,s1.ce.bayes=sl1.ce.bayes))
sl1.x<-matrix(c(0.85,0.95,1.05,1.15,1.85,1.95,2.05,2.15,2.85,2.95,3.05,3.15),nrow=3,ncol=4,byrow=T)

#Mu 
int.cb.bayes<-c(cbeta1.jags$summary[[1]][1,1],cbeta2.jags$summary[[1]][1,1],cbeta3.jags$summary[[1]][1,1])
int.cb.bayesl95<-c(cbeta1.jags$summary[[2]][1,1],cbeta2.jags$summary[[2]][1,1],cbeta3.jags$summary[[2]][1,1])
int.cb.bayesu95<-c(cbeta1.jags$summary[[2]][1,5],cbeta2.jags$summary[[2]][1,5],cbeta3.jags$summary[[2]][1,5])
int.ce.bayes<-c(ceps1.jags$summary[[1]][1,1],ceps2.jags$summary[[1]][1,1],ceps3.jags$summary[[1]][1,1])
int.ce.bayesl95<-c(ceps1.jags$summary[[2]][1,1],ceps2.jags$summary[[2]][1,1],ceps3.jags$summary[[2]][1,1])
int.ce.bayesu95<-c(ceps1.jags$summary[[2]][1,5],ceps2.jags$summary[[2]][1,5],ceps3.jags$summary[[2]][1,5])
int.c.lCI<-as.matrix(data.frame(int.cb.l95=as.numeric(cbCI[1,8:10]),int.ce.l95=as.numeric(ceCI[1,8:10]),int.cb.bayes.l95=int.cb.bayesl95,int.ce.bayes.l95=int.ce.bayesl95))
int.c.uCI<-as.matrix(data.frame(int.cb.u95=as.numeric(cbCI[2,8:10]),int.ce.u95=as.numeric(ceCI[2,8:10]),int.cb.bayes.u95=int.cb.bayesu95,int.ce.bayes.u95=int.ce.bayesu95 ))
mu1.c<-as.matrix(data.frame(int.cb=cbmeans$interceptmean,int.ce=cemeans$interceptmean,int.cb.bayes=int.cb.bayes,s1.ce.bayes=int.ce.bayes))

#Theta Pop 
sd.cb.bayes<-c(cbeta1.jags$summary[[1]][5,1],cbeta2.jags$summary[[1]][5,1],cbeta3.jags$summary[[1]][5,1])
sd.cb.bayesl95<-c(cbeta1.jags$summary[[2]][5,1],cbeta2.jags$summary[[2]][5,1],cbeta3.jags$summary[[2]][5,1])
sd.cb.bayesu95<-c(cbeta1.jags$summary[[2]][5,5],cbeta2.jags$summary[[2]][5,5],cbeta3.jags$summary[[2]][5,5])
sd.ce.bayes<-c(ceps1.jags$summary[[1]][5,1],ceps2.jags$summary[[1]][5,1],ceps3.jags$summary[[1]][5,1])
sd.ce.bayesl95<-c(ceps1.jags$summary[[2]][5,1],ceps2.jags$summary[[2]][5,1],ceps3.jags$summary[[2]][5,1])
sd.ce.bayesu95<-c(ceps1.jags$summary[[2]][5,5],ceps2.jags$summary[[2]][5,5],ceps3.jags$summary[[2]][5,5])
sd.c.lCI<-as.matrix(data.frame(sd.cb.l95=as.numeric(cbCI[1,17:19]),sd.ce.l95=as.numeric(ceCI[1,17:19]),sd.cb.bayes.l95=sd.cb.bayesl95,sd.ce.bayes.l95=sd.ce.bayesl95))
sd.c.uCI<-as.matrix(data.frame(sd.cb.u95=as.numeric(cbCI[2,17:19]),sd.ce.u95=as.numeric(ceCI[2,17:19]),sd.cb.bayes.u95=sd.cb.bayesu95,sd.ce.bayes.u95=sd.ce.bayesu95 ))
sd1.c<-as.matrix(data.frame(sd.cb=cbmeans$popsd,sd.ce=cemeans$popsd,sd.cb.bayes=sd.cb.bayes,s1.ce.bayes=sd.ce.bayes))





par(las=1,cex.axis=2.4,lwd=1.5,mar=c(6,6,2,2),mfrow=c(1,3),oma=c(2,5,2,2))
colmat<-matrix(rep(c("cornsilk","lightblue"),each=3,times=2),ncol=4,nrow=3)
shapemat<-matrix(rep(c(21,23),each=6,times=2),nrow=3,ncol=4)

#Slope1
plot(sl1.x,sl1.c,xlab="",xaxt="n",type="n",ylim=c(-0.2,5),ylab="")
abline(h=0.6,col="blue",lty=8)
arrows(sl1.x,sl1.c,sl1.x,sl1.c.lCI,angle=90,length=0.05)
arrows(sl1.x,sl1.c,sl1.x,sl1.c.uCI,angle=90,length=0.05)
points(sl1.x,sl1.c,pch=shapemat,bg=colmat,cex=2.5)
text(2,4.8,expression(beta[prey]),cex=4)
axis(1,at=seq(3),labels= simclutch,tick=F,lty=0)
#mtext("n Populations for Random Intercept",side=1,line=5,cex= 1.5)
mtext("Parameter Estimate",side=2,line=5,cex=2.5,las=0)


#Mu
plot(sl1.x,mu1.c,xlab="",xaxt="n",type="n",ylim=c(-25,10),ylab="")
abline(h=-1,col="blue",lty=8)
arrows(sl1.x,mu1.c,sl1.x,int.c.lCI,angle=90,length=0.05)
arrows(sl1.x,mu1.c,sl1.x,int.c.uCI,angle=90,length=0.05)
points(sl1.x,mu1.c,pch=shapemat,bg=colmat,cex=2.5)
text(2,9.2,expression(mu[pop]),cex=4)
axis(1,at=seq(3),labels= simclutch,tick=F,lty=0)
mtext("Binomial Sample (Clutch) Size",side=1,line=5,cex= 2.5)
legend("bottomright",c("BetaBin. / OLRE","Overdispersed Bin. / OLRE", "BetaBin. / BetaBin.", "Overdispersed Bin. / BetaBin."),pch=c(21,21,23,23),pt.bg=colmat[1,],cex=1.8,pt.cex=3,title="Overdispersion / Model",bty="n",title.col="blue")

#SD
plot(sl1.x,sd1.c,xlab="",xaxt="n",type="n",ylim=c(0,7),ylab="")
abline(h=0.5,col="blue",lty=8)
arrows(sl1.x,sd1.c,sl1.x,sd.c.lCI,angle=90,length=0.05)
arrows(sl1.x,sd1.c,sl1.x,sd.c.uCI,angle=90,length=0.05)
points(sl1.x,sd1.c,pch=shapemat,bg=colmat,cex=2.5)
text(2,6.8,expression(sigma[pop]),cex=4)
axis(1,at=seq(3),labels=simclutch,tick=F,lty=0)
#mtext("n Populations for Random Intercept",side=1,line=5,cex= 1.5)


##########################################
# Proportion Negative Graph
##########################################

pmatrix<-read.csv('SLOPE 2 PROPORTION NEGATIVE MEANS.CSV')
pmatrixl95<-read.csv('SLOPE 2 PROPORTION NEGATIVE LOWERCI.CSV')
pmatrixu95<-read.csv('SLOPE 2 PROPORTION NEGATIVE UPPER CI.CSV')
pmatrix<-as.matrix(pmatrix[,-1]); pmatrixl95<-as.matrix(pmatrixl95[,-1]);pmatrixu95<-as.matrix(pmatrixu95[,-1])


	par(mfrow=c(1,3),mar=c(2,6,1,1),oma=c(2,2,1,1),las=1,cex.axis=2.5)	
	barcol1<-"cornsilk" ; barcol2<-"lightblue"
	
	
	x<-barplot(pmatrix[,1:3],beside=T,ylim=c(0.2,0.6),xpd=F,col=c(barcol1,barcol2),names.arg=c("0.1/0.1","1 / 1.5","2/3"))
	arrows(x,pmatrix[,1:3],x,pmatrixl95[,1:3],angle=90,length=0.05)
	arrows(x,pmatrix[,1:3],x,pmatrixu95[,1:3],angle=90,length=0.05)
	mtext("Proportion Slope Values >0",side=2,las=0,line=5,cex=1.8)
	mtext(expression(paste("Overdispersion(", phi,"/",sigma[epsilon],")")),side=1,line=3,cex=1.24)
	text(1.5,0.6,"A",cex=3,xpd=T)
	legend(1.5,0.55,c("Beta-Binomial Data","Overdispersed Binomial Data"),pch=22,pt.bg=c(barcol1,barcol2),cex=2,bty="n")
	
	x2<-barplot(pmatrix[,4:6],beside=T,ylim=c(0.2,0.6),xpd=F,col=c(barcol1,barcol2),names.arg=simpops)
	arrows(x2,pmatrix[,4:6],x2,pmatrixl95[,4:6],angle=90,length=0.05)
	arrows(x2,pmatrix[,4:6],x2,pmatrixu95[,4:6],angle=90,length=0.05)
	mtext("Number of Populations",side=1,line=2.6,cex=1.24)
	text(1.5,0.6,"B",cex=3,xpd=T)
	
	x3<-barplot(pmatrix[,7:9],beside=T,ylim=c(0.2,0.6),xpd=F,col=c(barcol1,barcol2),names.arg=simclutch)
	arrows(x3,pmatrix[,7:9],x3,pmatrixl95[,7:9],angle=90,length=0.05)
	arrows(x3,pmatrix[,7:9],x3,pmatrixu95[,7:9],angle=90,length=0.05)
	mtext("Binomial Sample (Clutch) Size",side=1,line=2.6,cex=1.24)
	text(1.5,0.6,"C",cex=3,xpd=T)	





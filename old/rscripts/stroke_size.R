mice <- read.csv("/Users/micheldelange/Documents/mice/data/stroke_size2.csv",header=TRUE)
dim(mice_bl)

female <- grepl('female',mice$animal_type)
PLT <- grepl('PLT',mice$animal_type)
Pound <- grepl('Pound',mice$animal_type)
length(female)


type <- rep("aWT",dim(mice)[1])
type[grepl('PLT',mice$animal_type)] <- "PLT"
type[grepl('Pound',mice$animal_type)] <- "Pound"


### some analysis of stroke size, unscaled ###

ss <- mice$stroke_size
mean(ss)
range(ss)
plot(density(ss))

length(ss)
df <- data.frame(ss,female,type)
table(type,female)

par(mfrow=c(2,1))
boxplot(ss~type,main="stroke size by type")
boxplot(ss~female,main="stroke size by sex",xaxt='n')
axis(1, at=c(1,2), labels=c("male","female"))
grid()

t.test(ss~female)

summary(aov(ss~female))
summary(aov(ss~type))


m <- lm(ss~type*female)
AIC(m)
m2 <- lm(ss~type+female)
AIC(m2)

summary(m)




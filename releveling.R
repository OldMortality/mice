mice.df$typeperiod <- factor(paste(mice.df$type,mice.df$period,sep='-'),ordered = F)
m2 <- lmer(logity ~ typeperiod + (1|mouse.name) , 
                data=mice.df)
anova(m2)

w <- waldInterval(model=m2,Z=Z,FUN=exp)
w

mice.df$typeperiod <- relevel(mice.df$typeperiod,ref='aWT-4')
m2 <- lmer(logity ~ typeperiod + (1|mouse.name) , 
           data=mice.df)
anova(m2)

w <- waldInterval(model=m2,Z=Z,FUN=exp)
w

phi <- samples$BUGSoutput$sims.list$phi
alpha0 <- samples$BUGSoutput$sims.list$alpha0
alpha1 <- samples$BUGSoutput$sims.list$alpha1
alpha2 <- samples$BUGSoutput$sims.list$alpha2
alpha3 <- samples$BUGSoutput$sims.list$alpha3
alpha4 <- samples$BUGSoutput$sims.list$alpha4
alpha5 <- samples$BUGSoutput$sims.list$alpha5
alpha6 <- samples$BUGSoutput$sims.list$alpha6
alpha7 <- samples$BUGSoutput$sims.list$alpha7
alpha8 <- samples$BUGSoutput$sims.list$alpha8
alpha9 <- samples$BUGSoutput$sims.list$alpha9
beta1 <- samples$BUGSoutput$sims.list$beta1
beta2 <- samples$BUGSoutput$sims.list$beta2
gamma1 <- samples$BUGSoutput$sims.list$gamma1
gamma2 <- samples$BUGSoutput$sims.list$gamma2
gamma3 <- samples$BUGSoutput$sims.list$gamma3
gamma4 <- samples$BUGSoutput$sims.list$gamma4
gamma5 <- samples$BUGSoutput$sims.list$gamma5
delta1 <- samples$BUGSoutput$sims.list$delta1
delta2 <- samples$BUGSoutput$sims.list$delta2
delta3 <- samples$BUGSoutput$sims.list$delta3
delta4 <- samples$BUGSoutput$sims.list$delta4
delta5 <- samples$BUGSoutput$sims.list$delta5
delta6 <- samples$BUGSoutput$sims.list$delta6
delta7 <- samples$BUGSoutput$sims.list$delta7
delta8 <- samples$BUGSoutput$sims.list$delta8
delta9 <- samples$BUGSoutput$sims.list$delta9
delta10 <- samples$BUGSoutput$sims.list$delta10
eta_hat <- samples$BUGSoutput$sims.list$eta_hat
q <- expit(eta_hat)

p <- samples$BUGSoutput$sims.list$p
p[1,1,1]
b0 <- samples$BUGSoutput$sims.list$b0
b <- samples$BUGSoutput$sims.list$b
sigma <- samples$BUGSoutput$sims.list$sigma
sigma2 <- samples$BUGSoutput$sims.list$sigma2
p <- samples$BUGSoutput$sims.list$p
b0 <- samples$BUGSoutput$sims.list$b0
b <- samples$BUGSoutput$sims.list$b

DIC <- samples$BUGSoutput$DIC
chi <- samples$BUGSoutput$sims.list$chi
chirep <- samples$BUGSoutput$sims.list$chirep
print("got posteriors")
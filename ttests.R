# ### t.tests
# # per here is the column, so
# #  1 = time0
# #  2 = time3
# #  3 = time10
# #  4 = time17
# #  5 = time24
# do.ttest <- function(per) {
#   g1 <- getPropsByMouse(df.props,malePound.mice)[,per]
#   g2 <- getPropsByMouse(df.props,maleWT.mice)[,per]
#   t.test(g1,g2)
# }
# 
# 
# for (i in 1:5) {
#   print(i)
#   t <- do.ttest(i)
#   print(t)
# }
# 
# par(mfrow=c(1,1))
# main = paste(celltype,' t-tests: CIs mn(Pound) - mn(WT')
# plot('',xlim=c(0,24),ylim=c(-1,1),xaxt='n',
#      main=main,xlab='day')
# xbreaks <- c(0,3,10,17,
#              24)
# abline(v=xbreaks,lty=3,col='grey')
# abline(h=seq(-1,1,0.1),lty=3,col='grey')
# axis(side = 1, at=xbreaks)
# for (i in 1:5) {
#   t <- do.ttest(i)
#   ci <- t$conf.int
#   segments(xbreaks[i],ci[1],xbreaks[i],ci[2])
# }
# abline(h=0,col='red')
# 
# 
# # paired t-test for one type of mouse, period compared to period 0
# doPairedttest <- function(per,group) {
#   g1 <- getPropsByMouse(df.props,group)[,per]
#   g2 <- getPropsByMouse(df.props,group)[,1]
#   t.test(g1,g2,paired=T)
# }
# 
# group = maleWT.mice
# doPairedttest(2,group)$p.value
# doPairedttest(3,group)$p.value
# doPairedttest(4,group)$p.value
# doPairedttest(5,group)$p.value
# group = malePound.mice
# doPairedttest(2,group)$p.value
# doPairedttest(3,group)$p.value
# doPairedttest(4,group)$p.value
# doPairedttest(5,group)$p.value
# 

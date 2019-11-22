

```r
mice <- read.csv("/Users/micheldelange/Documents/mice/data/mice.csv",header=TRUE)[1:40,]


# number of mice
N <- dim(mice)[1]
# number of time periods (0,3,10,17,24,31)
periods <- 6

######### FIRST LOT ######################
cell_type <- "CD19pos_B220pos"
N_0 <- mice$Live_cells_0
R_0 <- mice$CD19pos_B220pos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD19pos_B220pos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD19pos_B220pos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD19pos_B220pos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD19pos_B220pos_24
N_31 <- mice$Live_cells_31
R_31 <- mice$CD19pos_B220pos_31

source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice5.r")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 240
##    Unobserved stochastic nodes: 533
##    Total graph size: 3839
## 
## Initializing model
## 
##   |                                                          |                                                  |   0%  |                                                          |+++++                                             |  10%  |                                                          |++++++++++                                        |  20%  |                                                          |+++++++++++++++                                   |  30%  |                                                          |++++++++++++++++++++                              |  40%  |                                                          |+++++++++++++++++++++++++                         |  50%  |                                                          |++++++++++++++++++++++++++++++                    |  60%  |                                                          |+++++++++++++++++++++++++++++++++++               |  70%  |                                                          |++++++++++++++++++++++++++++++++++++++++          |  80%  |                                                          |+++++++++++++++++++++++++++++++++++++++++++++     |  90%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
##   |                                                          |                                                  |   0%  |                                                          |*                                                 |   2%  |                                                          |**                                                |   5%  |                                                          |****                                              |   8%  |                                                          |*****                                             |  10%  |                                                          |******                                            |  12%  |                                                          |********                                          |  15%  |                                                          |*********                                         |  18%  |                                                          |**********                                        |  20%  |                                                          |***********                                       |  22%  |                                                          |************                                      |  25%  |                                                          |**************                                    |  28%  |                                                          |***************                                   |  30%  |                                                          |****************                                  |  32%  |                                                          |******************                                |  35%  |                                                          |*******************                               |  38%  |                                                          |********************                              |  40%  |                                                          |*********************                             |  42%  |                                                          |**********************                            |  45%  |                                                          |************************                          |  48%  |                                                          |*************************                         |  50%  |                                                          |**************************                        |  52%  |                                                          |****************************                      |  55%  |                                                          |*****************************                     |  58%  |                                                          |******************************                    |  60%  |                                                          |*******************************                   |  62%  |                                                          |********************************                  |  65%  |                                                          |**********************************                |  68%  |                                                          |***********************************               |  70%  |                                                          |************************************              |  72%  |                                                          |**************************************            |  75%  |                                                          |***************************************           |  78%  |                                                          |****************************************          |  80%  |                                                          |*****************************************         |  82%  |                                                          |******************************************        |  85%  |                                                          |********************************************      |  88%  |                                                          |*********************************************     |  90%  |                                                          |**********************************************    |  92%  |                                                          |************************************************  |  95%  |                                                          |************************************************* |  98%  |                                                          |**************************************************| 100%
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-3.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-4.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-5.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-6.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-7.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-8.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-9.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-10.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-11.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-12.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-13.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-14.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-15.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-16.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-17.png)

```
##      2.5%     97.5% 
## 0.4901716 1.0987714
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-18.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-19.png)

```
##     2.5%    97.5% 
## 1.166207 2.465271
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-20.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-21.png)

```
##      2.5%     97.5% 
## 0.2933059 0.6428816
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-22.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-23.png)

```
##      2.5%     97.5% 
## 0.5035617 0.9393074
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-24.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-25.png)

```
##      2.5%     97.5% 
## 0.7018151 1.3310102
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-26.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-27.png)

```
##      2.5%     97.5% 
## 0.6080282 1.1477795
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-28.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-29.png)

```
##      2.5%     97.5% 
## 0.6748537 1.2732475
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-30.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-31.png)

```
##      2.5%     97.5% 
## 0.5180600 0.9704089
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-32.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-33.png)

```
##      2.5%     97.5% 
## 0.3862126 1.1221160
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-34.png)

```
##     2.5%    97.5% 
## 1.001913 2.788167
```

```
## No id variables; using all as measure variables
## No id variables; using all as measure variables
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-35.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-36.png)

```
## [1] "WT 0 -0.124254986019966 0.015754601377978"
```

```
## [1] "WT 3 -0.0998220917378181 0.0120650272877392"
```

```
## [1] "WT 10 -0.122133099892044 0.014993868633609"
```

```
## [1] "WT 17 -0.112251513800279 0.0138100267134841"
```

```
## [1] "WT 24 -0.119715219190574 0.014710037669208"
```

```
## [1] "WT 31 -0.101740417348166 0.012529303845825"
```

```
## [1] "PLT\t 0 -0.211029216897827 -0.0753155391076899  ***"
```

```
## [1] "PLT\t 3 -0.17613464886358 -0.0624418946091731  ***"
```

```
## [1] "PLT\t 10 -0.208373792132388 -0.0743956708348427  ***"
```

```
## [1] "PLT\t 17 -0.193935774988096 -0.0688455305170521  ***"
```

```
## [1] "PLT\t 24 -0.204037202149014 -0.073222946658791  ***"
```

```
## [1] "PLT\t 31 -0.1781879041812 -0.0635246197592503  ***"
```

```
## [1] "Pound 0 -0.012335865682841 0.0601845491321347"
```

```
## [1] "Pound 3 -0.00931355798107563 0.0454666847393959"
```

```
## [1] "Pound 10 -0.0119332627825172 0.0588651398283344"
```

```
## [1] "Pound 17 -0.0106276611671813 0.0530352597625363"
```

```
## [1] "Pound 24 -0.0115632900720824 0.0569062925329707"
```

```
## [1] "Pound 31 -0.00943163222876367 0.0467710797055554"
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-37.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-38.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-39.png)

```r
######### SECOND LOT ####### N_0 <- mice$CD19pos_0
cell_type = "CD19pos_B220pos_MHCII"
R_0 <- mice$CD19pos_B220pos_MHCIIpos_0
R_3 <- mice$CD19pos_B220pos_MHCIIpos_3
R_10 <- mice$CD19pos_B220pos_MHCIIpos_10
R_17 <- mice$CD19pos_B220pos_MHCIIpos_17
R_24 <- mice$CD19pos_B220pos_MHCIIpos_24
R_31 <- mice$CD19pos_B220pos_MHCIIpos_31

######### THIRD LOT ######################

cell_type <- "CD11bpos" 
N_0 <- mice$Live_cells_0
R_0 <- mice$CD11bpos_0
N_3 <- mice$Live_cells_3
R_3 <- mice$CD11bpos_3
N_10 <- mice$Live_cells_10
R_10 <- mice$CD11bpos_10
N_17 <- mice$Live_cells_17
R_17 <- mice$CD11bpos_17
N_24 <- mice$Live_cells_24
R_24 <- mice$CD11bpos_24
N_31 <- mice$Live_cells_31
R_31 <- mice$CD11bpos_31
source("/Users/micheldelange/Documents/mice/rscripts/bayesian_mice5.r")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-40.png)

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 240
##    Unobserved stochastic nodes: 533
##    Total graph size: 3839
## 
## Initializing model
## 
##   |                                                          |                                                  |   0%  |                                                          |+++++                                             |  10%  |                                                          |++++++++++                                        |  20%  |                                                          |+++++++++++++++                                   |  30%  |                                                          |++++++++++++++++++++                              |  40%  |                                                          |+++++++++++++++++++++++++                         |  50%  |                                                          |++++++++++++++++++++++++++++++                    |  60%  |                                                          |+++++++++++++++++++++++++++++++++++               |  70%  |                                                          |++++++++++++++++++++++++++++++++++++++++          |  80%  |                                                          |+++++++++++++++++++++++++++++++++++++++++++++     |  90%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
##   |                                                          |                                                  |   0%  |                                                          |*                                                 |   2%  |                                                          |**                                                |   5%  |                                                          |****                                              |   8%  |                                                          |*****                                             |  10%  |                                                          |******                                            |  12%  |                                                          |********                                          |  15%  |                                                          |*********                                         |  18%  |                                                          |**********                                        |  20%  |                                                          |***********                                       |  22%  |                                                          |************                                      |  25%  |                                                          |**************                                    |  28%  |                                                          |***************                                   |  30%  |                                                          |****************                                  |  32%  |                                                          |******************                                |  35%  |                                                          |*******************                               |  38%  |                                                          |********************                              |  40%  |                                                          |*********************                             |  42%  |                                                          |**********************                            |  45%  |                                                          |************************                          |  48%  |                                                          |*************************                         |  50%  |                                                          |**************************                        |  52%  |                                                          |****************************                      |  55%  |                                                          |*****************************                     |  58%  |                                                          |******************************                    |  60%  |                                                          |*******************************                   |  62%  |                                                          |********************************                  |  65%  |                                                          |**********************************                |  68%  |                                                          |***********************************               |  70%  |                                                          |************************************              |  72%  |                                                          |**************************************            |  75%  |                                                          |***************************************           |  78%  |                                                          |****************************************          |  80%  |                                                          |*****************************************         |  82%  |                                                          |******************************************        |  85%  |                                                          |********************************************      |  88%  |                                                          |*********************************************     |  90%  |                                                          |**********************************************    |  92%  |                                                          |************************************************  |  95%  |                                                          |************************************************* |  98%  |                                                          |**************************************************| 100%
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-41.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-42.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-43.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-44.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-45.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-46.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-47.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-48.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-49.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-50.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-51.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-52.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-53.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-54.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-55.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-56.png)

```
##      2.5%     97.5% 
## 0.2971058 0.7123241
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-57.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-58.png)

```
##      2.5%     97.5% 
## 0.4261387 0.9622195
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-59.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-60.png)

```
##     2.5%    97.5% 
## 1.251598 2.916352
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-61.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-62.png)

```
##      2.5%     97.5% 
## 0.9606323 1.6598856
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-63.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-64.png)

```
##      2.5%     97.5% 
## 0.9864734 1.7063090
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-65.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-66.png)

```
##     2.5%    97.5% 
## 1.471666 2.536435
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-67.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-68.png)

```
##      2.5%     97.5% 
## 0.9592123 1.6648667
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-69.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-70.png)

```
##      2.5%     97.5% 
## 0.7978264 1.3774222
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-71.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-72.png)

```
##     2.5%    97.5% 
## 1.008260 3.171723
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-73.png)

```
##      2.5%     97.5% 
## 0.9109331 2.7430892
```

```
## No id variables; using all as measure variables
## No id variables; using all as measure variables
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-74.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-75.png)

```
## [1] "WT\t 0 -0.201930401418283 -0.0514493526727626  ***"
```

```
## [1] "WT\t 3 -0.226149779835856 -0.0593071736661388  ***"
```

```
## [1] "WT\t 10 -0.22806390726222 -0.0597856099263279  ***"
```

```
## [1] "WT\t 17 -0.264750234103687 -0.0721758837725117  ***"
```

```
## [1] "WT\t 24 -0.227500967185826 -0.0598723063477132  ***"
```

```
## [1] "WT\t 31 -0.206867795737684 -0.0533372894176992  ***"
```

```
## [1] "PLT 0 -0.0802413479443685 0.0280239142702539"
```

```
## [1] "PLT 3 -0.0924758921951069 0.0324268476639807"
```

```
## [1] "PLT 10 -0.094170895190998 0.0325060917112603"
```

```
## [1] "PLT 17 -0.114074865293592 0.0395236593302959"
```

```
## [1] "PLT 24 -0.0927851095522536 0.0323718840793977"
```

```
## [1] "PLT 31 -0.0827985916285236 0.0288965956903965"
```

```
## [1] "Pound 0 -0.152025752119462 0.00596169647267286"
```

```
## [1] "Pound 3 -0.158480691328762 0.0063614049193974"
```

```
## [1] "Pound 10 -0.159112621225489 0.00635785082237698"
```

```
## [1] "Pound 17 -0.159969083210115 0.00656981255246542"
```

```
## [1] "Pound 24 -0.159137602730228 0.00639540074409369"
```

```
## [1] "Pound 31 -0.153857436932764 0.00623424738128076"
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-76.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-77.png)![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-78.png)


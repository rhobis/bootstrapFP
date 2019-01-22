
library(bootstrapFP)

### Generate population data ---
N   <- 100; n <- 50
x   <- rgamma(N, scale=10, shape=5)
y   <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
pik <- sampling::inclusionprobabilities(x, n)

### Draw a dummy sample ---
s  <- sample(N, n)

### Bootstrap paramters ---
B <- 1000
D <- 10


### Gross
bootstrapFP(y = y[s], B=B, method = "ppGross")  #should give an errors
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppGross")  


### Booth
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppBooth")  


### Chao and Lo (1985)
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppChaoLo85")  


### Chao and Lo (1994)
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppChaoLo94")  


### Bickel and Freedman
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppBi")  



### Sitter (pseudo-population)
bootstrapFP(y = y[s], pik = n/N, B=B, method = "ppSitter") 



### Holmberg
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHolmberg") #should give error
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHolmberg", design = 'brewer')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'brewer')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'randomSystematic')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'tille')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'sampford')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'maxEntropy')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'poisson')
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppHolmberg", design = 'systematic')


### Chauvet
bootstrapFP(y = y[s], pik = pik[s], B=B, D=D, method = "ppChauvet")



### HotDeck
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'brewer', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'randomSystematic', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'tille', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'sampford', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'maxEntropy', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'poisson', x=x, s=(1:N %in% s))
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "ppHotDeck", design = 'systematic', x=x, s=(1:N %in% s))


### Efron
bootstrapFP(y = y[s], B=B, method = "dEfron", design = 'brewer') #should give an error
bootstrapFP(y = y[s], B=B, pik = n/N, method = "dEfron", design = 'brewer')


### McCarthy and Snowden
bootstrapFP(y = y[s], B=B, pik = n/N, method = "dMcCarthySnowden", design = 'brewer')


### Rao and Wu
bootstrapFP(y = y[s], pik = n/N, B=B, method = "dRaoWu")




### Sitter (direct)
bootstrapFP(y = y[s], pik = n/N, B=B, method = "dSitter")



### Antal and Tillé (2011)
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='brewer')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='randomSystematic')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='tille')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='maxEntropy')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='sampford')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='poisson')
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "dAntalTille_UPS", design='systematic')



### Rao, Wu and Yue
bootstrapFP(y = y[s], pik = n/N, B=B, method = "wRaoWuYue") 



### Chipperfield and Preston
bootstrapFP(y = y[s], pik = n/N, B=B, method = "wChipperfieldPreston")


### Generalised
bootstrapFP(y = y[s], pik = pik[s], B=B, method = "wGeneralised", distribution = 'normal')





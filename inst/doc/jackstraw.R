## ----setup,include=FALSE,echo=FALSE-------------
opts_chunk$set(fig.align='center', fig.width=6, fig.height=6, tidy=TRUE, cache=FALSE, warning=FALSE, message=TRUE)
options(keep.source = TRUE, width=50)
desc <- packageDescription("jackstraw")

## ----sim_pca_data-------------------------------
library(jackstraw)
library(corpcor)

set.seed(1)
B = c(runif(100, min=0.1, max=1), rep(0,900))
L = c(rep(1, 10), rep(-1, 10))
L = L / sd(L)
E = matrix(rnorm(1000*20), nrow=1000)
Y = B %*% t(L) + E

dim(Y)
Y[1:5,1:5]

## ----sim_pca_PA, dependson="sim_pca_data"-------
PA = permutationPA(Y, B=10, threshold=0.05)

plot(PA$p, pch=20, main="Permutation Parallel Analysis P-values", ylab="P-values", xlab="Principal Component")

## ----sim_pca_1pc, dependson="sim_pca_data"------
svd.out = fast.svd(Y)

par(mfrow=c(2,1))
plot(svd.out$d^2/sum(svd.out$d^2), pch=20, main="The scree plot", xlab="PC", ylab="Percent Variance Explained")
plot(svd.out$d[1] * svd.out$v[,1], pch=20, main="1st PC", xlab="Observation", ylab="Magnitude")

## ----sim_pca_jackstraw, dependson="sim_pca_data", cache=TRUE----
js.pca = jackstraw(Y, r = 1, method = "PCA", s = 100, B = 100, verbose = FALSE)

hist(js.pca$p.value, 10, col="black")

## ----sim_pca_pvalues, dependson="sim_pca_jackstraw"----
par(mfrow=c(1,2))
hist(js.pca$p.value[1:100], 10, col="black", main="Alternative P-values")
hist(js.pca$p.value[101:1000], 10, col="black", main="Null P-values")

## ----sim_lfa_data, cache=TRUE-------------------
library(jackstraw)
library(lfa)

set.seed(2)
m=5000; n=100; pi0=.5
m0 = round(m*pi0)
m1 = m-round(m*pi0)
B = matrix(0, nrow=m, ncol=1)
B[1:m1,] = matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=1)
L = matrix(rnorm(n), nrow=1, ncol=n)
BL = B %*% L
prob = exp(BL)/(1+exp(BL))

dat = list(
  Y = matrix(rbinom(m*n, 2, as.numeric(prob)), m, n),
  H = c(rep(1, m1), rep(0, m0)))

dim(dat$Y)
dat$Y[1:5,1:5]

## ----sim_lfa_jackstraw, dependson="sim_lfa_data", cache=TRUE----
js.lfa = jackstraw(dat$Y, r = 2, method = "LFA-R", s = 200, B = 10)

hist(js.lfa$p.value, 10, col="black")

## ----sim_lfa_pvalues, dependson="sim_lfa_jackstraw"----
par(mfrow=c(1,2))
hist(js.lfa$p.value[which(dat$H==1)], 10, col="black", main="Alternative P-values")
hist(js.lfa$p.value[which(dat$H==0)], 10, col="black", main="Null P-values")


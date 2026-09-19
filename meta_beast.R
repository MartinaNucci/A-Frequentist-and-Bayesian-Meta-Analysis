#1
library(foreign)
library(metafor)


BF<-read.dta("C:\\Users\\marti\\OneDrive\\Desktop\\Biostat\\0.ES\\meta_breast\\BF.dta")

#2 and 3
res.fe<-rma(BF$logor,sei=BF$selogor,method="FE",level=90,slab=BF$study)
forest(res.fe)

#4
res.re<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,slab=BF$study)
res.re
forest(res.re)

#5
plot(weights(res.fe),weights(res.re),xlim=c(0,25),ylim=c(0,25))
lines(weights(res.fe),weights(res.fe),type="l")
par(mfrow=c(1,2))
forest(res.fe,showweights=TRUE)
title("FE")
forest(res.re,showweights=TRUE)
title("RE")

#6

par(mfrow=c(1,2))
funnel(res.re)
trimfill(res.re)
funnel(trimfill(res.re))
regtest(res.re)
#7
# cumulative analysis
cumul.re <- cumul(res.re, order = BF$year)
forest(cumul.re)
# leave-one-out analysis
loo <- leave1out(res.re)

forest(
  x = loo$estimate,
  ci.lb = loo$ci.lb,
  ci.ub = loo$ci.ub,
  slab = gsub("^-", "", loo$slab),
  refline = 0,
 #xlab = "Pooled effect after excluding the study",
  xlim = c(-0.5, 0.3),
  alim = c(-0.5, 0.3),
  
  # I-squared column
  ilab = round(loo$I2, 1),
  ilab.xpos = 0.1
)
# column header
text(
  x = 0.1,
  y = length(loo$estimate) + 2,
  labels = expression(I^2),
  font = 2
)

#8
res.re.PC<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,slab=BF$study,subset=BF$des=="PC")
res.re.HC<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,slab=BF$study,subset=BF$des=="HC")
res.re.CS<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,slab=BF$study,subset=BF$des=="CS")
forest(res.re,order=order(BF$des),slab=BF$des)
res.meta.design<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,
                     slab=BF$study,mods=~factor(BF$des))
res.meta.year<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,
                   slab=BF$study,mods=~BF$year)
forest(res.re,order=order(BF$year),slab=BF$year)
forest(res.re,order=order(BF$out),slab=BF$out)
res.meta.outcome<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,
                      slab=BF$study,mods=~factor(BF$out))
forest(res.re,order=order(BF$country),slab=BF$country)
res.meta.country<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,
                      slab=BF$study,mods=~factor(BF$country))
res.meta.reg<-rma(BF$logor,sei=BF$selogor,method="DL",level=90,
                      slab=BF$study,mods=~factor(BF$reg))
# plot the meta-regression line
regplot(res.meta.year,
                xlab = "year",
               ylab = "Estimated effect (log-OR)",
        ci=F)

# I-squared comparison table: FE, RE, and meta-regression models
tab_I2 <- data.frame(
  Model = c("FE","RE","design","year","outcome","country"  ),
  I2 = c(res.fe$I2,res.re$I2,res.meta.design$I2,res.meta.year$I2,res.meta.outcome$I2,res.meta.country$I2))
tab_I2$I2 <- round(tab_I2$I2, 2)


# BAYESIAN META-ANALYSIS
library(R2jags)

N<-dim(BF)[1]
beta.hat<-BF$logor
se.hat<-BF$selogor


par(mfrow=c(1,2))
# FIXED-EFFECT MODEL
model.file<- ("C:\\Users\\marti\\OneDrive\\Desktop\\Biostat\\0.ES\\meta_breast\\model_FE.txt")

info <- list ("N", "beta.hat", "se.hat")

inits <- function() {list (b=rnorm(1,0,1))}
parameters <- c("b")
res<-jags(data=info, inits=inits, parameters.to.save=parameters,digits=5, model.file=model.file, n.iter=15000, n.burnin=5000, n.thin = 1, n.chains=3) # save them in a matrix

names(res)
res$BUGSoutput

sims <- res$BUGSoutput$sims.matrix
bpost.fe <- sims[,1]
# posterior samples on the OR scale
OR.post <- exp(bpost.fe)
# quantiles
CI <- exp(quantile(bpost.fe, c(0.025, 0.5, 0.975)))
# histogram
hist(OR.post, breaks = 100)
# posterior mean
abline(v = exp(mean(bpost.fe)), col = 2, lwd = 3)
# 95% credible interval
abline(v = CI[c(1,3)], col = "blue", lwd = 2, lty = 2)

exp(mean(bpost.fe))
CI

# RANDOM-EFFECTS MODEL
model.file<- ("C:\\Users\\marti\\OneDrive\\Desktop\\Biostat\\0.ES\\meta_breast\\model_RE.txt")
info <- list ("N", "beta.hat", "se.hat")
inits <- function() {list (b=rnorm(1,0,1), beta=rnorm(N,0,1), prec=runif(1,0,30))}
parameters <- c("b", "beta", "tau2", "beta.pred")
res<-jags(data=info, inits=inits, parameters.to.save=parameters, model.file=model.file, n.iter=15000, n.burnin=5000, n.thin = 1, n.chains=3)


# check the chains
res$BUGSoutput
sims<-res$BUGSoutput$sims.matrix
#head(sims)
dim(sims)
bpost.re<-(sims[,1]) # random-effects coefficient
tau2post<-sims[,29] #tau^2
bpred<-(sims[,27]) # predicted beta
I2<-tau2post/(tau2post+mean(BF$selogor^2))



sum(exp(bpost.re)<1)/length(bpost.re)
sum(exp(bpred)<1)/length(bpred)
sum(bpred<0)/length(bpost.re)


CI_tau<-(quantile(tau2post, c(0.025,  0.975)))
CI_re <- exp(quantile(bpost.re, c(0.025, 0.975)))
CI_pred <- exp(quantile(bpred, c(0.025,  0.975)))
CI_I2<-(quantile(I2, c(0.025,  0.975)))

c(exp(mean(bpost.re)), CI_re)
c((mean(tau2post)),CI_tau)
c(exp(mean(bpred)),CI_pred)
c((mean(I2)),CI_I2)

par(mfrow=c(1,2))
# histogram of beta_RE
hist(exp(bpost.re) ,100)
abline(v=exp((mean(bpost.re))),col=2, lwd=3)
abline(v = CI_re[c(1,2)], col = "blue", lwd = 2, lty = 2)
# histogram of I-squared
hist((I2) ,100)
abline(v=((mean(I2))),col=2, lwd=3)
abline(v = CI_I2[c(1,2)], col = "blue", lwd = 2, lty = 2)


# histogram of beta_pred
hist(exp(bpred) ,100)
abline(v=exp((mean(bpred))),col=2, lwd=3)
abline(v = CI_pred[c(1,2)], col = "blue", lwd = 2, lty = 2)



par(mfrow=c(1,1))
plot(density(exp(bpred))$x,density(exp(bpred))$y,type="l",ylim=c(0,15),xlab="beta",ylab="density",col=2,xlim=c(0.3,1.2))
lines(density(exp(bpost.re))$x,density(exp(bpost.re))$y,type="l",col=1)
legend("topleft",legend = c("posterior distr.", "posterior predictive distr."), col = c(1,2), lty = 1, cex = 0.5)
abline(v=0, lty=2)

# forest plot
forest(c(BF$logor, mean(bpost.fe), mean(bpost.re)),
       ci.lb = c(BF$logor - 1.96 * BF$selogor,
                 quantile(bpost.fe, 0.025),
                 quantile(bpost.re, 0.025)),
       ci.ub = c(BF$logor + 1.96 * BF$selogor,
                 quantile(bpost.fe, 0.975),
                 quantile(bpost.re, 0.975)),
       annotate = TRUE,
       refline = 1,
       #transf = exp,
       slab = c(as.character(BF$study), "FE", "RE"),
       psize = 0.5)

# example histograms for study-specific effects
par(mfrow = c(2,3))

# studies to plot
studies <- c("beta[2]", "beta[3]", "beta[4]",
             "beta[17]", "beta[18]", "beta[23]")

for(i in studies){
  # transform to the OR scale
  OR <- exp(sims[, i])
  # 95% credible interval
  IC <- quantile(OR, c(0.025, 0.975))
  # histogram
  hist(OR,breaks = 50,main = i,xlab = "OR",col = "lightgray",border = "white")
  # posterior mean
  abline(v = mean(OR), col = "red",lwd = 3)
  # credible interval
  abline(v = IC,col = "blue",lwd = 2,lty = 2)
  # line of no effect
  abline(v = 1,lty = 3)
}

#### BAYESIAN META-REGRESSION
model.file<- ("C:\\Users\\marti\\OneDrive\\Desktop\\Biostat\\0.ES\\meta_breast\\model_REG.txt")

inits <- function() {list (b=rnorm(1,0,1/100), gamma=rnorm(1,0,1/100), beta=rnorm(N,0,1/100), prec=runif(1,0,100))}
parameters <- c("b", "gamma", "beta", "tau2")

x<-as.numeric(as.factor(BF$year))

info <- list ("N", "beta.hat", "se.hat", "x")
res<-jags(data=info, inits=inits, parameters.to.save=parameters, model.file=model.file, n.iter=150000, n.burnin=50000, n.thin = 1, n.chains=3)

res$BUGSoutput
sims<-res$BUGSoutput$sims.matrix
head(sims)
gammapost<-sims[,28]
tau2post<-sims[,29]
bpost<-sims[,1]
I2res<-tau2post/(tau2post+mean(BF$selogor^2))

IC_m_b<-exp(quantile(bpost,c(0.025,0.975)))
IC_m_gamma<-(quantile(gammapost,c(0.025,0.975)))
IC_m_tau2<-(quantile(tau2post,c(0.025,0.975)))
IC_m_I2<-(quantile(I2res,c(0.025,0.975)))

c(mean(exp(bpost)),IC_m_b)
c(mean(gammapost),IC_m_gamma)
c(mean(tau2post),IC_m_tau2)
c(mean(I2res),IC_m_I2)


par(mfrow=c(1,2))
# histogram of beta
hist(exp(bpost) ,100)
abline(v=exp((mean(bpost))),col=2, lwd=3)
abline(v = IC_m_b[c(1,2)], col = "blue", lwd = 2, lty = 2)
# histogram of gamma
hist(gammapost ,100) 
abline(v=((mean(gammapost))),col=2, lwd=3)
abline(v = IC_m_gamma[c(1,2)], col = "blue", lwd = 2, lty = 2)



hist(gammapost,100)
abline(v=((mean(gammapost))),col=2, lwd=3)
hist(I2res,100)
abline(v=((mean(I2res))),col=2, lwd=3)
mean(gammapost)
quantile(gammapost,c(0.025,0.5,0.975))
quantile(I2res,c(0.025,0.5,0.975))
mean(I2res)
(mean(bpost))
(quantile(bpost,c(0.025,0.5,0.975)))

# plot residual I-squared and overall I-squared
par(mfrow=c(1,2))
plot(density(gammapost)$x,density(gammapost)$y,xlab="meta-regressor",ylab="density",type="l")
abline(v=0,col="red")
plot(density(I2res)$x,density(I2res)$y,type="l",lty=2, xlab="I2",ylab="density")
lines(density(I2)$x,density(I2)$y,lty=1)
legend(0.4,3,legend=c("I2","residual I2"),lty=1:2,cex=0.5)


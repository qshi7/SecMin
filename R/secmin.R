# NOTE:
# 1. modify output
# 2. dataset dimension match with the users' raw data
# 3. make the cheating cases output as a matrix for both LZT and KLD
# 4. make seperate arguments/extract.function for calling different plots (make vars definable? in plot())
# 5. make plots a seperate function from secmin, which is mainly for estimate results

#' secmin
#'
#' secmin fit the lognormal model and return parameter summaries and plots
#'
#' @param RT data frame containing RT information
#' @param Resp data frame containing item response
#'
#' @import runjags
#' @import coda
#' @importFrom graphics par abline hist
#' @import PerFit
#' @importFrom stats qchisq quantile
#' @import ltm
#'
#' @return a list
#' @export

secmin = function(RT,Resp) {

  model <- "
  Poisson model...

model {
	for (i in 1:N) {
		for (j in 1:m) {
			RT[i,j] ~ dlnorm(beta[j] - speed[i], alpha[j]^2)
		}
		speed[i] ~ dnorm(0, prec.speed)
	}
		prec.speed ~ dgamma(.001,.001)
		sd.speed <- sqrt(1/prec.speed)

	for (j in 1:m) {
		alpha[j] <- ipar[j,1]
		beta[j] <- ipar[j,2]
		ipar[j,1:2] ~ dmnorm(mu.ipar[1:2], prec.ipar[1:2,1:2])
	}

  mu <- c(1,3.5)
  V[1,1] <- 5
  V[2,1] <- 0
  V[1,2] <- 0
  V[2,2] <- 5

	mu.ipar[1:2] ~ dmnorm(mu[1:2], prec.ipar[1:2,1:2])
	prec.ipar[1:2,1:2] ~ dwish(V[1:2,1:2], 2)
	var.ipar[1:2,1:2] <- inverse(prec.ipar[1:2,1:2])
}

Data{
list(RT = structure(.Data = RT,
.Dim = c(nrow(RT), ncol(RT))),N = nrow(RT), m = ncol(RT))
}

Inits{
list( speed=speed_hat, prec.speed=prec.speed,ipar=cbind(alpha_hat,beta_hat), mu.ipar=mu.ipar, prec.ipar=prec.ipar)
}

Inits{
list( speed=speed_hat, prec.speed=prec.speed,ipar=cbind(alpha_hat,beta_hat), mu.ipar=mu.ipar, prec.ipar=prec.ipar)

}

"
out <- run.jags(model, n.chains = 2, burnin = 1000,
                sample = 8000, adapt = 2000, monitor=c("alpha","beta","mu.ipar","var.ipar","speed","sd.speed"),
                )

out.summary <- extract.runjags(add.summary(out),what = 'summary')
summary.table <- out.summary[c(1:40,47:146),c(1,3:5,9,11)]
plots = plot(out, vars='alpha[1]',plot.type = c('trace','histogram','autocorr','ecdf')) #make vars definable for users


### Convert output to mcmc objects
#n.burn <- 1000
#thin.int <- 1
#s <- n.burn+thin.int
#sample = 8000

### Create mcmc objects

#up.out.speed <- mcmc(out2$BUGSoutput$sims.list$speed, start=s, thin=thin.int)
#up.out.sd.speed <- mcmc(out$BUGSoutput$sims.list$sd.speed, start=s, thin=thin.int)

#up.out.alpha <- mcmc(out$BUGSoutput$sims.list$alpha, start=s, thin=thin.int)
#up.out.beta <- mcmc(out$BUGSoutput$sims.list$beta, start=s, thin=thin.int)
#up.out.mu <- mcmc(out$BUGSoutput$sims.list$mu.ipar, start=s, thin=thin.int)
#up.out.Var <- out$BUGSoutput$sims.list$var.ipar
#temp <- cbind(up.out.Var[,1,1], up.out.Var[,2,2], up.out.Var[,1,2])
#up.out.var <- mcmc(temp, start=s, thin=thin.int)

### Create summary objects
#sum.speed <- summary(up.out.speed)
#sum.sd.speed <- summary(up.out.sd.speed)

#sum.alpha <- summary(up.out.alpha)
#sum.beta <- summary(up.out.beta)
#sum.mu <- summary(up.out.mu)
#sum.var <- summary(up.out.var)

#####################################################
### Compute/save parameter estimates
### Person parameters
#ppars<- as.data.frame(round(sum.speed$statistics[,1:2],3)) # function output
#names(ppars) <- c('tau.hat','tau.sd')

### Item parameters
#ipars <- as.data.frame(round(cbind(sum.alpha$statistics[,1:2], sum.beta$statistics[,1:2]),3)) # function output
#names(ipars) <- c('alpha.hat','sd.alpha','beta.hat','sd.beta')

#####################################################
#####################################################
### MCMC convergence stuff --- will be updated in ver. 2

### Geweke's convergence diagnostic --- function output
#round(geweke.diag(up.out.sd.speed, frac1=.1, frac2=.5)$z, digits=2)
#round(geweke.diag(up.out.mu, frac1=.1, frac2=.5)$z, digits=2)
#round(geweke.diag(up.out.var, frac1=.1, frac2=.5)$z, digits=2)

#gew <- cbind(
#  geweke.diag(up.out.alpha, frac1=.1, frac2=.5)$z,
#  geweke.diag(up.out.beta, frac1=.1, frac2=.5)$z)
#round(ifelse(abs(gew) > 1.9,gew,NA),digits=2)

##################################################### PFS calculation               ####
##################### IRT based PFS calculation      ##########################
###### HT index
Ht.psf.A <- PerFit::Ht(Resp)
Ht.cut.A <- PerFit::cutoff(Ht.psf.A,Blvl=.05,Nrep=1000)
HT.plot <- plot(Ht.psf.A, cutoff.obj=Ht.cut.A,Type="Both",Blv1=0.05,CIlv1=0.9,Xcex=0.8,col.hist="grey",col.area="NA",col.ticks="NA", title="",Xlabel="")

###### Lzstar  index

Lzstar.psf.A <- lzstar(Resp,IRT.PModel = "1PL")
Lzstar.cut.A <- cutoff(Lzstar.psf.A,Blvl=.05,Nrep=1000)
Lzstar.plot <- plot(Lzstar.psf.A, cutoff.obj=Lzstar.cut.A,Type="Both",Blv1=0.05,CIlv1=0.9,Xcex=0.8,col.hist="grey",col.area="NA",col.ticks="NA",title="")


###### NCI index
NCI.psf.A <- NCI(Resp)
NCI.cut.A <- cutoff(NCI.psf.A ,Blvl=.05,Nrep=1000)
NCI.plot <- plot(NCI.psf.A, Type="Both",Blv1=0.05,CIlv1=0.9,Xcex=0.8,col.hist="grey",col.area="NA",col.ticks="NA",title="")


##################### RT based PFS calculation       ##########################
#### Lt index Marianti & Fox
rt.alpha.A <- as.vector(summary.table[c(1:20,3)]) # need to modify according to data
rt.beta.A <- as.vector(summary.table[c(21:40,3)]) # need to modify according to data
rt.tau.A <- as.matrix(summary.table[c(41:140),3]) # need to modify according to data
Form.A.RT <- as.matrix(RT)
I <- ncol(Form.A.RT)
J <- nrow(Form.A.RT)
Z.A=matrix(0,nrow=J,ncol=I)
for(i in 1:I){
  for(j in 1:J){
    Z.A[j,i]<- ((log(Form.A.RT[j,i])-(rt.beta.A[i]-rt.tau.A[j]))*rt.alpha.A[i])^2
  }
}
Z_sum.A <- apply(Z.A,1,sum)
hist(Z_sum.A,col="grey",main = "",xlab="Lz",breaks=20) ## add plot
abline(v=qchisq(0.95,20)) # 170 should be the number of items


# KL index programming

Cheating=Form.A.RT
Cheating[is.na(Cheating)] <- 0
average_cheating <- (apply(Cheating,2,sum)/J)
average_cheating
w2 <- sum(average_cheating)
w2
w_average_cheat <- average_cheating*(1/w2)
w_average_cheat
w1 <- apply(Cheating,1,sum)

KLdata <- diag(1/w1) %*% Cheating

summary(KLdata)
start.time <- Sys.time()
KLD=rep(0,J)
for (i in 1:J){
  KLD[i]= sum(w_average_cheat %*% log(w_average_cheat/KLdata[i,]));
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
Critical <-quantile(KLD,probs=c(0.5))+1.5*(quantile(KLD,probs=c(0.75))-quantile(KLD,probs=c(0.5)))
Critical
KLD.plot = c(hist(KLD),abline(v=Critical)) # function output

LZT.cheating.cases <- which(Z_sum.A>Critical)
KLD.cheating.cases <- which(KLD>Critical) # function output



return(list(summary = summary.table, LZT.cheating.cases,KLD.cheating.cases))

}

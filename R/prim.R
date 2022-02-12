#' prim
#'
#' prim fit the lognormal model and return parameter estimates
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
#plots = plot(out, vars='alpha[1]',plot.type = c('trace','histogram','autocorr','ecdf')) #make vars definable for users




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



return(summary =(list(unlist(out),summary.table, LZT.cheating.cases,KLD.cheating.cases)))

}

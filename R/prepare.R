#' prepare
#'
#' prepare prepare the datasets and create values in the Global Environment to be used in secmin later
#'
#' @param RT data frame containing RT information
#' @param Resp data frame containing item response
#'
#' @export

prepare <- function(RT, Resp) {
  RT[RT== 0] <- NA
  num.RT <<- RT
  RT <<- as.matrix(RT)
  num.RT[is.na(num.RT)] <<-1
  logRTs <<- log(num.RT)
  alpha_hat <<- 1/apply(logRTs,2,sd)
  beta_hat <<- apply(logRTs,2,mean)
  speed_hat <<- mean(beta_hat) - apply(logRTs,1,mean)
  Resp <<- as.data.frame(Resp)

  prec.speed <<- 1
  mu.ipar <<- c(mean(alpha_hat),mean(beta_hat))
  prec.ipar <<- diag(1,2,2)
}




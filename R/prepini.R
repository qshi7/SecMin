#' prepini
#'
#' prepini prepares the initial values to be used in the lognormal model later
#' and saves the RT dataset as \code{matrix} and the Resp dataset as a \code{data.frame}
#'
#' @param RT A required \eqn{N × J} \code{matrix} or \code{data.frame}
#' containing the response times of \eqn{N} individuals to \eqn{J} items. Missing
#' values need to be coded as \code{NA}.
#' @param Resp A required \eqn{N × J} \code{matrix} or \code{data.frame}
#' containing the item responses of \eqn{N} individuals to \eqn{J} items. The
#' correct responses are marked by 1 and incorrect responses by 0. Missing values
#' need to be coded as \code{NA}.
#'
#' @return \code{prepini} returns two datasets and several initial values
#'   \item{RT}{The RT data is saved as a matrix with zero values marked as \code{NA}.}
#'   \item{Resp}{The Resp data is saved as a \code{data.frame}.}
#'   \item{ini.alpha}{A vector containing initial alpha values for each item.}
#'   \item{ini.beta}{A vector containing initial beta values for each item.}
#'   \item{ini.speed}{A vector containing initial person speed values for each individual.}
#'   \item{ini.prec.speed}{The initial value for the precision of speed}
#'   \item{ini.mu.ipar}{The initial values for the means of the item parameters.}
#'   \item{ini.prec.ipar}{The initial values for the variance-covariance matrix of the item parameters.}
#'
#' @export

prepini <- function(RT, Resp) {
  RT.new <- as.matrix(RT)
  RT.new[RT.new== 0] <- NA
  num.RT <- RT.new
  num.RT[is.na(num.RT)] <-1
  logRTs <- log(num.RT)

  ini<- list(RT <- RT.new,
             Resp <- as.data.frame(Resp),
             al_hat <- 1/apply(logRTs,2,sd),
             be_hat <- apply(logRTs,2,mean),
             sp_hat <- mean(be_hat) - apply(logRTs,1,mean),
             prec.sp <- 1,
             mu.ipar <- c(mean(al_hat),mean(be_hat)),
             prec.ipar <- diag(1,2,2)
             )
  names(ini) <- c("RT", "Resp", "ini.alpha", "ini.beta", "ini.speed",
                  "ini.prec.speed", "ini.mu.ipar", "ini.prec.ipar")

  class(ini) <- "ini"
  invisible(ini)

}



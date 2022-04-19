#' dv
#'
#' dv gives histogram of RT distribution as well as the proportion of correctness
#'
#' @param RT A required \eqn{N × J} \code{matrix} containing the response
#' times of \eqn{N} individuals to \eqn{J} items. Zero and missing values need
#' to be coded as \code{NA}.
#' @param Resp A required \eqn{N × J} \code{matrix} or \code{data.frame}
#' containing the item responses of \eqn{N} individuals to \eqn{J} items. The
#' correct responses are marked by 1 and incorrect responses by 0. Missing values
#' need to be coded as \code{NA}.
#' @param color A non-required \code{character string} concerning the color to be used to fill the bars.
#'
#' @import graphics
#'
#' @return a 3*4 matrix of plots containing the RT distribution of each item as well as the proportion of correctness.
#' @export
#'
#' @examples
#' x <- as.data.frame(matrix(c(1,2,3,4,5,6), nrow = 2))
#' dv(Resp,RT, color="grey")


dv<- function(RT, Resp, color ="lightgray"){
  Resp<- as.data.frame(Resp)
  RT<- as.data.frame(RT)
  nresp<-nrow(Resp)
  nrt<- ncol(RT)
  loop <- 1: nrt
  p <- colSums(Resp)
  p1<- round(p/nresp,2)

  old.par <- par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1),oma =c(0,0,0,0))

  par(mfrow = c(3,4), oma =c(0.5,0.5,0.5,0.5), mar=c(2,2,3,2))

  for (i in loop){
    x<-RT[,i]
    hist(x,
         main = paste("RT Distr. for Item",i),col= color,
         breaks = 15)
    mtext(paste("Correct Proportion:", p1[i]), cex= 0.5, col=2)

  }
  par(old.par)
}


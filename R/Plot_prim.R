#' plot_prim
#'
#' plot_prim plot the four graphics for any assigned individual
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
#' @return plot
#' @export

plots = plot(out, vars='alpha[3]',plot.type = c('trace','histogram','autocorr','ecdf'))


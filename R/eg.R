#' eg
#'
#' eg gives description of the dataset
#'
#' @param object a matrix or vector
#'
#' @import psych
#'
#' @return an output
#' @export
#'
#' @examples
#' x <- matrix(c(1,2,3,4,5,6), nrow = 2)
#' summary(x)

eg <- function(object) {
  describe(object)
}

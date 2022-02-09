#' histall
#'
#' histall gives histogram of all numeric variables
#'
#' @param data a data frame
#'
#' @import purrr
#' @import ggplot2
#' @import tidyr
#'
#' @return histogram
#' @export
#'
#' @examples
#' x <- as.data.frame(matrix(c(1,2,3,4,5,6), nrow = 2))
#' histall(x)

histall <- function(data) {
  data %>%
    keep(is.numeric) %>%
    gather() %>%
    ggplot(aes(value)) +
    facet_wrap(~key,scales = "free") +
    geom_histogram(bins=30)
}

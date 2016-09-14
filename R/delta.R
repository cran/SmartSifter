#' delta
#'
#' the delta function
#'
#'
#' @name delta
#' @param a A number.
#' @param b A number.
#' @return The function delta(\code{a},\code{b}), if \code{a} == \code{b}, return 1; else, return 0.


delta <- function (a, b)
{
  if (a == b) {result = 1}
  else {result = 0}
  result
}


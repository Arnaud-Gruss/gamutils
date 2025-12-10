#' @title Mode of a categorical variable
#'
#' @description This function computes the mode of a categorical variable
#'
#' @param x categorical variable of interest.
#'
#' @return the mode of a categorical variable
#'
#' @export
modefunc <- function( x ) {

  a <- table( x )
  names( a[a == max( a )] )

}

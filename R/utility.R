#' @title Mode of a categorical variable
#'
#' @description This function computes the mode of a categorical variable.
#'
#' @param x categorical variable of interest.
#'
#' @return The mode of a categorical variable.
#'
#' @export
mode_factor <- function( x ) {

  a <- table( x )
  names( a[a == max( a )] )

}

#' @title Maximum element from a vector
#'
#' @description This function finds the maximum element from a vector.
#'
#' @param x vector of interest
#'
#' @return The maximum element from the vector of interest.
#'
#' @export
get_max <- function( x ) {

  if ( all( is.na( x ) ) ) {
    ( x[1] )
  } else {
    x[which( x == max( x, na.rm = T ) )][1]
  }

}

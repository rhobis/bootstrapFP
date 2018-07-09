#' Check if a number is integer
#'
#' Check if \code{x} is an integer number, differently from \code{is.integer},
#' which checks the type of the object \code{x}
#'
#' @param x a scalar or a numeric vector
#' @param tol a scalar, indicating the tolerance
#'
#'
#' @note From the help page of function \code{\link[base]{is.integer}}
#'
#' @keywords internal

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
}

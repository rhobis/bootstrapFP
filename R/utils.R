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




#' Select the random part of a pseudo-population
#' 
#' Helper function that generates the fixed part of a pseudo-population in 
#' function /code{ppBS_srs()}.
#' 
#' @param ... parameters of the function, depending on the bootstap method chosen.
#' @param method string indicating the bootstrap method
#' 
#' @keywords internal


select_Uc <- function(..., method){
    
    # Initialisation
    method <- match.arg(method, c('Gross'
                                  ,'Booth'
                                  ,'ChaoLo85'
                                  ,'ChaoLo94'
                                  ,'BickelFreedman'
                                  ,'Sitter'
    )
    )
    
    arguments <- list(...)
    y <- ifelse(is.null(arguments$y), NULL, arguments$y)
    n <- ifelse(is.null(arguments$n), NULL, arguments$n)
    N <- ifelse(is.null(arguments$N), NULL, arguments$N)
    k <- ifelse(is.null(arguments$k), NULL, arguments$k)
    
    tol <- .Machine$double.eps^0.5 #tolerance to use with floor()
    
    # Generate random part of bootstrap pseudo-population
    if( identical(method, 'Gross') ){
        
        return( NULL )
        
    } else if( identical(method, 'Booth') ){
        
        # k  <- floor( N/n + tol)
        s <- as.logical( sampling::srswor( (N-n*k), n) )
        return(y[s]) 
        
    } else if( identical(method, 'ChaoLo94') ){
        
        # k <- floor( N/n + tol)
        if(N-n*k > 0){
            s <- as.logical( sampling::srswr( (N-n*k), n) )    
        }else return(NULL)
        
        return( y[s] )
        
    } else if( identical(method, 'ChaoLo85') ){
        # k <- floor( N/n + tol)
        G <- function(x) return( (1 - n/x) * x*(n-1) / ((x-1)*n) )
        q <- (G(N) - G(n*(k+1)) ) / (G(n*k) - G( n*(k+1)) )
        if( q <= runif(1) ){
            return( NULL )
        } else return( y )
        
    } else if( identical(method, 'BickelFreedman') ){
        
        # k <- floor( N/n + tol)
        q <- ( 1 - (N-n*k)/n ) * ( 1 - (N-n*k)/(N-1) )
        if( q <= runif(1) ){
            return( NULL )
        } else return(y)
        
    } else if( identical(method, 'Sitter') ){
        
        f  <- n/N
        # k  <- floor( (N/n)*(1 - (1-f)/n) + tol)
        a1 <- ( n*k - n + 1 ) / ( n*(n-1)*(n*k - 1) )
        a2 <- k / ( n * ( n*(k+1) - 1 ) )
        q  <- (1-f)/(n*(n-1)) - a2
        q  <- q / ( a1-a2 )
        
        
        if( q <= runif(1) ){
            return( NULL )
        } else return( y )
        
        
    } 
    
}

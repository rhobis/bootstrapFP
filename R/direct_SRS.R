#' Direct bootstrap methods for simple random sampling
#'
#' @param y vector of sample values
#' @param N scalar, representing the population size
#' @param B scalar, number of bootstrap replications
#' @param method a string indicating the bootstrap method to be used, available
#'     methods are: 'Efron', 'McCarthySnowden', 'RaoWu', 'Sitter'.
#'
#' @details
#' See Mashreghi et al. (2016) for details about the algorithm.
#' 
#' @references
#' Mashreghi Z.; Haziza D.; Léger C., 2016. A survey of bootstrap methods in 
#' finite population sampling. Statistics Surveys 10 1-52.
#' 
#' 
#' @importFrom stats rbinom



directBS_srs <- function(y, N, B, method){
    
    ### Check input ---
    method <- match.arg(method,
                        c('Efron',
                          'McCarthySnowden',
                          'RaoWu',
                          'Sitter')
    )
    
    n <- length(y)
    f <- n/N
    
    if( identical(method, 'Efron') ){
        CC <- 1
        nn <- 1
        k  <- n
    }else if( identical(method, 'McCarthySnowden')){
        CC  <- 1
        nn  <- 1
        k   <- round( (n-1)/(1-f) )
    }else if( identical(method, 'RaoWu')){
        k  <- n-1
        CC <- sqrt( (k*(1-f)/(n-1)) )
        nn <- 1
    }else if( identical(method, 'Sitter')){
        CC <- 1
        nn <- n/(2-f)
        
        ff <- nn/n
        k  <- n*(1-ff)/(nn*(1-f)) 
        if( is_wholenumber(k)) k <- round(k)
        kf <- floor(k)
        kc <- kf + 1
        q  <- (1/kf - 1/k) / (1/kf - 1/kc)
        u  <- rbinom(1,1,q)
        k  <- ceiling( n*(1-ff) / (nn*(1/f)) ) + u
    }
    n1 <- k*nn
    
    ym <- mean(y)
    y1 <- ym + CC*(ym - y)
    
    tot <- vector('numeric', length = B)
    for(b in seq_len(B)){
        s  <- sapply(seq_len(k), function(i) sample(y1, size=nn, replace=FALSE ) )
        
        tot[b] <- N/n1 * sum(s)
    }
    
    return( var(tot) )
}

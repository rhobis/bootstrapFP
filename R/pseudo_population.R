#' Pseudo-population bootstrap for simple random sampling
#'
#' @param y vector of sample values
#' @param N scalar, represents the population size
#' @param B scalar, number of bootstrap replications
#' @param D scalar, number of replications for the double bootstrap (when applicable)
#' @param method a string indicating the bootstrap method to be used, available
#'     methods are: 'Gross', 'Booth', 'ChaoLo85', 'ChaoLo94', 'BickelFreedman', 'Sitter'
#'
#' @details
#' See Mashreghi et al. (2016) for details about these bootstrap methods.
#' 
#' @references
#' Mashreghi Z.; Haziza D.; Léger C., 2016. A survey of bootstrap methods in 
#' finite population sampling. Statistics Surveys 10 1-52.
#' 
#' @keywords internal

ppBS_srs <- function(y, N, B, D=1, method) {
    
    ### Check input ---
    method <- match.arg(method,
                        c('Gross',
                          'Booth',
                          'ChaoLo85',
                          'ChaoLo94',
                          'BickelFreedman',
                          'Sitter')
    )
    
    if( identical(method, 'Gross') & !is_wholenumber(N/n)) 
        stop("Gross method can be used only when N/n is an integer. Please, choose another method!")
    
    
    ### Initialise quantities ---
    tol <- .Machine$double.eps^0.5  
    n <- length(y)
    f <- n/N
    k <- switch(method,
                Gross    = N/n,
                Booth    = floor(N/n + tol),
                ChaoLo85 = floor(N/n) + tol,
                ChaoLo94 = floor(N/n + tol),
                BickelFreedman = floor(N/n + tol),
                Sitter = floor( (N/n)*(1 - (1-f)/n) + tol)
                )

    # finite part of pseudo-population
    Uf  <- rep(y, each = k )
    
    ### Bootstrap replicates ---
    Vb <- vector('numeric', length = D)
    for(d in seq_len(D) ){
        
        Uc <- select_Uc(y=y, N=N, n=n, k=k, method = method)
        U  <- c( Uf, Uc)
        n1 <- ifelse(is.null(Uc) & identical(method, 'Sitter'), n-1, n)
        
        tot <- vector('numeric', length = B)
        ### Bootstrap procedure
        for(b in seq_len(B)){
            #select resample
            sb <- U[ sample(length(U), n1, replace = FALSE) ]
            tot[b] <- N * mean(sb)
        }
        Vb[d] <- var(tot)
    }
    
    ### Return results ---
    return( mean(Vb) )
}





#' Pseudo-population bootstrap for simple random sampling
#'
#' @param y vector of sample values
#' @param pik vector of sample first-order inclusion probabilities
#' @param x vector of length N with values of the auxiliary variable for all population units,
#'     only required if method "HotDeck" is chosen
#' @param s logical vector of length N, TRUE for units in the sample, FALSE otherwise. 
#'     Alternatively, a vector of length n with the indices of the sample units.
#'     Only required for "HotDeck" method.
#' @param B scalar, number of bootstrap replications
#' @param D scalar, number of replications for the double bootstrap
#' @param method a string indicating the bootstrap method to be used, available
#'     methods are: 'Gross', 'Booth', 'ChaoLo85', 'ChaoLo94', 'BickelFreedman', 'Sitter'
#' @param smplFUN a function that takes as input a vector of length N of 
#' inclusion probabilities and return a vector of length N, either logical or a 
#' vector of 0s and 1s,  where \code{TRUE} or \code{1} indicate sampled
#' units and \code{FALSE} or \code{0} indicate non-sample units.
#' 
#' 
#' @references
#' Mashreghi Z.; Haziza D.; Léger C., 2016. A survey of bootstrap methods in 
#' finite population sampling. Statistics Surveys 10 1-52.
#' 
#' 
#' @keywords internal

ppBS_ups <- function(y, pik, B, D=1, method, smplFUN, x = NULL, s = NULL) {
    
    ### Check input ---
    method <- match.arg(method, c('Holmberg', 'Chauvet', 'HotDeck') )
    
    
    ### Initialisation ---
    ## Create finite part of pseudo-population
    n <- length(y)
    if( identical(method, 'HotDeck') ){
        D  <- 1
        U  <- vector('numeric', length = length(x) )
        xs <- x[s]
        j  <- sapply(x[!s], function(xi) which.min( abs(xi-xs) ) )
        U[s]  <- y
        U[!s] <- y[j]
        xb <- x
        xb[!s] <- xs[j]
        pkb <- sampling::inclusionprobabilities(xb, n)
    } else {
        r    <- 1/pik
        rint <- as.integer(r)
        rp   <- ( r - rint )
        Uf   <- rep(y, times = rint )
    }
    
    ### Bootstrap replicates ---
    Vb <- vector('numeric', length = D)
    for(d in seq_len(D) ){
        tot <- vector('numeric', length = B)
        for(b in seq_len(B)){
            ## random part of pseudo-population
            if( method %in% c('Holmberg', 'Chauvet') ){
                Uc  <- as.logical( sampling::UPpoisson(rp) )
                U  <- c(Uf, y[Uc])
                pkb <- c( rep(pik, times=rint), pik[Uc] )
                if( identical(method, 'Holmberg') )
                    pkb <- drop( sampling::inclusionprobabilities(pkb, n) )
            }
            
            #select resample
            sb <- as.logical(smplFUN( pkb ))
            tot[b] <- drop( sampling::HTestimator(U[sb], pkb[sb]) )
        }
        Vb[d] <- var(tot)
    }
    
    ### Return results ---
    return( mean(Vb) )
    
}



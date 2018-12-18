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
    n <- length(y)
    f <- n/N
    k <- switch(method,
                Gross    = N/n,
                Booth    = floor(N/n),
                ChaoLo85 = floor(N/n),
                ChaoLo94 = floor(N/n),
                BickelFreedman = floor(N/n),
                Sitter = floor( (N/n)*(1 - (1-f)/n) )
                )
    # # quantities for the random population, Uc
    # if( identical(method, 'Gross') ){
    #     
    #     k <- N/n 
    #     select_Uc <- function( ... ) return( NULL )
    #     
    # } else if( identical(method, 'Booth') ){
    #     
    #     k  <- floor( N/n )
    #     select_Uc <- function( ... ){
    #         s <- as.logical( sampling::srswor( (N-n*k), n) )
    #         return( y[s] )
    #     }
    #     
    # } else if( identical(method, 'ChaoLo94') ){
    #     
    #     k <- floor( N/n )
    #     select_Uc <- function( ... ){
    #         s <- as.logical( sampling::srswr( (N-n*k), n) )
    #         return( y[s] )
    #     }
    #     
    # } else if( identical(method, 'ChaoLo85') ){
    #     
    #     k <- floor( N/n )
    #     G <- function(x) return( (1 - n/x) * x*(n-1) / ((x-1)*n) )
    #     q <- (G(N) - G(n*(k+1)) ) / (G(n*k) - G( n*(k+1)) )
    #     
    #     select_Uc <- function( ... ){
    #         if( q <= runif(1) ){
    #             return( NULL )
    #         } else return( y )
    #     }
    #     
    # } else if( identical(method, 'BickelFreedman') ){
    #     
    #     k <- floor( N/n )
    #     q <- ( 1 - (N-n*k)/n ) * ( 1 - (N-n*k)/(N-1) )
    #     select_Uc <- function(){
    #         if( q <= runif(1) ){
    #             return( NULL )
    #         } else return(y)
    #     }
    #     
    # } else if( identical(method, 'Sitter') ){
    #     
    #     f  <- n/N
    #     k  <- floor( (N/n)*(1 - (1-f)/n) )
    #     a1 <- ( n*k - n + 1 ) / ( n*(n-1)*(n*k - 1) )
    #     a2 <- k / ( n * ( n*(k+1) - 1 ) )
    #     q  <- (1-f)/(n*(n-1)) - a2
    #     q  <- q / ( a1-a2 )
    #     
    #     select_Uc <- function(){
    #         if( q <= runif(1) ){
    #             return( NULL )
    #         } else return( y )
    #     }
    #     
    # } 
    
    
    # finite part of pseudo-population
    Uf  <- rep(y, each = k )
    
    ### Bootstrap replicates ---
    Vb <- vector('numeric', length = D)
    for(d in seq_len(D) ){
        
        Uc <- select_Uc(y=y, N=N, n=n, method = method)
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
#' @param design sampling procedure to be used for sample selection.
#'        Either a string indicating the name of the sampling design or a function;
#'        see section "Details" for more information.
#'
#' @details
#' Argument \code{design} accepts either a string indicating the sampling design
#' to use to draw samples or a function.
#' Accepted designs are "brewer", "tille", "maxEntropy", "poisson",
#' "sampford", "systematic", "randomSystematic".
#' The user may also pass a function as argument; such function should take as input
#' a vector of length N of inclusion probabilities and return a vector of length N,
#' either logical or a vector of 0s and 1s,  where \code{TRUE} or \code{1} indicate sampled
#' units and \code{FALSE} or \code{0} indicate non-sample units.
#' 
#' 
#' @references
#' Mashreghi Z.; Haziza D.; Léger C., 2016. A survey of bootstrap methods in 
#' finite population sampling. Statistics Surveys 10 1-52.
#' 

ppBS_ups <- function(y, pik, B, D=1, method, design, x = NULL, s = NULL) {
    
    ### Check input ---
    method <- match.arg(method, c('Holmberg', 'Chauvet', 'HotDeck') )
    
    
    if( !is.function(design) & !identical(method, 'Chauvet') ){
        design <- match.arg(design, c('brewer',
                                      'tille',
                                      'maxEntropy',
                                      'randomSystematic',
                                      'sampford',
                                      'poisson',
                                      'systematic')
        )
    }else if( identical(method, 'Chauvet') & !identical(design, 'poisson') ){
        design <- 'poisson'
        message( paste0("Sampling design set to 'Poisson', if your sample has been drawn with ",
                        "a different design, please choose a different bootstrap method!") )
    }
    
    
    ### Initialisation ---
    if( is.character(design)){
        # sampling function
        smplFUN <- switch(EXPR=design,
                          'brewer' = sampling::UPbrewer,
                          'tille' = sampling::UPtille,
                          'maxEntropy' = sampling::UPmaxentropy,
                          'randomSystematic' = sampling::UPrandomsystematic,
                          'sampford' = sampling::UPsampford,
                          'poisson' = 
                              function(pik){
                                  ss <- 0
                                  while(ss < 2){
                                      s  <- sampling::UPpoisson( pik )
                                      ss <- sum(s)
                                  }
                                  return( s )
                              },
                          'systematic' = sampling::UPsystematic
        )
    }else if( is.function(design) ){
        smplFUN <- design
    }else stop("Argument design is not well-specified: it should be either a string representing ",
               "one of the available sampling designs or an object of class function!")
    
    ## Create finite part of pseudo-population
    n <- length(y)
    if( identical(method, 'HotDeck') ){
        U <- vector('numeric', length = length(x) )
        xs <- x[s]
        j  <- sapply(x[!s], function(xi) which.min( abs(xi-xs) ) )
        U[s] <- y
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
            sb <- smplFUN( pkb )
            tot[b] <- drop( sampling::HTestimator(U[sb], pkb[sb]) )
        }
        Vb[d] <- var(tot)
    }
    
    ### Return results ---
    return( mean(Vb) )
    
}



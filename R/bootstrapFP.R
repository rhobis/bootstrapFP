#' Bootstrap algorithms for Finite Population sampling
#' 
#' Bootstrap variance estimation for finite population sampling.
#'
#'
#'
#' @param y vector of sample values
#' @param pik vector of sample first-order inclusion probabilities
#' @param B scalar, number of bootstrap replications
#' @param D scalar, number of replications for the double bootstrap (when applicable)
#' @param method a string indicating the bootstrap method to be used, see Details for more
#' @param design sampling procedure to be used for sample selection.
#'        Either a string indicating the name of the sampling design or a function;
#'        see section "Details" for more information.
#' @param x vector of length N with values of the auxiliary variable for all population units,
#'     only required if method "ppHotDeck" is chosen
#' @param s logical vector of length N, TRUE for units in the sample, FALSE otherwise. 
#'     Alternatively, a vector of length n with the indices of the sample units.
#'     Only required for "ppHotDeck" method.
#' @param distribution required only for \code{method='generalised'}, a string
#'     indicating the distribution to use for the Generalised bootstrap. 
#'     Available options are "uniform", "normal", "exponential" and "lognormal"
#'        
#'
#'
#' @details
#' Argument \code{design} accepts either a string indicating the sampling design
#' to use to draw samples or a function.
#' Accepted designs are "brewer", "tille", "maxEntropy", "poisson",
#' "sampford", "systematic", "randomSystematic", "srs", chao", "sunter", "srs", "fpdust", "fpdust_pps".
#' The user may also pass a function as argument; such function should take as input
#' the parameters passed to argument \code{design_pars} and return either a logical
#' vector or a vector of 0 and 1,  where \code{TRUE} or \code{1} indicate sampled
#' units and \code{FALSE} or \code{0} indicate non-sample units.
#' The length of such vector must be equal to the length of \code{x}
#' if \code{units} is not specified, otherwise it must have the same length of \code{units}.
#'
#'
#' \code{method} must be a string indicating the bootstrap method to use.
#' A list of the currently available methods follows, the sampling design they
#' they should be used with is indicated in square brackets.
#' The prefix "pp" indicates a pseudo-population method, the prefix "d"
#' represents a direct method, and the prefix "w" inicates a weights method.
#' For more details on these methods see Mashreghi et al. (2016).
#' 
#' \itemize{
#'     \item "ppGross" [SRSWOR]
#'     \item "ppBooth" [SRSWOR]
#'     \item "ppChaoLo85" [SRSWOR]
#'     \item "ppChaoLo94" [SRSWOR]
#'     \item "ppBickelFreedman" [SRSWOR]
#'     \item "ppSitter" [SRSWOR]
#'     
#'     \item "ppHolmberg" [UPSWOR]
#'     \item "ppChauvet"  [UPSWOR]
#'     \item "ppHotDeck"  [UPSWOR]
#'     
#'     \item "dEfron" [SRSWOR]
#'     \item "dMcCarthySnowden" [SRSWOR]
#'     \item "dRaoWu" [SRSWOR]
#'     \item "dSitter" [SRSWOR]
#'     \item "dAntalTille_ups" [UPSWOR]
#'     
#'     \item "wRaoWuYue"    [SRSWOR]
#'     \item "wChipperfieldPreston"    [SRSWOR]
#'     \item "wGeneralised" [any]
#' 
#' } 
#'
#'
#' @return 
#' The bootstrap variance of the Horvitz-Thompson estimator.
#'
#'
#' @examples
#'
#' library(bootstrapFP)
#' 
#' ### Generate population data ---
#' N   <- 20; n <- 5
#' x   <- rgamma(500, scale=10, shape=5)
#' y   <- abs( 2*x + 3.7*sqrt(x) * rnorm(N) )
#' pik <- n * x/sum(x)
#' 
#' ### Draw a dummy sample ---
#' s  <- sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(.1, .9) )
#' n  <- sum(s)
#' 
#' ### Estimate bootstrap variance ---
#' bootstrapFP(y = y[s], pik = n/N, B=100, method = "ppSitter")
#' bootstrapFP(y = y[s], pik = pik[s], B=100, method = "ppHolmberg", design = 'brewer')
#' bootstrapFP(y = y[s], pik = pik[s], B=100, D=10, method = "ppChauvet")
#' bootstrapFP(y = y[s], pik = pik[s], B=100, method = "wGeneralised", distribution = 'normal')
#' 
#' 
#' 
#'
#'
#'
#' @references
#' Mashreghi Z.; Haziza D.; Léger C., 2016. A survey of bootstrap methods in 
#' finite population sampling. Statistics Surveys 10 1-52.
#' 
#' 
#' 
#' 
#' @export
#' 
#' @import sampling 
#'




bootstrapFP <- function(y, pik, B, D=1, method, design, x=NULL, s=NULL, distribution='uniform' ){
    
    
    ### Check input ---
    
    method <- match.arg(method, 
                        c( 'ppGross', 
                           'ppBooth', 
                           'ppChaoLo85', 
                           'ppChaoLo94',
                           'ppBickelFreedman',
                           'ppSitter',
                           'ppHolmberg',
                           'ppChauvet',
                           'ppHotDeck',
                           'dEfron',
                           'dMcCarthySnowden',
                           'dRaoWu',
                           'dSitter',
                           'dAntalTille_ups',
                           'wRaoWuYue',
                           'wChipperfieldPreston',
                           'wGeneralised'
                        )
    )
    
    
    if( identical(method, 'wGeneralised') ){
        distribution <- match.arg(distribution, 
                                  c("uniform", 
                                    "normal", 
                                    "exponential", 
                                    "lognormal"))
    }
    
    
    
    
    # if( identical(method, 'Gross') & !is_whole(N/n)) 
    #     stop("Gross method can be used only when N/n is an integer. Please, choose another method!")
    
    
    
    
    
    
    n <- length(y)
    lp <- length(pik)
    
    if( !identical( class(pik), "numeric" ) ){
        stop( "The argument 'pik' should be a numeric vector!")
    }else if( lp < 2 & !(method %in% c('ppGross', 'ppBooth', 'ppChaoLo85', 
                                       'ppChaoLo94', 'ppBickelFreedman',
                                       'ppSitter') )){
        stop( "The 'pik' vector is too short!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some 'pik' values are outside the interval [0, 1]")
    }else if( any(pik %in% c(NA, NaN, Inf)) ){
        stop( "The 'pik' vector contains invalid values (NA, NaN, Inf)" )
    }
    
    if( !(class(y) %in% c("numeric", "integer")) ){
        stop( "The argument 'y' should be a numeric vector!")
    }else if( n < 2 ){
        stop( "The 'y' vector is too short!" )
    }else if( any(y %in% c(NA, NaN, Inf)) ){
        stop( "The 'y' vector contains invalid values (NA, NaN, Inf)" )
    }
    
    if( any(y<0) ){
        message( "Some 'y' values are negative, continuing anyway...")
    }
    
    
    
    ## pseudo-population, srs
    if( method %in% c('ppGross',
                      'ppBooth',
                      'ppChaoLo85', 
                      'ppChaoLo94',
                      'ppBickelFreedman',
                      'ppSitter') ){
        if( length(unique(pik)) > 1 ) stop("pik values should be all equal!")
        N <- (1/pik) * n
    }
    
    ## pseudo-population, ups
    
    
    
    
    ### Initialise variables ---
    n <- length(y)
    
    ### Bootstrap ---
    
    out <- switch(method, 
                  'ppGross' = ppBS_srs(y, N, B, D, method = 'Gross'), 
                  'ppBooth' = ppBS_srs(y, N, B, D, method = 'Booth'), 
                  'ppChaoLo85' = ppBS_srs(y, N, B, D, method = 'ChaoLo85'), 
                  'ppChaoLo94' = ppBS_srs(y, N, B, D, method = 'Chaolo94'), 
                  'ppBickelFreedman' = ppBS_srs(y, N, B, D, method = 'BickelFreedman'), 
                  'ppSitter' = ppBS_srs(y, N, B, D, method = 'Sitter'),
                  'ppHolmberg' = ppBS_ups(y, pik, B, D, method = 'Holmberg', design),
                  'ppChauvet' =  ppBS_ups(y, pik, B, D, method = 'Chauvet', design = 'poisson'),
                  'ppHotDeck' =  ppBS_ups(y, pik, B, D, method = 'HotDeck', design, x=x, s=s),
                  # 'dEfron',
                  # 'dMcCarthySnowden',
                  # 'dRaoWu',
                  # 'dSitter',
                  # 'dAntalTille_ups',
                  # 'wRaoWuYue',
                  # 'wChipperfieldPreston',
                  'wGeneralised' = generalised(y, pik, B, distribution = distribution)
    )
    
    ### Return
    return( out )
    
}
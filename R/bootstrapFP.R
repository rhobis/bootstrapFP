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
#'
#'
#'
#'
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




bootstrapFP <- function(y, pik, B, D=1, method, design ){
    
    
    ### Check input ---
    
    method <- match.arg(method, 
                        c( 'ppGross', 
                           'ppBooth', 
                           'ppChaoLo85', 
                           'ppChaoLo94', 
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
    
    
    
    ly <- length(y)
    lp <- length(pik)
    
    if( !identical( class(pik), "numeric" ) ){
        stop( "The argument 'pik' should be a numeric vector!")
    }else if( lp < 2 ){
        stop( "The 'pik' vector is too short!" )
    }else if( any(pik<0)  | any(pik>1) ){
        stop( "Some 'pik' values are outside the interval [0, 1]")
    }else if( any(pik %in% c(NA, NaN, Inf)) ){
        stop( "The 'pik' vector contains invalid values (NA, NaN, Inf)" )
    }
    
    if( !(class(y) %in% c("numeric", "integer")) ){
        stop( "The argument 'y' should be a numeric vector!")
    }else if( ly < 2 ){
        stop( "The 'y' vector is too short!" )
    }else if( any(y %in% c(NA, NaN, Inf)) ){
        stop( "The 'y' vector contains invalid values (NA, NaN, Inf)" )
    }
    
    
    if( any(y<0) ){
        message( "Some 'y' values are negative, continuing anyway...")
    }
    
    
    ### Bootstrap ---
    
    out <- switch(method, 
                  'ppGross', 
                  'ppBooth', 
                  'ppChaoLo85', 
                  'ppChaoLo94', 
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
    
    ### Return
    return( out )
    
}
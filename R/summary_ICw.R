#' @importFrom jtools md_table
#' @S3method summary ICw
#' @export summary.ICw
#' @export

summary.ICw <- function(x, digits = NULL)
{
  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
    sig.digits <- TRUE
  }else{
    NrDigits <- digits
    sig.digits <- FALSE
  }

  DF <- data.frame(IC = x$IC,
                   IC.weights = x$IC.weights,
                   rel.IC.weights = x$rel.IC.weights)

  md_table(DF, digits = NrDigits, sig.digits = sig.digits, align = 'c')

}

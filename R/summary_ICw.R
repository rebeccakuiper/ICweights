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
    align <- NULL
  }else{
    NrDigits <- digits
    sig.digits <- FALSE
    align <- 'c'
  }

  DF <- data.frame(IC = x$IC,
                   IC.weights = x$IC.weights,
                   ratio.IC.weights = x$ratio.IC.weights)

  cat("\n")
  cat("Per hypothesis/model, the information criterion value (IC), its weight (IC.weights), and its ratio / relative support (ratio.IC.weights) versus the other hypotheses/models: \n")
  cat("\n")
  print(md_table(DF, digits = NrDigits, sig.digits = sig.digits, align = align))

}

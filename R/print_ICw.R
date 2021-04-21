#' @S3method print ICw
#' @export print.ICw
#' @export

print.ICw <- function(x, digits = NULL)
{

  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
  }else{
    NrDigits <- digits
  }

  DF <- data.frame(IC = x$IC,
                   IC.weights = x$IC.weights,
                   ratio.IC.weights = x$ratio.IC.weights)

  cat("\n")
  cat("Per hypothesis/model, the information criterion value (IC), its weight (IC.weights), and its ratio / relative support (ratio.IC.weights) versus the other hypotheses/models: \n")
  cat("\n")
  print(DF, digits = NrDigits, right = F)

  return(invisible(DF))
  #return(invisible(x))

}

# IC <- myIC; ICw <- IC.weights(IC); print(ICw)

#' @S3method print ICw
#' @export print.ICw
#' @export

print.ICw <- function(x)
{

  x <- as.list(x)
  DF <- data.frame(IC = x$IC,
                   IC.weights = x$IC.weights,
                   rel.IC.weights = x$rel.IC.weights)
  print(DF, digits = 4, right = F)

  return(invisible(DF))

}

# IC <- myIC; ICw <- IC.weights(IC); print(ICw)
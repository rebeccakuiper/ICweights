#' @importFrom jtools md_table
#' @S3method summary ICw
#' @export summary.ICw
#' @export

summary.ICw <- function(x)
{
  x <- as.list(x)
  DF <- data.frame(IC = x$IC,
                   IC.weights = x$IC.weights,
                   rel.IC.weights = x$rel.IC.weights)

  md_table(DF, digits = 2, align = 'c') # TO DO set nr digit as input value

}

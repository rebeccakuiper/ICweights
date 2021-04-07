#'  @importFrom jtools md_table
#'  @export

summary.ICw <- function(x)
{
  x <- as.data.frame(x)

  MD_DF <- md_table(x, digits = 2, align = 'c') # TO DO set nr digit as input value

  class(MD_DF) <- "summary.ICw"
  MD_DF

}

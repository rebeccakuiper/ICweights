
#' Calculating IC weights based on IC values (AIC, ORIC, GORIC(A), BIC, SIC, ...)
#'
#' This function transforms IC values into IC weights: IC values denote the ordering of hypotheses/models, while IC weights quantify the relative strength of hypotheses/models.
#'
#' @param IC A vector or one-column matrix with information criteria (AIC, ORIC, GORIC(A), BIC, SIC, ...) values of length 'NrHypos', where 'NrHypos' stands for the number of hypotheses/models.
#' @param Name_Hypo Optional. Vector containing 'NrHypos' characters which will be used for labeling the hypothesis. Default: H1, H2, ....
#'
#' @return IC weights, which quantify the relative strength of hypotheses/models.
#' @export print.ICw
#' @export summary.ICw
#' @export
#' @examples
#'
#' # library(ICweights)
#'
#' IC <- myIC # Example based on 3 hypotheses.
#' IC.weights(IC)
#'
#'
#' # Change labels of hypotheses #
#' #
#' # Example 1: Let us say that we evaluated H1: beta1 > beta2; H2: beta1 = beta2; and H3: beta1 < beta2.
#' Name_Hypo <- c("Higher", "Equal", "Smaller")
#' IC.weights(IC, Name_Hypo)
#' #
#' # Example 2: Let us say that we evaluated a linear, quadratic, and cubic model.
#' # Notably, this can also be represented by
#' # H1: beta_linear, beta_quadratic = 0, beta_cubic = 0;
#' # H2: beta_linear, beta_quadratic, beta_cubic = 0; and
#' # Hunc: beta_linear, beta_quadratic, beta_cubic (i.e., no restrictions on the beta's).
#' Name_Hypo <- c("Linear", "Quadratic", "Cubic")
#' IC.weights(IC, Name_Hypo)
#'
#'
#' # Output options #
#' #
#' # Example 3: This examples shows some output options
#' ICw <- IC.weights(IC)
#' print(ICw)
#' summary(ICw)
#' print(ICw, digits = 4)
#' summary(ICw, digits = 4)
#' # In Rstudio, use 'ICw$' to see what output there is. The following is available:
#' #' ICw$IC
#' ICw$IC.weights
#' ICw$ratio.IC.weights
#'
#'
#' # PT weights #
#' #
#' # Example 4: This examples shows how to calculate PT weights.
#' # Notably, one is interested in PT weights when the log likelihood for two or more hypotheses are (approximately) equal.
#' # Then, the comparison between those hypotheses is solely based on the PT values.
#' # The IC weights will then equal the PT weights.
#' # In that case, there is support for the overlap (boundary) of these hypotheses.
#' # Thus, when the IC weights equal the PT weights for a (sub)set of hypotheses,
#' # then there is support for the overlap (boundary) of these hypotheses.
#' #
#' # Example code:
#' # est <- coef(fit.object)
#' # VCOV_est <- vcov(fit.object)
#' # results <- goric(est, VCOV = VCOV_est, H1, comparison = "complement", type = "gorica")
#' # IC.weights(2*results$result[,3])$IC.weights
#' ## Note that the penalty is 2*`penalty'
#' ## So, do not use: IC.weights(results$result[,3])$IC.weights
#'




IC.weights <- function(IC, Name_Hypo = NULL) {

  # Checks on input
  #
  if(length(dim(IC)) != 2 & !is.vector(IC)){
    print(paste0("The argument IC should either be a vector or a matrix (with one column)."))
    stop()
  }
  if(length(dim(IC)) == 2){
    if(dim(IC)[2] != 1){
      print(paste0("The argument IC can be a matrix but should have one column then (or it should be a vector)."))
      stop()
    }
  }
  NrHypos <- length(IC)
  #
  if(is.null(Name_Hypo)){
    Name_Hypo <- paste0("H", 1:NrHypos)
  }
  #
  if(length(Name_Hypo) != NrHypos){
    print(paste("The argument 'Name_Hypo' should consist of NrHypos = ", NrHypos, " elements (all characters)."))
    stop()
    if(!all(is.numeric(Name_studies)) & !all(is.character(Name_studies))){
      print(paste("The argument 'Name_Hypo' should consist of solely characters (NrHypos = ", NrHypos, " characters)."))
      stop()
    }
  }

  names(IC) <- Name_Hypo
  #
  minIC <- min(IC)
  weight_m <- exp(-0.5*(IC-minIC)) / sum(exp(-0.5*(IC-minIC)))
  names(weight_m) <- Name_Hypo
  #
  ratio.IC.weights <- weight_m %*% t(1/weight_m)
  rownames(ratio.IC.weights) <- Name_Hypo
  colnames(ratio.IC.weights) <- paste0("vs ", Name_Hypo)


  # Ouput
  #DF <- data.frame(IC = IC,
  #           IC.weights = weight_m,
  #           ratio.IC.weights = ratio.IC.weights)
  #class(DF) <- c("ICw", "data.frame")
  #DF

  final <- list(IC = IC,
              IC.weights = weight_m,
              ratio.IC.weights = ratio.IC.weights)
  class(final) <- c("ICw", "list")
  final
  #print.ICw(final)
  #return(invisible(final))

}




#' Calculating IC weights based on IC values (AIC, ORIC, GORIC(A), BIC, SIC, ...)
#'
#' This function transforms IC values into IC weights: IC values denote the ordering of hypotheses/models, while IC weights quantify the relative strength of hypotheses/models.
#'
#' @param IC A vector or one-column matrix with information criteria (AIC, ORIC, GORIC(A), BIC, SIC, ...) values of length 'NrHypos', where 'NrHypos' stands for the number of hypotheses/models.
#' @param Name_Hypo Optional. Vector containing 'NrHypos' characters which will be used for labelling the hypothesis. Default: H1, H2, ....
#'
#' @return IC weights, which quantify the relative strength of hypotheses/models.
#' @export
#' @examples
#'
#' IC <- myIC # Example based on 3 hypotheses.
#' IC.weights(IC)
#'
#' # Change labels of hypotheses #
#' # For example, let us say that we tested a linear model vs quadratic vs cubic.
#' Name_Hypo <- c("Linear", "Quadratic", "Cubic")
#' IC.weights(IC, Name_Hypo)


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


  weight_m <- matrix(NA, nrow = 1, ncol = NrHypos)
  minIC <- min(IC)
  weight_m <- exp(-0.5*(IC-minIC)) / sum(exp(-0.5*(IC-minIC)))
  rel.IC.weights <- weight_m %*% t(1/weight_m)
  rownames(rel.IC.weights) <- Name_Hypo
  colnames(rel.IC.weights) <- paste0("vs ", Name_Hypo)


  # Ouput
  final <- list(IC = IC,
                IC.weights = weight_m,
                rel.IC.weights = rel.IC.weights)
  return(final)

}

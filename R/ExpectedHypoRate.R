
#' Calculating the expected (true) hypothesis rate(s).
#'
#' This function calculates the expected (true) hypothesis rate for the hypotheses in the set under a specific hypothesis (say, H0 or H1 or even another hypothesis).
#'
#' @param mu A vector denoting the population value of the (standardized) parameters. It is assumed that the (standardized) parameter estimates come from a multivariate normal distribution with mean mu.
#' @param VCOV A square matrix denoting the covariance matrix of the (standardized) parameters. It is assumed that the (standardized) parameter estimates come from a multivariate normal distribution with covariance matrix VCOV.
#' @param NrHypos The number of theory-based hypotheses in the set of candidate hypotheses (is a scalar with an integer value).
#' @param Hypos A vector of strings containing the NrHypos theory-based hypotheses.
#' @param Safeguard Indicator of which safeguard-hypothesis should be used: "unconstrained" (default; i.e., all possible theories including the one specified), "none" (only advised when set of hyptheses cover all theories), or (only when 'NrHypos = 1') "complement" (i.e., the remaining theories).
#' @param nsim Optional. The number of iterations used in the calculation. By default, nsim = 1000. The higher nsim, the higher the precision but als the higher the computation time.
#' @param seed Optional. The seed value used in the calculation. By default, seed = 123. This can be used for a sensitivty check. If for another seed value the results unnegligibly differ, then nsim has to be increased.
#'
#' @return This function renders 3 types of output.
#' The first (mle_Hypos) is the percentage of times the sampled (maximum likelihood) estimates are in agreement with the hypotheses in the set, when sampled from a multivariate normal distribution with mean mu and covariance matrix VCOV. This is independent from the GORICA, in a way purely based on the log likelhood. Bear in mind that the estimates are always in agreement with the unconstrained hypothesis.
#' The second (HR_Hypos) is the percentage of times each hypothesis in the set is preferred under that distribution; rendering the Hypothesis Rate (HR) for H1. This is based on the GORICA and thus both the log likelihood and the penalty are taken into account.
#' The third are the average GORICA weight values for the hypotheses in the set; thus either for H1 and its complement or H1 and the unconstrained hypothesis.
#' @importFrom restriktor goric
#' @importFrom MASS mvrnorm
#' @export
#' @examples
#'
#' # Example based on 6 standardized parameters of interest.
#' p <- 6
#' mu <- rep(0, p) # mu is in accordance with H0: X1 = X2 = X3 = X4 = X5 = X6 (or H0 <- "X1 == X2; X3 == X4; X5 == X6")
#' VCOV <- diag(p) # Here, the variances of the parameters are 1 and the covariances are 0
#' NrHypos <- 1
#' H1 <- "X1 < X2; X3 < X4; X5 < X6"
#' Hypos <- c(H1)
#' Safeguard <- "complement"
#' ExpectedHypoRate(mu, VCOV, NrHypos, Hypos, Safeguard)
#'
#' # Other covariance matrix of the standardized parameters #
#' # If we assume a correlation between the standardized parameters of rho = 0.25 and a sample size of 100, then
#' rho <- 0.25
#' n <- 100
#' sigma <- matrix(rho, ncol = p, nrow = p)
#' for(i in 1:p){
#'  sigma[i,i] <- 1
#' }
#' VCOV <- sigma / n
#' ExpectedHypoRate(mu, VCOV, NrHypos, Hypos, Safeguard)
#'
#' # Other mean vector of the standardized parameters #
#' # In the example above, it is assumed that the truth is a null hypothesis (mu = vector of 0's).
#' # One can also assume another distribution, e.g. stating that the standardized estimates have a specific effect:
#' mu <- c(-0.1, 0.1, -0.1, 0.1, -0.1, 0.1) # mu is now in accordance with H1 <- "X1 < X2; X3 < X4; X5 < X6"
#' # One can do this for multiple values corresponding to specific effect sizes for example.
#' ExpectedHypoRate(mu, VCOV, NrHypos, Hypos, Safeguard)
#'
#' # Multiple hypotheses of interest #
#' NrHypos <- 2
#' H1 <- "X1 < X2; X3 < X4; X5 < X6"
#' H2 <- "X1 > X2; X3 > X4; X5 > X6"
#' Hypos <- c(H1, H2) # Currently, in case of multiple theory-based hypotheses in the set, one cannot use the complement of these as safeguard, so:
#' Safeguard <- "unconstrained" # which is the default, so could be left out
#' ExpectedHypoRate(mu, VCOV, NrHypos, Hypos, Safeguard)
#'

ExpectedHypoRate <- function(ES, VCOV, NrHypos, Hypos, Safeguard = "unconstrained", nsim = 1000, seed = 123) {

  if(length(VCOV) == 1 & length(mu) != 1){
    print(paste("The covariance matrix VCOV is 1 times 1, while the mean vector mu does not consist of 1 element. Change one of these such that VCOV is p times p matrix and mu is a p-vector."))
    stop()
  }else  if(length(VCOV) > 1){
    if(length(dim(VCOV)) > 2){
      print(paste("The covariance matrix VCOV should be an p times p matrix (i.e., a square matrix) and not an array with more than 2 dimensions."))
      stop()
    }else if(dim(VCOV)[1] != dim(VCOV)[2]){
      print(paste("The covariance matrix VCOV should be a square matrix, that is, an p times p matrix."))
      stop()
    }else if(dim(VCOV)[1] != length(mu)){
      print(paste("The covariance matrix VCOV is p times p, while the mean vector mu does not consist of p elements. Change one of these such that VCOV is p times p matrix and mu is a p-vector."))
      stop()
    }
  }
  #
  if(length(NrHypos) != 1){
    print(paste("The number of hypotheses (NrHypos) should be a scalar (integer)."))
    stop()
  }
  if(length(NrHypos) == 1){
    if(NrHypos %% 1 != 0){
      print(paste("The number of hypotheses (NrHypos) should be an integer value."))
      stop()
    }
  }
  #
  if(length(Hypos) != NrHypos){
    print(paste("The argument (Hypos) should consist of NrHypos  = ", NrHypos, " elements (which are character strings / text elements)."))
    stop()
  }
  if(!is.character(Hypos)){
    print(paste("The argument (Hypos) does not contain character strings / text elements."))
    stop()
  }
  #
  if(NrHypos == 1){
    if(length(Safeguard) != 1){
      print(paste("The type of safeguard-hypothesis (Safeguard) should be one word ('unconstrained', 'none', or (if NrHypos = 1) 'complement')."))
      stop()
    }
    if(Safeguard != "unconstrained" & Safeguard != "none" & Safeguard != "complement"){
      print(paste("The type of safeguard-hypothesis (Safeguard) should be 'unconstrained', 'none', or (if NrHypos = 1) 'complement'."))
      stop()
    }
  }


  for(HypoTeller in 1:NrHypos){
    eval(parse(text = paste0("H", HypoTeller, " <<- Hypos[(HypoTeller)]")))
  }
  HypoSet <- noquote(paste0("H", 1:NrHypos, collapse = ", "))

  nrhypos <- NrHypos + (Safeguard != "none")


  # Generate data from multivariate normal distribution with mean mu and covariance matrix VCOV
  set.seed(seed)
  X <- mvrnorm(nsim, mu, VCOV, empirical=FALSE)

  LL <- matrix(NA, nrow = nsim, ncol = nrhypos)
  PT <- matrix(NA, nrow = nsim, ncol = nrhypos)
  weight <- matrix(NA, nrow = nsim, ncol = nrhypos)
  for (i in 1:nsim) {
    #i <- 1
    est <- X[i,]
    names(est) <- paste0("X", 1:p)
    #result <- goric(est, VCOV = VCOV, HypoSet, comparison = Safeguard, type = 'gorica')
    eval(parse(text = paste0("result <- restriktor:::goric(est, VCOV = VCOV, ",
                             HypoSet,
                             ", comparison = Safeguard, type = 'gorica')")))
    #
    LL[i,] <- result$result[,2]
    PT[i,] <- result$result[,3]
    weight[i,] <- result$result[,5]

  }

  index <- 1:nsim
  mle_Hypos <- matrix(NA, nrow = 1, ncol = nrhypos)
  HR_Hypos <- matrix(NA, nrow = 1, ncol = nrhypos)
  if(nrhypos > 2){
    for(HypoTeller in 1:nrhypos){
      #HypoTeller = 1
      mle_Hypos[HypoTeller] <- mean(apply(LL[index, HypoTeller] >= LL[index, -HypoTeller], 1, all)) # % of times mle's in agreement with H1
      HR_Hypos[HypoTeller] <- mean(apply(weight[index, HypoTeller] > weight[index, -HypoTeller], 1, all)) # % of times H1 preferred
    }
  }else{
    for(HypoTeller in 1:nrhypos){
      #HypoTeller = 1
      mle_Hypos[HypoTeller] <- mean(LL[index, HypoTeller] >= LL[index, -HypoTeller]) # % of times mle's in agreement with H1
      HR_Hypos[HypoTeller] <- mean(weight[index, HypoTeller] > weight[index, -HypoTeller]) # % of times H1 preferred
    }
  }
  av.GORICA.weights <- matrix(colMeans(weight[index,]), nrow = 1, ncol = nrhypos)  # average weight values for the hypotheses in the set

  if(Safeguard == "none"){
    colnames(mle_Hypos) <- colnames(HR_Hypos) <- colnames(av.GORICA.weights) <- paste0("H", 1:nrhypos)
  }else{
    colnames(mle_Hypos) <- colnames(HR_Hypos) <- colnames(av.GORICA.weights) <- c(paste0("H", 1:(nrhypos-1)), Safeguard)
  }
  rownames(mle_Hypos) <- rownames(HR_Hypos) <- rownames(av.GORICA.weights) <- ""

  # Ouput
  final <- list(mle_Hypos = mle_Hypos,
                HR_Hypos = HR_Hypos,
                av.GORICA.weights = av.GORICA.weights)
  return(final)

}


## --------------------------------------------- ##
#                                                 #
#       sMM                                       #
#       rMM                                       #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#                       sMM(Km, Vm, ...)                                            #
## ------------------------------------------------------------------------------- ##
#' Kinetics for An Enzyme of Known Parameters
#' @description Simulates the kinetics of an michaelian enzyme of unknown parameters.
#' @usage sMM(Km, Vm, unit_S = 'mM', unit_v = 'au', replicates = 3, error = a, sd = 0.05*Vm, outlier = 0, plot = TRUE)
#' @param Km Michaelis contant.
#' @param Vm maximal velocity.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param replicates number of replicates for the dependent variable.
#' @param error it should be one among c('fixed', 'absolute', 'relative').
#' @param sd standard deviation of the error.
#' @param outlier an integer, when different from zero, a value of v choseen at random will be affected by outlier-fold its corresponding error.
#' @param plot, logical. If TRUE data are plotted.
#' @detail The error type 'fixed' assumes that all velocity values are affected by the same fixed error. The option 'absolute' implies normally distributed error of constante absolute magnitud. For this obtion, random numbers of mean zero and the SD set by the user are added to the computed values of v. The parameter 'relative' carries an error that is proportional to the value of v. For that, each v is multiplied by a number of mean unity and the SD set.
#' @return Returns a dataframe where the two first column are the substrate concentration and the velocity without error. The last two columns are the mean velocity and its standard error. In between, when required, the different velocity replicates affected of error are given.
#' @author Juan Carlos Aledo
#' @examples sMM(Km = 1, Vm = 1)
#' @seealso rMM()
#' @export

sMM <- function(Km, Vm, unit_S = 'mM', unit_v = 'au', replicates = 3,
                error = 'a', sd = 0.05*Vm, outlier = 0, plot = TRUE, ...){

  set.seed(123)

  ## --------- Choosing values for independent variable, [S] ---------- ##
  z <- list(...)
  if (!is.null(z$S)){
    S <- z$S
  } else {
    S <- c(seq(from = Km/10, to = 2*Km, by = 0.4*Km), 5*Km, 10*Km)
  }
  ## ---------------- Formating the output dataframe ------------------ ##
  output <- as.data.frame(matrix(rep(NA, length(S)*(replicates + 2)),
                                 ncol = replicates + 2))
  names(output) <- c('S', 'v', LETTERS[1:replicates])
  output$S <- S

  ## ----------- Computing the dependent variable, v ------------------ ##
  v <- Vm*S/(Km + S)
  output$v <- round(v, 2)

  for (i in 1:nrow(output)){
    for (j in 1:replicates){
      if (error == 'a' | error == 'absolute'){
        ve <- v[i] + rnorm(1, mean = 0, sd = sd)
      } else if (error == 'r' | error == 'relative'){
        ve <- v[i] * rnorm(1, mean = 1, sd = sd)
      } else if (error == 'f' | error == 'fixed'){
        ve <- v[i] + ((-1)^j) * sd
      } else {
        stop ("A proper error type must be provided")
      }
      output[i, j+2] <- round(ve, 2) # v for the replicate j-th
    }
    if (error == 'f' | error == 'fixed'){
      if (replicates %% 2 == 1){ output[i, j+2] <- v[i] }
    }
  }

  ## -------------------- Dealing with outliers ----------------------- ##
  if (outlier != 0){
    i <- sample(1:nrow(output), 1)
    j <- sample(1:replicates, 1)
    if (error == 'a'){
      output[i, j+2] <- v[i] + rnorm(n = 1, mean = 0, sd = outlier*sd)
    } else if (error == 'r'){
      output[i, j+2] <- v[i] * rnorm(n = 1, mean = 1, sd = outlier*sd)
    } else if (error == 'f'){
      output[i, j+2] <- v[i] + ((-1)^sample(c(1,2), 1)) * outlier*sd
    }
  }
  ol <- c(i,j)

  ## --------------- Computing mean and sd if required ---------------- ##
  output[output<0] <- 0
  if (ncol(output) > 3){
    v <- output[,-c(1,2)]
    output$v_mean <- round(apply(v, MARGIN = 1, mean), 2)
    output$v_sd <- round(apply(v, MARGIN = 1, sd), 2)
  } else if (ncol(output) == 3){
    output$v_mean <- output$A
    output$v_sd <- 0
  } else {
    output$v_mean <- output$v
    output$v_sd <- 0
  }

  ## ---------------- Returning and showing results ------------------ ##
  if (plot){
    plot(c(0,output$S), c(0,output$v_mean), ty = 'b',
         ylim = c(0, Vm + 0.1*Vm),
         xlab = paste('[S] (', unit_S, ')', sep = ""),
         ylab = paste('v (', unit_v, ')', sep = ""))
    abline(h = Vm, lty = 2)
    abline(v = Km, lty = 2)

    if (replicates > 1){
      arrows(output$S, output$v_mean-output$v_sd,
             output$S, output$v_mean+output$v_sd,
             length=0.05, angle=90, code=3)
    }
  }

  attr(output, 'outlier') <- ol
  return(output)
}


## ------------------------------------------------------------------------------- ##
#           rMM(Km_range, Vm_range, ...)                                            #
## ------------------------------------------------------------------------------- ##
#' Kinetics for An Enzyme of Unknown Parameters
#' @description Simulates the kinetics of an michaelian enzyme of unknown parameters.
#' @usage rMM(Km_range = c(1,10), Vm_range = c(10,100), unit_S = 'mM', unit_t = 'au', replicates = 1, error = 'a', sd = 0.1*Vm, outlier = 0, plot = TRUE)
#' @param Km_range range of possible values for Km.
#' @param Vm_range range of possible values for Vm.
#' @param Km Michaelis contant.
#' @param Vm maximal velocity.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param replicates number of replicates for the dependent variable.
#' @param error it should be one among c('fixed', 'absolute', 'relative').
#' @param sd standard deviation of the error.
#' @param outlier an integer, when different from zero, a value of v choseen at random will be affected by outlier-fold its corresponding error.
#' @param plot, logical. If TRUE data are plotted.
#' @detail The error type 'fixed' assumes that all velocity values are affected by the same fixed error. The option 'absolute' implies normally distributed error of constante absolute magnitud. For this obtion, random numbers of mean zero and the SD set by the user are added to the computed values of v. The parameter 'relative' carries an error that is proportional to the value of v. For that, each v is multiplied by a number of mean unity and the SD set.
#' @return Returns a dataframe where the two first column are the substrate concentration and the velocity without error. The last two columns are the mean velocity and its standard error. In between, when required, the different velocity replicates affected of error are given.
#' @author Juan Carlos Aledo
#' @examples
#' @seealso
#' @export

rMM <- function(Km_range = c(1,10),
                Vm_range = c(10,100),
                unit_S = 'mM',
                unit_v = 'au',
                replicates = 1,
                error = 'a',
                sd = 0.1*Vm,
                outlier = 0,
                plot = TRUE){

  set.seed(123)

  Km_population <- seq(from = Km_range[1],
                       to = Km_range[2],
                       by = (Km_range[2] - Km_range[1])/20)

  Vm_population <- seq(from = Vm_range[1],
                       to = Vm_range[2],
                       by = (Vm_range[2] - Vm_range[1])/20)

  Km <- sample(Km_population, size = 1)
  Vm <- sample(Vm_population, size = 1)

  ## --------- Choosing values for independent variable, [S] ---------- ##
  S <- c(seq(from = Km/10, to = 2*Km, by = 0.4*Km), 5*Km, 10*Km)

  ## ---------------- Formating the output dataframe ------------------ ##
  output <- as.data.frame(matrix(rep(NA, length(S)*(replicates + 2)),
                                 ncol = replicates + 2))
  names(output) <- c('S', 'v', LETTERS[1:replicates])
  output$S <- S

  ## ----------- Computing the dependent variable, v ------------------ ##
  v <- Vm*S/(Km + S)
  output$v <- round(v, 2)

  for (i in 1:nrow(output)){
    for (j in 1:replicates){
      if (error == 'a' | error == 'absolute'){
        ve <- v[i] + rnorm(1, mean = 0, sd = sd)
      } else if (error == 'r' | error == 'relative'){
        ve <- v[i] * rnorm(1, mean = 1, sd = sd)
      } else if (error == 'f' | error == 'fixed'){
        ve <- v[i] + ((-1)^j) * sd
      }
      output[i, j+2] <- ve
    }
    if (replicates %% 2 == 1){
      output[i, j+2] <- v[i]
    }
  }

  ## -------------------- Dealing with outliers ----------------------- ##
  if (outlier != 0){
    i <- sample(1:nrow(output), 1)
    j <- sample(1:replicates, 1)
    if (error == 'a'){
      output[i, j+2] <- v[i] + rnorm(n = 1, mean = 0, sd = outlier*sd)
    } else if (error == 'r'){
      output[i, j+2] <- v[i] * rnorm(n = 1, mean = 1, sd = outlier*sd)
    } else if (error == 'f'){
      output[i, j+2] <- v[i] + ((-1)^sample(c(1,2), 1)) * outlier*sd
    }
  }
  ol <- c(i,j)

  ## --------------- Computing mean and sd if required ---------------- ##
  output[output<0] <- 0
  if (ncol(output) > 3){
    v <- output[,-c(1,2)]
    output$v_mean <- round(apply(v, MARGIN = 1, mean), 2)
    output$v_sd <- apply(v, MARGIN = 1, sd)
  } else if (ncol(output) == 3){
    output$v_mean <- round(output$A, 2)
    output$v_sd <- 0
  } else {
    output$v_mean <- round(output$v, 2)
    output$v_sd <- 0
  }

  ## ---------------- Returning and showing results ------------------ ##
  if (plot){
    plot(output$S, output$v_mean, ty = 'b',
         ylim = c(0, Vm+(0.1*Vm)),
         xlab = paste('[S] (', unit_S, ')', sep = ""),
         ylab = paste('v (', unit_v, ')', sep = ""))
    abline(h = Vm, lty = 2)
    abline(v = Km, lty = 2)

    if (replicates > 1){
      arrows(output$S, output$v_mean-output$v_sd,
             output$S, output$v_mean+output$v_sd,
             length=0.05, angle=90, code=3)
    }

  }

  attr(output, 'Km') <- Km
  attr(output, 'Vm') <- Vm
  return(output)
}

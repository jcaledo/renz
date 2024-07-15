## --------------------------------------------- ##
#                                                 #
#                    int.MM                       #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#                            int.MM(data)                                           #
## ------------------------------------------------------------------------------- ##
#' Linearization of The Integrated Michaelis-Menten Equation
#' @description Estimates the kinetic parameters using an linearized form of the integrated Michaelis-Menten equation.
#' @usage int.MM(data, unit_S = 'mM', unit_t = 'min')
#' @param data a dataframe with two columns. The first column contains the values of the independent variable time, t, and the second column contains the substrate concentrations.
#' @param unit_S concentration unit.
#' @param unit_t time unit.
#' @details The r-squared value of the model can be checked using attributes().
#' @return A list of two elements. The first element is named vector containing the Km and Vm. The second element is a dataframe where the first two columns are the original data and the last two columns are the transformed variables. Also a linear plot of the transformed variables together with the parameters values are provided.
#' @examples int.MM(data = sE.progress(So = 10, time = 5, Km = 4, Vm = 50)[, c(1,3)])
#' @importFrom stats lm
#' @importFrom graphics abline
#' @export

int.MM <- function(data, unit_S = 'mM', unit_t = 'min'){

  So <- data[1,2]
  ## ---------------------- Transformed independent variable --------------------- ##
  # x = (So - St)/t
  data$x <- (So - data[,2])/data[,1]

  ## ----------------------- Transformed dependent variable ---------------------- ##
  # y = (1/t)ln(So/St)
  data$y <- log(So/data[,2])/data[,1]

  ## ----------------------------- Model fitting --------------------------------- ##

  data <- data[is.finite(rowSums(data)),]
  # data$weights <- 1/(data[,3])^weighting # inverse of the variance
  model <- lm(data$y ~ data$x)
  # if (weighting){
  #   model <- lm(data$y ~ data$x, weights = data$weights)
  # } else {
  #   model <- lm(data$y ~ data$x)
  # }

  ## ------------------------ Computuing parameters ------------------------------ ##
  Km <- round(unname(-1/model$coefficients[2]), 3)
  Vm <- round(unname(Km * model$coefficients[1]), 3)

  ## ------------------- Plotting the transformed variables ---------------------- ##
  parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
  plot(data$x, data$y,
       xlab = paste("(So - St)/t (", unit_S, "/", unit_t, ")", sep = ""),
       ylab = paste("(1/t) log(So/St) (1/", unit_t, ")", sep = ""),
       main = parameters)
  abline(model)

  ## --------------------------------- Output ------------------------------------ ##
  KmVm <- c(Km, Vm)
  names(KmVm) <- c("Km", "Vm")

  output <- list(KmVm, data)
  names(output) <- c('parameters', 'data')

  attr(output, 'x') <- '(So - St)/t'
  attr(output, 'y') <- '(1/t) log(So/St)'
  attr(output, 'r-squared') <- summary(model)$r.squared

  return(output)
}


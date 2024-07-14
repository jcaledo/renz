## --------------------------------------------- ##
#                                                 #
#                  dir.MM                         #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#             dir.MM(data, unit_S = 'mM', unit_v = 'au', plot = T)                  #
## ------------------------------------------------------------------------------- ##
#' Non-linear Least-squares Fitting of the MM equation
#' @description Non-linear least-squares fitting of the Michaelis-Menten equation.
#' @usage dir.MM(data, unit_S = 'mM', unit_v = 'au', plot = TRUE)
#' @param data a dataframe with two columns. The first column contains the values of the independent variable (substrate concentration), and the second column contains the initial rates.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param plot logical. If true, the data and fitted curve are plotted.
#' @details This function invokes nls() to carry out the fitting.
#' @return A list of two elements. The first one is a vector containing the enzyme kinetic parameters. The second one is a dataframe with the original data plus the fitted value of v.
#' @examples dir.MM(ONPG[, c(1,2)])
#' @importFrom stats complete.cases
#' @importFrom stats nls
#' @importFrom graphics points
#' @export

dir.MM <- function(data, unit_S = 'mM', unit_v = 'au', plot = TRUE){

  ## ---------------------------- Removing incomplete data ----------------------- ##
  data <- data[complete.cases(data), ]

  ## ---------------------------- Estimating the seed ---------------------------- ##
  t <- ecb(data, plot = FALSE)
  K <- t$fitted_parameters[1]
  V <- t$fitted_parameters[2]
  seed = list(Km = K, Vm = V)

  ## ----------------------------- Fitting the curve ----------------------------- ##
  names(data) <- c('S', 'v')
  model <- nls(data$v ~ (Vm * data$S)/(Km + data$S), data =  data, start = seed )

  ## --------------------------- Computuing parameters --------------------------- ##
  Km <- round(summary(model)$coefficient[1,1], 3)
  sd_Km <- round(summary(model)$coefficient[1,2], 3)

  Vm <- round(summary(model)$coefficient[2,1], 3)
  sd_Vm <- summary(model)$coefficient[2,2]

  ## --------------------------- Fitted velocity lues ---------------------------- ##
  mm.eq <- function(x) {(Vm * x)/(Km + x)}
  data$fitted_v <- mm.eq(data$S)

  ## ------------------- Plotting the transformed variables ---------------------- ##
  if (plot){
    parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
    plot(data$S, data$v,
         xlab = paste("[S] (", unit_S, ")", sep = ""),
         ylab = paste("v (", unit_v, ")", sep = ""),
         main = parameters)

    x <- seq(from = 0, to = max(data$S), by = max(data$S)/1000)
    y <- mm.eq(x)
    points(x, y, ty = 'l')

  }

  ## --------------------------------- Output ------------------------------------ ##
  KmVm <- c(Km, Vm)
  names(KmVm) <- c("Km", "Vm")

  output <- list(KmVm, data)
  names(output) <- c('parameters', 'data')

  return(output)
}


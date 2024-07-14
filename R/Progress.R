## --------------------------------------------- ##
#                                                 #
#             sE.progress                         #
#             fE.progress                         #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#    sE.progress(So, time, Km, Vm, unit_S, unit_t, I, Kic, Kiu, replicates, ...)    #
## ------------------------------------------------------------------------------- ##
#' Progress Curve for Enzyme-Catalyzed Reaction
#' @description Simulates the evolution of the substrate concentration along time.
#' @usage sE.progress(So, time, Km, Vm, unit_S = 'mM', unit_t = 'min',
#'                    I = 0, Kic = Inf, Kiu = Inf, replicates = 3,
#'                    error = 'a', sd = 0.005, plot = TRUE)
#' @param So initial substrate concentration.
#' @param time reaction timespan.
#' @param Km Michaelis constant.
#' @param Vm maximal velocity.
#' @param unit_S concentration unit.
#' @param unit_t time unit.
#' @param I inhibitor concentration.
#' @param Kic competitive inhibition constant.
#' @param Kiu uncompetitive inhibition constant.
#' @param replicates number of replicates for the dependent variable
#' @param error it should be one among c('absolute', 'relative').
#' @param sd standard deviation of the error.
#' @param plot logical. If TRUE, the progress curve is plotted.
#' @details When sd is different to 0, then an absolute error normally distributed is added to the variable St.
#' @return Returns a dataframe where the two first columns are time and St (without error). The two last columns are the mean and sd of the variable St.
#' @examples sE.progress(So = 10, time = 5, Km = 4, Vm = 50, plot = FALSE)
#' @seealso fE.progress()
#' @importFrom VGAM lambertW
#' @importFrom stats rnorm
#' @export

sE.progress <- function(So, time, Km, Vm,
                       unit_S = 'mM', unit_t = 'min',
                       I = 0, Kic = Inf, Kiu = Inf,
                       replicates = 3,
                       error = 'a',
                       sd = 0.005,
                       plot = TRUE){

  ## -------------- Km and Vm apparent when I is present -------------- ##
  Km_a <- Km*( (1 + I/Kic)/(1 + I/Kiu) )
  Vm_a <- Vm/(1 + I/Kic)
  time <- seq(from = 0, to = time, by = (time/100))

  ## ---------------- Formating the output dataframe ------------------ ##
  output <- as.data.frame(matrix(rep(NA, length(time)*(replicates + 2)),
                                 ncol = replicates + 2))
  names(output) <- c('t', 'St', LETTERS[1:replicates])
  output$t <- time

  ## -------- Computing the variable substrate as function of t ------- ##
  set.seed(123)
  counter <- 0
  for (t in time){

    counter <- counter + 1
    argument <- (So/Km_a)*exp((-Vm_a*t + So)/Km_a)
    w <- VGAM::lambertW(argument)
    St <- Km_a*w # Substrate at time t
    output$St[counter] <-  St

    for (j in 1:replicates){
      if (error == 'r' | error == 'relative'){
        Se <- St * rnorm(1, mean = 1, sd = sd)
        output[counter, j+2] <- Se
      } else if (error == 'a' | error == 'absolute'){
        Se <- St + rnorm(1, mean = 0, sd = sd)
        output[counter, j+2] <- Se
      }
    }
  }

  ## -------------- Stop when St drops below a threshold -------------- ##
  output[1, -1] <- So
  if (So < 2*sd) {stop ('So lower than twice the SD')}
  output <- output[output$St > 2*sd, ]
  output[output < 0] <- 0

  ## --------------- Computing mean and sd if required ---------------- ##
  if (ncol(output) > 3){
    Substrate <- output[,-c(1,2)]
    output$S_mean <- apply(Substrate, MARGIN = 1, mean)
    output$S_sd <- apply(Substrate, MARGIN = 1, sd)
  } else if (ncol(output) == 3){
    output$S_mean <- output$A
    output$S_sd <- 0
  } else {
    output$S_mean <- output$St
    output$S_sd <- 0
  }

  ## -------- Plotting the results ------------- ##
  if (plot){
    plot(output$t, output$S_mean, ty = 'l', col = 'blue',
         xlab = paste("Time(", unit_t, ")", sep = ""), ylab = paste('[S]', unit_S))
    # arrows(output$t, output$S_mean-output$S_sd,
    #        output$t, output$S_mean-output$S_sd, length=0.05, angle=90, code=3)
  }

  return(output)
}

## ------------------------------------------------------------------------------- ##
#                       fE.progress(data)                                            #
## ------------------------------------------------------------------------------- ##
#' Fitted Progress Curve for Enzyme-Catalyzed Reaction
#' @description Fits the progress curve of an enzyme-catalyzed reaction.
#' @usage fE.progress(data, unit_S = 'mM', unit_t = 'min')
#' @param data a dataframe where the first column is the time and the second column is the substrate concentration.
#' @param unit_S concentration unit.
#' @param unit_t time unit.
#' @return Returns a list with two elements. The first one contains the fitted kinetic parameters, the second one is a dataframe giving the fitted substrate concentration time course.
#' @examples data <- sE.progress(So = 10, time = 5, Km = 4, Vm = 50, plot = FALSE)
#' @examples fE.progress(data[, c(1,3)])
#' @references Biochem Mol Biol Educ.39:117-25 (10.1002/bmb.20479).
#' @seealso sEprogress(), int.MM()
#' @importFrom VGAM lambertW
#' @importFrom stats nls
#' @importFrom graphics points
#' @export

fE.progress <- function(data, unit_S = 'mM', unit_t = 'min'){

  names(data) <- c('t', 'St')
  So <- data$St[1]

  ## ----------------------- Estimating the seed ------------------------ ##
  t <- int.MM(data)
  seed = list(Km = unname(t$parameters[1]), Vm = unname(t$parameters[2]))


  ## ------------------------ Fitting the curve ---.--------------------- ##
  model <- nls(St ~ (Km * VGAM::lambertW((So/Km)*exp((-Vm*t + So)/Km))),
               data = data, start = seed, trace = TRUE)

  Km <- round(summary(model)$coefficient[1,1], 3)
  sd_Km <- round(summary(model)$coefficient[1,2], 3)

  Vm <- round(summary(model)$coefficient[2,1], 3)
  sd_Vm <- summary(model)$coefficient[2,2]

  ## ------------------------ Fitted St values ------------------------- ##
  argument <- (So/Km)*exp((-Vm*data$t + So)/Km)
  w <- VGAM::lambertW(argument)
  fitted_St <- Km*w # Substrate at time t according to the fitted curve
  data$fitted_St <- fitted_St

  ## --------------------------- Plotting data ------------------------- ##
  parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")

  plot(data$t, data$St, ty = 'p', col = 'red', pch = 20,
       xlab = paste("time (", unit_t, ")", sep = ""),
       ylab = paste("[S] (", unit_S, ")", sep = ""),
       main = parameters)

  points(data$t, data$fitted_St, ty = 'l', col = 'blue')

  ## ------------------------------- Output ---------------------------- ##
  KmVm <- c(Km, Vm)
  names(KmVm) <- c('Km', 'Vm')

  output = list(KmVm, data)
  names(output) <- c('parameters', 'data')

  return(output)

}




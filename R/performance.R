## --------------------------------------------- ##
#                                                 #
#             performance.vo                      #
#             performance.intMM                   #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#    performance.vo(Km, Vm, error, sd)    #
## ------------------------------------------------------------------------------- ##
#' Performance of Methods Based on Initial Velocities
#' @description Evaluates the performance of the methods to estimate Km and Vm based on initial velocities.
#' @usage performance.vo(Km, Vm, unit_S = 'mM', unit_v = 'au', error, sd)
#' @param Km Michaelis contant.
#' @param Vm maximal velocity.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param error it should be one among c('absolute', 'relative').
#' @param sd standard deviation of the error.
#' @details
#' @return Returns a dataframe containing the fitted parameters under different conditions.
#' @author Juan Carlos Aledo
#' @examples
#' @seealso sE.progress(), int.MM()
#' @export

performance.vo <- function(Km, Vm, unit_S = 'mM', unit_v = 'au', error = 'a', sd){

  if (error == 'a' | error == 'absolute'){
    data <- sMM(Km, Vm, error = 'a', sd = sd)[, c(1,3:5)]
  } else if (error == 'r' | error == 'relative'){

  } else {
    stop("A proper error type should be provided")
  }

  output <- data.frame(method = c('LB', 'LBw','HW', 'EH', 'ECB', 'Non-linear'),
                       error = rep(error, 6),
                       sd = rep(sd, 6),
                       Km = rep(NA, 6), Vm = rep(NA,6))

  lb <- lb(data,weighting = FALSE)$fitted_parameters
  lbw <- lb(data,weighting = TRUE)$fitted_parameters
  hw <- hw(data)$fitted_parameters
  eh <- eh(data)$fitted_parameters
  ecb <- ecb(data)$fitted_parameters

  ## -- Compute velocity mean for non-linear fitting
  ddata <- data.frame(S = data$S, v = rep(NA, nrow(data)))
  ddata$v <- apply(data[,-1], MARGIN = 1, mean)
  nl <- dir.MM(ddata)$parameters

  output$Km[1] <- lb[1]
  output$Km[2] <- lbw[1]
  output$Km[3] <- hw[1]
  output$Km[4] <- eh[1]
  output$Km[5] <- ecb[1]
  output$Km[6] <- nl[1]

  output$Vm[1] <- lb[2]
  output$Vm[2] <- lbw[2]
  output$Vm[3] <- hw[2]
  output$Vm[4] <- eh[2]
  output$Vm[5] <- ecb[2]
  output$Vm[6] <- nl[2]

  rp <- c(Km, Vm)
  names(rp) <- c('Km', 'Vm')
  attr(output, 'real_parameters') <- rp

  return(output)
}

## ------------------------------------------------------------------------------- ##
#       performance.intMM(Km, Vm, time, error, sd)    #
## ------------------------------------------------------------------------------- ##
#' Performance of the Linearized Integrated MM Fitting
#' @description Evaluates the performance of the int.MM() function when estimating the kinetic parameters.
#' @usage performance.intMM(Km, Vm, time, error, sd)
#' @param Km Michaelis contant.
#' @param Vm maximal velocity.
#' @param unit_S concentration unit.
#' @param unit_t time unit.
#' @param error it should be one among c('absolute', 'relative').
#' @param sd standard deviation of the error.
#' @details The time interval must be carefully chosen, otherwise artefacts are possible. When sd is different to 0, then an abosolute error normally distributed is added to the variable St.
#' @return Returns a dataframe containing the fitted parameters under different conditions (So/Km ratios).
#' @author Juan Carlos Aledo
#' @examples
#' @seealso sE.progress(), int.MM()
#' @export

performance.intMM <- function(Km, Vm, time, error, sd){

  output <- data.frame(So = c(rep(0.1*Km,3), rep(Km,3), rep(3*Km,3), rep(10*Km,3)),
                       So_Km = c(rep(0.1,3), rep(1,3), rep(3,3), rep(10,3)),
                       R2 = rep(NA, 12), Km = rep(NA, 12), Vm = rep(NA, 12),
                       n = rep(NA, 12))

  So <- c(0.1*Km, Km, 3*Km, 10*Km)

  ns <- Kms <- Vms <- r2 <- c()
  for (i in 1:length(So)){
    d <- sE.progress(So[i], time = time, Km = Km, Vm = Vm, error = error, sd = sd)
    a <- int.MM(data.frame(t = d$t, St = d$A))
    b <- int.MM(data.frame(t = d$t, St = d$B))
    c <- int.MM(data.frame(t = d$t, St = d$C))
    ns <- c(ns, rep(nrow(d),3))
    r2 <- c(r2, attributes(a)[4], attributes(b)[4], attributes(b)[4])
    Kms <- c(Kms, a$parameters[1], b$parameters[1], c$parameters[1])
    Vms <- c(Vms, a$parameters[2], b$parameters[2], c$parameters[2])
  }
  r2 <- round(as.numeric(r2), 3)
  output$n <- ns
  output$R2 <- r2
  output$Km <- Kms
  output$Vm <- Vms

  parameters <- c(Km, Vm)
  names(parameters) <- c('Km', 'Vm')
  attr(output, 'real.parameters') <- parameters

  return(output)
}

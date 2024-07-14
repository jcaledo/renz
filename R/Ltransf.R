## --------------------------------------------- ##
#                                                 #
#        lb  (Lineweaver-Burk)                    #
#        hw (Hanes-Woolf)                         #
#        eh (Eadie-Hofstee)                       #
#        ecb (Eisenthal & Cornish-Bowden)         #
#                                                 #
## --------------------------------------------- ##

## ------------------------------------------------------------------------------- ##
#           lb(data, unit_S, unit_v, weighting = F, plot = T)                       #
## ------------------------------------------------------------------------------- ##
#' Lineweaver-Burk Transformation
#' @description Obtains Km and Vm using double reciprocal transformation
#' @usage lb(data, unit_S = 'mM', unit_v = 'au', weighting = FALSE, plot = TRUE)
#' @param data a dataframe where the first column is the independent variable, [S], and the remaining columns (as many as experiment replicates) correspond to the dependent variable, v.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param weighting logical. When TRUE the weight v^4 is employed.
#' @param plot logical. If TRUE the data and fitted line are plotted.
#' @return A double reciprocal plot and the Km and Vm computed using averaged 1/v (when more than one replicate is provided). In addition, this function returns a list of five elements. The first and second ones are vectors with the Km and Vm, respectively, computed individually for each replicate. The third one provides the R-squared values of the fits. The fourth element of the list gives the fitted Km and Vm. The last element of the list is a dataframe with the values of the transformed variables.
#' @examples lb(ONPG[, c(1,2)], weighting = TRUE)
#' @references J. Am. Chem. Soc.1934, 56, 3,658-666 (doi.org/10.1021/ja01318a036)
#' @seealso hw(), eh(), ecb()
#' @importFrom stats complete.cases
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom graphics arrows
#' @export

lb <- function(data, unit_S = 'mM', unit_v = 'au', weighting = FALSE, plot = TRUE){

  ## ---------------------------- Removing incomplete data ----------------------- ##
  data <- data[complete.cases(data), ]
  colnames(data)[1] <- 'S'

  ## -------------------------- Taking double reciprocal ----------------------------- ##
  tdata <- apply(data, MARGIN = 2, function(x)  1/x)
  tdata <- as.data.frame(round(tdata, 4))
  names(tdata)[1] <- 'inv_S'


  ## ----------------------- Replicates model fitting ----------------------------- ##
  Kms <- Vms <- R2 <- c()
  for (j in 2:ncol(tdata)){
    model <- lm(tdata[,j] ~ tdata$inv_S )
    Vm <- round(1/model$coefficients[1], 2)
    Km <- round(Vm * model$coefficients[2], 2)
    Vms <- c(Vms, Vm)
    Kms <- c(Kms, Km)
    R2 <- c(R2,  summary(model)$r.squared)
  }

  ## -------- Computing mean and sd if required of the transformed data ---------- ##
  if (ncol(tdata) > 2){
    tdata$inv_v <- round(apply(tdata[,-1], MARGIN = 1, mean), 4)
    tdata$sd <- round(apply(tdata[,-c(1, ncol(tdata))], MARGIN = 1, sd), 4)
  } else {
    names(tdata)[2] <- 'inv_v'
  }
  ## -------------------------- Mean model fitting -------------------------------- ##
  if (weighting){
    model <- lm(tdata$inv_v ~ tdata$inv_S, weights = 1/(tdata$inv_v)^4)
  } else {
    model <- lm(tdata$inv_v ~ tdata$inv_S)
  }

  Vm <- round(1/model$coefficients[1], 2)
  Km <- round(Vm * model$coefficients[2], 2)

  ## --------------------- Plotting the transformed data ------------------------- ##
  if (plot == TRUE){
    parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
    plot(tdata$inv_S, tdata$inv_v, ty = 'p',
         ylim = c(0, max(tdata[,-1]) + 0.1*max(tdata[,-1])),
         xlab = paste('1/[S] (1/', unit_S, ')', sep = ""),
         ylab = paste('1/v (1/', unit_v, ')', sep = ""), main = parameters)
    abline(model)

    if (ncol(tdata) > 2){
      arrows(tdata$inv_S, tdata$inv_v - tdata$sd,
             tdata$inv_S, tdata$inv_v + tdata$sd,
             length=0.05, angle=90, code=3)
    }
  }


  fitted_parameters <- c(Km, Vm)
  names(fitted_parameters) <- c('Km', 'Vm')
  output <- list(unname(Kms), unname(Vms), R2, fitted_parameters, tdata)
  names(output) <- c('Kms', 'Vms', 'R2s', 'fitted_parameters', 'inverse_data')
  return(output)
}

## ------------------------------------------------------------------------------- ##
#                hw(data, unit_S = 'mM', unit_v = 'au', plot = T)                   #                   #
## ------------------------------------------------------------------------------- ##
#' Hanes-Woolf Transformation
#' @description Obtains Km and Vm using the Hanes-Woolf transformation.
#' @usage hw(data, unit_S = 'mM', unit_v = 'au', plot = TRUE)
#' @param data a dataframe where the first column is the independent variable, [S], and the remaining columns (as many as experiment replicates) correspond to the dependent variable, v.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param plot logical. If TRUE the data and fitted line are plotted.
#' @return A dataframe with the values of the transformed variables is returned. The fitted Km and Vm are given as attributes of this dataframe.
#' @examples hw(ONPG[, c(1,2)])
#' @seealso lb(), eh(), ecb()
#' @importFrom stats complete.cases
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @export

hw <- function(data, unit_S = 'mM', unit_v = 'au', plot = TRUE){

  ## ---------------------------- Removing incomplete data ----------------------- ##
  data <- data[complete.cases(data), ]
  colnames(data)[1] <- 'S'

  ## ---------------------- Transforming original data  -------------------------- ##
  tdata <- apply(data, MARGIN = 2, function(x)  1/x)
  tdata <- as.data.frame(round(tdata, 4))
  tdata[,1] <- 1/tdata[,1]
  if (ncol(tdata) > 2){
    tdata[,2:ncol(tdata)] <- apply(tdata[,-1], MARGIN = 2, function(x)  x*tdata[,1])
  } else {
    tdata[,2] <- tdata[,1]*tdata[,2]
  }

  ## ----------------------- Replicates model fitting ----------------------------- ##
  Kms <- Vms <- R2 <- c()
  for (j in 2:ncol(tdata)){
    model <- lm(tdata[,j] ~ tdata$S )
    Vm <- round(1/model$coefficients[2], 2)
    Km <- round(Vm * model$coefficients[1], 2)
    Vms <- c(Vms, Vm)
    Kms <- c(Kms, Km)
    R2 <- c(R2,  summary(model)$r.squared)
  }

  ## -------- Computing mean and sd if required of the transformed data ---------- ##
  if (ncol(tdata) > 2){
    tdata$S_v <- round(apply(tdata[,-1], MARGIN = 1, mean), 4)
    tdata$sd <- round(apply(tdata[,-c(1, ncol(tdata))], MARGIN = 1, sd), 4)
  } else {
    names(tdata)[2] <- 'S_v'
  }

  ## ---------------------------- Model fitting ---------------------------------- ##
  model <- lm(tdata$S_v ~ tdata$S)
  Vm <- round(1/model$coefficients[2], 2)
  Km <- round(Vm*model$coefficients[1], 2)

  ## --------------------- Plotting the transformed data ------------------------- ##
  parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
  if (plot){
    plot(tdata$S, tdata$S_v, ty = 'p',
         ylim = c(0, max(tdata[,-1]) + 0.1*max(tdata[,-1])),
         xlab = paste('[S] (', unit_S, ')', sep = ""),
         ylab = paste('[S]/v (', unit_S, '/', unit_v, ')', sep = ""), main = parameters)
    abline(model)

    if (ncol(tdata) > 2){
      arrows(tdata$S, tdata$S_v - tdata$sd,
             tdata$S, tdata$S_v + tdata$sd,
             length=0.05, angle=90, code=3)
    }
  }

  fitted_parameters <- c(Km, Vm)
  names(fitted_parameters) <- c('Km', 'Vm')
  output <- list(unname(Kms), unname(Vms), R2, fitted_parameters, tdata)
  names(output) <- c('Kms', 'Vms', 'R2s', 'fitted_parameters','transformed_data')

  return(output)
}

## ------------------------------------------------------------------------------- ##
#                  eh(data, unit_S = 'mM', unit_v = 'au', plot = T)                 #                                       #
## ------------------------------------------------------------------------------- ##
#' Eadie-Hofstee Transformation
#' @description Obtain Km and Vm using the Eadie-Hofstee transformation.
#' @usage eh(data, unit_S = 'mM', unit_v = 'au', plot = TRUE)
#' @param data a dataframe where the first column is the independent variable, [S], and the remaining columns (as many as experiment replicates) correspond to the dependent variable, v.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param plot logical. If TRUE the data and fitted line are plotted.
#' @return A dataframe with the values of the transformed variables is returned. The fitted Km and Vm are given as attributes of this dataframe.
#' @examples eh(ONPG[, c(1,2)])
#' @seealso lb(), hw(), ecb()
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom graphics abline
#' @importFrom graphics arrows
#' @export

eh <- function(data, unit_S = 'mM', unit_v = 'au', plot = TRUE){

  ## ---------------------------- Removing incomplete data ----------------------- ##
  data <- data[complete.cases(data), ]
  colnames(data)[1] <- 'S'

  replicates <- ncol(data) - 1
  ## -------------------- Transformed independent variable ----------------------- ##
  if (replicates > 1){
    tdata <- apply(data[,-1], MARGIN = 2, function(x)  x/data[,1])
    tdata <- as.data.frame(tdata)
    tdata$v_S <- apply(tdata, MARGIN = 1, mean)
    tdata$x_sd <- apply(tdata, MARGIN = 1, sd)
  } else {
    tdata <- data.frame(v_S = data[,2]/data[,1])
  }
  x_max <- max(tdata)

  ## --------------------- Transformed dependent variable ------------------------ ##
  if (replicates > 1){
    tdata[, (ncol(tdata) + 1):(ncol(tdata) + ncol(data) - 1)] <- data[,-1]
    tdata$v <- apply(data[,-1], MARGIN = 1, mean)
    tdata$y_sd <- apply(data[,-1], MARGIN = 1, sd)
  } else {
    tdata$v <- data[,2]
  }
  y_max <- max(data)

  ## ----------------------- Replicates model fitting ----------------------------- ##
  Kms <- Vms <- R2s <- c()
  if (replicates > 1){
    for (j in 1:replicates){
      model <- lm(tdata[,(replicates + 2 + j)] ~ tdata[,j] )
      Km <- round(-model$coefficients[2], 2)
      Vm <- round(model$coefficients[1], 2)
      Vms <- c(Vms, Vm)
      Kms <- c(Kms, Km)
      R2s<- c(R2s,  summary(model)$r.squared)
    }
  }

  ## ---------------------------- Model fitting ---------------------------------- ##
  model <- lm(tdata$v ~ tdata$v_S)
  Km <- round(-model$coefficients[2], 2)
  Vm <- round(model$coefficients[1], 2)
  if (replicates == 1){
    Vms <- Vm
    Kms <- Km
    R2s <- summary(model)$r.squared
  }

  ## --------------------- Plotting the transformed data ------------------------- ##
  parameters <- paste('Km: ', Km, '     Vm: ', Vm, sep = "")
  if (plot){
    plot(tdata$v_S, tdata$v, ty = 'p',
         xlim = c(0, x_max + 0.1 * x_max),
         ylim = c(0, y_max + 0.1 * y_max),
         xlab = paste('v/[S] (', unit_v, '/', unit_S, ')', sep = ""),
         ylab = paste('v (',  unit_v, ')', sep = ""), main = parameters)
    abline(model)

    if (ncol(tdata) > 2){
      arrows(tdata$v_S, tdata$v - tdata$y_sd,
             tdata$v_S, tdata$v + tdata$y_sd,
             length=0.05, angle=90, code=3)

      arrows(tdata$v_S - tdata$x_sd, tdata$v,
             tdata$v_S + tdata$x_sd, tdata$v,
             length=0.05, angle=90, code=3)
    }
  }

  fitted_parameters <- c(Km, Vm)
  names(fitted_parameters) <- c('Km', 'Vm')
  output <- list(unname(Kms), unname(Vms), R2s, fitted_parameters,tdata)
  names(output) <- c('Kms', 'Vms', 'R2s', 'fitted_parameters','transformed_data')

  return(output)
}

## ------------------------------------------------------------------------------- ##
#                     ecb(data, unit_S = 'mM', unit_v = 'au', plot = T)             #                                   #
## ------------------------------------------------------------------------------- ##
#' Eisenthal & Cornish-Bowden
#' @description Obtains Km and Vm using the Eisenthal & Cornish-Bowden method.
#' @usage ecb(data, unit_S = 'mM', unit_v = 'au', plot = TRUE)
#' @param data a dataframe where the first column is the independent variable, [S], and the remaining columns (as many as experiment replicates) correspond to the dependent variable, v.
#' @param unit_S concentration unit.
#' @param unit_v time unit.
#' @param plot logical. If TRUE data are plotted.
#' @details For each experimental replicate the observations (S, v) are plotted as lines in the Km-Vm parameter space, instead of points in observation space. Afterwards, the lines tend to intersect at a common point, whose coordinates provide the kinetic parameters. Nevertheless, since the observations are subject to error, there is no unique intersection point for all the lines. In this case, the method computes all the pair-wise intersections. Then, the median value from each series is taken to be the best estimate of Km and Vm. This procedure is repeated as many times as replicates and finally the mean and sd is returned.
#' @return Returns a list with the estimated values of Km and Vm.
#' @examples ecb(ONPG[, c(1,2)])
#' @references Biochem.J.(1974) 139:715-720 (10.1042/bj1390715)
#' @seealso lb(), hw(), eh()
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats complete.cases
#' @importFrom graphics points
#' @importFrom graphics abline
#' @export

ecb <- function(data, unit_S = 'mM', unit_v = 'au', plot = TRUE){

  ## -------------------- Variables' names ---------------------- ##
  colnames(data)[1] <- 'S'
  colnames(data)[2:length(colnames(data))] <- LETTERS[1:length(colnames(data)) - 1]
  data <- data[complete.cases(data), ]

  ## -------------------- Computing v/[S] ----------------------- ##
  if (ncol(data) > 2){
    tdata <- apply(data[,-1], MARGIN = 2, function(x)  x/data[,1])
    mean_v <- apply(data[,-1], MARGIN = 1, mean)
    mean_v_S <- apply(tdata, MARGIN = 1, mean)
  } else {
    tdata <- data[,2]/data[,1]
    mean_v <- data[,2]
    mean_v_S <- tdata
  }

  tdata <- as.data.frame(tdata)
  names(tdata)[1] <- 'v_S'
  tdata <- cbind(data, tdata)
  names(tdata)[2] <- 'v'

  nc <- ncol(data) + 1 # number of the first colum for v/S data


  ## ------------- Pair-wise intersection and median ----------- ##
  line.intersect <- function(v, v_S){
    Kms <- Vms <- c()
    for (i in 1:(length(v)-1)){
      for (j in (i+1):length(v)){
        Km <- (v[j] - v[i])/(v_S[i] - v_S[j])
        Kms <- c(Kms, Km)
        Vm <- (v_S[i]*(v[j] - v[i])/(v_S[i] - v_S[j])) + v[j]
        Vms <- c(Vms, Vm)
      }
    }
    Km <- median(Kms, na.rm = TRUE)
    Vm <- median(Vms, na.rm = TRUE)
    return(c(Km, Vm))
  }
  ## ------------------------------------------------------------ ##


  ## --------------------- Estimating Km and Vm ----------------- ##
  Kms <- Vms <- c()
  c <- 0
  for (i in 2:(nc - 1)){
    p <- line.intersect(tdata[,i], tdata[,nc + c])
    Kms <- c(p[1], Kms)
    Vms <- c(p[2], Vms)
    c <- c + 1
  }
  Km <- paste('Km: ', round(mean(Kms), 3), ' \u00b1 ', round(sd(Kms), 3), ' ', unit_S)
  Vm <- paste('Vm: ', round(mean(Vms), 3), ' \u00b1 ', round(sd(Vms), 3), ' ', unit_v)
  fitted_parameters <- c(round(mean(Kms, na.rm = TRUE), 3),  round(mean(Vms, na.rm = TRUE), 3))
  names(fitted_parameters) <- c('Km', 'Vm')

  ## --------------------- Plotting results --------------------- ##
  if (plot){
    plot(0, 0, ty = 'n',
         xlim = c(-tdata$S[nrow(tdata)], tdata$S[nrow(tdata)]),
         ylim = c(0, 3*max(tdata$v, na.rm = TRUE)),
         xlab = paste("Km (", unit_S, ")", sep = ""),
         ylab = paste("Vm (", unit_v, ")", sep = ""))

    for (i in 1:nrow(tdata)){
      points(c(0, -tdata$S[i]), c(mean_v[i], 0), pch = 20)
      abline(mean_v[i], mean_v_S[i])
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
  }

  output <- list(fitted_parameters, Km, Vm)
  names(output) <- c('fitted_parameters', 'Km', 'Vm')
  return(output)
}



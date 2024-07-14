## --------------------------------------------- ##
#                                                 #
#                    bibi                         #
#                                                 #
## --------------------------------------------- ##


## ------------------------------------------------------------------------------- ##
#                            bibi(data)                                             #
## ------------------------------------------------------------------------------- ##
#' Kinetic Mechanisms and Parameters for Bi-Bi Reactions
#' @description Discriminates between sequential and ping-pong mechanisms and estimates the kinetic parameters.
#' @usage bibi(data, unit_a = "mM", unit_b = "mM", unit_v = "ua", vice_versa = FALSE)
#' @param data either a dataframe or the path to a text file containing the data (see details).
#' @param unit_a concentration unit for substrate A.
#' @param unit_b concentration unit for substrate B.
#' @param unit_v velocity unit.
#' @param vice_versa logical. When FALSE the variable substrate is A. If TRUE, then the variable substrate is B.
#' @details Either the txt file or the dataframe containing the data must conform to the following format: a table with three columns and as many rows as conditions were assessed. The first and second columns are named 'a' and 'b' and they give the concentrations for substrate A and B, respectively. The third column, named 'rate', provides the assessed rates.
#' @return A list with three elements: (i) a character vector giving the kinetic parameters Vmax, KiA, Km_A and Km_B values; (ii) a numeric vector giving the apparent inverse of Vmax for each concentration of substrate B (intercepts of primary representation); and (iii) a numeric vector giving the apparent specificity constant for each concentration of substrate B (slopes from primary representations).
#' @examples bibi(data = hk)
#' @importFrom utils read.table
#' @export

bibi <- function(data, unit_a = "mM", unit_b = "mM", unit_v = "ua", vice_versa = FALSE){

   if (!is.data.frame(data)){
     data <- read.table(data, header = TRUE)
     if (gregexpr("M", data$a[1])[[1]] != -1){
       data <- data[-1, ]
     }
   }
   ## Make sure that variables are numeric ones
   data$a <- as.numeric(as.character(data$a))
   data$b <- as.numeric(as.character(data$b))
   data$rate <- as.numeric(as.character(data$rate))

   if (vice_versa){
   ## ----------- Vice versa -------------------------------##
   # We make the substrate 'a' to be 'b' and  vice versa
   data.viceversa <- data[,c(2,1,3)]
   data.viceversa <- data.viceversa[order(data$b),]
   colnames(data.viceversa) <- c("a", "b", "rate")
   data <- data.viceversa
   ## ----------------------------------------------------##
   }

   ## Substrate A as variable and substrate B as parameter
   data$inv_a <- 1/data$a
   data$inv_rate <- 1/data$rate

   b <- unique(data$b) # Different concentrations of substrate B
   slopes <- c()
   intercepts <- c()
   r2 <- c() # r.squared

     for (i in 1:length(b)){
       a <- c() # Concentrations of substrate A
       v <- c() # Initial rates

       for (j in 1:nrow(data)){
         if (data$b[j] == b[i]){
           a <- c(a, data$a[j])
           v <- c(v, data$rate[j])
         }
       }

       inv_a <- 1/a
       inv_v <- 1/v
       fit <- lm(inv_v~inv_a)
       slopes <- c(slopes, fit[[1]][2])
       intercepts <- c(intercepts, fit[[1]][1])
       r2 <- c(r2, summary(fit)[8])
     }

     ## Primary representation
     plot(data$inv_a,
          data$inv_rate, xlim = c(-min(data$inv_a), max(data$inv_a)),
          ylim = c(0, max(data$inv_rate)),
          xlab = "1/a", ylab = "1/v")
     abline(h = 0, col = 'red')
     abline(v = 0, col = 'red')
     for (i in 1:length(slopes)){
       abline(intercepts[i], slopes[i])
     }

    ## Vmax and Km_B
    inv_b <- 1/b
    sec1 <- lm(intercepts~inv_b)

    Vm <- round(1/sec1[[1]][1], 2) # Numeric value
    Vmax <- paste(Vm, unit_v) # Result with units

    Km.b <- round(sec1[[1]][2] * Vm, 2)
    Km_B <- paste(Km.b, unit_b)

    ## Km_A and KiA
    sec2 <- lm(slopes~inv_b)

    Km.a <- round(sec2[[1]][1] * Vm, 2)
    Km_A <- paste(Km.a, unit_a)

    KiA <- round(sec2[[1]][2]*Vm/Km.b, 2)

    kinetic_parameters <- c(Vmax, KiA, Km_A, Km_B)
    names(kinetic_parameters) <- c("Vmax", "KiA", "Km_A", "Km_B")

    output = list(kinetic_parameters, intercepts, slopes)
    names(output) <- c("kinetic_parameters", "inv_Vmax_app", "specificity_const_app")
    return(output)
}

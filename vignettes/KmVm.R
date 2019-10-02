## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = FALSE--------------------------------------------------------
library(renz)
library(knitr)
E <- rMM(plot = FALSE)
Km <- attributes(E)$Km
Vm <- attributes(E)$Vm
S = c(0.5, 1, 2, 4, 7, 10, 13, 20)
v <- round(Vm * S / (Km + S), 2)
data <- data.frame(S,v)
names(data) <- c('[S] (mM)', 'v (au)')
kable(data)


## ---- fig.align='center',out.extra='angle=90', echo=FALSE----------------
plot(S, v, ty = 'b', ylim = c(0, Vm + 0.25*Vm), 
     ylab = 'v (au)', xlab = '[S] (mM)')
abline(h = 70, lty = 2, col = 'blue')
abline(v = (35*Km)/(Vm - 35), lty = 2, col = 'blue')
abline(h = Vm, lty = 2, col = 'red')
abline(v = Km, lty = 2, col = 'red')


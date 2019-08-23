#MaferBravoG
#Ubicación

#Librerias
library(deSolve)
library(scatterplot3d)
## Modelo SEIR sin demografía
SEIRsd <- function(t, estado, parametros) {
  with(as.list(c(estado, parametros)), {
    dS <-  -beta*I*S/sum(estado)
    dE <-  beta*I*S/sum(estado) - theta*E
    dI <-  theta*E-gama*I
    dR <-  gama*I
    list(c(dS, dE, dI, dR))
  })
}
#beta: S -> E
#theta: E -> I
#gama: I -> R
parametros <- c(beta =25, gama = 10, theta = 40)
##Número de suceptibles, latencia, infectados y recuperadosen la población
estado <- c(S = 50, E=25, I = 15, R = 10 )
#Periodo de tiempo (del tiempo 0 al 10 y calcular las derivadas cada 0.01 tiempo)
tiempo <- seq(0, 10, by = 0.01)
#Objeto a graficar
fSEIRsd <- ode(y = estado, times = tiempo, func = SEIRsd, parms = parametros)
#Archivo pdf para la gráfica
pdf("Grafica_SEIR_sd.pdf")
#Gráfica SEIR con demografía
plot(fSEIRsd,col="orange2")
dev.off()
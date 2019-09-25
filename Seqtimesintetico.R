#María Fernanda Bravo García
#Manual Seq.time

#Librerias 
library(seqtime)
library(ggplot2)
library(reshape2)
library(RCy3)
cytoscapePing()

#To generate a simulated dataset with 50 species, 40 samples. 

#1) Interaction matrix: Klemm-Eguiluz algorithm, as this produces modular and scale-free networks.
N = 50 
S = 40
#pep: porcentanje de edges positivos  c:connectance 
A = generateA(N,"klemm", pep = 10, c = 0.05)
save(A, file = "Matriz.RData")
#Adjusting connectance to 0.05"
#Initial edge number 520"
#Initial connectance 0.191836734693878"
#Number of edges removed 348"
#Final connectance 0.0497959183673469"
#Final connectance: 0.0497959183673469"
#Initial edge number 172"
#Initial connectance 0.0497959183673469"
#Number of negative edges already present: 50"
#Converting 105 edges into negative edges"
#Final connectance: 0.0497959183673469"
#Final arc number (excluding self-arcs) 122"
#Final negative arc number (excluding self-arcs) 105"
#PEP: 13.9344262295082"

rownames(A)=c(1:N)
colnames(A)=rownames(A)
matrizKlemm=plotA(A, header = "Klemm-Eguiluz interaction matrix")
save(matrizKlemm, file = "MatrizKlemm.RData")

# "Largest value: 0.487130493974374"
# "Smallest value: -0.5"

#2) Abundancias iniciales con gLV (diferentes a las generadas por una distribución de Poisson)
dataset = generateDataSet(S, A)
dataset = seqtime::normalize(dataset)
dataset = melt(dataset)
colnames(dataset) = c("Species", "Sample", "Abundance")
ggplot(data=dataset, aes(x=dataset$Sample, y=dataset$Abundance, width=1)) + geom_bar(aes(y = dataset$Abundance, x= dataset$Sample, fill=dataset$Species), data=dataset, stat="identity", show.legend=F) + theme(aspect.ratio=.4) + theme_classic()+ ylab("Relative abundance") + xlab("Sample")
save(dataset, file = "abundanciasgLV.RData")
#3) Para la red a partir de la matriz de interacción 
plotA(A, method = "network")

#[1] "Largest value: 0.466089489392224"
#[1] "Smallest value: -0.5"
#[1] "Initial edge number 172"
#[1] "Initial connectance 0.0497959183673469"
#[1] "Final connectance: 0.0518707482993197"

#4)Simulate a test time series with:
#SIN PERTURBACIÓN 
##the Ricker community model
out.rickers=ricker(N,A=A, tend=250)
tsplot(out.rickers,main="Ricker sin ")
save(out.rickers, file= "TimeRickerSin.RData")
## the gLV
out.glvs = glv (N, A=A, tend=25)
tsplot(out.glvs, main="gLV sin")
save(out.glvs, file = "TimegLVSin.RData" )
## Hubbell
out.hubbells=simHubbell(N=N, M=N, I=1500, d=N, m=0.1, tskip=250, tend=500)
tsplot(out.hubbells, main="Hubbell sin")
save(out.hubbells, file = "TimeHubbellSin.RData" )
## SOI (I= individuos? )
out.sois= soi(N, A=A, I=1500, tend=250)
tsplot(out.sois, main= "SOI sin")
save(out.sois, file = "TimeSOISin.RData" )

#CON PERTURBACIÓN
gc=runif(50, min =-1, max=1)
pert=perturbation(times = c(100), durations = (50), growthchanges = gc)
##the Ricker community model
out.rickerc=ricker(N,A=A, tend=250, perturb = pert)
tsplot(out.rickerc,main="Ricker con ")
save(out.rickerc, file = "TimeRickerCon.RData")
## the gLV ( Please provide as many growth changes as species!)
out.glvc = glv (N, A=A, tend=250, perturb = pert)
tsplot(out.glvc, main="gLV con")
save(out.glvc, file = "TimegLVCon.RData")
## Hubbell
out.hubbellc=simHubbell(N=N, M=N, I=1500, d=N, m=0.1, tskip=250, tend=500, perturb = pert)
tsplot(out.hubbellc, main="Hubbell con")
save(out.hubbellc, file = "TimeHubbellCon.RData")
## SOI (I= individuos? )
out.soic= soi(N, A=A, I=1500, tend=250, perturb = pert)
tsplot(out.soic, main= "SOI con")
save(out.soic, file = "TimeSOICon.RData")
#5) LIMITS to test how well it can infer the known interaction matrix from the time series
mfrow=c(1,1)
rs=limits(out.rickers)$Aest
save(rs, file = "RedRickerSin.RData")

rc=limits(out.rickerc)$Aest
save(rc, file="RedRickerCon.RData")

gs=limits(out.glvs)$Aest
save(gs, file= "RedgLVSin.RData")

gc=limits(out.glvc)$Aest
save(gc, file="RedgLVCon.RData")

hs=limits(out.hubbells)$Aest
save(hs, file="RedHubbellSin.RData")

hc=limits(out.hubbellc)$Aest
save(hc, file="RedHubbellCon.RData")

ss=limits(out.sois)$Aest
save(ss, file="RedSOISin.RData")

sc=limits(out.soic)$Aest
save(sc, file="RedSOICon.RData")

#Para la comparativa entre matrices
par(mfrow=c(1,3))
plotA(rs,header="RickerSin")
plotA(A,header="Conocida")
plotA(rc,header="RickerCon")

plotA(gs,header="gLVSin")
plotA(A,header="Conocida")
plotA(gc,header="gLVCon")

plotA(hs,header="HubbellSin")
plotA(A,header="Conocida")
plotA(hc,header="HubbellCon")

plotA(ss,header="SOISin")
plotA(A,header="Conocida")
plotA(sc,header="SOICon")
par(mfrow=c(1,1))

#REDES 
na=graph_from_adjacency_matrix(A, mode="directed", weighted=TRUE)
plot(na, main="Conocida")

nrs=graph_from_adjacency_matrix(rs, mode="directed", weighted=TRUE)
plot(nrs, main="Ricker Sin")

nrc=graph_from_adjacency_matrix(rc, mode="directed", weighted=TRUE)
plot(nrc, main="Ricker Con")

ngs=graph_from_adjacency_matrix(gs, mode="directed", weighted=TRUE)
plot(ngs, main="gLV Sin")

ngc=graph_from_adjacency_matrix(gc, mode="directed", weighted=TRUE)
plot(ngc, main="gLV Con")

nhs=graph_from_adjacency_matrix(hs, mode="directed", weighted=TRUE)
plot(nhs, main="Hubbell Sin")

nhc=graph_from_adjacency_matrix(hc, mode="directed", weighted=TRUE)
plot(nhc, main="Hubbell Con")

nss=graph_from_adjacency_matrix(ss, mode="directed", weighted=TRUE)
plot(nss, main="SOI Sin")

nsc=graph_from_adjacency_matrix(sc, mode="directed", weighted=TRUE)
plot(nsc, main="SOI Con")

#Para ver lo de diferentes pep (%edges positivos)
N = 50 
S = 40
#pep: porcentanje de edges positivos  c:connectance 
A0 = generateA(N,"klemm", pep = 0, c = 0.05)
save(A0, file = "MatrizA0.RData")
A25 = generateA(N,"klemm", pep = 25, c = 0.05)
save(A25, file = "MatrizA25.RData")
A50 = generateA(N,"klemm", pep = 50, c = 0.05)
save(A50, file = "MatrizA50.RData")
A75 = generateA(N,"klemm", pep = 75, c = 0.05)
save(A75, file = "MatrizA75.RData")
#Warning message para 75
#In modifyA(A = A, perc = (100 - pep), symmetric = negedge.symm,  :
#The matrix has more negative edges than are required to reach the desired negative edge percentage!
A100 = generateA(N,"klemm", pep = 100, c = 0.05)
save(A100, file = "MatrizA100.RData")

#Simulate a test time series SIN PERTURBACIÓN CON LAS DIFERENTES MATRICES
##the Ricker community model
out.rickersA0=ricker(N,A=A0, tend=250)
tsplot(out.rickersA0,main="Ricker sin A0")
save(out.rickersA0, file= "TimeRickerSinA0.RData")

out.rickersA25=ricker(N,A=A25, tend=250)
tsplot(out.rickersA25,main="Ricker sin A25")
save(out.rickersA25, file= "TimeRickerSinA25.RData")
#EXPLOSION: 
#Error in seq.default(0, 1, 1/nrow(x)) : 'by' must be of length 1
out.rickersA50=ricker(N,A=A50, tend=250)
tsplot(out.rickersA50,main="Ricker sin A50 ")
save(out.rickersA50, file= "TimeRickerSinA50.RData")
#EXPLOSION: 
#Error in seq.default(0, 1, 1/nrow(x)) : 'by' must be of length 1
out.rickersA75=ricker(N,A=A75, tend=250)
tsplot(out.rickersA75,main="Ricker sin A75")
save(out.rickersA75, file= "TimeRickerSinA75.RData")
#EXPLOSION: 
#Error in seq.default(0, 1, 1/nrow(x)) : 'by' must be of length 1
out.rickersA100=ricker(N,A=A100, tend=250)
tsplot(out.rickersA100,main="Ricker sin A100")
save(out.rickersA100, file= "TimeRickerSinA100.RData")
## the gLV
out.glvsA0 = glv (N, A=A0, tend=25)
tsplot(out.glvsA0, main="gLV sin A0")
save(out.glvsA0, file = "TimegLVSinA0.RData" )

out.glvsA25 = glv (N, A=A25, tend=25)
tsplot(out.glvsA25, main="gLV sin A25")
save(out.glvsA25, file = "TimegLVSinA25.RData")

out.glvsA50 = glv (N, A=A50, tend=25)
tsplot(out.glvsA50, main="gLV sinA50")
save(out.glvsA50, file = "TimegLVSinA50.RData" )

out.glvsA75 = glv (N, A=A75, tend=25)
tsplot(out.glvsA75, main="gLV sin A75")
save(out.glvsA75, file = "TimegLVSinA75.RData" )

out.glvsA100 = glv (N, A=A100, tend=25)
tsplot(out.glvsA100, main="gLV sin A100")
save(out.glvsA100, file = "TimegLVSinA100.RData" )

## SOI (I= individuos? )
out.soisA0= soi(N, A=A0, I=1500, tend=250)
tsplot(out.soisA0, main= "SOI sin A0")
save(out.soisA0, file = "TimeSOISinA0.RData" )

out.soisA25= soi(N, A=A25, I=1500, tend=250)
tsplot(out.soisA25, main= "SOI sin A25")
save(out.soisA25, file = "TimeSOISinA25.RData" )

out.soisA50= soi(N, A=A50, I=1500, tend=250)
tsplot(out.soisA50, main= "SOI sin A50")
save(out.soisA50, file = "TimeSOISinA50.RData" )

#Error in params[[5]][[l]] : subscript out of bounds
out.soisA75= soi(N, A=A75, I=1500, tend=250)
tsplot(out.soisA75, main= "SOI sin A75")
save(out.soisA75, file = "TimeSOISinA75.RData" )
#Error in params[[5]][[l]] : subscript out of bounds
out.soisA100= soi(N, A=A100, I=1500, tend=250)
tsplot(out.soisA100, main= "SOI sin A100")
save(out.soisA100, file = "TimeSOISinA100.RData" )

#LIMITS de las que si crecieron
mfrow=c(1,1)
rsA0=limits(out.rickersA0)$Aest
save(rsA0, file = "RedRickerSinA0.RData")

rsA25=limits(out.rickersA25)$Aest
save(rsA25, file = "RedRickerSinA25.RData")

gsA0=limits(out.glvsA0)$Aest
save(gsA0, file= "RedgLVSinA0.RData")

gsA25=limits(out.glvsA25)$Aest
save(gsA25, file= "RedgLVSinA25.RData")

gsA50=limits(out.glvsA50)$Aest
save(gsA50, file= "RedgLVSinA50.RData")

gsA75=limits(out.glvsA75)$Aest
save(gsA75, file= "RedgLVSinA75.RData")

gsA100=limits(out.glvsA100)$Aest
save(gsA100, file= "RedgLVSinA100.RData")

ssA0=limits(out.soisA0)$Aest
save(ssA0, file="RedSOISinA0.RData")

ssA25=limits(out.soisA25)$Aest
save(ssA25, file="RedSOISinA25.RData")

ssA50=limits(out.soisA50)$Aest
save(ssA50, file="RedSOISinA50.RData")

#Redes de diferentes matrices
naA0=graph_from_adjacency_matrix(A0, mode="directed", weighted=TRUE)
createNetworkFromIgraph(naA0, "Conocida A0")

naA25=graph_from_adjacency_matrix(A25, mode="directed", weighted=TRUE)
createNetworkFromIgraph(naA25, "Conocida A25")

naA50=graph_from_adjacency_matrix(A50, mode="directed", weighted=TRUE)
createNetworkFromIgraph(naA50, "Conocida A50")

naA75=graph_from_adjacency_matrix(A75, mode="directed", weighted=TRUE)
createNetworkFromIgraph(naA75, "Conocida A75")

naA100=graph_from_adjacency_matrix(A100, mode="directed", weighted=TRUE)
createNetworkFromIgraph(naA100, "Conocida A100")

nrsA0=graph_from_adjacency_matrix(rsA0, mode="directed", weighted=TRUE)
createNetworkFromIgraph(nrsA0, "Ricker Sin A0")
#plot(nrsA0, main="Ricker Sin A0")

nrsA25=graph_from_adjacency_matrix(rsA25, mode="directed", weighted=TRUE)
createNetworkFromIgraph(nrsA25, "Ricker Sin A25")
#plot(nrsA25, main="Ricker Sin A25")

ngsA0=graph_from_adjacency_matrix(gsA0, mode="directed", weighted=TRUE)
createNetworkFromIgraph(ngsA0, "gLV Sin A0")
#plot(ngsA0, main="gLV Sin A0")

ngsA25=graph_from_adjacency_matrix(gsA25, mode="directed", weighted=TRUE)
createNetworkFromIgraph(ngsA25, "gLV Sin A25")
#plot(ngsA25, main="gLV Sin A25")

ngsA50=graph_from_adjacency_matrix(gsA50, mode="directed", weighted=TRUE)
createNetworkFromIgraph(ngsA50, "gLV Sin A50")
#plot(ngsA50, main="gLV Sin A50")

ngsA75=graph_from_adjacency_matrix(gsA75, mode="directed", weighted=TRUE)
createNetworkFromIgraph(ngsA75, "gLV Sin A75")
#plot(ngsA75, main="gLV Sin A75")

ngsA100=graph_from_adjacency_matrix(gsA100, mode="directed", weighted=TRUE)
createNetworkFromIgraph(ngsA100, "gLV Sin A100")
#plot(ngsA100, main="gLV Sin A100")

nssA0=graph_from_adjacency_matrix(ssA0, mode="directed", weighted=TRUE)
createNetworkFromIgraph(nssA0, "SOI Sin A0")
#plot(nssA0, main="SOI Sin A0")

nssA25=graph_from_adjacency_matrix(ssA25, mode="directed", weighted=TRUE)
createNetworkFromIgraph(nssA25, "SOI Sin A25")
#plot(nssA25, main="SOI Sin A25")

nssA50=graph_from_adjacency_matrix(ssA50, mode="directed", weighted=TRUE)
createNetworkFromIgraph(nssA50, "SOI Sin A50")
#plot(nssA50, main="SOI Sin A50")

#Analisis de cytoscape: análisis de la red
aConocida=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/Conocida.csv", header = TRUE, sep=",")
aConocidaA0=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/ConocidaA0.csv", header = TRUE, sep=",")
aConocidaA25=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/ConocidaA25.csv", header = TRUE, sep=",")
aConocidaA50=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/ConocidaA50.csv", header = TRUE, sep=",")
aConocidaA75=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/ConocidaA75.csv", header = TRUE, sep=",")
aConocidaA100=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/ConocidaA100.csv", header = TRUE, sep=",")


aRickerCon=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/Ricker_Con.csv", header = TRUE, sep=",")
aRickerSin=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/Ricker_Sin.csv", header = TRUE, sep=",")
aRickerSinA0=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/Ricker_SinA0.csv", header = TRUE, sep=",")
aRickerSinA25=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/Ricker_SinA25.csv", header = TRUE, sep=",")

agLVCon=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_Con.csv", header = TRUE, sep=",")
agLVSin=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_Sin.csv", header = TRUE, sep=",")
agLVSinA0=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_SinA0.csv", header = TRUE, sep=",")
agLVSinA25=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_SinA25.csv", header = TRUE, sep=",")
agLVSinA50=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_SinA50.csv", header = TRUE, sep=",")
agLVSinA75=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_SinA75.csv", header = TRUE, sep=",")
agLVSinA100=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/gLV_SinA100.csv", header = TRUE, sep=",")

aSOICon=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/SOI_Con.csv", header = TRUE, sep=",")
aSOISin=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/SOI_Sin.csv", header = TRUE, sep=",")
aSOISinA0=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/SOI_SinA0.csv", header = TRUE, sep=",")
aSOISinA25=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/SOI_SinA25.csv", header = TRUE, sep=",")
aSOISinA50=read.table("~/Microbiología/Tesis/Tesis2/Tesis2/SOI_SinA50.csv", header = TRUE, sep=",")

#Transitivity (Tiende a disminuir pero en naA100 aumenta casi igual a naA0)

transitivity(naA0)
transitivity(na)
transitivity(naA25)
transitivity(naA50)
transitivity(naA75)
transitivity(naA100)
transitivity(nrc)
transitivity(nrs)
transitivity(nrsA0)
transitivity(nrsA25)
transitivity(ngc)
transitivity(ngs)
transitivity(ngsA0)
transitivity(ngsA25)
transitivity(ngsA50)
transitivity(ngsA75)
transitivity(ngsA100)
transitivity(nhc)
transitivity(nhs)
transitivity(nsc)
transitivity(nss)
transitivity(nssA0)
transitivity(nssA25)
transitivity(nssA50)


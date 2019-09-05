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
N = 50 
S = 40
rownames(A)=c(1:N)
colnames(A)=rownames(A)
rA=plotA(A, method = "network")
#####ESTO ES LO QUE NO CORRE
createNetworkFromGraph(rA, "Red Conocida")

rownames(rs)=c(1:N)
colnames(rs)=rownames(rs)
plotA(rs, method = "network")

rownames(rc)=c(1:N)
colnames(rc)=rownames(rc)
plotA(rc, method = "network")

rownames(gs)=c(1:N)
colnames(gs)=rownames(gs)
plotA(gs, method = "network")

rownames(gc)=c(1:N)
colnames(gc)=rownames(gc)
plotA(gc, method = "network")

rownames(hs)=c(1:N)
colnames(hs)=rownames(hs)
plotA(hs, method = "network")

rownames(hc)=c(1:N)
colnames(hc)=rownames(hc)
plotA(hc, method = "network")

rownames(ss)=c(1:N)
colnames(ss)=rownames(ss)
plotA(ss, method = "network")

rownames(sc)=c(1:N)
colnames(sc)=rownames(sc)
plotA(sc, method = "network")



#NO CORRER ESTO
 setwd("~/Microbiología/Tesis/Tesis2/Tesis2")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/Matriz.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/MatrizKlemm.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/abundanciasgLV.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimegLVCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimegLVSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeHubbellCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeHubbellSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeRickerCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeRickerSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeSOICon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/TimeSOISin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedgLVCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedgLVSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedHubbellCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedHubbellSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedRickerCon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedRickerSin.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedSOICon.RData")
 load("~/Microbiología/Tesis/Tesis2/Tesis2/RedSOISin.RData")

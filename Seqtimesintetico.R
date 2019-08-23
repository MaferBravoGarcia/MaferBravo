#María Fernanda Bravo García
#Manual Seq.time
#Librerias 
library(seqtime)
library(ggplot2)
library(reshape2)

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
save(out.rickers, file= "TimeRickerSin")
## the gLV
out.glvs = glv (N, A=A, tend=25)
tsplot(out.glvs, main="gLV sin")
save(out.glvs, file = "TimegLVSin" )
## Hubbell
out.hubbells=simHubbell(N=N, M=N, I=1500, d=N, m=0.1, tskip=250, tend=500)
tsplot(out.hubbells, main="Hubbell sin")

## SOI (I= individuos? )
out.sois= soi(N, A=A, I=1500, tend=250)
tsplot(out.sois, main= "SOI sin")


#CON PERTURBACIÓN
pert=perturbation(times = c(100), durations = (50))
##the Ricker community model
out.rickerc=ricker(N,A=A, tend=250, perturb = pert)
tsplot(out.rickerc,main="Ricker con ")

## the gLV ( Please provide as many growth changes as species!)
out.glvc = glv (N, A=A, tend=250, perturb = pert)
tsplot(out.glvc, main="gLV con")

## Hubbell
out.hubbellc=simHubbell(N=N, M=N, I=1500, d=N, m=0.1, tskip=250, tend=500, perturb = pert)
tsplot(out.hubbellc, main="Hubbell con")

## SOI (I= individuos? )
out.soic= soi(N, A=A, I=1500, tend=250, perturb = pert)
tsplot(out.soic, main= "SOI con")

#5) LIMITS to test how well it can infer the known interaction matrix from the time series
mfrow=c(1,1)
rs=limits(out.rickers)$rs
save(rs, file = "Red Ricker sin")

rc=limits(out.rickerc)$rc
save(rc, file="Red Ricker con")

gs=limits(out.glvs)$gs
save(gs, file= "Red gLV sin")

gc=limits(out.glvc)$gc
save(gc, file="Red gLV con")

hs=limits(out.hubbells)$hs
save(hs, file="Red Hubbell sin")

hc=limits(out.hubbellc)$hc
save(hc, file="Red Hubbell con")

ss=limits(out.sois)$ss
save(ss, "Red SOI sin")

sc=limits(out.soic)$sc
save(sc, file="Red SOI con")



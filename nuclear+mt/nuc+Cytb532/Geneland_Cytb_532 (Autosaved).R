library(Geneland)
setwd("~/MegaSync/Doutorado/Dados/@Filogeografia/Geneland/nuclear+mt/nuc+Cytb532")
#### Files to be used

geno.temp <- read.table("input_nuc_geneland.txt",na.string="-999")
geno <- as.matrix(geno.temp)
geno

coord.temp <- read.table("coordinates.txt",na.string='-999')
coord <- as.matrix(coord.temp)
coord

mt.temp <- read.table('Cytb_geneland_532.txt', na.string='-999')  ###desmarcar quando usar mitocondrial
mt <- as.matrix(mt.temp)
mt

###se tiver dado morfologicos:

#morf.temp <- read.table('morfologia.txt', na.string='NA')
#morf <- as.matrix(morf.temp)

all.freq.model='Uncorrelated' #### ucorrelated or correlated allele frequencies



############Loop for multiple runs#######################################################################

maindir <- "~/MegaSync/Doutorado/Dados/@Filogeografia/Geneland/nuclear+mt/nuc+Cytb532/Result" ##só nuclear
#maindir <- "/Users/felipemedeiros/Desktop/Ischnocnema/@Geneland/results10numtDNA" ##nuclear e mit


nrun <- 3        ### number of runs ### obs: aumentar pra no minimo 10
burnin <- 200     ### burnin ####aumentar pra 200
chain <- 1000000  ### number of steps in chain ###aumentar pra 5000000
freq <-      1000 ### sampling frequency 5000
pop <- 8          ### max number of pops 

for(irun in 1:nrun)
{
  
  path.mcmc <- paste(maindir,'/',irun,"/",sep="")
  
  
  ### create directory
  if (file.exists(path.mcmc)){
    setwd(file.path(path.mcmc))
  } else {
    dir.create(file.path(path.mcmc))
    setwd(file.path(path.mcmc))
  }
  
  ####
  
  system(paste(path.mcmc))
  MCMC(coordinates=coord,
       geno.dip.codom=geno,
       geno.hap=mt, #desmarcar quando usar mitocondrial
       rate.max=172, #número de individuos
       nb.nuclei.max=500, #até 3x o numero de indivíduos
       varnpop=TRUE,
       delta.coord=0.5,
       npopmax=pop,
       spatial=TRUE,
       freq.model="Uncorrelated",
       nit=chain,
       thinning=freq,
       path.mcmc=path.mcmc)
  
  ## MCMC postprocessing
  PostProcessChain(coordinates=coord,
                   path.mcmc=path.mcmc,
                   nxdom=100, #antes:200
                   nydom=100, #antes:200
                   burnin=burnin)
}


############ log posterior density calculation ###
lpd <- rep(NA,nrun)

for(irun in 1:nrun)
{
  path.mcmc <- paste(maindir,'/',irun,"/",sep="")
  path.lpd <- paste(path.mcmc,"log.posterior.density.txt",sep="")
  lpd[irun] <- mean(scan(path.lpd)[-(1:burnin)])
}


order(lpd,decreasing=TRUE)
lpd


##### plot number of pops


for(irun in 1:nrun)
{
  path.mcmc <- paste(maindir,'/',irun,"/",sep="")
  setwd(file.path(path.mcmc))
  Plotnpop(path.mcmc=path.mcmc, burnin=200,printit=TRUE,
           file="Number_of_Clusters.pdf",format="pdf")
}


for(irun in 1:nrun)
{
  path.mcmc <- paste(maindir,'/',irun,"/",sep="")
  setwd(file.path(path.mcmc))
  PosteriorMode(coordinates=coord, path.mcmc="./", file="map.pdf",printit=TRUE,format="pdf")
}



#### map of cluster membership probability
x11()
PlotTessellation(coord, path.mcmc= './',printit=TRUE,path='./')


library(dplyr)


library(dartR)


gl.install.vanilla.dartR(flavour="dev_bernd")

packageVersion("dartR")
library(ggplot2)






#number of populations

n.pop=3

#number of individuals

n.ind=50

#number of loci (snps)

n.loc=500



#number of generations to run

n.gen=20



#mutation rate

m.rate = 1e-7



#start with a list of population based on pstart

p <- list()

pstart <- glSim(n.ind = n.ind, n.snp.nonstruc = n.loc, ploidy = 2 )


for (i in 1:n.pop)
  
{
  
  p[[i]] <- gl.sim.ind(pstart, n = n.ind, popname = paste0("P",i))
  
  indNames(p[[i]])<- paste0("p",i,"_",1:n.ind)
  
}


names(p) <- paste0("P",1:n.pop)






px <- c(1,2)

py <- c(2,2)



plot(px, py, cex=7)

points(px, py, cex=3, pch=as.character(1:length(px)))



#emigration (constant #individuals)

#emi.table <- matrix(c(0,0,5,0), 2,2)

emi.table <- matrix(c(0,5,0,0,0,5,0,0,0), 3,3)

emi.table



#generation zero

pp <- do.call("rbind", p)







for(i in 1:n.gen)
  
{
  
  ### reproduce and mutate each population
  
  for (ii in 1:n.pop)
    
  {
    
    p[[ii]] <- gl.sim.ind(p[[ii]],n = n.ind, popname =  paste0("P",ii))
    
    p[[ii]] <- gl.sim.mutate(p[[ii]], m.rate)
    
  }
  
  
  ###disperse between population according to emi.table
  
  p <- gl.sim.emigration(p, 0.1, emi.table = emi.table)
  
  
  
  if (!(i %% 5)) {
    
    cat(paste("Generation:", i,"/" ,n.gen, "\n"))
    
    flush.console()
    
  }
  
  #combine the populations
  
  pp <- do.call("rbind", p)
  
}



#prettify output and create a full genlight object

ploidy(pp)<- rep(2, nrow(pp))

locNames(pp)<- 1:nLoc(pp)

indNames(pp) <- 1:nrow(pp)

pp <- gl.compliance.check(pp, verbose = 0)

# create coordinates

tt <- table(pop(pp))

pp@other$latlon <- data.frame(lon =rep(px,tt)+runif(nrow(pp)*0.2), lat=rep(py,tt)+runif(nrow(pp)*0.2))





pas <- gl.report.pa(pp, verbose=0)

pas

mm <- matrix(0, nPop(pp),nPop(pp))

for (i in 1:nrow(pas)) mm[pas[i,1], pas[i,2]] <- pas$priv2[i]

for (i in 1:nrow(pas)) mm[pas[i,2], pas[i,1]] <- pas$priv1[i]

mm

gl.map.interactive(pp, matrix = mm, symmetric = FALSE)


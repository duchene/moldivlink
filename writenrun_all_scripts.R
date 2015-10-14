settings.data <- read.table("Settings.csv", head = T, sep = ",", stringsAsFactors = F, row.names = 1)

for(j in c(1:150)){
    file.name <- settings.data$runName[j]
    if(file.name %in% dir()) next
    system(paste("mkdir", file.name))
    setwd(file.name)

### l0 loads the beast programs and xml file writer.
    
    l0 <- 'log.ann.path <- "~/PortableDucheneHo/BEAST1.7.5/bin/loganalyser"
tree.ann.path <- "~/PortableDucheneHo/BEAST1.7.5/bin/treeannotator"
code.path <- "~/Portabletraitdivmolproject/apoint_scripts/"

require(phangorn)
require(geiger)
require(phytools)

dir.cur <- getwd()
setwd(code.path)
for(i in dir()) source(i)
setwd(dir.cur)'

# l1 and l2 simulate the data, and there is a while loop in case there is an error and it needs to repeat.

	l1 <- paste0("tree.sim <- try(tr.mu.sp(traitstart = ", settings.data$traitstart[j], ", trait.r = ", settings.data$trait.r[j], ", direct = ", settings.data$direct[j], ", sprerror = ", settings.data$sprerror[j], ", regcoefspr = ", settings.data$regcoefspr[j], ", regcoefmu = ", settings.data$regcoefmu[j], "))")
    
    l2 <- paste0("while(class(tree.sim) == 'try-error')", "tree.sim <- try(tr.mu.sp(traitstart = ", settings.data$traitstart[j], ", trait.r = ", settings.data$trait.r[j], ", direct = ", settings.data$direct[j], ", sprerror = ", settings.data$sprerror[j], ", regcoefspr = ", settings.data$regcoefspr[j], ", regcoefmu = ", settings.data$regcoefmu[j], "))")
    
# l3 runs beast and R8S.

    l3 <- c(paste0('get.xml2(list(list(max(branching.times(tree.sim[[1]])), tree.sim[[1]]$tip.label)), list(tree.sim[[1]], as.matrix(tree.sim[[2]])))'),
    'system("beast sim.xml")',
    paste0('write.tree(tree.sim[[1]], file = paste0("', file.name, '", "_sim.tre"))'),
    paste0('save(tree.sim, file = paste0("', file.name, '", ".Rdata"))'),
    paste0('loganalyser.command <- paste(log.ann.path, "-burnin 7500 sim.log", paste0("', file.name,'", "_log.csv"))'),
    'system(loganalyser.command)',
    'treeannotator.command <- paste(tree.ann.path, "-burnin 7500 sim.trees", "BEASTestimated.tree")',
    'system(treeannotator.command)',
    'erase.files <- grep("[.]log|[.]trees|[.]txt|[.]state|[.]e|[.]o", dir(), value=T)',
    'sapply(erase.files, function(x) system(paste("rm", x)))',
    'beastchron <- read.nexus("BEASTestimated.tree")',
    'substchron <- optim.pml(pml(NJ(dist.dna(tree.sim[[2]], model = "GG95")), as.phyDat(tree.sim[[2]]), model = "GTR"), optNni = T)$tree',
    'roottime <- max(branching.times(tree.sim[[1]]))',
    'l <- 10^(-1:6)',
    'plchrons <- list()',
    'cv <- vector()',
    'for (i in 1:8){',
    'plchrons[[i]] <- chronopl(substchron, lambda = l[i], age.min = roottime*0.95, age.max = roottime*1.05, CV= T)',
    'cv <- c(cv, sum(attr(plchrons[[i]], "D2")))',
    '}',
    'plchron <- plchrons[[which(cv == min(cv))]]',
    'write.tree(plchron, "NPRSestimated.tree")'
    )
    
    runfile <- c(l0, l1, l2, l3)
    writeLines(runfile, con = paste0(file.name, "_script.R"))
    
# s1 and s2 create the bash code to start the run.

    s1 <- "#!/bin/csh
#PBS -l wd 
#PBS -q normal
#PBS -l walltime=20:00:00,jobfs=4000Mb
#PBS -l mem=20144MB
#PBS -l software=R
setenv TMPDIR $PBS_JOBFS
module load java
module load beast/2.1.1
module load R/3.0.2"

	s2 <- paste("R --vanilla <", paste0(file.name, "_script.R"))
    
    runscript <- c(s1, s2)
    writeLines(runscript, con = paste0(file.name, "_job.sh"))
    submit.command <- paste("qsub", paste0(file.name, "_job.sh"))
    system(submit.command)
    setwd("..")
    print(paste("Completed scripts for", file.name))

}



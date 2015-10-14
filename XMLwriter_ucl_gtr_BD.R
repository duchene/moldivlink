get.xml2 <- function(cal.data, sim.data, file.name = "sim"){

# BLOCK1 - Header

b1 <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">
"

# BLOCK2 - Alignment

b2 <- c("<data", paste0("id=\"simulationduchene\""), "name=\"alignment\">")

for(i in 1:nrow(sim.data[[2]])){
    sec <- paste0("<sequence id=\"seq_", rownames(sim.data[[2]])[i], "\" taxon=\"", rownames(sim.data[[2]])[i], "\" totalcount=\"4\" value=\"", paste(toupper(sim.data[[2]][i, ]), collapse = ""), "\"/>")
    b2 <- c(b2, sec)
}
b2 <- c(b2, "</data>")

# BLOCK3 - Mathematical distributions and tree

b3 <- "<map name=\"Beta\">beast.math.distributions.Beta</map>
<map name=\"Exponential\">beast.math.distributions.Exponential</map>
<map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>
<map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>
<map name=\"Gamma\">beast.math.distributions.Gamma</map>
<map name=\"Uniform\">beast.math.distributions.Uniform</map>
<map name=\"prior\">beast.math.distributions.Prior</map>
<map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>
<map name=\"OneOnX\">beast.math.distributions.OneOnX</map>
<map name=\"Normal\">beast.math.distributions.Normal</map>
"

tr.topo <- sim.data[[1]]
tr.topo$edge.length <- abs(rnorm(n = length(tr.topo$edge.length)))

b3 <- c(b3, paste0("<tree id=\"Tree.t:simulationduchene\" spec=\'beast.util.TreeParser\' newick=\"", write.tree(tr.topo), "\""), "taxa=\"@simulationduchene\"/>")

# BLOCK4 - State and population model

b4 <- "<run chainLength=\"50000000\" id=\"mcmc\" spec=\"MCMC\">
    <state id=\"state\" storeEvery=\"5000\">
        <input name=\'stateNode\' idref=\'Tree.t:simulationduchene\'/>
	<parameter dimension=\"4\" id=\"freqParameter.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.25</parameter>
	<parameter id=\"rateAC.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
	<parameter id=\"rateAG.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
	<parameter id=\"rateAT.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
	<parameter id=\"rateCG.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
	<parameter id=\"rateGT.s:simulationduchene\" lower=\"0.0\" name=\"stateNode\">1.0</parameter>
	<parameter id=\"gammaShape.s:simulationduchene\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"ucldStdev.c:simulationduchene\" lower=\"0.0\" name=\"stateNode\" upper=\"5.0\">0.5</parameter>
        <stateNode dimension=\"162\" id=\"rateCategories.c:simulationduchene\" spec=\"parameter.IntegerParameter\">1</stateNode>
	<parameter id=\"birthRate2.t:simulationduchene\" lower=\"0.0\" name=\"stateNode\" upper=\"10000.0\">1.0</parameter>
	<parameter id=\"relativeDeathRate2.t:simulationduchene\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">0.5</parameter>
        <parameter id=\"ucldMean.c:simulationduchene\" name=\"stateNode\">1.0</parameter>
    </state>
"

# BLOCK5 - Prior 1

b5 <- "<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
            <distribution birthDiffRate=\"@birthRate2.t:simulationduchene\" id=\"BirthDeath.t:simulationduchene\" relativeDeathRate=\"@relativeDeathRate2.t:simulationduchene\" spec=\"beast.evolution.speciation.BirthDeathGernhard08Model\" tree=\"@Tree.t:simulationduchene\"/>
            <prior id=\"BirthRatePrior.t:simulationduchene\" name=\"distribution\" x=\"@birthRate2.t:simulationduchene\">
                <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"3000.0\"/>
            </prior>
	    <prior id=\"DeathRatePrior.t:simulationduchene\" name=\"distribution\" x=\"@relativeDeathRate2.t:simulationduchene\">
                <Uniform id=\"Uniform.011\" name=\"distr\"/>
            </prior>
            <prior id=\"GammaShapePrior.s:simulationduchene\" name=\"distribution\" x=\"@gammaShape.s:simulationduchene\">
	    	<Exponential id=\"Exponential.0\" name=\"distr\">
		    <parameter id=\"RealParameter.0\" lower=\"0.0\" name=\"mean\" upper=\"0.0\">1.0</parameter>
		</Exponential>
	    </prior>
	    <prior id=\"RateACPrior.s:simulationduchene\" name=\"distribution\" x=\"@rateAC.s:simulationduchene\">
	    	<Gamma id=\"Gamma.0\" name=\"distr\">
		    <parameter id=\"RealParameter.01\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
		    <parameter id=\"RealParameter.02\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
		</Gamma>
	    </prior>
	    <prior id=\"RateAGPrior.s:simulationduchene\" name=\"distribution\" x=\"@rateAG.s:simulationduchene\">
                <Gamma id=\"Gamma.01\" name=\"distr\">
		           <parameter id=\"RealParameter.03\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.04\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">20.0</parameter>
                </Gamma>
	    </prior>
	    <prior id=\"RateATPrior.s:simulationduchene\" name=\"distribution\" x=\"@rateAT.s:simulationduchene\">
                <Gamma id=\"Gamma.02\" name=\"distr\">
		           <parameter id=\"RealParameter.05\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.06\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
	    </prior>
	    <prior id=\"RateCGPrior.s:simulationduchene\" name=\"distribution\" x=\"@rateCG.s:simulationduchene\">
                <Gamma id=\"Gamma.03\" name=\"distr\">
		           <parameter id=\"RealParameter.07\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.08\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
	    </prior>
	    <prior id=\"RateGTPrior.s:simulationduchene\" name=\"distribution\" x=\"@rateGT.s:simulationduchene\">
                <Gamma id=\"Gamma.04\" name=\"distr\">
		           <parameter id=\"RealParameter.09\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\">0.05</parameter>
                    <parameter id=\"RealParameter.010\" lower=\"0.0\" name=\"beta\" upper=\"0.0\">10.0</parameter>
                </Gamma>
	    </prior>
	    <prior id=\"MeanRatePrior.c:simulationduchene\" name=\"distribution\" x=\"@ucldMean.c:simulationduchene\">
                <Uniform id=\"Uniform.01\" name=\"distr\" upper=\"Infinity\"/>
            </prior>
            <prior id=\"ucldStdevPrior.c:simulationduchene\" name=\"distribution\" x=\"@ucldStdev.c:simulationduchene\">
                <Exponential id=\"Exponential.01\" name=\"distr\">
                    <parameter estimate=\"true\" id=\"RealParameter.011\" name=\"mean\">1.0</parameter>
                </Exponential>
            </prior>
"

# BLOCK6 - Prior 2 - Calibrations

b6 <- vector()

for(i in 1:length(cal.data)){
    cal.name <- paste0("cal", i)
    cal.block1 <- c(paste0("<distribution id=\"", cal.name, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:simulationduchene\">"), paste0("<taxonset id=\"", cal.name, "\" spec=\"TaxonSet\">")) 
    prevcals <- vector()
    if(i > 1){
	    for(j in 1:(i-1)){
    		prevcals <- c(prevcals, cal.data[[j]][[2]])
    	}
    }
    cal.block2 <- vector()
    for(j in 1:length(cal.data[[i]][[2]])){
    	if(cal.data[[i]][[2]][j] %in% prevcals){
    		cal.block2 <- c(cal.block2, paste0("<taxon idref=\"", cal.data[[i]][[2]][j], "\"/>"))
    	} else {
    		cal.block2 <- c(cal.block2, paste0("<taxon id=\"", cal.data[[i]][[2]][j], "\" spec=\"Taxon\"/>"))
    	}
    }
    cal.block3 <- c("</taxonset>", paste0("<Normal id=\"Normal.0", i, "\" name=\"distr\">"), paste0("<parameter estimate=\"true\" id=\"RealParameter.a", i, "\" name=\"mean\">", round(cal.data[[i]][[1]], 3), "</parameter>"), paste0("<parameter estimate=\"true\" id=\"RealParameter.b", i, "\" name=\"sigma\">", round((cal.data[[i]][[1]] / 10000), 5), "</parameter>"), "</Normal>", "</distribution>")
    b6 <- c(b6, cal.block1, cal.block2, cal.block3)
}

b6 <- c(b6, "</distribution>")

# BLOCK7

b7 <- "<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">
            <distribution data=\"@simulationduchene\" id=\"treeLikelihood.simulationduchene\" spec=\"TreeLikelihood\" tree=\"@Tree.t:simulationduchene\">
                <siteModel gammaCategoryCount=\"4\" id=\"SiteModel.s:simulationduchene\" shape=\"@gammaShape.s:simulationduchene\" spec=\"SiteModel\">
		    <parameter estimate=\"false\" id=\"mutationRate.s:simulationduchene\" name=\"mutationRate\">1.0</parameter>
		    <parameter estimate=\"false\" id=\"proportionInvariant.s:simulationduchene\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\">0.0</parameter>
		    <substModel id=\"gtr.s:simulationduchene\" rateAC=\"@rateAC.s:simulationduchene\" rateAG=\"@rateAG.s:simulationduchene\" rateAT=\"@rateAT.s:simulationduchene\" rateCG=\"@rateCG.s:simulationduchene\" rateGT=\"@rateGT.s:simulationduchene\" spec=\"GTR\">
                        <parameter estimate=\"false\" id=\"rateCT.s:simulationduchene\" lower=\"0.0\" name=\"rateCT\">1.0</parameter>
			<frequencies frequencies=\"@freqParameter.s:simulationduchene\" id=\"estimatedFreqs.s:simulationduchene\" spec=\"Frequencies\"/>
                    </substModel>
                </siteModel>
                <branchRateModel clock.rate=\"@ucldMean.c:simulationduchene\" id=\"RelaxedClock.c:simulationduchene\" rateCategories=\"@rateCategories.c:simulationduchene\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"@Tree.t:simulationduchene\">
                    <LogNormal S=\"@ucldStdev.c:simulationduchene\" id=\"LogNormalDistributionModel.c:simulationduchene\" meanInRealSpace=\"true\" name=\"distr\">
                        <parameter estimate=\"true\" id=\"RealParameter.012\" lower=\"0.0\" name=\"M\" upper=\"1.0\">1.0</parameter>
                    </LogNormal>
                </branchRateModel>
            </distribution>
        </distribution>
    </distribution>
"

# BLOCK8 - Operators

b8 <- "<operator id=\"BirthRateScaler.t:simulationduchene\" parameter=\"@birthRate2.t:simulationduchene\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"DeathRateScaler.t:simulationduchene\" parameter=\"@relativeDeathRate2.t:simulationduchene\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"treeScaler.t:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:simulationduchene\" weight=\"3.0\"/>

    <operator id=\"treeRootScaler.t:simulationduchene\" rootOnly=\"true\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:simulationduchene\" weight=\"3.0\"/>

    <operator id=\"UniformOperator.t:simulationduchene\" spec=\"Uniform\" tree=\"@Tree.t:simulationduchene\" weight=\"30.0\"/>

    <operator delta=\"0.01\" id=\"FrequenciesExchanger.s:simulationduchene\" spec=\"DeltaExchangeOperator\" weight=\"0.1\">
    	<parameter idref=\"freqParameter.s:simulationduchene\"/>
    </operator>

    <operator id=\"RateACScaler.s:simulationduchene\" parameter=\"@rateAC.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>
    
    <operator id=\"RateAGScaler.s:simulationduchene\" parameter=\"@rateAG.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateATScaler.s:simulationduchene\" parameter=\"@rateAT.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateCGScaler.s:simulationduchene\" parameter=\"@rateCG.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"RateGTScaler.s:simulationduchene\" parameter=\"@rateGT.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"gammaShapeScaler.s:simulationduchene\" parameter=\"@gammaShape.s:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>

    <operator id=\"ucldStdevScaler.c:simulationduchene\" parameter=\"@ucldStdev.c:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"CategoriesRandomWalk.c:simulationduchene\" parameter=\"@rateCategories.c:simulationduchene\" spec=\"IntRandomWalkOperator\" weight=\"10.0\" windowSize=\"1\"/>

    <operator id=\"CategoriesSwapOperator.c:simulationduchene\" intparameter=\"@rateCategories.c:simulationduchene\" spec=\"SwapOperator\" weight=\"10.0\"/>

    <operator id=\"CategoriesUniform.c:simulationduchene\" parameter=\"@rateCategories.c:simulationduchene\" spec=\"UniformOperator\" weight=\"10.0\"/>

    <operator id=\"ucldMeanScaler.c:simulationduchene\" parameter=\"@ucldMean.c:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"1.0\"/>

    <operator id=\"relaxedUpDownOperator.c:simulationduchene\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
        <parameter idref=\"ucldMean.c:simulationduchene\" name=\"up\"/>
        <tree idref=\"Tree.t:simulationduchene\" name=\"down\"/>
    </operator>
"

# BLOCK 9 - Loggers

b9 <- paste0("<logger fileName=\"", file.name, ".log\" id=\"tracelog\" logEvery=\"5000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">")

b9 <- c(b9, "<log idref=\"posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
        <log idref=\"treeLikelihood.simulationduchene\"/>
        <log id=\"TreeHeight.t:simulationduchene\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:simulationduchene\"/>
        <log idref=\"BirthDeath.t:simulationduchene\"/>
	<parameter idref=\"freqParameter.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"rateAC.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"rateAG.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"rateAT.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"rateCG.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"rateGT.s:simulationduchene\" name=\"log\"/>
	<parameter idref=\"gammaShape.s:simulationduchene\" name=\"log\"/>
        <parameter idref=\"ucldStdev.c:simulationduchene\" name=\"log\"/>
        <log branchratemodel=\"@RelaxedClock.c:simulationduchene\" id=\"rate.c:simulationduchene\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"@Tree.t:simulationduchene\"/>
        <parameter idref=\"ucldMean.c:simulationduchene\" name=\"log\"/>
")

cals <- vector()
for(i in 1:length(cal.data)){
	cals <- c(cals, paste0("<log idref=\"cal", i, ".prior\"/>"))
}

b9 <- c(b9, cals, "</logger>")

b10 <- paste0("<logger id=\"screenlog\" logEvery=\"5000\">
        <log idref=\"posterior\"/>
        <log arg=\"@posterior\" id=\"ESS.0\" spec=\"util.ESS\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
    </logger>

    <logger fileName=\"", file.name, ".trees\" id=\"treelog.t:simulationduchene\" logEvery=\"5000\" mode=\"tree\">
        <log branchratemodel=\"@RelaxedClock.c:simulationduchene\" id=\"TreeWithMetaDataLogger.t:simulationduchene\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:simulationduchene\"/>
    </logger>

</run>

</beast>
")

complete.xml <- c(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)

writeLines(complete.xml, con = paste0(file.name, ".xml"))

}

#' Phylogenetic tree to xml for beast2
#'
#' This function creates a xml document from a phylogenetic tree in newick format. 
#' The xml document is to be used for phylogeographic analysis with beast2.
#' @import ape
#' @import XML
#' @param save_path Where to save the document 
#' @param tree_path Where the tree is stored 
#' @param label_sep Symbol separating information in tip labels 
#' @param location_pos Position of location in tip label 
#' @param Name Name of the output files produced by beast2 
#' @param chain_length Mcmc 
#' @param chain_store Mcmc 
#' @param state_store Mcmc 
#' @param tracelog_logevery Mcmc 
#' @param screenlog_logevery Mcmc 
#' @param treelog_logevery Mcmc 
#' @param treeWithTraitlog_logevery Mcmc
#'
#' @return Xml document.
#'
#' @export
beast2.xml.parser <- function(save_path, tree_path, label_sep, location_pos, Name, chain_length, chain_store, state_store, tracelog_logevery, screenlog_logevery, treelog_logevery, treeWithTraitlog_logevery){
# example 
# beast2.xml.parser(save_path='~/Desktop/testnew_2.xml', tree_path='~/Desktop/Imperial College/Urop/Data/Subtype_01_AE_0.newick', label_sep='_', location_pos=1, Name='01_AE_01', chain_length=3000000, chain_store=300, state_store=5000, tracelog_logevery=2000, screenlog_logevery=10000, treelog_logevery=2000, treeWithTraitlog_logevery=10000)	

#disable scientific notation
options(scipen=999)

tree <- ape::read.tree(paste(tree_path))
tree_newick <- ape::write.tree(tree, file='');

bxml <- XML::newXMLDoc()

root <- XML::newXMLNode('beast', 
attrs= c(beautitemplate='Standard', beautistatus='noAutoSetClockRate', namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood', required='', version='2.6'), doc=bxml)

# create data frame with information for alignment
df <- data.frame( 'id_df' = paste('seq_', tree$tip.label, sep=''), 'spec_df' = 	'Sequence', 'taxon_df'= tree$tip.label)
	
	# filter out location (and date) 	
	data_split <- strsplit(as.character(df$taxon_df), split=label_sep)
	df$location <- c(unlist(lapply(data_split, function(x) x[[location_pos]])))
	#df$date <- lapply(data_split, function(x) x[[date_pos]]) 
	
#####################################################################

#add alignment/data node to XML

	data_node <- XML::newXMLNode('data',
		attrs= c(id=paste(Name), spec='Alignment',name='alignment'), 
		parent=root, 
		doc=bxml)
		
	for(i in seq( 1, length(df$taxon_df))) {
		dummy <- XML::newXMLNode('sequence', 
		attrs = c(id=paste(df$id_df[i]), spec=paste(df$spec_df[i]), taxon=paste(df$taxon_df[i]), 
		totalcount='4', value='-'), doc=bxml, parent=data_node)
	}


#####################################################################

# add map nodes


dummy <- XML::newXMLNode('map', attrs=c(name='Uniform'),
	'beast.math.distributions.Uniform', 
	doc=bxml, 
	parent=root);	
dummy <- XML::newXMLNode('map', attrs=c(name='Exponential'),
	'beast.math.distributions.Exponential', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='LogNormal'),
	'beast.math.distributions.LogNormalDistributionModel', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='Normal'),
	'beast.math.distributions.Normal', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='Beta'),
	'beast.math.distributions.Beta', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='Gamma'),
	'beast.math.distributions.Gamma', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='LaplaceDistribution'),
	'beast.math.distributions.LaplaceDistribution', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='prior'),
	'beast.math.distributions.Prior', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='InverseGamma'),
	'beast.math.distributions.InverseGamma', 
	doc=bxml, 
	parent=root);
dummy <- XML::newXMLNode('map', attrs=c(name='OneOnX'),
	'beast.math.distributions.OneOnX', 
	doc=bxml, 
	parent=root);
	
#####################################################################

#add run node all mcmc specifics are captured in here
	
	
	run_node <- XML::newXMLNode('run', attrs=c(id='mcmc', spec='MCMC',
			chainLength=paste(chain_length), storeEvery=paste(chain_store)),
			doc=bxml, parent=root);
	
		state_node <- XML::newXMLNode('state', 
			attrs=c(id='state', spec='State', storeEvery=paste(state_store)),
			parent=run_node, doc=bxml);
			
			tree_node <- XML::newXMLNode('tree', 
				attrs=c(id=paste('Tree.t:', Name, sep=''),
				spec='beast.evolution.tree.Tree', name='stateNode'), 
				parent=state_node, doc=bxml);
			
				taxonset_node <- XML::newXMLNode('taxonset',
					attrs=c(id=paste('TaxonSet.',Name, sep=''), spec='TaxonSet'),
					parent=tree_node, doc=bxml);
					
					dummy <- XML::newXMLNode('alignment', attrs=c(idref=paste(Name)),
					parent=taxonset_node, doc=bxml);
			dummy <- XML::newXMLNode ('stateNode',
				attrs=c(id='rateIndicator.s:location',
				spec='parameter.BooleanParameter', dimension=paste(length(unique(df$location, incomparables='false'))*(length(unique(df$location, incomparables='false'))-1)/2)), 
				'true',parent=state_node, doc=bxml);
			dummy <- XML::newXMLNode('parameter', 
				attrs=c(id='relativeGeoRates.s:location',
				spec='parameter.RealParameter', dimension=paste(length(unique(df$location, incomparables='false'))*(length(unique(df$location, incomparables='false'))-1)/2), name='stateNode'), 
				'1.0', parent=state_node, doc=bxml);
			dummy <- XML::newXMLNode('parameter', 
				attrs=c(id='traitClockRate.c:location',
				spec='parameter.RealParameter', name='stateNode'), 
				'1.0', parent=state_node, doc=bxml);

			dummy <<- XML::newXMLNode('init', 
				attrs=c(spec='beast.util.TreeParser', id=paste('NewickTree.t:',
				Name, sep=''), IsLabelledNewick='true', adjustTipHeights='false',
				initial=paste('@Tree.t:', Name, sep=''), 
				taxa=paste('@', Name, sep=''), newick=paste(tree_newick)),parent=run_node, doc=bxml);
				
				
		distribution_node1 <- XML::newXMLNode('distribution', 
			attrs=c(id='posterior', spec='util.CompoundDistribution'), 
			parent=run_node, doc=bxml);
			
			distribution_node2 <- XML::newXMLNode('distribution', 
				attrs=c(id='prior', spec='util.CompoundDistribution'), 
				parent=distribution_node1, doc=bxml);
				
				
				prior_node1 <- XML::newXMLNode('prior', 
					attrs=c(id='geoclockPrior.c:location', name='distribution', 
					x='@traitClockRate.c:location'), 
					parent=distribution_node2, doc=bxml);
					
					gamma_node1 <- XML::newXMLNode('Gamma', 
						attrs=c(id='Gamma.8', name='distr'), 
						parent=prior_node1, doc=bxml);
						
						dummy <- XML::newXMLNode('parameter', 
							attrs=c(id='RealParameter.43', 
								spec='parameter.RealParameter', 
								estimate='false', name='alpha'), 
							'0.001', parent=gamma_node1, doc=bxml);
							
						dummy <- XML::newXMLNode('parameter', 
							attrs=c(id='RealParameter.44', 
								spec='parameter.RealParameter', 
								estimate='false', name='beta'), 
							'1000.0', parent=gamma_node1, doc=bxml);
							
							
				prior_node2 <- XML::newXMLNode('prior', 
					attrs=c(id='relativeGeoRatesPrior.s:location', 
						name='distribution', 
						x='@relativeGeoRates.s:location'), 
					parent=distribution_node2, doc=bxml);
					
					gamma_node2 <- XML::newXMLNode('Gamma', 
						attrs=c(id='Gamma.9', name='distr'), 
						parent=prior_node2, doc=bxml);
						
						dummy <- XML::newXMLNode('parameter', 
							attrs=c(id='RealParameter.45', 
								spec='parameter.RealParameter', 
								estimate='false', name='alpha'), 
							'1.0', parent=gamma_node2, doc=bxml);
							
						dummy <- XML::newXMLNode('parameter', 
							attrs=c(id='RealParameter.46', 
								spec='parameter.RealParameter', 
								estimate='false', name='beta'), 
							'1.0', parent=gamma_node2, doc=bxml);
							
							
				
			prior_node3 <- XML::newXMLNode('prior', 
					attrs=c(id='nonZeroRatePrior.s:location', 
						name='distribution'), 
					parent=distribution_node2, doc=bxml);
					
					x_node <- XML::newXMLNode('x', 
						attrs=c(id='nonZeroRates.s:location', 
							spec='util.Sum'), 
						parent=prior_node3, doc=bxml);
						
						dummy <- XML::newXMLNode('arg', 
							attrs=c(idref='rateIndicator.s:location'),
								parent=x_node, doc=bxml);
							
					distr_node <- XML::newXMLNode('distr', 
							attrs=c(id='Poisson.2', 
								spec='beast.math.distributions.Poisson', 
								offset='4.0'), 
							parent=prior_node3, doc=bxml)
							
							dummy <- XML::newXMLNode('parameter', 
								attrs=c(id='RealParameter.47', 
									spec='parameter.RealParameter', 
									estimate='false', name='lambda'),
								'0.693', parent=distr_node, doc=bxml);
								
			distribution_node3 <- XML::newXMLNode('distribution', attrs=c(id='likelihood', spec='util.CompoundDistribution', useThreads='true'), parent=distribution_node1, doc=bxml);
			
				distribution_node4 <- XML::newXMLNode('distribution', attrs=c(id='traitedtreeLikelihood.location', spec='AncestralStateTreeLikelihood', tag='location', tree= paste('@Tree.t:',Name,sep='')), parent=distribution_node3, doc=bxml);
					
					data_node2 <- XML::newXMLNode('data', attrs=c(id='location', spec='AlignmentFromTrait'), parent=distribution_node4, doc=bxml);
					
						dummy <- XML::newXMLNode('traitSet', attrs=c(id='traitSet.location', spec='beast.evolution.tree.TraitSet', taxa=paste('@TaxonSet.', Name, sep=''), traitname='discrete'), paste(paste(df$taxon_df,df$location,sep='='), collapse=','), parent=data_node2, doc=bxml);
						
						dummy <- XML::newXMLNode('userDataType', attrs=c(id='traitDataType.location', spec='beast.evolution.datatype.UserDataType', codeMap=paste(paste(lapply(sort(unique(df$location, incomparables='false')), function(x) paste(x,which(x==sort(unique(df$location, incomparables='false')))[[1]]-1, sep='=')), collapse=','), paste('? = ', paste(seq(0,length(sort(unique(df$location, incomparables='false')))-1,1), collapse=' '), sep=''), sep=','), codelength='-1', states=paste(length(sort(unique(df$location, incomparables='false'))))), parent=data_node2, doc=bxml);
						
						
						siteModel_node <- XML::newXMLNode('siteModel', attrs=c(id='geoSiteModel.s:location', spec='SiteModel', gammaCategoryCount='1'), parent=distribution_node4, doc=bxml);
						
							dummy <- XML::newXMLNode('parameter', attrs=c(id='mutationRate.s:location', spec='parameter.RealParameter', estimate='false', name='mutationRate'),'1.0', parent=siteModel_node, doc=bxml);
							 
							dummy <- XML::newXMLNode('parameter', attrs=c(id='gammaShape.s:location', spec='parameter.RealParameter', estimate='false', name='shape'),'1.0', parent=siteModel_node, doc=bxml);
							
							dummy <- XML::newXMLNode('parameter', attrs=c(id='proportionInvariant.s:location', spec='parameter.RealParameter', estimate='false', lower='0.0', name='proportionInvariant', upper='1.0'),'0.0', parent=siteModel_node, doc=bxml);
							
							substModel_node <- XML::newXMLNode('substModel', attrs=c(id='svs.s:location', spec='SVSGeneralSubstitutionModel', rateIndicator='@rateIndicator.s:location', rates='@relativeGeoRates.s:location'), parent=siteModel_node, doc=bxml);
							
								freq_node <- XML::newXMLNode('frequencies', attrs=c(id='traitfreqs.s:location', spec='Frequencies'), parent=substModel_node, doc=bxml);
									dummy <- XML::newXMLNode('parameter', attrs=c(id='traitfrequencies.s:location', spec='parameter.RealParameter', dimension=paste(length(sort(unique(df$location, incomparables='false')))), name='frequencies'), paste(1/length(unique(df$location, incomparables='false'))), parent=freq_node, doc=bxml);
									
									
						dummy <- XML::newXMLNode('branchRateModel', attrs=c(id='StrictClockModel.c:location', spec='beast.evolution.branchratemodel.StrictClockModel', clock.rate='@traitClockRate.c:location'), parent=distribution_node4, doc=bxml);
						
			dummy <- XML::newXMLNode('operator', attrs=c(id='georateScaler.s:location', spec='ScaleOperator', parameter='@relativeGeoRates.s:location', scaleAllIndependently='true', scaleFactor='0.99', weight='30.0'), parent=run_node, doc=bxml);
			
			dummy <- XML::newXMLNode('operator', attrs=c(id='indicatorFlip.s:location', spec='BitFlipOperator', parameter='@rateIndicator.s:location', weight='30.0'), parent=run_node, doc=bxml);
			
			dummy <- XML::newXMLNode('operator', attrs=c(id='geoMuScaler.c:location', spec='ScaleOperator', parameter='@traitClockRate.c:location', scaleFactor='0.9', weight='3.0'), parent=run_node, doc=bxml);
			
			dummy <- XML::newXMLNode('operator', attrs=c(id='BSSVSoperator.c:location', spec='BitFlipBSSVSOperator', indicator='@rateIndicator.s:location', mu='@traitClockRate.c:location', weight='30.0'), parent=run_node, doc=bxml);
			
			
			
			logger_node1 <- XML::newXMLNode('logger', attrs=c(id='tracelog', spec='Logger', fileName=paste(Name,'log', sep='.'), logEvery=paste(tracelog_logevery), model='@posterior', sanitiseHeaders='true', sort='smart'), parent=run_node, doc=bxml);
			
				dummy <- XML::newXMLNode('log', attrs=c(idref='posterior'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='likelihood'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='prior'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='rateIndicator.s:location'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='relativeGeoRates.s:location'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='traitClockRate.c:location'), parent=logger_node1, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(id='geoSubstModelLogger.s:location', spec='SVSGeneralSubstitutionModelLogger', dataType='@traitDataType.location', model='@svs.s:location'), parent=logger_node1, doc=bxml);
						
				
				
			logger_node2 <- XML::newXMLNode('logger', attrs=c(id='screenlog', spec='Logger', logEvery=paste(screenlog_logevery), sanitiseHeaders='true'), parent=run_node, doc=bxml);
			
				dummy <- XML::newXMLNode('log', attrs=c(idref='posterior'), parent=logger_node2, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='likelihood'), parent=logger_node2, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(idref='prior'), parent=logger_node2, doc=bxml);
				
				
			logger_node3 <- XML::newXMLNode('logger', attrs=c(id=paste('treelog.t',Name, sep=':'), spec='Logger', fileName='$(tree).trees', logEvery=paste(treelog_logevery), mode='tree', sanitiseHeaders='true'), parent=run_node, doc=bxml);
				
				dummy <- XML::newXMLNode('log', attrs=c(id=paste('TreeWithTraitLogger.t', Name, sep=':'), spec='beast.evolution.tree.TreeWithMetaDataLogger', tree=paste('@Tree.t',Name, sep=':')), parent=logger_node3, doc=bxml);
				
			
			logger_node4 <- XML::newXMLNode('logger', attrs=c(id='treeWithTraitLogger.location', spec='Logger', fileName='location_tree_with_trait.trees', logEvery=paste(treeWithTraitlog_logevery), mode='tree'), parent=run_node, doc=bxml);
			
				log_node4 <- XML::newXMLNode('log', attrs=c(id='TreeWithTraitLogger.0', spec='beast.evolution.tree.TreeWithTraitLogger', tree=paste('@Tree.t', Name, sep=':')), parent=logger_node4, doc=bxml);
				
					dummy <- XML::newXMLNode('metadata', attrs=c(idref='posterior'), parent=log_node4, doc=bxml);
					
					dummy <- XML::newXMLNode('metadata', attrs=c(idref='traitedtreeLikelihood.location'), parent=log_node4, doc=bxml);

################################################################################

#save the final bxml
XML::saveXML(bxml, save_path, encoding="UTF-8", standalone='no')
}

################################################################################


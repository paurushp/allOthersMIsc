
# Noisy-OR model for knowledge integration \citation{Axiomating Noisy-OR}

Noisy.OR = function(data, testGenes){
	N.OR = matrix(1, nrow = length(tsetGenes), ncol = length(testGenes))
	for(i in 1:ncol(data)){
		mat = 1-(matrix(data[,i], nrow = length(tsetGenes), ncol = length(testGenes)))
		N.OR=N.OR*(mat)
	}
	N.ORprob=1-N.OR
	return(N.ORprob)
}
Indep.prior = function(data, testGenes){
	IPr = matrix(1, nrow = length(tsetGenes), ncol = length(testGenes))
	for(i in 1:ncol(data)){
		mat = matrix(data[,i], nrow = length(tsetGenes), ncol = length(testGenes))
		IPr=IPr*(mat)
	}
	return(IPr)
}
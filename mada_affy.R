library(affy)
prostrate.data=ReadAffy() # Read the data (all the cel files in the working directory)
eset.mas5 = mas5(prostrate.data) # normalize with MAS5
exprSet.nologs = exprs(eset.mas5) # raw expression data
colnames(exprSet.nologs) # List the column (chip) names
exprSet = log(exprSet.nologs, 2) # log transformation
write.table(exprSet, file="Prostrate_cancer_mas5_matrix.csv", quote=F, sep="\t") # Exposrt the table as a csv file
data.mas5calls = mas5calls(prostrate.data)# Run the Affy A/P call algorithm on the CEL files we processed above, Performs the Wilcoxon signed rank-based gene expression presence/absence detection algorithm
data.mas5calls.calls = exprs(data.mas5calls)# Get the actual A/P calls
write.table(data.mas5calls.calls, file="Prostrate_cancer_mas5calls.csv", quote=F, sep="\t") # Print the calls as a matrix


##############################
# Normalization for raw data #
##############################
# with trimmed mean normalization

exprSetRaw = read.delim("file") # Read an expression matrix
trmean = apply(exprSetRaw, 2, mean, trim=0.02) # Calculate trimmed mean for all columns (for rows change 2 to 1, trim value gives the fraction to be trimmed)
sd = apply(exprSetRaw, 2, sd) # variation between the gene on each chip
median = apply(exprSetRaw, 2, median) # compare the medians for each chip
mean.of.trmeans = mean(trmean)
exprSet.trmean = (exprSetRaw/trmean)*mean.of.trmeans
write.table("exprSet.trmean", file="data_norm.csv", quote=F, sep="\t")
exprSet = exprSet.trmean

# quantile normalization using limma

library(limma)
exprSet.quantile = normalizeQuantiles(exprSet)
exprSet = exprSet.quantile

###################################
# Differentially expressed genes  #LV
###################################
exprSet=read.delim("file.csv", head=T, sep="\t") #read file if needed
colnames(exprSet) # see the sample names (chip names)
dataset1=exprSet[1, c(colnames(exprSet)[1:2])] # can use the column numbers instead with exprSet[1,c(1,2)]
dataset2=exprSet[1, c(colnames(exprSet)[3:4])]
t.test.genes=t.test(dataset1,dataset2, "two.sided")
t.test.genes$p.value
ttestAll=apply(exprSet, 1, function(x) {t.test(x[1:2], x[3:4]) $p.value})

cand.net.mat = as(cand.net, "matrix")
	diag(cand.net.mat)=0
	candNets[[i]]=c(cand.net.mat)
	geneSet = unlist(dimnames(cand.net.mat)[1])
	names = list()
	for(j in 1:length(geneSet)){
		name = unlist(strsplit(geneSet[j], ":"))
		names[j] = name[2]
	}
	for(i in 1:ncol(cand.net.mat)){
		rownames(cand.net.mat)[i]=names[[i]]
		colnames(cand.net.mat)[i]=names[[i]]
	}
	Original = cand.net.mat
	genes = unlist(names)
	GO = goSimMatrix(genes, "mean", "Lin") # GO similarity of Genes
	DOMAIN = ProtDomSim(featuresA, genes) # interpro domain similarity
	DOMAIN.INT = getDomainInt(domD, genes) # number of interacting domain pairs across the proteins
	noiseKegg2 = hideInfo(cand.net.mat, genes)
	noiseKegg2 = addInfo(noiseKegg2, genes)
	noiseGraph = graph.adjacency(noiseKegg, mode="directed", weighted = NULL, diag = TRUE, add.colnames = NULL, add.rownames = NA)
	PATHWAY = s.path(noiseGraph, genes) # Shortest path length from the KEGG
	DBASE = subgraphs(gr, genes) # Shortest path length from PathwayCommons
	PCOM=s.path(G.sel, genes)
	DBASE2=s.path(gr, genes)

res=createLFM(datas)
	
	nw.prior = PKMCMC(Sample=500000, Burnin=500000, datas, res, genes, timelag=50, threshold=0.001 , interval=1,tria=1)


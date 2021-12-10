###### change gene format ######
###### ID~~gene symbol #########


Clean_of_many_tags <- function(x){
	###########
	pattern = " (1 of many)."
	patterns = paste(pattern,1:100,sep='')
	patterns = c(" (1 of many)",patterns)
	class(x)
	for(p in patterns){
		x = gsub(p,'',x,fixed=T)
	}
	return(x)
}


Process_gene_features <- function(gene_features){
	ID = gene_features$V1
	symbol = gene_features$V2
	#### remove some tags #######
	symbol_cl = Clean_of_many_tags(symbol)
	####
	out = paste(ID,symbol,sep='~~')
	return(out)
}

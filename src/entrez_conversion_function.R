
exper<-read.delim("disease_experiments_filtered.tsv",header=F)
know<-read.delim("disease_knowledge_filtered.tsv",header=F)
text<-read.delim("disease_textmining_filtered.tsv",header=F)

get.entrez.ids.disease<-function(textmining,knowledge,experiment) {
 	###libraries and prepwork###
	require(DBI)
 	require(org.Hs.eg.db)
 	require(goseq)
 	require(foreach)
 	require(doMC)
 	registerDoMC(cores=1)
	SYMBOLS<-unlist(as.list(org.Hs.egSYMBOL))
	sy.eid<-names(SYMBOLS)
	names(sy.eid)<-SYMBOLS
	xx <- org.Hs.egENSEMBLPROT
	mapped_genes <- mappedkeys(xx)
	xx <- as.list(xx[mapped_genes])
	ep.eid<-goseq:::reversemapping(xx)
	###following sections all work the same, basically seperates the ids##
	###into either ensemble ids or gene symbols, then converts to entrez##
	###ids and creates a new table with pertinent info##
	#####handles textmining table###
	text<-textmining
	ensp<-text[grep("ENSP",text[,1]),]
	ensp<-if(length(which(is.na(names(ep.eid[as.character(ensp[,1])]))))>0) {
	ensp<-ensp[-(which(is.na(ep.eid[as.character(ensp[,1])]))),]} else {ensp}
	other<-text[-(grep("ENSP",text[,1])),]
	other<-other[-(grep("rRNA",other[,1])),]
	other<-other[-(grep("\\.",other[,1])),]
	other<-other[-(which(is.na(sy.eid[as.character(other[,2])]))),]
	if(dim(other)[1]==0) {
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		text.fin<-data.frame(as.character(b),y)
		} else{ if(dim(ensp)[1]==0) {
		text.fin<-data.frame(sy.eid[as.character(other[,2])],other)} else{
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		ensp<-data.frame(as.character(b),y)
		other<-data.frame(sy.eid[as.character(other[,2])],other)
		text.fin<-rbind(other,ensp)}}
	text.fin<-data.frame(text.fin[,1],text.fin[,3:4],text.fin[,6], 
		rep("Literature",dim(text.fin)[1]), text.fin[,5], 
		rep("Textmining",dim(text.fin)[1]))
	colnames(text.fin)<-1:7
	
	###handles knowledge table####
	know<-knowledge
	ensp<-know[grep("ENSP",know[,1]),]
	ensp<-if(length(which(is.na(names(ep.eid[as.character(ensp[,1])]))))>0) {
	ensp<-ensp[-(which(is.na(names(ep.eid[as.character(ensp[,1])])))),]} else 			{ensp}
	other<-know[-(grep("ENSP",know[,1])),]
	other<-other[-(grep("rRNA",other[,1])),]
	other<-other[-(grep("\\.",other[,1])),]
	other<-other[-(which(is.na(sy.eid[as.character(other[,2])]))),]
	if(dim(other)[1]==0) {
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		know.fin<-data.frame(as.character(b),y)
		} else{if(dim(ensp)[1]==0) {
		know.fin<-data.frame(sy.eid[as.character(other[,2])],other)} else{
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		ensp<-data.frame(as.character(b),y)
		other<-data.frame(sy.eid[as.character(other[,2])],other)
		know.fin<-rbind(other,ensp)}}
	know.fin<-data.frame(know.fin[,1],know.fin[,3:4],know.fin[,8], 
		know.fin[,6:7],rep("Knowledge",dim(know.fin)[1]))
	colnames(know.fin)<-1:7
	
	####handles experiment table#####
	exper<-experiment
	ensp<-exper[grep("ENSP",exper[,1]),]
	ensp<-if(length(which(is.na(ep.eid[as.character(ensp[,1])])))>0) {
	ensp<-ensp[-(which(is.na(ep.eid[as.character(ensp[,1])]))),]} else {ensp}
	other<-exper[-(grep("ENSP",exper[,1])),]
	other<-other[-(grep("rRNA",other[,1])),]
	other<-other[-(grep("\\.",other[,1])),]
	other<-other[-(which(is.na(sy.eid[as.character(other[,2])]))),]
	if(dim(other)[1]==0) {
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		exper.fin<-data.frame(as.character(b),y)} else{if(dim(ensp)[1]==0) {
		exper.fin<-data.frame(sy.eid[as.character(other[,2])],other)} else{
		x<-ep.eid[as.character(ensp[,1])]
		z<-which(sapply(x,length)>1)
		y<-ensp[-as.vector(z),]
		sub.m<-foreach(i=1:length(z),.combine='rbind') %do% {
			b<-x[[z[i]]]
			m<-ensp[as.vector(z[i]),]
			m<-m[rep(1,time=length(b)),]
			data.frame(as.character(b),m)}
		b<-x[-z]
		ensp<-data.frame(as.character(b),y)
		other<-data.frame(sy.eid[other[,1]],other)
		exper.fin<-rbind(other,ensp)}}
	exper.fin<-data.frame(exper.fin[,1],exper.fin[,3:4],exper.fin[,8], 
		exper.fin[,6:7],rep("Experimental",dim(exper.fin)[1]))
	colnames(exper.fin)<-1:7

	fin<-rbind(text.fin,know.fin,exper.fin)
	colnames(fin)<-c("Entrez_ID","Gene_Symbol","Disease_Ontology_ID",
		"Confidence_Score","Source","Evidence","Type")
	fin
###final table has 7 columns. the 'condfidence score' is from the Jensen database###
###'Source' is which database the info was dervied from. Evidence is a z-score###
###for textmining, knowledge is type, e.g. curated, etc., experimental is either###
###sample number or a p-value. Finally, type is which table the row is derived from##
### e.g. knowledge disease table####
}

####example of running function###

fin<-get.entrez.ids.disease(text,know,exper)







library(tm)
list_file<-dir()

group_mining<-function(x){
  contig_species<-x[5]
  contig_species<-tolower(contig_species)
  contig_species<-strsplit(contig_species,' ')[[1]]
  contig_species<-contig_species[!grepl('>',contig_species)]
  contig_species<-removePunctuation(contig_species)
  contig_species_tokenize<-table(factor(contig_species))
  contig_species_tokenize<-contig_species_tokenize[order(contig_species_tokenize,decreasing = T)]
  return(head(contig_species_tokenize,10))
}

for(i in 1:length(list_file)){
  temp<-read.table(list_file[i],sep='\t',stringsAsFactors = F,quote = "",header=T,fill = T)
  temp_contig<-temp[grepl('human|homo',temp$contig_simplspecies,ignore.case = T)&!grepl('virus',temp$contig_simplspecies,ignore.case = T),1]
  temp_contig<-temp[grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus',temp$contig_simplspecies,ignore.case = T),1]
  temp_contig<-factor(temp_contig)
  list_homo_contigs<-levels(temp_contig)
  temp2<-temp[!(temp$contig_name %in% list_homo_contigs),]
  temp2<-temp2[order(temp2$contig_length,decreasing = T),]
  text_eval<-paste0('data',i,'<-temp2')
  eval(parse(text=text_eval))
  print(i)
}


for(i in 1:length(list_file)){
text_eval<-paste0('dat<-data',i)
eval(parse(text=text_eval))  

temp<-as.list(by(dat,INDICES = dat$contig_name,
                 function(x) paste(names(group_mining(x)),collapse = ' ')))
temp_dat<-data.frame(contig_names=names(temp),contig_species=unname(unlist(temp)))
temp_list<-dat[ which(dat$contig_name %in% temp_dat$contig_names) & !duplicated(dat$contig_name),c(1,2)]
temp_dat$contig_length<-temp_list[ match(temp_dat$contig_names,temp_list$contig_name),2]
temp_dat<-temp_dat[order(temp_dat$contig_length,decreasing = T),c(1,3,2)]
text_eval<-paste0('clean_data',i,'<-temp_dat')
eval(parse(text=text_eval)) 
print(i)
}



rm(list=ls()[!grepl('data',ls())])




#name_f<-unname(sapply(unname(sapply(data1$full_name,function(x) strsplit(x,' ')[[1]][1])),
#                      function(x) substr(x,2,nchar(x))[1]))
#name_g<-unname(sapply(data1$contig_simplspecies,function(x) strsplit(x,' ')[[1]][1]))
#sum(name_f!=name_g)


blast_result_categorize<-function(data_path,output_path){

###Make_data_frame
Main_data<-data.frame()
Current_data<-data.frame()

###Index preparation
file_result<-file(data_path)
blast_result<-readLines(file_result)
close(file_result)
index_Cstart<-grepl('^Query=',blast_result)
index_Length<-grepl('Length=',blast_result)
index_Glist<-grepl('(?!\\d+\\s+) \\d+\\s+(\\d+\\.\\d+|\\d+e-\\d+)\\s*$', blast_result,perl = T)

###Subgroup data for 1st loop
index_Step1<-which(index_Cstart|index_Length|index_Glist)
state_Length<-'Unknown'
state_list<-'Unrecord'

print('Step 1')
###Progress bar
pb <- txtProgressBar(min = 0, max = length(index_Step1), style = 3)

### Step1. Generalized Information
for(i in 1:length(index_Step1)){
  ### Progress report
  setTxtProgressBar(pb, i)
  
  ### Current line mark
  lmark_decider<-c(New_start=index_Cstart[index_Step1[i]],
                   Length=index_Length[index_Step1[i]],
                   Align_list=index_Glist[index_Step1[i]])
  lmark<-names(which(lmark_decider))
  #print(lmark)
  
  switch(lmark,
         New_start={
           ### Arrange previous data
           if(state_list=='Record'){
             Current_data<-data.frame(contig_name,contig_length,contig_score,contig_simplspecies)}
           ### Joint prvious data
           Main_data<-rbind(Main_data,Current_data)
           Current_data<-data.frame()
           contig_name<-strsplit(blast_result[index_Step1[i]],'[\\ =]')[[1]][3]
           ### State define
           state_Length<-'Unknown'
           state_list<-'Unrecord'
         },
         Length={
           if(state_Length=='Unknown'){
             contig_length<-strsplit(blast_result[index_Step1[i]],'[\\ =]')[[1]][2]
             state_Length<-'Known'}
         },
         Align_list={
           switch(state_list,
                  Unrecord={
                    vector_list<-strsplit(blast_result[index_Step1[i]],'[ ]+')[[1]]
                    contig_score<-as.numeric(vector_list[length(vector_list)-1])
                    contig_simplspecies<-paste(vector_list[1:(length(vector_list)-2)],collapse = ' ')
                    state_list<-'Record'
                  },
                  Record={
                    vector_list<-strsplit(blast_result[index_Step1[i]],'[ ]+')[[1]]
                    current_score<-as.numeric(vector_list[length(vector_list)-1])
                    if(current_score==contig_score){
                      contig_simplspecies<-c(contig_simplspecies,
                                             paste(vector_list[1:(length(vector_list)-2)],collapse = ' '))}
                  })
           
         })
}


### Step 2. Filled up information

Main_data$full_name<-NA
Main_data$transcript_length<-NA
list_contigs<-unname(sapply(blast_result[index_Cstart],
                            function(x) strsplit(x,'[\\ =]')[[1]][3]))
table_contig_index<-data.frame(row.names = list_contigs,index=which(index_Cstart))
index_Detail<-grepl('^>',blast_result)
index_Step2<-which(index_Cstart|index_Length|index_Detail)

###Progress bar
pb <- txtProgressBar(min = 0, max = dim(Main_data)[1], style = 3)
print('Step 2')
current_contig<-'NONE'
for(i in 1:dim(Main_data)[1]){
  setTxtProgressBar(pb, i)
  if(Main_data[i,1]!=current_contig){
    current_contig<-as.character(Main_data[i,1])
    current_index<-which(index_Step2==table_contig_index[current_contig,])
    Main_data[i,'full_name']<-paste(blast_result[index_Step2[current_index+2]:
                                                (index_Step2[current_index+3]-1)],collapse =' ')
    Main_data[i,'transcript_length']<-strsplit(blast_result[index_Step2[current_index+3]],'[\\ =]')[[1]][2]    
    current_index<-current_index+4
  }else{
    Main_data[i,'full_name']<-paste(blast_result[index_Step2[current_index]: 
                                                   (index_Step2[current_index+1]-1)],collapse = ' ')
    Main_data[i,'transcript_length']<-strsplit(blast_result[index_Step2[current_index+1]],'[\\ =]')[[1]][2]    
    current_index<-current_index+2
  }
}



write.table(Main_data,output_path,quote = F,col.names = T,row.names = F,sep = '\t')

}

### Loop start
output_path_header<-'./Blast_result/'
if (!file.exists(output_path_header)){
  dir.create(output_path_header)
}

file_all<-list.files('./',recursive=T)
file_blastresult<-grep('nr_result.out$',file_all,value=T)

for(i in 1:length(file_blastresult)){
  temp_id<-strsplit(file_blastresult[i],'/')[[1]]
  data_id<-temp_id[grep('trinity_out',temp_id)-1]
  print(paste('Start working on',data_id))
  blast_result_categorize(file_blastresult[i],paste0(output_path_header,data_id,'_out.txt'))
}


##\\d+\\s+(\\d+\\.\\d+|\\d+e-?\\d+)\\s*$

#RSEM_original<-read.table('data/RSEM.genes.results',header=T)
#Main_data$name<-as.character(Main_data$name)
#Main_data$gene_id<-sapply(Main_data$name,function(x) strsplit(x,split = '_seq')[[1]][1])
#RSEM_original$species_name<-Main_data$species[ match(RSEM_original$gene_id,Main_data$gene_id) ]
#RSEM_original<-RSEM_original[order(RSEM_original$TPM,decreasing = T),]

#Main_data$len<-as.numeric(Main_data$len)
#Main_data<-Main_data[order(Main_data$len,decreasing = T),]
#temp<-Main_data[!(grepl('homo',Main_data$species) | grepl('human',Main_data$species)),]


#temp$True_positive<-grepl('Pseudo',temp$species)
#accuracy<-(sum(temp$True_positive[temp$len>583])+sum(!temp$True_positive[temp$len<583]))/length(temp$len)
#ggplot(temp)+geom_dotplot(aes(x=len,colour=True_positive))+ggtitle('300 Virus reads')+geom_vline(x=583)+xlab('Contigs length')

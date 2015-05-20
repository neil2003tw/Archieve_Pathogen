library(tm)


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
  return(Main_data)
}



blast_result_filter<-function(temp_raw){
  temp_contig<-temp_raw[(grepl('human|homo',temp_raw$contig_simplspecies,ignore.case = T)&!grepl('virus',temp_raw$contig_simplspecies,ignore.case = T))
                        |grepl('pan troglodytes|gorilla|Pongo abelii|pan paniscus',temp_raw$contig_simplspecies,ignore.case = T),1]
  temp_contig<-factor(temp_contig)
  list_homo_contigs<-levels(temp_contig)
  temp2<-temp_raw[!(temp_raw$contig_name %in% list_homo_contigs),]
  temp_filtered<-temp2[order(as.numeric(as.character(temp2$contig_length)),decreasing = T),]
  temp_filtered$contig_name<-factor(as.character(temp_filtered$contig_name))
  return(temp_filtered)
}



blast_result_duplicate_remover<-function(temp_filtered){
  temp<-as.list(by(temp_filtered,INDICES = temp_filtered$contig_name,
                   function(x) paste(names(group_mining(x)),collapse = ' ')))
  temp_dat<-data.frame(contig_names=names(temp),contig_species=unname(unlist(temp)),stringsAsFactors = F)
  temp_dat$contig_length<-as.numeric(as.character(temp_filtered[ match(temp_dat$contig_names,
                                                                       temp_filtered$contig_name),2]))
  temp_dat<-temp_dat[order(temp_dat$contig_length,decreasing = T),c(1,3,2)]
  return(temp_dat)
}

Trinity_subgrouper<-function(Trinity_path,lookup_table,counts,reads_header,current_counts){
  con<-file(Trinity_path)
  Trinity_file<-readLines(con)
  close(con)
  contig_index<-grep('^>',Trinity_file)
  contig_name<-unname(sapply(Trinity_file[contig_index], function(x) strsplit(x,'>|\\s',perl = T)[[1]][2]))
  contig_list<-data.frame(contig_name,contig_index,stringsAsFactors = F)
  lookup_table_inuse<-head(lookup_table,counts)
  new_fasta_file<-c()
  for(i in 1:counts){
    temp_index<-which(contig_list$contig_name %in% lookup_table_inuse$contig_names[i])
    contig_sub_name<-strsplit(lookup_table_inuse$contig_names[i],'_')[[1]][1]
    new_fasta_file<-c(new_fasta_file,paste0('>comp',current_counts,'_c0_seq1'))
    new_fasta_file<-c(new_fasta_file,
                      Trinity_file[ (contig_list$contig_index[temp_index]+1)
                                    :(contig_list$contig_index[temp_index+1]-1)])
    current_counts<<-current_counts+1
  }
  sub_annotate<<-data.frame(new_id=paste0('comp',(current_counts-counts):current_counts,'_c0_seq1'),
                            sample_id=reads_id,
                            species=head(lookup_table$contig_species,counts))
  return(new_fasta_file)
}


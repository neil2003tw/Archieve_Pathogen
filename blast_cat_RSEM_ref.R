source('/Users/NeilWu/Google Drive/Research_Pathogen/code/tool.R')
HDfas_path<-'/Users/NeilWu/Google Drive/Research_Pathogen/data_source/House_keeping.fasta'

output_path_header<-'./Blast_result/'
if (!file.exists(output_path_header)){
  dir.create(output_path_header)
}

file_all<-list.files('./',recursive=T)
file_blastresult<-grep('nr_result.out$',file_all,value=T)
file_Trinity<-grep('Trinity.fasta$',file_all,value=T)
new_fasta_file<-c()
new_annotate_file<-data.frame()

for(i in 1:length(file_blastresult)){
  temp_id<-strsplit(file_blastresult[i],'/')[[1]]
  data_id<-temp_id[grep('trinity_out',temp_id)-1]
  print(paste('Start working on',data_id))
  temp_raw<-blast_result_categorize(file_blastresult[i],paste0(output_path_header,data_id,'_out.txt'))
  temp_filtered<-blast_result_filter(temp_raw)
  lookup_table<-blast_result_duplicate_remover(temp_filtered) 
  sub_fasta<-Trinity_subgrouper(file_Trinity[i],lookup_table,10,data_id)
  new_fasta_file<-c(new_fasta_file,sub_fasta)
  new_annotate_file<-rbind(new_annotate_file,sub_annotate)
}

con<-file(HDfas_path)
HKfasta_file<-readLines(con)
close(con)

new_fasta_file<-c(new_fasta_file,HKfasta_file)
write.table(new_fasta_file,'RSEM_ref.fasta',quote = F,row.names = F,col.names = F)
write.table(new_annotate_file,'RSEM_ref_table',quote = F,row.names = F,col.names = F)



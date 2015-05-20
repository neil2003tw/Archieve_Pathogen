## Havn't add gene situation, archieved
#fasta_generator<-function(reference,readlength=75,readnumber,seeds=1){
#  set.seed(seeds)
#  genome<-as.vector(read.table(reference,sep='')[1,1])
#  reads_file<-c(rep('',readnumber*2))
#  reads_start<-as.integer(runif(readnumber,1,nchar(genome)-readlength))
#  for(i in 1:readnumber){
#    reads_file[2*i-1]<-paste0('>',strsplit(reference,split = '\\.')[[1]][1],'_startfrom_',reads_start[i],'_',i)
#    reads_file[2*i]<-substr(genome,reads_start[i],reads_start[i]+74)
#  }
#  write.table(reads_file,file = paste0(strsplit(reference,split = '\\.')[[1]][1],'_reads.fasta'),quote = F,row.names=F,col.names = F)
#  return(reads_file)
#}

reference_generator<-function(num,length,setseeds=T,seeds=1){
  if(setseeds){set.seed(seeds)}
  all<-c()
  for(i in 1:num){
  temp<-factor(as.integer(runif(length,1,5)))
  levels(temp)<-c('A','T','C','G')
  all<-c(all,Reduce(function(x,y) paste0(x,y),temp))
  }
  return(all)
}

fasta_generator_pseudo<-function(reference,readlength=75,readnumber,name){
  genome<-reference
  reads_file<-c(rep('',readnumber*2))
  gene_length<-runif(1, min = 1000, max = 2000)
  print(paste0('GLe',gene_length))
  gene_locus<-as.integer(runif(1, min=1, max=nchar(genome)-gene_length))
  print(paste0('GLo',gene_locus))
  reads_start<-as.integer(runif(readnumber, min = gene_locus, max = gene_locus+gene_length-readlength))
  for(i in 1:readnumber){
    reads_file[2*i-1]<-paste0('>',name,'_startfrom_',reads_start[i],'_',i)
    reads_file[2*i]<-substr(genome,reads_start[i],reads_start[i]+readlength-1)
  }
  write.table(reads_file,file = paste0(name,'_reads.fasta'),quote = F,row.names=F,col.names = F)
  return(reads_file)
}


fake_ref<-c()
for(i in 1:20){
  fake_ref[i]<-reference_generator(1,rnorm(1,6000,1000),setseeds = F)
}

### Make ref to reads pipe
for(i in 1:20){
  temp<-c()
  temp[1]<-paste0('>Pseudo_',i,'_length_',nchar(fake_ref[i]))
  temp[2]<-fake_ref[i]
  t=2
  for(s in seq(1,nchar(fake_ref[i]),by = 80)){
    temp[t]<-substr(fake_ref[i],s,s+79)
    t=t+1
  }
  write.table(temp,paste0('Pseudo_',i,'.fasta'),quote = F,row.names=F,col.names = F)
}

###
for(i in 1:20){
  fasta_generator_pseudo(fake_ref[i],readlength = 150,readnumber = 30000,name = refs[i])
}


### Change older fasta to greater fasta
refs<-dir('code/Pseudo_reference/')
refs<-refs[!grepl('reads',refs)]
dir.create('code/Pseudo_reference/corrected/')
for(s in 1:length(refs)){
P1<-read.table(paste0('code/Pseudo_reference/',refs[s]),stringsAsFactors = F)
print(paste0('Working on ',refs[s]))
P1_v2<-c()
P1_v2[1]<-P1[1,1]
t=2
for(i in seq(1,nchar(P1[2,1]),by = 80)){
  P1_v2[t]<-substr(P1[2,1],i,i+79)
  t=t+1
}
write.table(P1_v2,paste0('code/Pseudo_reference/corrected/',refs[s]),quote = F,row.names=F,col.names = F)
}


### Reload reference to lists
refs<-dir('data_source/Pseudo_reference/corrected_ref/')
refs<-refs[!grepl('reads',refs)]
fake_ref<-c()
for(i in 1:length(refs)){
  temp_expr<-parse(text=paste0('read.table(\'data_source/Pseudo_reference/corrected_ref/',
                               refs[i],'\',stringsAsFactors = F)'))
  temp<-eval(temp_expr)
  temp_seq<-paste0(temp[-1,1],collapse = '')
  fake_ref[i]<-temp_seq
}

### See gene locus length
reads<-dir('.')
reads<-reads[grepl('reads',reads)]
for(i in 1:length(refs)){
  temp_expr<-parse(text=paste0('read.table(\'./',reads[i],'\',stringsAsFactors = F)'))
  temp<-eval(temp_expr)
  temp<-temp[grepl('>Pseudo_',temp[,1]),]
  temp_len<-as.integer(unname(sapply(temp,function(x) strsplit(x,'_')[[1]][4])))
  print(max(temp_len)-min(temp_len))
}

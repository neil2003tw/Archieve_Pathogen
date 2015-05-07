library(tm)

###word construct
useless_words<-c('clone','sequence','dna','component','noncoding','long','rna','tpa','gpa','bac','protein','predicted','synthetic','construct','mitochondrion','transcrip','cdna','isolate')

###Begin
result_nr<-readLines(file('data/nr_clean_result_20v0.5kr.out'))
query<-grepl(']$',result_nr)|grepl('^Query',result_nr)|grepl('[1-9]:[1-9]',result_nr)
state=0
score=0
count=0
organized_data<-data.frame()
current_data<-data.frame()
for(i in 1:length(query)){
  ###Progress
  print(paste0((i*100)/length(query),'%'))
  
  if(xor(state,query[i])){
    if(!state){
      
      ###Save text mining result
      if(i != 1){
        text_tokenize<-table(factor(text_tokenize))
        if(max(text_tokenize)==min(text_tokenize)){
          current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize))
        }else if(length(text_tokenize)==0){
          current_data$species<-NA
        }else if(max(text_tokenize)==2){
          current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize[text_tokenize>1]))
        }else{
          current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize[text_tokenize>2]))
        }
      }
      ### Bind data
      organized_data<-rbind(organized_data,current_data)
      ### Create new data
      current_data<-data.frame(name=unname(strsplit(result_nr[i],'[\\ =]')[[1]][3]))
      current_data$len<-unname(strsplit(result_nr[i],'[\\ =]')[[1]][5])
    }else if(state){
      text_tokenize<-c()
      score=0
      count=0
    }
    state=!state
  }else if(!state){
    ##if it is blast mapping line
    ## Older for proto blast format'[a-z]+\\|[A-Z]+'
    if(grepl('\\d+\\s+(\\d+\\.\\d+|\\d+e-?\\d+)\\s*$', result_nr[i],perl = T) & count<4){
      score_grid<-strsplit(result_nr[i],'[ ]+')[[1]]
      ### if its score is highest and is numeric
      if( as.numeric(score_grid[length(score_grid)-1])>=score){ #| !is.numeric(score_grid[length(score_grid)-1])){
        print(paste(score,' ',i))
      ## if it belongs to older blast format
        if(grepl('\\|',result_nr[i])){
          a<-strsplit(result_nr[i],'\\| ')[[1]][2]
          a<-stripWhitespace(a)
          a<-removeNumbers(a)
          a<-removePunctuation(a)
          a<-tolower(a)
          a<-removeWords(a,stopwords('english'))
          a<-removeWords(a,useless_words)
          a<-strsplit(a,split = ' ')[[1]]
          a<-a[nchar(a)>2]
          text_tokenize<-c(text_tokenize,a)
          count=count+1
          score=as.numeric(score_grid[length(score_grid)-1])
        }else{
          ## Rseudo Reference
          a<-strsplit(result_nr[i],'  ')[[1]][2]
          text_tokenize<-c(text_tokenize,a)
          count=count+1
        }
      }
    }
  }
}


text_tokenize<-table(factor(text_tokenize))
if(max(text_tokenize)==min(text_tokenize)){
  current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize))
}else if(length(text_tokenize)==0){
  current_data$species<-NA
}else if(max(text_tokenize)==2){
  current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize[text_tokenize>1]))
}else{
  current_data$species<-Reduce(function(x,y) paste(x,y,sep = ' '),names(text_tokenize[text_tokenize>2]))
}
organized_data<-rbind(organized_data,current_data)

#rm(list = ls()[ls()!='organized_data'])

RSEM_original<-read.table('data/RSEM.genes.results',header=T)
organized_data$name<-as.character(organized_data$name)
organized_data$gene_id<-sapply(organized_data$name,function(x) strsplit(x,split = '_seq')[[1]][1])
RSEM_original$species_name<-organized_data$species[ match(RSEM_original$gene_id,organized_data$gene_id) ]
RSEM_original<-RSEM_original[order(RSEM_original$TPM,decreasing = T),]

organized_data$len<-as.numeric(organized_data$len)
organized_data<-organized_data[order(organized_data$len,decreasing = T),]
temp<-organized_data[!(grepl('homo',organized_data$species) | grepl('human',organized_data$species)),]


temp$True_positive<-grepl('Pseudo',temp$species)
accuracy<-(sum(temp$True_positive[temp$len>583])+sum(!temp$True_positive[temp$len<583]))/length(temp$len)
ggplot(temp)+geom_dotplot(aes(x=len,colour=True_positive))+ggtitle('300 Virus reads')+geom_vline(x=583)+xlab('Contigs length')

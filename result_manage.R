temp<-read.csv('../result_2_out.csv',stringsAsFactors = F,fill=T,header=T)
temp<-read.table('../result_Hep_out.txt',stringsAsFactors = F,fill=T,sep='\t')
temp<-temp[,-5]

temp_ordered<-temp[order(temp[,1],decreasing = T),]
split_grep<-function(x){
  line<-strsplit(x,split=' ')
  line<-line[[1]]
  final<-grep('^[A-Z][a-z]+',line,value = T)
  if(any(c('Homo','Human','Sapiens') %in% final)){
    final<-'Homo'
  }
  if(any(c('Pongo','Pan') %in% final)){
    final<-'Homo'
  }
  
  
  return(final)
}
trial<-sapply(temp_ordered[,4],split_grep,USE.NAMES = F)
temp_ordered[,4]<-sapply(trial,function(x) x[1])
temp_ordered[,4]<-factor(temp_ordered[,4])
levels(temp_ordered[,4])[21]<-'Others'
temp_sub<-temp_ordered[temp_ordered[,4]!='Homo',]


ggplot(temp_sub,aes(x=V1,fill=V4))+geom_dotplot(stackgroups = T,method = 'histodot',binwidth = 50)+guides(fill=guide_legend(title='Race'))+xlab('Contig length')
ggplot(temp_ordered,aes(x=V1,fill=V4))+geom_dotplot(stackgroups = T,method = 'histodot',binwidth = 20)+guides(fill=guide_legend(title='Race'))+xlab('Contig length')

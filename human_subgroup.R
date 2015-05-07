set.seed(1234)
reads_index_p1<-sample(1:(462442128/2),30000000)*2
reads_index_p2<-reads_index_p1-1
reads_index<-sort(c(reads_index_p1,reads_index_p2))
print('index done!')

read1<-read.table('flux_homo_R1.fasta')
read1_sub<-read1[reads_index,]
write.table(read1_sub,'flux_homo_sub_R1.fasta',quote=F,row.names=F,col.names = F)
rm(read1)

read2<-read.table('flux_homo_R2.fasta')
read2_sub<-read2[reads_index,]
write.table(read2_sub,'flux_homo_sub_R2.fasta',quote=F,row.names=F,col.names = F)
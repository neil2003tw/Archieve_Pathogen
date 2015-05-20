from sys import argv, exit
from Bio.Blast import NCBIXML
import os.path
import re

def main():
	if len(argv) != 3:
		print("Usage: {0} dataset".format(argv[0]))
	
	dataset = argv[1]

	if not os.path.exists(dataset):
		print('dataset {0} not fount'.format(dataset))
		exit(1)
	
	
	result=open(dataset)
	blast_records=NCBIXML.parse(result)

	total=len(list(blast_records))+1	
	
	query_len=[0]*total
	query_score=[0]*total
	query_e=[0]*total
	query_hit=[0]*total	
	
	result.seek(0)
	blast_records=NCBIXML.parse(result)
	i=0
	for blast_record in blast_records:
		if len(blast_record.descriptions) is not 0:
			query_len[i]=blast_record.alignments[0].hsps[0].identities
			query_score[i]=blast_record.alignments[0].hsps[0].score
			query_e[i]=blast_record.alignments[0].hsps[0].expect
			query_hit[i]=blast_record.alignments[0].title
		i+=1	

	f = open(argv[2],'w')

	for i in range(total):
		line_str=str(query_len[i])+'\t'+str(query_score[i])+'\t'+str(query_e[i])+'\t'+str(query_hit[i])+'\n'
		f.write(line_str)


if __name__=="__main__":
	exit(main())

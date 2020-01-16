#!/usr/bin/python3

#make a dictionary of keys: common gene name, collection: exons keys, values: corresponding locus info

import sys
#sys.path.append("/Users/fss/codeProjects/pipeline/")

from collections import defaultdict

#input: read in the conserved gtf
filePath = sys.argv[1]
fileOut = sys.argv[2]

#initialize dictionary of dictionaries
d=defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

#define chromosome list
chrom_list = ['chr'+str(i) for i in range(1,23)] #modify this for mouse (should be 20 instead of 23)
chrom_list.extend(['chrX','chrY'])
print(chrom_list)

gene_order=[]

#read input gtf and get rid of anything that isn't autosome or X/Y; get rid of NR noncoding transcripts, get rid of poorly annotated LOC transcripts, make sure the feature is an exon (since we're interested in introns)
count=0
for line in open(filePath,'r').readlines():

	count+=1
	if count%10000:
		print(count)

	line = line.rstrip().split('\t')
	chrom = line[0]
	source = line[1]
	feature = line[2]
	start = int(line[3])
	end = int(line[4])
	strand = line[6]
	info = line[-1]

	gene_id = line[-1].split(' ')[1][1:-2]
	NM_id = line[-1].split(' ')[-3][1:-2]

	if gene_id not in gene_order:
		gene_order.append(gene_id)

	duple = [start,end]

	#initialize things that don't change from transcript to transcript per gene
	d[gene_id]['chrom'] = chrom
	d[gene_id]['strand'] = strand
	d[gene_id]['NM_id'][NM_id].extend(duple)

#sort the lists of start/end locations by tss_id if not already sorted (due to nature of extend and the format of gtfs, these lists should already be sorted, but better to be safe)
for gene_id in d:
	longest = 0
	for NM_id in d[gene_id]['NM_id']:
		d[gene_id]['NM_id'][NM_id] = sorted(d[gene_id]['NM_id'][NM_id],key=int)
		length = d[gene_id]['NM_id'][NM_id][-1]-d[gene_id]['NM_id'][NM_id][0]
		if length > longest:
			longest = length
			d[gene_id]['NM_longest'] = {NM_id:[]}
		d[gene_id]['NM_id'][NM_id] = list(zip(*[iter(d[gene_id]['NM_id'][NM_id][1:])]*2))

print(d)
print(len(gene_order))

f = open(fileOut,'w')

for gene_id in gene_order:
	x=0
	NM_longest = list(d[gene_id]['NM_longest'].keys())[0]
	for duple in d[gene_id]['NM_id'][NM_longest]:
		x = ['TOSS' if duple not in d[gene_id]['NM_id'][key] else 'KEEP' for key in d[gene_id]['NM_id'].keys()]
		if 'TOSS' not in x:
			chrom = d[gene_id]['chrom']
			start = str(duple[0])
			end = str(duple[1])
			strand = d[gene_id]['strand']
			string = 'gene_id %s; gene_name %s; transcript_id %s;' % (gene_id,gene_id,NM_longest)
			f.write('%s\tjarey\tintron\t%s\t%s\t.\t%s\t.\t%s\n' % (chrom,start,end,strand,string))

f.close()

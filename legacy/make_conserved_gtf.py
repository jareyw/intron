#!usr/bin/python

#make a dictionary of keys: common gene name, collection: exons keys, values: corresponding locus info

import sys
sys.path.append("/Users/fss/codeProjects/pipeline/")

from collections import defaultdict

#input: read in the full gtf (AceView has ~3e6 lines, RefSeq has ~1e6 lines)
filePath = sys.argv[1]
fileOut = sys.argv[2]

#line by line, create dictionary of all exons by gene name
d=defaultdict(lambda: defaultdict(list))
chrom_list = ['chr'+str(i) for i in range(1,23)] #modify this for mouse (should be 20 instead of 23)
chrom_list.extend(['chrX','chrY'])
print chrom_list

#read input gtf and get rid of anything that isn't autosome or X/Y; get rid of NR noncoding transcripts, get rid of poorly annotated LOC transcripts, make sure the feature is an exon (since we're interested in introns)
for line in open(filePath,'r').readlines():
	line = line.split('\t')
	gene_id = line[-1].split(' ')[1][1:-2]
	transcript_id = line[-1].split(' ')[-3][1:-2]
	chrom = line[0]
	feature = line[2]
	start = int(line[3])
	end = int(line[4])
	strand = line[6]
	duple = [start,end]

	if chrom not in chrom_list:
		continue
	if transcript_id.startswith('NR'):
		continue
	if gene_id.startswith('LOC') or feature != "exon":
		continue

	key = gene_id
	d[key]['duple'].extend(duple)
	d[key]['chrom'] = line[0]
	d[key]['strand'] = strand

count=0
#for each gene name, obtain the broadest coordinates as a duple (list)
for gene in d:
	mincoord = min(d[gene]['duple'])
	maxcoord = max(d[gene]['duple'])
	if maxcoord - mincoord > 1000000:
		count+=1
		print [gene,mincoord,maxcoord]
		d.pop(gene,None)
	d[gene]['locus'] = [mincoord,maxcoord]
print count

'''
#for each gene name, compare the coordinates to ask: does any other gene have coordinates that overlap?

#if so, discard both genes from the dictionary [note: there could be a stacking issue with deletion, since if you iteratively delete while the loop is running you could not be catching overlap with deleted items] best to store all bad genes in a list and then dict.pop() all of them at once

#don't forget about chromosomes!
'''
count=0
badList = []
for gene in d:
	x,y = d[gene]['locus']

	count+=1
	if count%1000==0:
		print count

	for gene2 in d:
		if d[gene]['chrom'] != d[gene2]['chrom']:
			continue
		# if d[gene]['strand'] != d[gene2]['strand']:
		# 	continue
		if gene == gene2:
			continue

		x2,y2 = d[gene2]['locus']

		if x2 < y and y2 > x:
			if gene not in badList:
				badList.append(gene)
			if gene2 not in badList:
				badList.append(gene2)
			break

print badList
print len(badList)
print len(d)

for gene in badList:
	d.pop(gene,None)

print len(d)

f = open(fileOut,'w')

for line in open(filePath,'r').readlines():
	newline = line.split('\t')
	gene_id = newline[-1].split(' ')[1][1:-2]
	transcript_id = newline[-1].split(' ')[-3][1:-2]
	if gene_id in d and newline[2] == "exon" and transcript_id.startswith('NM'):
		f.write(line)

f.close()


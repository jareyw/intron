#!/usr/bin/python

import pysam
import sys
import re
import getopt
import multiprocessing

#current iteration assumes .bam inputs, reads in gtf file, takes in output location

#gtffile should be limited to only introns

try:
	opts, args = getopt.getopt(sys.argv[1:],'h:b:g:o:',['bamfolderpaths=','gtffile=','output='])
except getopt.GetoptError as err:
	print (err)
	sys.exit(2)
	
for opt, arg in opts:
	if opt == '-h':
		sys.exit()
	elif opt in ('-b','--bamfolderpaths'):
		id_list = []
		for line in open(arg,'r').readlines():
			line = line.split('\t')
			id_list.append(line)
	elif opt in ('-g','--gtffile'):
		intron_list=[]
		for line in open(arg,'r').readlines():
			intron_list.append(line[0:-1])
		gtf = arg
	elif opt in ('-o','--output'):
		output = arg

def intronCounter(file,short_id,intron_list,output,count):
	fileout = open(output + '/%s.gtf' % (short_id.rstrip()),'w')
	#write initial row of intron ids

	count+=1
	print count

	print file
	print pysam.idxstats(file)
	
	bamfile = pysam.AlignmentFile(file,'rb')

	#calculate mappable reads for tcga (star) bams
	#mapped_reads = reduce(lambda x, y: x + y, [int(l.split('\t')[2]) for l in pysam.idxstats(file)])

	#calculate mappable reads for hisat generated bams
	#mapped_reads = reduce(lambda x, y: x + y, [int(l.split('\t')[2]) for l in pysam.idxstats(file).rstrip('\n').split('\n')])
	#mapped_frags = float(mapped_reads/2) #this is over-conservative, since some frags will only have 1 mapped read (not 2)

	mapped_reads = reduce(lambda x, y: x + y, [int(l.rstrip('\n').split('\t')[2]) for l in pysam.idxstats(file)])
	print mapped_reads

	subcount=0
	for newline in intron_list:
		line = newline.split('\t')

		subcount+=1
		if subcount %10000 == 0:
			print subcount

		chrom = line[0]
		start = int(line[3])
		end = int(line[4])
		strand = line[6]

		fileout.write(newline + '\t')

		if end - start <= 10:
			fileout.write('0' + '\n')
			continue

		intron_count=0
		span_count=0
		for read in bamfile.fetch(chrom,start+5,end-5):
			if not read.is_duplicate and not read.is_qcfail and not read.is_secondary and read.is_proper_pair and read.get_tag('XS') == strand: #include XS tag for hisat2 stranded aligned bams (and read.get_tag('XS') == strand)
				span_count+=1
				if 'N' not in read.cigarstring:
					intron_count+=1

		#calculate rpkm per genomic locus for every intron
		intron_length = end - start
		reads_per_mil = float(mapped_reads/1000000)
		rpm = float(intron_count)/reads_per_mil
		rpk = float(intron_length)/1000
		rpkm = float(rpm)/rpk

		fileout.write(str(intron_count)+'\t'+str(span_count)+'\t'+str(rpkm)+'\n')
	fileout.close()

if __name__ == '__main__':
	jobs=[]
	for file, short_id in id_list:
		count=0
		p = multiprocessing.Process(target=intronCounter,args=(file,short_id,intron_list,output,count))
		jobs.append(p)
		p.start()
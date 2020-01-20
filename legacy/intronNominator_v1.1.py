#!/usr/bin/python3

import sys
import re
import glob
import heapq
import numpy
import getopt
import os

#for filtering
from collections import Counter

'''
Arguments
'''

try:
	opts, args = getopt.getopt(sys.argv[1:],'h:g:m:l:d:r:f:o:',['gtffolderpath=','metric=','readlength=','datatype=','ranktype=','filter','output='])
except getopt.GetoptError as err:
	print (err)
	print("'gtffolderpath=','metric=','readlength','datatype=','ranktype=','filter','output='")
	sys.exit(2)

print("'gtffolderpath=','metric=','readlength','datatype=','ranktype=','filter','output='")
cons_dict = {} #initialize to empty dictionary, which is default

for opt, arg in opts:
	if opt == '-h':
		sys.exit()
		print("Documentation not complete. See administrator.")
	elif opt in ('-g','--gtffolderpath'):
		folderpath = arg
	elif opt in ('-m','--metric'):
		metric = arg #median or mean
	elif opt in ('-l','--readlength'):
		readlength = int(arg) #need to figure this out empirically by using samtools view
	elif opt in ('-d','--datatype'):
		datatype = str(arg)
	elif opt in ('-r','--ranktype'):
		ranktype = str(arg)
	elif opt in ('-f','--filter'):
		spkm_min = int(arg)
	elif opt in ('-o','--output'):
		output = arg

print(output)
fileout = open('%s.txt' % (output),'w')

'''
Functions
'''
#overall, this version of the code 160619 is just for TCGA data
def file_is_empty(path):
    return os.stat(path).st_size==0

#fast counting solution to determine most common element in list
def Most_Common(lst):
	data = Counter(lst)
	return data.most_common(1)[0][0]

'''
Main Code
'''
d={}
i = 0
for file in glob.glob("%s/*.gtf" % (folderpath)):
	print(file)
	count=0
	subcount=0

	#get millions of counts retroactively (there may be a way to implement this more intelligently)
	for line in open(file,'r').readlines():
		line = line.split('\t')
		intron_count = int(line[-3])
		start=int(line[3])
		end=int(line[4])
		length = end-start
		rpkm = float(line[-1].rstrip())
		try:
			mr = intron_count*1000/float(length*rpkm)
			if mr:
				print(mr)
				break
		except:
			continue

	for line in open(file,'r').readlines():
		line = line.split('\t')
		info = ','.join(map(str,line[:-3]))
		start=int(line[3])
		end=int(line[4])
		sample_id = file.split('/')[-1][:-4] #modify depending on what you want your sample names to be

		intron_id = line[-4].split(' ')[-1][:-1]+'.'+str(start) #custom intron id

		length = end-start
		added_length = length+2*(readlength-1)

		# subcount+=1
		# if subcount %10000 == 0:
		# 	print subcount
		# 	break

		try:
			intron_count = int(line[-3])
			sum_count = int(line[-2])
			rpkm = float(line[-1].rstrip())
			rpm = float(intron_count)/mr
		except:
			intron_count=0
			sum_count=0
			rpm=0
			rpkm=0

		span_count = sum_count-intron_count

		#can filter at this point since the results you're getting are wacky if the constitutive option is used
		if cons_dict and intron_id not in cons_dict:
			continue

		if span_count != 0:
			spkm = span_count*1000.0/float(2*(readlength-1)*mr)
		else:
			spkm = 0.0

		try:
			ratio = float(rpkm/spkm)
		except:
			ratio = 0.0

		if intron_id not in d:
			d[intron_id] = {'span_count':[span_count],'rpkm':[rpkm], 'rpm':[rpm], 'spkm':[spkm], 'ratio':[ratio], 'info':info, 'intron_id':intron_id, 'sample_id':[sample_id]}
		else:
			d[intron_id]['span_count'].append(span_count)
			d[intron_id]['rpkm'].append(rpkm)
			d[intron_id]['rpm'].append(rpm)
			d[intron_id]['spkm'].append(spkm)
			d[intron_id]['ratio'].append(ratio)
			d[intron_id]['sample_id'].append(sample_id)

'''
Computing Intron Metrics, Sorting, and Filtering
'''

intron_list=[] #compute metrics for introns to be ranked by
for intron in d:
	d[intron]['%s_mean' % ranktype] = sum((d[intron][ranktype]))/float(len(d[intron][ranktype]))
	d[intron]['%s_median' % ranktype] = numpy.median((d[intron][ranktype]))

	#should filter by ratio
	d[intron]['ratio_mean'] = sum((d[intron]['ratio']))/float(len(d[intron]['ratio']))

	#get rid of the following lines before appending to the intron_list
	#filter for minimum spkm > 20 (not sure how much this will impact the output) - will have to test run
	if any(i<spkm_min for i in d[intron]['span_count']):
		print(d[intron]['sample_id'])
		print(d[intron]['span_count'])
		continue

	#filter for reasonable ratio of true intron reads/spanning reads
	if d[intron]['ratio_mean'] > 1: #weakness of applying this filter is that many of the "IR" features it cuts are cassette exons
		continue

	intron_list.append(d[intron])

print(len(intron_list))

#option for metric used (so far, sum and median)
#sorts by whichever metric you decide on
high_values = heapq.nlargest(len(intron_list),intron_list,key=lambda s: s['%s_%s' % (ranktype, metric)])

'''
Writing Table
'''
#write the header row
sample_id_list = high_values[0]['sample_id']
fileout.writelines('intron_id'+'\t'+'info'+'\t')
sample_id_string = '\t'.join(sample_id_list)
fileout.write(sample_id_string)
fileout.write('\n')
print(sample_id_list)

for item in high_values:
	fileout.write(item['intron_id']+'\t')
	fileout.write(item['info']+'\t')
	data_string = '\t'.join(map(str,item[datatype])) #here you decide what to output ranked by rpkm
	fileout.write(data_string)
	fileout.write('\n')

fileout.close()






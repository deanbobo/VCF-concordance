#!/usr/bin/env python2

import argparse
import itertools
from argparse import RawTextHelpFormatter
import vcf
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='This script checks two VCFs against each other. One should be a truth dataset, the other is being iterrogated. The sanger data can be converted with fasta-to-vcf.py', formatter_class=RawTextHelpFormatter)
parser.add_argument('-ngs', metavar='<ngs>' , help='bgzipped & tabix indexed VCF file (this file is being interrogated)', required=True)
parser.add_argument('-truth', metavar='<truth>', help='bgzipped & tabix indexed VCF file (this is the truth dataset)', required=True)
parser.add_argument('-sample', metavar='<sample>', help='sample name (must match the sample name in both VCFs)', required=True)
parser.add_argument('-bed', metavar='<bed>', help='BED file (regions to be interrogated)', required=True)
parser.add_argument('-strict', action='store_true', help='penalizes FNs even if there is no coverage at site')
parser.add_argument('-verbose', action='store_true', help='prints more information (verbose mode)')
parser.add_argument('-FNreport', action='store_true', help='write a report of FNs which is used by FN-coverage-alleleFreq-lookup.py.')
parser.add_argument('-FPreport', action='store_true', help='write a report of FPs')
args = parser.parse_args()
truth_file = args.truth
ngs_file = args.ngs
sample = args.sample
bedfile = args.bed

#for each variant in the truth dataset, lookup the variant in the NGS dataset
#if exists, good. If it doesn't -> false negative

#set counters to zero
fn_counter = 0
truth_calls = 0
ngs_calls = 0
fn_rate=None


print "\n\nAnalysis for sample: ",sample


if args.FNreport:
	
	if args.strict:
		FNfile = open(ngs_file+'.'+sample+'.strict.FNs','w')
	else:
		FNfile = open(ngs_file+'.'+sample+'.lenient.FNs','w')
	FNfile.write("#Sample\tChrom\tPosition\tFN_count\tCoverage\n")	
	
	
if args.FPreport:
	if args.strict:
		FPfile = open(ngs_file+'.'+sample+'.strict.FPs','w')
		
	else:
		FPfile = open(ngs_file+'.'+sample+'.lenient.FPs','w')
	FPfile.write("#Sample\tChrom\tPosition\tFP_count\tCoverage\n")

def FNout(sample, chrom, position, FNcount, coverage):
	if coverage == '' or coverage is None:
		coverage='NA'
	FNfile.write(sample + "\t" + str(chrom) + "\t" + str(position) + "\t" + str(FNcount) + "\t" + str(coverage) + "\n")
	return
	
def FPout(sample, chrom, position, FPcount, coverage):
	if coverage == '' or coverage is None:
		coverage='NA'
	FPfile.write(sample + "\t" + str(chrom) + "\t" + str(position) + "\t" + str(FPcount) + "\t" + str(coverage) + "\n")
	return
	
	
	

#########################################
####### FALSE NEGATIVE CALCULATOR #######
#########################################


#open the truth VCF tabix indexed gzipped
truth_reader = vcf.Reader(open(truth_file, 'r'))

#open the NGS dataset tabix index gzipped
ngs_reader = vcf.Reader(open(ngs_file, 'r'))

if args.verbose:
	print "\nFalse Negative Analysis:"


for region in BedTool(bedfile):
	chrom = region[0]
	start = int(region[1])
	stop = int(region[2])
	truth_iterator=None
	truth_test_iterator=None

 	try:
		truth_iterator = truth_reader.fetch(chrom, start, stop)
		truth_test_iterator, truth_iterator = itertools.tee(truth_iterator)
	except:
		if args.verbose:
			print "No truth variants found within bed region:",str(chrom),":",str(start),"-",str(stop)
		pass
	
	#if iterator is not empty
	if truth_test_iterator and any(truth_test_iterator) is True:
	
		#go through truth iterator
		for truth_variant in truth_iterator:

			#get truth variants state:   0 = ref/ref  1=alt/ref  2=alt/alt
			truth_state = None
			truth_chrom=truth_variant.CHROM
			truth_pos=truth_variant.POS
			truth_genotype=truth_variant.genotype(sample)['GT']
			
			try:
				truth_state = truth_variant.genotype(sample).gt_type # 0=hom ref; 1=het; 2=hom alt
			except:
				if args.verbose:
					print "cannot get truth variant state at position: ",truth_chrom, truth_pos
				pass
			#add to truth counter (use number from truth_state)
			if truth_state:
				truth_calls = truth_calls + truth_state  # add number of calls to truth calls
				
				#print "print0 truth: ",truth_chrom, truth_pos, truth_state


				ngs_iterator=None
				ngs_test_iterator=None
				ngs_depth=None
	
				try:
					ngs_iterator = ngs_reader.fetch(str(truth_chrom),truth_pos, truth_pos)
					ngs_test_iterator, ngs_iterator = itertools.tee(ngs_iterator)
					#print "there is an ngs_iterator" 
				except:					
					#print "print entered except because ngs_iterator cannot be defined"
					# if truth_state > 0:
# 						
# 						fn_counter = fn_counter + truth_state
# 						print "fn_counter: ",fn_counter
# 						if args.FNreport:
# 							FNout(sample,truth_chrom,truth_pos,truth_state,ngs_depth)						
# 						if args.verbose:
# 							print "FN: truth:", truth_chrom, truth_pos, truth_genotype, "\t NGS: NOT CALLED"
					pass

				if ngs_test_iterator and any(ngs_test_iterator) is True:

					for ngs_variant in ngs_iterator:
					
						ngs_chrom = None
						ngs_pos = None
						ngs_genotype=None
						ngs_state=None
						ngs_filter=None
						ngs_depth=None
				
						try:
							ngs_chrom = ngs_variant.CHROM
							ngs_pos = ngs_variant.POS
						except:
							if args.verbose:
								print "Cannot find NGS variant at position: ", truth_chrom, truth_pos
							pass
				
						try:
							ngs_state = ngs_variant.genotype(sample).gt_type
							ngs_genotype = ngs_variant.genotype(sample)['GT']
						except:
							if args.verbose:
								print "Cannot get NGS state or genotype at position: ", truth_chrom, truth_pos
							pass
					
						try: 
							ngs_depth = ngs_variant.genotype(sample)['DP']
						except:
							if args.verbose:
								print "Cannot get NGS depth at position: ", truth_chrom, truth_pos
							pass
							
						try:
							ngs_filter = ngs_variant.FILTER[0]
						except:
							pass
						#print "print1 ngs:",ngs_chrom, ngs_pos, ngs_state

						if ngs_state is None:
							fn_counter = fn_counter + truth_state
							if args.verbose:
								print "FN: truth:", truth_chrom, truth_pos, truth_genotype, "\t NGS: NOT CALLED"
							if args.FNreport:
								FNout(sample,truth_chrom,truth_pos,truth_state,ngs_depth)


						elif args.strict and ngs_filter == "LowQual" and truth_state > ngs_state:
							fn_counter = fn_counter + truth_state  #penalize the entire truth state
							if args.verbose:
								print "FN:  truth:",truth_chrom, truth_pos, truth_genotype, "\t NGS: ", ngs_chrom, ngs_pos, ngs_genotype, "LowQual"
							if args.FNreport:
								FNout(sample,truth_chrom,truth_pos,truth_state-ngs_state,ngs_depth)

						
						elif truth_state > ngs_state:
							fn_counter = fn_counter + (truth_state - ngs_state) #penalize any missing variants
							if args.verbose:
								print "FN:  truth:",truth_chrom, truth_pos, truth_genotype, "\t NGS: ", ngs_chrom, ngs_pos, ngs_genotype
							if args.FNreport:
								FNout(sample,truth_chrom,truth_pos,truth_state-ngs_state,ngs_depth)
						#print "fn_counter: ",fn_counter
							
				else:
					if truth_state > 0:
						fn_counter = fn_counter + truth_state
						if args.FNreport:
							FNout(sample,truth_chrom,truth_pos,truth_state,ngs_depth)
							#printFN(1,1,1,1)
						if args.verbose:
							print "FN: truth:", truth_chrom, truth_pos, truth_genotype, "\t NGS: NOT CALLED" #ngs_chrom, ngs_pos, ngs_genotype, "NOT CALLED"
	 

if not truth_calls:
	if args.verbose:
		print "Cannot calculate FN rate. No truth variants found!"
else:
	fn_rate = float(fn_counter) / float(truth_calls) * 100
	#print "FN rate: ",fn_rate, "fn counter: ",fn_counter,"truth calls: ",truth_calls



#########################################
####### FALSE POSITIVE CALCULATOR #######
#########################################

#set counters
fp_counter = 0
fp_rate=None

if args.verbose:
	print "\nFalse Positive Analysis:"

for region in BedTool(bedfile):
	#print "chrom: ",chrom, ":", start, "-",stop
	chrom = region[0]
	start = int(region[1])
	stop = int(region[2])
	ngs_iterator2=None
	ngs_test_iterator2=None
 	try:
		ngs_iterator2 = ngs_reader.fetch(chrom, start, stop)
		ngs_iterator2, ngs_test_iterator2 = itertools.tee(ngs_iterator2)
	except:
		if args.verbose:
			print "No NGS variants found within bed region:",str(chrom),":",str(start),"-",str(stop)
		pass
	
	#go through ngs variants in iterator
	if ngs_test_iterator2 and any(ngs_test_iterator2) is True:


		ngs_variant = None
		
		for ngs_variant in ngs_iterator2:
			
			ngs_state = None
			ngs_genotype = None
			ngs_chrom = None
			ngs_pos = None
			ngs_chrom = None
			ngs_filter = None
			truth_chrom = None
			truth_pos = None
			truth_genotype = None
			truth_state = None
			ngs_depth=None
			
			ngs_chrom = ngs_variant.CHROM
			ngs_pos = ngs_variant.POS 	

			try:
				ngs_state = ngs_variant.genotype(sample).gt_type	#get ngs variants state:   0 = ref/ref  1=alt/ref  2=alt/alt	
				ngs_genotype = ngs_variant.genotype(sample)['GT']
			except:
				if args.verbose:
					print "Cannot get NGS state and/or genotpye at position: ", ngs_chrom, ngs_pos
				pass
				
			try:
				ngs_filter = ngs_variant.FILTER[0]
			except:
				pass 

			# add number of calls to truth calls
			if ngs_state:
				ngs_calls = ngs_calls + ngs_state  

				truth_iterator2=None
				truth_test_iterator2=None

				try:
					truth_iterator2 = truth_reader.fetch(str(ngs_chrom),ngs_pos, ngs_pos)
					truth_iterator2, truth_test_iterator2 = itertools.tee(truth_iterator2)
				except:

					# if ngs_state > 0 and args.strict and ngs_filter != 'LowQual':
# 						#fp_counter = fp_counter + ngs_state
# 						if args.verbose:
# 							print "FP: NGS: ",ngs_variant.CHROM, ngs_variant.POS, ngs_variant.genotype(sample)['GT'], "\t Truth: NOT CALLED\tSkipped because of LowQual filter"
# 				
# 					elif ngs_state > 0:
# 						fp_counter = fp_counter + ngs_state
# 						if args.verbose:
# 							print "FP: NGS: ",ngs_variant.CHROM, ngs_variant.POS, ngs_variant.genotype(sample)['GT'], "\t Truth: NOT CALLED"
					pass

				if truth_test_iterator2 and any(truth_test_iterator2) is True:
					
					truth_variant=None

					for truth_variant in truth_iterator2:
				
						truth_chrom = None
						truth_pos = None
						truth_genotype = None
						truth_state = None
						
						try:
							truth_chrom = truth_variant.CHROM
							truth_pos = truth_variant.POS
						except:
							if args.verbose:
								"Cannot find truth variant at position: ", ngs_chrom, ngs_pos
							pass
									
						try:
							truth_genotype = truth_variant.genotype(sample)['GT']
							truth_state = truth_variant.genotype(sample).gt_type
						except:
							if args.verbose:
								"Cannot get genotype or state at position: ", ngs_chrom, ngs_pos
							pass
					

						if truth_state is None:
							
							if args.strict and ngs_filter == 'LowQual' and ngs_state > 0:
								if args.verbose:
									print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype, "\t Truth: NOT CALLED\tskipped because LowQual filter in NGS"
							elif ngs_state > truth_state:
								fp_counter = fp_counter + ngs_state
								if args.FPreport:
									FPout(sample,ngs_chrom,ngs_pos,ngs_state,ngs_depth)
								if args.verbose:
									print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype, "\t Truth: NOT CALLED"

						
						
						else:
							if args.strict and ngs_filter == 'LowQual' and ngs_state > truth_state:
								if args.verbose:
									print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype, "\t Truth: ", truth_chrom, truth_pos, truth_genotype, "\tskipped because LowQual filter in NGS"	
							elif ngs_state > truth_state:
								fp_counter = fp_counter + (ngs_state - truth_state)
								if args.FPreport:
									FPout(sample,ngs_chrom,ngs_pos,ngs_state - truth_state,ngs_depth)
								if args.verbose:
									print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype,  "\t Truth: ", truth_chrom, truth_pos, truth_genotype
				
				else:
					if args.strict and ngs_filter == 'LowQual' and ngs_state > 0:
						if args.verbose:
							print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype, "\t Truth: NOT CALLED\tskipped because LowQual filter in NGS"	
					elif ngs_state > 0:
						fp_counter = fp_counter + ngs_state
						if args.FPreport:
							FPout(sample,ngs_chrom,ngs_pos,ngs_state,ngs_depth)
						if args.verbose:
							print "FP:  NGS:",ngs_chrom, ngs_pos, ngs_genotype,  "\t Truth: NOT CALLED" 	
				
						

 	else:
		continue



if not ngs_calls:
	if args.verbose:
		print "Cannot calculate FP rate. No NGS variants found!" 
else:
	fp_rate = float(fp_counter) / float(ngs_calls) * 100



print "\nSUMMARY:"
print "FN counter:\t", fn_counter
print "FP counter:\t", fp_counter
print "Expected Calls:\t", truth_calls
print "Actual Calls:\t", ngs_calls

if fn_rate is None:
	print "FN Rate:\t0.00 %"
else:
	print "FN Rate:\t", "%.2f" % fn_rate, "%"


if fp_rate is None:
	print "FP Rate:\t0.00 %"
else:
	print "FP Rate:\t", "%.2f" % fp_rate,  "%"







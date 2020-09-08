#!/usr/bin/env python3

import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import pysam
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


#check if truth file exists 
if not os.path.isfile(truth_file):
	sys.exit('The -truth file '+truth_file+' provided  does not exist.\nExiting.')
	
#check if ngs file exists 
if not os.path.isfile(ngs_file):
	sys.exit('The -ngs file '+ngs_file+' provided  does not exist.\nExiting.')

#check if truth file is indexed 
if not (os.path.isfile(truth_file+'.tbi') or os.path.isfile(truth_file+'.csi')):
	sys.exit('The -truth file '+truth_file+' provided is not indexed. No .tbi or .csi file found. Use tabix -p vcf <filename>\nExiting.')

#check if ngs file is indexed
if not (os.path.isfile(ngs_file+'.tbi') or os.path.isfile(ngs_file+'.csi')):
	sys.exit('The -ngs file '+ngs_file+' provided is not indexed. No .tbi or .csi file found. Use tabix -p vcf <filename>\nExiting.')

print("\n\nAnalysis for sample: ",sample)


if args.FNreport:
	FNfile = open(ngs_file+'.'+sample+'.FN.report.txt','w')
	FNfile.write("#Sample\tChrom\tPosition\tFN_count\tCoverage\n")	
	
	
if args.FPreport:
	FPfile = open(ngs_file+'.'+sample+'.FP.report.txt','w')
	FPfile.write("#Sample\tChrom\tPosition\tFP_count\tCoverage\n")

def FNout(sample, chrom, position, FNcount):
	FNfile.write(sample + "\t" + str(chrom) + "\t" + str(position) + "\t" + str(FNcount)  + "\n")
	
	
def FPout(sample, chrom, position, FPcount):
	FPfile.write(sample + "\t" + str(chrom) + "\t" + str(position) + "\t" + str(FPcount) + "\n")
	
	
# def FNchecker(truth_ref, truth_genotypes, ngs_genotype):
# 	# A C G T
# 	truth_alts = { 'A':0,  'C':0,  'G':0,  'T':0 }
# 	ngs_alts = { 'A':0,  'C':0,  'G':0,  'T':0 }
# 	for i in ['A','C','G','T']:
# 		ngs_alts[i] = ngs_genotype.count(i)
# 		truth_alts[i] = truth_genotype.count(i)
# 	
# 	#assign zero to reference alelles
# 	truth_alts[truth_ref]=0
# 	ngs_alts[truth_ref]=0
# 	
# 	#subtract the two dictionarys
# 	fn_dict = {}
# 	for i in ['A','C','G','T']:
# 		fn_dict[i] = truth_alts[i] - ngs_alts[i]
# 	
# 	fn_checker_sum = sum(fn_dict.values())
# 	return(fn_checker_sum)
# 	
# 	
# 	
# def FPchecker(truth_ref, truth_genotype, ngs_genotype):
# 	# A C G T
# 	truth_alts = { 'A':0,  'C':0,  'G':0,  'T':0 }
# 	ngs_alts = { 'A':0,  'C':0,  'G':0,  'T':0 }
# 	for i in ['A','C','G','T']:
# 		ngs_alts[i] = ngs_genotype.count(i)
# 		truth_alts[i] = truth_genotype.count(i)
# 	
# 	#assign zero to reference alelles
# 	truth_alts[truth_ref]=0
# 	ngs_alts[truth_ref]=0
# 	
# 	#subtract the two dictionarys
# 	fp_dict = {}
# 	for i in ['A','C','G','T']:
# 		fp_dict[i] = ngs_alts[i] - truth_alts[i]
# 	
# 	fp_checker_sum = sum(fp_dict.values())
# 	return(fp_checker_sum)
	

#########################################
####### FALSE NEGATIVE CALCULATOR #######
#########################################



#open the truth VCF tabix indexed gzipped
truth_in = pysam.VariantFile(truth_file)

#open the NGS dataset tabix index gzipped
ngs_in = pysam.VariantFile(ngs_file)

if args.verbose:
	print("\nFalse Negative Analysis:")

fn_counter = 0

for region in BedTool(bedfile):
	chrom = region[0]
	start = int(region[1])
	stop = int(region[2])
	
	#go through truth iterator
	for truth in truth_in.fetch(str(chrom), start, stop):

		truth_ref = truth.ref
		truth_chr=truth.contig
		truth_pos=truth.pos
		
		#set these to none. 
		truth_alt=None
		truth_state=None
		truth_genotype=None
		fn_checker_sum=None

		#die if truth is multiallelic
		if len(truth.alts) > 1:
			print("ERROR: The truth site: ", truth_chr, truth_pos, "is a multiallelic site. This will cause errors for this script. The VCFs are expected to be normalized. Use bcftools norm.")
			sys.exit('Truth VCF not normalized. Exiting.')  
		
		#assign alt allele
		truth_alt = str(truth.alts[0])
		
		# skip truth indels 
		if len(truth.ref) > 1 or len(truth.alts[0]) > 1:
			if args.verbose:
				print("WARNING: The truth site: ", truth_chr, truth_pos, "is an indel. Skipping.")
			continue #goes to next truth
		
		#assign truth state
		try:
			truth_state = sum(0 if v is None else v for v in list(truth.samples[sample]['GT'])) #0 for hom ref, 1 for het, 2 for hom alt
		except:
			if args.verbose:
				print("WARNING: Cannot get truth variant state at position: ",truth_chr, truth_pos)
			pass
		
# 		#assign truth genotype
# 		try:
# 			truth_genotype=str(truth.alleles[truth.samples[sample]['GT'][0]])+str(truth.alleles[truth.samples[sample]['GT'][1]])
# 		except:
# 			if args.verbose:
# 				print("WARNING: Cannot get truth variant genotype at position", truth_chr, truth_pos)
# 			pass
			
		#add to truth counter (use number from truth_state)
		if truth_state:
			truth_calls = truth_calls + truth_state  # add number of calls to truth calls

			found_ngs=None

			for ngs in ngs_in.fetch(str(truth_chr), truth_pos-1, truth_pos):
				
				ngs_state=None
				ngs_filter=None
				ngs_depth=None
				fn_checker_sum=None
				
				#get basic site-level information				
				try:
					ngs_chr = ngs.contig
					ngs_pos = ngs.pos
					ngs_ref = ngs.ref
					ngs_alt = str(ngs.alts[0])
				except:
					if args.verbose:
						print("WARNING: Cannot find NGS variant at position: ", truth_chr, truth_pos)
					pass
				
				#die if NGS is multialleic
				if len(ngs.alts) > 1:
					print("ERROR: The NGS site: ", ngs_chr, ngs_pos, "is a multiallelic site. This will cause errors for this script. The VCFs are expected to be normalized. Use bcftools norm.")
					sys.exit('NGS VCF not normalized. Exiting.')  
				
				#find the fetched NGS variant that matches the Sanger variant
				if truth_alt == ngs_alt:
					
					found_ngs=True #need this to leave the NGS searching loop 
					
					#assign NGS state
					try:
						ngs_state = sum(0 if v is None else v for v in list(ngs.samples[sample]['GT'])) #0 for hom ref, 1 for het, 2 for hom alt
					except:
						if args.verbose:
							print("WARNING: Cannot get NGS state at position:", truth_chr, truth_pos)
						pass	
	
					if ngs_state is None or ngs_state == 0:
						fn_counter = fn_counter + truth_state
						if args.verbose:
							print("FN: TRUTH_INFO: truth_chr:", truth_chr, "truth_pos:", truth_pos, "truth_ref:", truth_ref, "truth_alt:", truth_alt, "truth_state:", truth_state, "NGS_INFO: genotype_not_called")
						if args.FNreport:
							FNout(sample,truth_chr,truth_pos,truth_state)
					
					else:  #both variants are present in vcf			
						
						#check that reference alleles match (they should be len=1 because indels were skipped above)
						if truth_ref != ngs_ref: 
							sys.exit('ERROR: The truth and NGS reference alleles do not match. This is not possible. Exiting.') 
						
						fn_checker_sum = truth_state - ngs_state
					
						if fn_checker_sum > 0:
										
							fn_counter = fn_counter + fn_checker_sum
											
							if args.verbose:
								print("FN: TRUTH_INFO: truth_chr:", truth_chr, "truth_pos:", truth_pos, "truth_ref:", truth_ref, "truth_alt:", truth_alt, "truth_state:", truth_state, "NGS_INFO: ngs_chr", ngs_chr, "ngs_pos:", ngs_pos, "ngs_alt:", ngs_alt, "ngs_state:", ngs_state)

							if args.FNreport:
								FNout(sample,truth_chr,truth_pos,fn_checker_sum)
									
				else: #keep searching if not find a matched alt allele. This goes until out of NGS variants in fetch iterable.
					continue

			# if a matching alt allele is not found. 
			if found_ngs is False:
				if args.verbose:
					print("FN: TRUTH_INFO: truth_chr:", truth_chr, "truth_pos:", truth_pos, "truth_ref:", truth_ref, "truth_alt:", truth_alt, "truth_state:", truth_state, "NGS_INFO: genotype_not_called")
				if args.FNreport:
					FNout(sample, truth_chr, truth_pos, truth_state)
				fn_counter = fn_counter + truth_state  #add the ngs state to the FP counter.


if not truth_calls:
	if args.verbose:
		print("Cannot calculate FN rate. No truth variants found!")
else:
	fn_rate = float(fn_counter) / float(truth_calls)
	print("FN_rate: ","{:.2%}".format(fn_rate), "fn_counter: ",fn_counter,"truth_calls: ",truth_calls)



#########################################
####### FALSE POSITIVE CALCULATOR #######
#########################################

#open the ngs VCF tabix indexed gzipped
ngs_in = pysam.VariantFile(ngs_file)

#open the truth dataset tabix index gzipped
truth_in = pysam.VariantFile(truth_file)

if args.verbose:
	print("\nFalse Positive Analysis:")
	
fp_counter = 0
	
for region in BedTool(bedfile):
	chrom = region[0]
	start = int(region[1])
	stop = int(region[2])
	
	#go through ngs iterator
	for ngs in ngs_in.fetch(str(chrom), start, stop):
	
		ngs_ref = ngs.ref
		ngs_chr=ngs.contig
		ngs_pos=ngs.pos
		
		#set these to none. 
		ngs_alt=None
		ngs_state=None
		
		#die if ngs is multiallelic
		if len(ngs.alts) > 1:
			print("ERROR: The NGS site: ", ngs_chr, ngs_pos, "is a multiallelic site. This will cause errors for this script. The VCFs are expected to be normalized. Use bcftools norm.")
			sys.exit('NGS VCF not normalized. Exiting.')  
		
		# skip NGS indels 
		if len(ngs.ref) > 1 or len(ngs.alts[0]) > 1:
			if args.verbose:
				print("WARNING: The ngs site: ", ngs_chr, ngs_pos, "is an indel. Skipping.")
			continue #goes to next ngs
		
		
		#get NGS alt allele
		ngs_alt = str(ngs.alts[0])
		
		
		try:
			ngs_state = sum(0 if v is None else v for v in list(ngs.samples[sample]['GT'])) #0 for hom ref, 1 for het, 2 for hom alt
		except:
			if args.verbose:
				print("WARNING: Cannot get NGS variant state at position: ",ngs_chr, ngs_pos)
			pass
			
		#add to ngs counter (use number from ngs_state)
		if ngs_state:
			ngs_calls = ngs_calls + ngs_state  # add number of calls to ngs calls
			
			truth_depth=None
			found_truth=None

			truth_state=None

			for truth in truth_in.fetch(str(ngs_chr), ngs_pos-1, ngs_pos):
			
				# get truth site-level information
				try:
					truth_chr = truth.contig
					truth_pos = truth.pos
					truth_ref = truth.ref
					truth_alt = str(truth.alts[0])
				except:
					if args.verbose:
						print("WARNING: Cannot find truth variant at position: ", ngs_chr, ngs_pos)
					pass

				#die if truth is multialleic
				if len(truth.alts) > 1:
					print("ERROR: The truth site: ", truth_chr, truth_pos, "is a multiallelic site. This will cause errors for this script. The VCFs are expected to be normalized. Use bcftools norm.")
					sys.exit('Truth VCF not normalized. Exiting.')  


				#find the fetched truth variant that matches the Sanger variant
				if ngs_alt == truth_alt:
				
					found_truth=True #need this to leave the truth searching loop 
					
					# get truth sample state
					try:
						truth_state = sum(0 if v is None else v for v in list(truth.samples[sample]['GT'])) #0 for hom ref, 1 for het, 2 for hom alt
					except:
						if args.verbose:
							print("WARNING: Cannot get truth state at position: ", ngs_chr, ngs_pos)
						pass
					
# 					# get truth sample genotype
# 					try:
# 						truth_genotype=str(truth.alleles[truth.samples[sample]['GT'][0]])+str(truth.alleles[truth.samples[sample]['GT'][1]])
# 					except:
# 						if args.verbose:
# 							print("WARNING: Cannot get truth genotype at position:", ngs_chr, ngs_pos)
# 						pass

					# if the truth variant is not called
					if truth_state is None or truth_state == 0:
						fp_counter = fp_counter + ngs_state
						if args.verbose:
							print("FP: NGS_INFO: ngs_chr:", ngs_chr, "ngs_pos:", ngs_pos, "ngs_ref:", ngs_ref, "ngs_alt:", ngs_alt, "ngs_state:", ngs_state, "TRUTH_INFO: genotype_not_called")
						if args.FPreport:
							FPout(sample, ngs_chr, ngs_pos, ngs_state)
							
					# if a truth variant is called (both VCFs have variant called)
					else:  		
					
						#check that reference alleles match (they should be len=1 because indels were skipped above)
						if ngs_ref != truth_ref: 
							sys.exit('ERROR: The ngs and truth reference alleles do not match. This is not possible. Exiting.') 
							
						#submit genotypes to FP checker
						fp_checker_sum = ngs_state - truth_state
						
						if fp_checker_sum > 0:
						
							fp_counter = fp_counter + fp_checker_sum
							
							if args.verbose:
								print("FP: NGS_INFO: ngs_chr:", ngs_chr, "ngs_pos:", ngs_pos, "ngs_ref: ", ngs_ref, "ngs_alt:", ngs_alt, "ngs_state:" ,ngs_state, "TRUTH_INFO: ", "truth_chr:", truth_chr, "truth_pos:", truth_pos, "truth_state:", truth_state)
							if args.FPreport:
								FPout(sample,ngs_chr,ngs_pos,fp_checker_sum)
								
								
				else: #keep searching if not find a matched alt allele. This goes until out of truth variants in fetch iterable.
					continue
					
			# if a matching alt allele is not found. 
			if found_truth is False:
				if args.verbose:
					print("FP: NGS_INFO: ngs_chr:", ngs_chr, "ngs_pos:", ngs_pos, "ngs_ref:", ngs_ref, "ngs_alt:", ngs_alt, "ngs_state:", ngs_state, "TRUTH_INFO: genotype_not_called")
				if args.FPreport:
					FPout(sample, ngs_chr, ngs_pos, ngs_state)
				fp_counter = fp_counter + ngs_state  #add the ngs state to the FP counter.


if not ngs_calls:
	if args.verbose:
		print("Cannot calculate FP rate. No ngs variants found!")
else:
	fp_rate = float(fp_counter) / float(ngs_calls)
	print("FP_rate: ", "{:.2%}".format(fp_rate), "fp_counter:", fp_counter, "ngs_calls: ",ngs_calls)









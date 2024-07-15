#!/usr/bin/env python3

import sys
import argparse
import bisect
import time
import tracemalloc
import random as rand

def generate_genome(glen: int):
	"""
	Generates "genome" of length glen
	"""
	nucleotides = "ATCG"
	genome = "".join( [nucleotides[rand.randrange(4)] for i in range(glen) ] )
	return genome

def get_reads(rnum: int, rlen: int, genome: str, glen:int):
	"""
	Generates reads from genome sequence
	
	Parameters:
	- rnum: number of reads
	- rlen: read length
	- genome: genome sequence
	- glen: length of genome
	
	Returns:
	- list of reads
	"""
	if rlen > glen:
		sys.stderr.write("Error: Read length must be smaller than genome length.\n")
		exit(1)
	# the largest possible index of the beginning of any read
	last_index = glen - rlen
	reads = []
 
	# generate one read each cycle
	for i in range(rnum):
		rand_index = rand.randrange(last_index+1)
		reads.append(genome[rand_index : (rand_index+rlen)])
	
	return reads
	
def get_suffixarray(genome: str, glen: len):
	"""
	Generates suffix array of the genome
	
	Parameters:
	- genome: genome sequence
	- glen: length of genome
	
	Returns:
	- Suffix array of the genome
	"""
	class Suffix:
		def __init__(self, suffix, index):
			self.suffix = suffix
			self.index = index
		def __gt__(self, other):
			return self.suffix > other.suffix

	# generate suffixes and keep them in alphabetical order
	suffix_table = []
	for i in range(glen):
		bisect.insort(suffix_table, Suffix(genome[i:], i))
	
	# get the indices from the suffixes
	sufarray = []
	for suf in suffix_table:
		sufarray.append(suf.index)

	return sufarray

def pattern_matching_suffix_array(genome: str, read: str, sufarray: list[str]):
	"""
	Map reads to genome with suffix array
	
	Parameters:
	- genome: genome sequence
	- read: sequence of read
	- sufarray: suffix array of the genome
	
	Returns:
	- A tuple containing indices of first and last match 
	"""
	glen = len(genome)
	rlen = len(read)
	min_index = 0
	max_index = glen - 1
	mid_index = None
	first = None
	last = None
	genome_fract = ""

	# find first match
	while min_index <= max_index:
		mid_index = int((min_index + max_index) / 2)
		genome_fract =  genome[sufarray[mid_index] : min(glen, sufarray[mid_index]+rlen)]
		if read > genome_fract:
			min_index = mid_index + 1
			mid_index += 1
		
		elif mid_index == 0:
			# `read` smaller than first suffix in `sufarray`
			break
		else:
			max_index = mid_index - 1

	genome_fract =  genome[sufarray[mid_index] : min(glen, sufarray[mid_index]+rlen)]
	if read == genome_fract:
		# first match found
		first = min_index
	else:
		# no match found
		sys.stderr.write("Error: read not found in genome\n")
		exit(1)
		return None

	# find last match
	min_index = first
	max_index = glen - 1
	while min_index <= max_index:
		mid_index = int((min_index + max_index) / 2)
		if read == genome[sufarray[mid_index] : sufarray[mid_index]+rlen]:
			min_index = mid_index + 1
		else:
			max_index = mid_index - 1
	last = max_index

	return (first, last)

def get_lc_and_fo(genome: str, glen: int):
	"""
	Get last character after cyclic rotations
	
	Parameters:
	- genome: genome sequence
	- glen: length of genome
	
	Returns:
	- A tuple: (list of last chars, dict of first occurance)
	"""
	class Table:
		def __init__(self, text, index):
			self.before = text[index:] + text[:index]
			self.last = text[index-1]
		def __gt__(self, other):
			return self.before > other.before
	
	# get all cyclic rotations of `genome`
	tables = []
	for i in range(glen):
		bisect.insort(tables, Table(genome, i))
	
	# get last column
	last_column = []
	for table in tables:
		last_column.append(table.last)
	
	# get first occurence
	first_occur = {}
	cur_symbol = None
	for i, table in enumerate(tables):
		if table.before[0] != cur_symbol:
			cur_symbol = table.before[0]
			first_occur[cur_symbol] = i

	return (last_column, first_occur)

def get_count(last_column: list):
	"""
	Get count matrix of BWT matching
	
	Parameters:
	- last_column: list of last characters from BWT transform
	
	Returns:
	- count list[list]: count matrix
	"""
	# initialize the count matrix
	count = {}
	unique_symbols = set(last_column)
	for symbol in unique_symbols:
		count[symbol] = [0]
	
	# fill in the count matrix
	for i, target in enumerate(last_column):
		for symbol in unique_symbols:
			count[symbol].append(count[symbol][i])
		count[target][i + 1] =  1 + count[target][i + 1]

	return count

def better_bw_matching(first_occurence: dict, last_column: list[str], read: str, count: list[list]):
	"""
	Map reads to BWT transformed genome sequence
	
	Parametrs:
	- first_occurence: dict of first occu. from BWT of genome
	- last_column: list of last characters from BWT of genome
	- read: sequence of read
	- count: count matrix
	"""
	top = 0
	bottom = len(last_column) - 1
	symbol = None
	while top <= bottom:
		if read != "":
			symbol = read[-1]
			read = read[0:-1]
			if symbol in last_column[top:bottom+1]:
				top = first_occurence[symbol] + count[symbol][top]
				bottom = first_occurence[symbol] + count[symbol][bottom+1] - 1
			else:
				return 0
		else:
			return bottom - top + 1

	sys.stderr.write("Error: reached unreachable region in better_bw_matching()\n")
	exit(1)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description = "Align reads to genome with suffix-array and BWT\n0: PatternMatchingWithSuffixArray\n1: BetterBWMatching\n2: checking if the two algorithms's results agree")
	parser.add_argument("-n", "--num_of_reads", type = int, required = True, help = "Number of reads to be generated (int)")
	parser.add_argument("-l","--len_of_read", type = int, required = True, help = "Length of reads (int)")
	parser.add_argument("-g","--len_of_genome", type = int, required = True, help = "Length of genome to be generated (int)")
	parser.add_argument("-a","--algorithm", type = str, required = True, help = "Type of mapping algorithm (str). [0|1|2]")
	args = parser.parse_args()

	rnum = args.num_of_reads
	rlen = args.len_of_read
	glen = args.len_of_genome
	
	genome = generate_genome(glen)
	reads = get_reads(rnum, rlen, genome, glen)

	# PatternMatchingWithSuffixArray
	if args.algorithm == "0":
		tracemalloc.start()
		sufarray = get_suffixarray(genome, glen)
		start_time = time.time()
		num_match = 0
		for read in reads:
			res = pattern_matching_suffix_array(genome, read, sufarray)
			num_match += res[1] - res[0] + 1
		current, peak = tracemalloc.get_traced_memory()
		print("PatternMatchingWithSuffixArray:")
		print("- runtime:", str(time.time() - start_time), "seconds")
		print("- results:", str(num_match), "matches found")
		print(f"- current memory usage: {current / 10**6}MB; peak: {peak / 10**6}MB")

	# BetterBWMatching
	elif args.algorithm == "1":
		tracemalloc.start()
		last_column, first_occurence = get_lc_and_fo(genome, len(genome))
		count = get_count(last_column)
		start_time = time.time()
		num_match = 0
		for read in reads:
			res = better_bw_matching(first_occurence, last_column, read, count)
			num_match += res
		current, peak = tracemalloc.get_traced_memory()
		print("BetterBWMatching:")
		print("- runtime:", str(time.time() - start_time), "seconds")
		print("- results:", str(num_match), "matches found")
		print(f"- current memory usage: {current / 10**6}MB; peak: {peak / 10**6}MB")
		tracemalloc.stop()
	
	# check if the two algorithms agree
	elif args.algorithm == "2":
		genome = generate_genome(glen)
		reads = get_reads(rnum, rlen, genome, glen)

		# PatternMatchingWithSuffixArray
		sufarray = get_suffixarray(genome, glen)
		num_match_array = 0
		for read in reads:
			res = pattern_matching_suffix_array(genome, read, sufarray)
			num_match_array += res[1] - res[0] + 1

		# BetterBWMatching
		last_column, first_occurence = get_lc_and_fo(genome, len(genome))
		count = get_count(last_column)
		maxi = 1
		num_match_bwm = 0
		for read in reads:
			res = better_bw_matching(first_occurence, last_column, read, count)
			maxi = max(maxi, res)
			num_match_bwm += res
		if num_match_bwm == num_match_array:
			print("Success. Results from the two algorithms matches")

	else:
		sys.stderr.write(f"Error: Input {args.algorithm} is not valid. Choose from [0|1|2]")
		exit(1)

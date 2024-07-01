#!/usr/bin/python

# Script name: DamageMasker.py
# Version: v1.0
# Author: Thomas Lesley Sitter
# Contact: lesleysitter@hotmail.com

# Description: Used to mask thymine and adenine sites in ancient SAM/BAM files created from single stranded libraries, to prevent damage from being incorporated into genotyping.

# Usability: When setting to hardmasking, al T's on the forward strand and all A's on the reverse strand are masked regardless of anything else
# Example: python DamageMasker.py --input_file Sample.processed.bam --masking R --ref_file Reference.fasta --output_file Sample_RefGuidedmasked.bam
# Usability: When setting the script to Reference guided masking, all T's on the forward strand are masked if the reference has a C, and all A's on the reverse strand will be masked if the reference has a G on that position
# Example: python DamageMasker.py --input_file Sample.processed.bam --masking H --output_file Sample_Hardmasked.bam
# Usability: When using edge masking, only T's on the 5' and 3' of the forward read are masked, and only A's on 5' and 3' of the reverse strand will be masked, the user can specify how many bases into these edges the masking runs.
# Example: python DamageMasker.py --input_file Sample.processed.bam --masking E --edge_count 3 --output_file Sample_Edgemasked.bam
# Optional: The user can also remove reads that are too short (default 0bp), or are not mapping with a high enough MapQ score (default 0) using the 'Filtering' option (-m F)
# Example: python DamageMasker.py --input_file Sample.processed.bam --output_file Sample_OnlyFiltered.bam --masking F --mapq_cutoff 37 --len_cutoff 25

# If python > v3 is installed as python3, make sure to run the software as 'python3 DamageMasker.py --input_file etc'

# Comment: Some features are experimental, and have had limited validation.

import os
import sys
import math
import argparse

# Set scoring variables for reference guided masking
nuc_fwd_dmg = {}
nuc_rev_dmg = {}
nuc_fwd_ins = {}
nuc_rev_ins = {}
nuc_fwd_mm = {}
nuc_rev_mm = {}
nuc_fwd_tot = {}
nuc_rev_tot = {}
nuc_total = {"Total":0}

# Check if the script is running with Python 3 since there is some python3 specific print statement further down the script
if sys.version_info.major < 3 or int(sys.version[:1]) < 3:
	print("\nError: This script requires Python 3. Please run it with a Python 3 interpreter.\n\n")
	sys.exit(1)

# Check if the pysam library is installed
try:
	import pysam
except ImportError:
	print("\nError: This script requires pysam to be installed.\nTry \'pip install pysam\' or \'pip3 install pysam\' (depending on your setup) to install it.\n\n")
	sys.exit(1)

# Check if the pysam library is installed
try:
	from Bio import SeqIO
except ImportError:
	print("\nError: This script requires biopython to be installed.\nTry \'pip install biopython\' or \'pip3 install biopython\' (depending on your setup) to install it.\n\n")
	sys.exit(1)

# Generate an STDout with the settings chosen, for reproducibility, and troubleshooting purposes
def settings_summary_printer(args):
	print(f"\nDamageMasker will now process {args.input_file}, using the following settings;")
	if args.masking == "R" :
		print(f"Masking type: Reference guided")
		print(f"Reference used: {args.ref_file}")
	if args.masking == "E" :
		print(f"Masking type: Edge masking (default)")
		print(f"Number of edge bases masked: {args.edge_count}")
	if args.masking == "H" :
		print(f"Masking type: Hard masking")
	if not args.mapq_cutoff == 0 :
		print(f"Remove reads with MapQ score below: {args.mapq_cutoff}")
	if not args.mapq_cutoff == 0 :
		print(f"Remove reads reads with length below: {args.len_cutoff} bp")
	print(f"Masked SAM/BAM will be saved to: {args.output_file}")
	if args.masking == "R" :
		print(f"Additional edit distance analysis is stored to: {args.output_file}_stats.tsv")
	print() # Just a paragraph breaker


# Determine if a bam or sam is provided
def sam_bam_picker(input_file,ref_file,output_file,mapq_cutoff,len_cutoff,masking,edge_count):

	# Check if  is mapped and MAPQ is above the cutoff
	if '.sam' in input_file :
		with pysam.AlignmentFile(input_file, 'r') as sam, pysam.AlignmentFile(output_file, 'w', header=sam.header) as output_sam:
			process_sam_bam(sam, output_sam, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count)

	elif '.bam' in input_file :
		with pysam.AlignmentFile(input_file, 'rb') as bam, pysam.AlignmentFile(output_file, 'wb', header=bam.header) as output_bam:
			process_sam_bam(bam, output_bam, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count)

def process_sam_bam(in_file, out_file, input_file, ref_file, output_file, mapq_cutoff, len_cutoff, masking,edge_count):

	if masking == "R":
		# Load the reference genome
		reference_dict = SeqIO.index(ref_file, "fasta")

	# Go through the reads in the FASTA
	for read in in_file:
		if not read.is_unmapped and read.mapping_quality >= mapq_cutoff and read.query_length >= len_cutoff:
			nuc_total['Total'] += 1 # Just counting the number of processed reads
			read_sequence = read.query_sequence # Get read sequence
			modified_qualities = read.query_qualities # Obtain the original quality values

			if masking == "R":
				# Load the reference sequence
				reference_seq = reference_dict[read.reference_name][read.reference_start:read.reference_end] # Get reference sequence

				# Deal with indels (annoying!!!)
				if len(read.cigartuples)>1 :
					pos_slider=0 # The cigar tuples give match (0), insertion (1) and deletion (2) information in tuples
					# and example would be [(0, 16), (2, 1), (0, 60)]. Here there are 16 matches, 1 deletion (of 1bp) in the read, and 60 more matches
					# Or [(0, 1), (1, 1), (0, 74)]. Here there is one match, 1 insertion (of additional nucleotide), and 74 more matches
					for cigar in read.cigartuples:
						if cigar[0] == 0 : 
							pos_slider += cigar[1]
						elif cigar[0] == 1 : # Add spacer to reference if a insertion occured in the read
							toggle=1
							reference_seq.seq = reference_seq.seq[:pos_slider] + ("-"*cigar[1]) + reference_seq.seq[pos_slider:]
							pos_slider += cigar[1]
						elif cigar[0] == 2 : # Add spacer to read if a deletion occured in the read
							read_sequence = read_sequence[:pos_slider] + ("-"*cigar[1]) + read_sequence[pos_slider:]
							pos_slider += cigar[1]

				# Forward read
				if not read.is_reverse:
					modified_sequence = ''.join(['N' if ref == 'C' and read == 'T' else read for ref, read in zip(reference_seq, read_sequence)])
					calc_CT = modified_sequence.count('N') # Calculate number of C>T transversions in read for scoring
					nuc_fwd_dmg[str(calc_CT)] = nuc_fwd_dmg.get(str(calc_CT), 0) + 1

					calc_indel= read_sequence.count('-') # Calculate number of indels in read for scoring
					nuc_fwd_ins[str(calc_indel)] = nuc_fwd_ins.get(str(calc_indel), 0) + 1
					calc_MM = ''.join(['Z' if ref != read and not read in ['X','Y','-'] else read for ref, read in zip(reference_seq, read_sequence)]).count('Z')
					nuc_fwd_mm[str(calc_MM)] = nuc_fwd_mm.get(str(calc_MM), 0) + 1
					num_total = modified_sequence.count('N') + calc_indel + calc_MM
					nuc_fwd_tot[str(num_total)] = nuc_fwd_tot.get(str(num_total), 0) + 1
				# Reverse read
				else:
					modified_sequence = ''.join(['N' if ref == 'G' and read == 'A' else read for ref, read in zip(reference_seq, read_sequence)])
					calc_GA = modified_sequence.count('N') # Calculate number of A<G transversions in read for scoring
					nuc_rev_dmg[str(calc_GA)] = nuc_rev_dmg.get(str(calc_GA), 0) + 1

					calc_indel= read_sequence.count('-') # Calculate number of indels in read for scoring
					nuc_rev_ins[str(calc_indel)] = nuc_rev_ins.get(str(calc_indel), 0) + 1
					calc_MM = ''.join(['Z' if ref != read and not read in ['X','Y','-'] else read for ref, read in zip(reference_seq, read_sequence)]).count('Z')
					nuc_rev_mm[str(calc_MM)] = nuc_rev_mm.get(str(calc_MM), 0) + 1
					num_total = modified_sequence.count('N') + calc_indel + calc_MM
					nuc_rev_tot[str(num_total)] = nuc_rev_tot.get(str(num_total), 0) + 1

				# Remove the gap spacers before saving the read back to the newly made SAM/BAM file
				modified_sequence = modified_sequence.replace('-', '')

			elif masking == "H":
				# Forward read
				if not read.is_reverse:
					modified_sequence = ''.join(['N' if read == 'T' else read for read in read_sequence])
				# Reverse read
				else:
					modified_sequence = ''.join(['N' if read == 'A' else read for read in read_sequence])

			elif masking == "E":
				if edge_count >= math.floor(len(read_sequence)/2):
					mask_numb = math.floor(len(read_sequence)/2)-1
				else:
					mask_numb = edge_count
				# Forward read
				if not read.is_reverse:
					modified_seq_L = ''.join(['N' if read == 'T' else read for read in read_sequence[:mask_numb]])
					modified_seq_R = ''.join(['N' if read == 'T' else read for read in read_sequence[len(read_sequence)-mask_numb:]])
					modified_sequence = modified_seq_L + read_sequence[mask_numb:len(read_sequence)-mask_numb] + modified_seq_R
				# Reverse read
				else:
					modified_seq_L = ''.join(['N' if read == 'A' else read for read in read_sequence[:mask_numb]])
					modified_seq_R = ''.join(['N' if read == 'A' else read for read in read_sequence[len(read_sequence)-mask_numb:]])
					modified_sequence = modified_seq_L + read_sequence[mask_numb:len(read_sequence)-mask_numb] + modified_seq_R

			elif masking == "F": # If masking is F it won't need to do anything, filtering already happend upstream of the function
				modified_sequence=read.query_sequence

			else:
				print(f"\nError: The masking setting '{masking}' is not recognized. Please use \'S\' for SoftMasking, \'H\' for HardMasking, and \'E\' for EdgeMasking.\n\n")
				sys.exit(1)

			read.query_sequence = modified_sequence  # Update the sequence field
			read.query_qualities = modified_qualities # Preserve the original quality values
			
			out_file.write(read)

def main():
	parser = argparse.ArgumentParser(description='Mask a SAM/BAM file for deaminated bases based on reference genome. The script can softmask, hardmask and edgemask', add_help=False)
	parser.add_argument('-m', '--masking', default="H", metavar='', help='Change masking behaviour.\n\'R\' for Reference based Masking.\n\'H\' for HardMasking.\n\'E\' for EdgeMasking.\n\'F\' for Filtering. (default: Hardmasking)\n')
	parser.add_argument('-i', '--input_file', metavar='', help='Input BAM or SAM file (mandatory)')
	parser.add_argument('-r', '--ref_file', default="NA", metavar='', help='Input reference genome file in FASTA format (mandatory for --masking \'R\')')
	parser.add_argument('-e', '--edge_count', type=int, default=5, metavar='', help='Number of bases to be masked from 5\' and 3\' edges if --masking \'E\' is turned on (default: 5)')
	parser.add_argument('-q', '--mapq_cutoff', type=int, default=0, metavar='', help='Remove reads below MAPQ cutoff value from output (default: 0)')
	parser.add_argument('-l', '--len_cutoff', type=int, default=0, metavar='', help='Remove reads below a certain length from output (default: 0)')
	parser.add_argument('-o', '--output_file', metavar='', default='output_modified.sam', help='Output SAM/BAM file with modified reads (default: \'output_modified.sam/.bam\')')
	parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
	args = parser.parse_args()

	# Check if the input BAM file exists
	if not os.path.exists(args.input_file):
		print(f"\nError: Input SAM/BAM file '{args.input_file}' not found.\n\n")
		return
		# Check if the input BAM file exists
	extention_test = args.input_file.lower()
	if '.sam' in extention_test :
		if not args.output_file.endswith(".sam") :
			args.output_file += ".sam"
	elif '.bam' in extention_test :
		if not args.output_file.endswith(".bam") :
			args.output_file += ".bam"
	else:
		print(f"\nError: It could not be determined if '{input_file}' is a SAM or BAM formatted file. Make sure the file has a .sam or .bam extention.\n\n")
		sys.exit(1)

	settings_summary_printer(args)

	if args.masking == "R":
		if not args.ref_file =="NA":
			if not os.path.exists(args.ref_file):
				print(f"\nError: Input reference FASTA file '{args.ref_file}' not found.\n\n")
				return
		else:
			print(f"\nError: When using the \'R\' Reference guided masking, the user has to supply a reference genome with --ref_file.\nNo Reference was supplied.\nPlease supply the reference genome against which the reads were mapped,\nor choose the \"H\" or \"E\" masking option instead by setting these using --masking.\n\n")
			return
	
	sam_bam_picker(args.input_file, args.ref_file, args.output_file, args.mapq_cutoff, args.len_cutoff, args.masking, args.edge_count)
	print(f"\nFinished processing {args.input_file}, {nuc_total['Total']} reads processed\n\n")

	if args.masking == "R":
		with open(args.output_file+"_stats.tsv", "a") as o:	# Output an additional stats tsv file containing "edit distance"
			print(f"Analysis for {args.input_file}, total reads analyzed: {nuc_total['Total']}\n", file=o, end="")
			print(f"count\tfwd_reads_dmg\tfwd_reads_indels\tfwd_reads_other_mismatches\tfwd_reads_mismatches_total\t", file=o, end="")
			print(f"rev_reads_dmg\trev_reads_indels\trev_reads_other_mismatches\trev_reads_mismatches_total\t", file=o, end="")
			print(f"tot_reads_dmg\ttot_reads_indels\ttot_reads_other_mismatches\ttot_reads_mismatches_total", file=o)
			for i in range(0,10):
				print(f"{i}\t", file=o, end="")
				print(f"{nuc_fwd_dmg.get(str(i), 0)}\t{nuc_fwd_ins.get(str(i), 0)}\t{nuc_fwd_mm.get(str(i), 0)}\t{nuc_fwd_tot.get(str(i), 0)}\t", file=o, end="")
				print(f"{nuc_rev_dmg.get(str(i), 0)}\t{nuc_rev_ins.get(str(i), 0)}\t{nuc_rev_mm.get(str(i), 0)}\t{nuc_rev_tot.get(str(i), 0)}\t", file=o, end="")
				print(f"{(nuc_fwd_dmg.get(str(i), 0)+nuc_rev_dmg.get(str(i), 0))}\t{(nuc_fwd_ins.get(str(i), 0)+nuc_rev_ins.get(str(i), 0))}\t", file=o, end="")
				print(f"{(nuc_fwd_mm.get(str(i), 0)+nuc_rev_mm.get(str(i), 0))}\t{(nuc_fwd_tot.get(str(i), 0)+nuc_rev_tot.get(str(i), 0))}", file=o)
				
if __name__ == "__main__":
	main()


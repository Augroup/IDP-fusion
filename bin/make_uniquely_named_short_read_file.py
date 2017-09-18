#!/usr/bin/env python
import sys, os, re
import sequence_basics

# pre: a fastq or fasta file.  If a fasta file is used it must have the 
#      file extension of .fa or .fasta
# post: an output file in the same format (fasta or fastq) as the input with 
#       new names for all the entries, each unique.

def main():
  if len(sys.argv) != 3:
    print sys.argv[0] + ' <short reads file> <short read out file>'
    sys.exit()

  reads_filename = sys.argv[1]
  out_filename = sys.argv[2]

  if re.search('\.fa$|\.fasta$',reads_filename):
    sequence_basics.fasta_to_unique_name_fasta(reads_filename,out_filename)
  else:
    sequence_basics.fastq_to_unique_name_fastq(reads_filename,out_filename)

main()

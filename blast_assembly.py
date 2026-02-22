#importing requirec packages
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="<add description of what script does>")
    parser.add_argument("-i", "--input", 
    help="first input file",
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

with open(infile, "r") as handle:
    seqs= [] #initializing list to store sequences
    ids= []
    for record in SeqIO.parse(handle, "fasta"): #using SeqIO to parse file
        ids += [record.id]
        seqs += [record.seq] #pulling out each sequence and adding to the lsit
longest_contig = seqs[0] #first contig is longest based on how SPADES orders its output
longest_contig = SeqRecord(Seq(longest_contig), id = ids[0]) #making a SeqRecord with sequence and ID

with open(outfile, 'w') as o: #writing output to outfile
    SeqIO.write(longest_contig, o, "fasta")

#blastn -query assembles -db {blast_db} -max_hsps 1 -out {output.blast_out} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
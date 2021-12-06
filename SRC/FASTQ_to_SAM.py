# Import libraries 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# Open FASTQ file 
f = open("sequences.fastq", "r")
lines = f.readlines()

# Create output file
outfile = open("clean.txt", "w")

# Clean data and write outfile to have name of organism and sequence 
for line in lines:
  first_char = line[0]
  split_line = line.split()
  if first_char == '@':
    name = split_line[1]
    outfile.write(name)
    outfile.write("\n")
  elif (first_char == 'A') | (first_char == 'T') | (first_char == 'G') | (first_char == 'C'):
    seq = split_line[0]
    if seq.isalpha() == True:
      sequence = seq
      outfile.write(sequence)
      outfile.write("\n")
     
# Close output file
outfile.close()

# Read in cleaned data
f = open("clean.txt", "r")
clean_lines = f.readlines()

# Create list to hold name of organism and sequences 
names = []
seq_collection = []

# Add organism and sequences to respective lists 
for line in clean_lines:
  first_char = line[0]
  if (first_char == 'A') | (first_char == 'T') | (first_char == 'G') | (first_char == 'C'):
    split_line = line.split()
    str_convert = ""
    seq = str_convert.join(split_line)
    seq_collection.append(seq)
  else:
    names.append(line)

# Create list to hold RNA sequences
RNA_seq = []

# Convert DNA sequences to RNA sequences
for i in range(len(seq_collection)):
  coding_dna = Seq(seq_collection[i])
  messenger_rna = coding_dna.transcribe()
  RNA_seq.append(messenger_rna)

# Record respective sequences 
a = SeqRecord(Seq(RNA_seq[0]), id="Alpha")
b = SeqRecord(Seq(RNA_seq[1]), id="Beta")
c = SeqRecord(Seq(RNA_seq[2]), id="Gamma")

# Convert Reference DNA to RNA
coding_DNA = Seq("GGNTNGGTTACATCTGGTTACTGTCCTGGGTAAATCATTTTTATAGAGATGGCCTTCCAAGTGGTTTTAAAATTTACTGAAGTTTTTAGGTCAATTATGTAAAACAAGTTTATTTGTAAATTTAGTCAACATACATAATTGACCTAAAAACTTCAGTAAATTTTAAAACCACTTGGAAGGCCATCTCTATAAAAATGATTTT")
reference_mRNA = str(coding_DNA.transcribe())

# Align each sequence in a multiple sequence alignment 
align = MultipleSeqAlignment([a, b, c],
annotations={"tool": "demo"},
column_annotations={"stats": reference_mRNA})

# Create output file
outfile_two = open("align.sam", "w")

# Write out name of organism and sequence in outfile 
for i in range(len(seq_collection)):
  outfile_two.write(">" + names[i])
  outfile_two.write(str(align[i].seq))
  outfile_two.write("\n")
  
outfile_two.close()




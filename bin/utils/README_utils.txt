These scripts are intended to get and parse sequence files from genbank 
and make them compatible with POTION/OrthoMCL 1.4 input.

#####################
# get_seq_by_ids.pl #
#####################

Takes as input a single argument, namely a list containing genbank IDs (one per
line). Downloads each genbank file in current directory.

############################
# extract_aa_nt_from_gb.pl #
############################

Takes as input two arguments. First argument is a directory containing one or
more genbank files containing CDS data. Second argument is output directory.

This script parses genbank files, extract CDS and translation data and generates
three files for each genbank file: A fasta file containing CDS sequences, a 
fasta file containing protein (translation) sequences and a dictionary file 
linking the temporary IDs present in CDS and translation files (same ID on both
files) to the original gene/protein IDs from genbank. After creating these
files, users can use the protein fasta files as input to OrthoMCL 1.4 and,
since the temporary IDs for each gene/protein pair are the same, users can take
the OrthoMCL output directly and use it as input to POTION, together with CDS
data. The dictionary file can be used at the end of the analysis to link POTION
results to original protein/CDS IDs.

##############################
# split_multi_organism_gb.pl #
##############################

Takes as input a single argument, namely the path to a genbank file containing
two or more species. Creates as output a genbank file for each species within
original genbank file.

###################
# orthoXML2MCL.pl #
###################

Takes as input a single argument, namely a xml file describing homology
relationships of sequences (orthoXML format), producing as result an homology
description file compatible with POTION (ORTHOMCL V 1.4 format)

#####################
# get_seq_by_ids.pl #
#####################

Takes as input three arguments: a text file with a list of IDs (one per line),
the path to a fasta file and the name of an output file. Produces a fasta file
containing the sequences in the list of ids.

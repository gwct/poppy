#!/usr/bin/python3
############################################################
# For popgen, 10.19
# Takes a FASTA file alignment and transposes the matrix
# so each line corresponds to a site. This (hopefullly) will
# allow for faster reading and processing. Writes this 
# transposed FASTA file as a .tfa file. Also writes a .title
# file containing the FASTA headers for each sequence.
############################################################

import sys, os, argparse

############################################################

def fastaLists(infile):
    indata = open(infile, "r").read();

    titles, seqs = [],[];

    indata = indata.split("\n");
    indata = filter(None, indata);

    first = True;
    for line in indata:
        if line[0] == ">":
            titles.append(line);
            if not first:
                seqs.append(seq);
            else:
                first = False;
            
            seq = "";
            continue;

        seq += line;

    
    seqs.append(seq);
    return titles, seqs;

############################################################

print("# Program call: " + " ".join(sys.argv) + "\n");

parser = argparse.ArgumentParser(description="");
parser.add_argument("-i", dest="input", help="The path to a FASTA alignment file with sequences from different groups.", default=False);
#parser.add_argument("-o", dest="output", help="The output file name.", default=False);
args = parser.parse_args();
# Input options.

if not args.input or not os.path.isfile(args.input):
    sys.exit(" * ERROR 1: A valid input FASTA alignment file must be specified with -i.");

print("# Reading input...");
titles, seqs = fastaLists(args.input);
#inseqs = core.fastaReader(args.input);

# print(len(titles));
# print(len(seqs));
# print(titles);

tmpfilename = os.path.splitext(args.input)[0] + ".tmp298534";
tfafilename = os.path.splitext(args.input)[0] + ".tfa";
titlefilename = os.path.splitext(args.input)[0] + ".title";

# print(tfafilename);
# print(titlefilename);

numseqs = len(seqs);
seqlen = len(seqs[0]);

with open(tfafilename, 'w') as tfafile, open(titlefilename, 'w') as titlefile:
    titlefile.write("\n".join(titles));
    for i in range(seqlen):
        outline = "";
        for j in range(numseqs):
            outline += seqs[j][i];
        tfafile.write(outline + "\n");




# tmpfile.write("\n".join(seqs));

# transpose = """
# awk '
# {
#     for (i=1; i<=NF; i++) {
#         a[NR,i] = $i
#     }
# }
# NF>p { p = NF }
# END {
#     for(j=1; j<=p; j++) {
#         str=a[1,j]
#         for(i=2; i<=NR; i++) {
#             str=str" "a[i,j];
#         }
#         print str
#     }
# }' """ + tmpfilename + " > " + tfafilename;

# print(transpose);
# os.system(transpose);
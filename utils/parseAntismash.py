import sys
import os.path
from optparse import OptionParser
from Bio import SeqIO

#### Collect Input ####
#######################
usage="""Takes a antismash GenBank, and outputs a txt file.
usage: %prog -i FILE [options]"""

parser = OptionParser(usage=usage, version="%prog 1.0")


parser.add_option("-i", "--in_file", metavar="FILE", dest="in_file", default=None,
                help="Specify the input FILE")


parser.add_option("-o", "--out_file", metavar="FILE", dest="out_file", default=None,
                help="Specify the path and name of the output fasta file you wish to create. "
                "Default will be the same as the in_file, but with a 'txt' suffix.")

(options, args) = parser.parse_args()


#### Variables and Names ####
#############################
#Figure out some names and paths
if options.in_file:
    in_file = os.path.abspath(options.in_file)
else:
    print ("You must specify an in_file. Use '-h' for help.")
    sys.exit()

#Figure out what our out_file is.
if options.out_file:
    out_file = os.path.abspath(options.out_file)
else:
    print ("You must specify an in_file. Use '-h' for help.")
    sys.exit()


def parseAntismashGBK(genbank, out):
    ##parse genbank file
    for seq_record in SeqIO.parse(genbank, "genbank"):
        for seq_feat in seq_record.features:
            if seq_feat.type == "protocluster":
               number = seq_feat.qualifiers["protocluster_number"][0]
               cluster = seq_feat.qualifiers["product"][0]
               out.write( number + "\t" + cluster +"\n")


#### Main ####
##############
in_file_handle = open (in_file, 'r', newline='')
out_file_handle = open (out_file, 'w+')

parseAntismashGBK(in_file_handle, out_file_handle)

in_file_handle.close()
out_file_handle.close()

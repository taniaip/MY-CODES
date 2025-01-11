import argparse
import gffutils
from Bio import SeqIO


# Set up argument parser
parser = argparse.ArgumentParser()


parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")


# Parse the arguments
args = parser.parse_args()

db = gffutils.create_db(args.gff3_file, dbfn=':memory:')

list_tuples=[]
tuple=()

def list_extraction_motif(gff3_file):
    db = gffutils.create_db(args.gff3_file, dbfn=':memory:')
    for motif in db.all_features():
        tuple=(motif.id, motif.start, motif.end)
        list_tuples.append(tuple)
    return list_tuples


def main(args):
    result = list_extraction_motif(args.gff3_file)
    print(result)


if __name__ == '__main__':
    main(args)
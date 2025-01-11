"""#!/usr/bin/env python

# standard stuff
import argparse

# stuff you need to install
import gffutils
import regex as re
# import re
from Bio import SeqIO

# argument parser
parser = argparse.ArgumentParser(
            description=
        '''
        Find T....TGTTTGTT motifs in a Blastocystis genome,
        and report them as gff3 features
        '''
        )

parser.add_argument(
        "-f", "--fasta",
        dest='fasta_file',
        type=str,
        required=True,
        help="Input FASTA file")
args = parser.parse_args()


def main(args) -> None:

    # load fasta
    fasta = SeqIO.index(args.fasta_file, 'fasta')

    # for seq_record in fasta.values():
    for seqid, seqrecord in fasta.items():

        # motif counter
        n: int = 0

        # print(seqid)
        seq: str = str( seqrecord.seq ).upper()
        # print(seq)

        fwd_pattern = r'(T|TA|TG)[ACTG]{4}((?:TGTTTGTT){s<=2})'
        for match in re.finditer(fwd_pattern, seq):
            n += 1
            start: int = match.start() + 1
            end: int = match.end()
            # print()
            # print(start, end)
            # print(match.group(0))
            # print(f'{match.group(1)}    {match.group(2)}')
            
            motif = gffutils.feature.Feature(
                    seqid = seqid,
                    source = 'blasto_motif.py',
                    featuretype = 'motif',
                    start = start,
                    end = end,
                    score = '.',
                    strand = '+',
                    frame = '.',
                    attributes = f'ID={seqid}_STOP_codon_motif_{n:04d}'
                    )

            print(motif)

        rev_pattern = r'((?:AACAAACA){s<=2})[ACTG]{4}(A|TA|CA)'
        for match in re.finditer(rev_pattern, seq):
            n += 1
            start: int = match.start() + 1
            end: int = match.end()
            # print()
            # print(start, end)
            # print(match.group(0))
            # print(f'{match.group(1)}    {match.group(2)}')
            
            motif = gffutils.feature.Feature(
                    seqid = seqid,
                    source = 'blasto_motif.py',
                    featuretype = 'motif',
                    start = start,
                    end = end,
                    score = '.',
                    strand = '-',
                    frame = '.',
                    attributes = f'ID={seqid}_STOP_codon_motif_{n:04d}'
                    )

            print(motif)

if __name__ == '__main__':
    main(args)
"""

import argparse

# stuff you need to install
import gffutils
import regex as re
from Bio import SeqIO

# argument parser
parser = argparse.ArgumentParser(
    description='''
        Find T....TGTTTGTT motifs in a Blastocystis genome,
        and report them as gff3 features
    '''
)

parser.add_argument(
    "-f", "--fasta",
    dest='fasta_file',
    type=str,
    required=True,
    help="Input FASTA file"
)
args = parser.parse_args()

# Function to count the number of substitutions
def count_substitutions(motif, canonical_motif):
    sum=0
    for a, b in zip(motif, canonical_motif) if a != b:
    sum +=1
    return sum

def main(args) -> None:
    # load fasta
    fasta = SeqIO.index(args.fasta_file, 'fasta')

    # for seq_record in fasta.values():
    for seqid, seqrecord in fasta.items():

        # motif counter
        n: int = 0

        # print(seqid)
        seq: str = str(seqrecord.seq).upper()
        # print(seq)

        # Canonical forward motif for comparison
        canonical_fwd = "TGTTTGTT"
        canonical_rev = "AACAAACA"

        # Forward pattern with up to 2 substitutions allowed
        fwd_pattern = r'(T|TA|TG)[ACTG]{4}((?:TGTTTGTT){s<=2})'
        for match in re.finditer(fwd_pattern, seq):
            n += 1
            start: int = match.start() + 1
            end: int = match.end()

            matched_motif = match.group(2)  # The actual motif part (TGTTTGTT)
            num_substitutions = count_substitutions(matched_motif, canonical_fwd)

            motif = gffutils.feature.Feature(
                seqid=seqid,
                source='blasto_motif.py',
                featuretype='motif',
                start=start,
                end=end,
                score='.',
                strand='+',
                frame='.',
                attributes=f'ID={seqid}_STOP_codon_motif_{n:04d};Substitutions={num_substitutions}'
            )

            print(motif)

        # Reverse pattern with up to 2 substitutions allowed
        rev_pattern = r'((?:AACAAACA){s<=2})[ACTG]{4}(A|TA|CA)'
        for match in re.finditer(rev_pattern, seq):
            n += 1
            start: int = match.start() + 1
            end: int = match.end()

            matched_motif = match.group(1)  # The actual motif part (AACAAACA)
            num_substitutions = count_substitutions(matched_motif, canonical_rev)

            motif = gffutils.feature.Feature(
                seqid=seqid,
                source='blasto_motif.py',
                featuretype='motif',
                start=start,
                end=end,
                score='.',
                strand='-',
                frame='.',
                attributes=f'ID={seqid}_STOP_codon_motif_{n:04d};Substitutions={num_substitutions}'
            )

            print(motif)

if __name__ == '__main__':
    main(args)

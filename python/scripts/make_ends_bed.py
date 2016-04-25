"""Use the scaffold fasta to produce a bed file representing the ends of each scaffold.

Scaffolds are only included if they are at least as 2x as long as the supplied window.
"""
import logging

logging.basicConfig(filename=snakemake.log.path,level=logging.INFO)


from pyfaidx import Fasta

import pandas as pd




# Settings
window = snakemake.params.window

# input
fasta_path = snakemake.input.fasta_path

# output
bed_path = snakemake.output.bed_path

# logging



def make_ends_bed_recs(fasta, window, excluded, included):
    """Use the scaffold fasta to produce a bed file representing the ends of each scaffold.

    Scaffolds are only included if they are at least as 2x as long as the supplied window.
    """
    for name, seq in fasta.records.items():
        seq_range = len(seq)

        chrom = name
        startL = 0
        endL = window
        startR = seq_range - window
        endR = seq_range

        if seq_range < (window * 2) + 1:
            excluded.append(chrom)
            logging.info(" {chrom} with length {range} was excluded from bed.".format(chrom=chrom,
                                                                                     range=seq_range))
            continue

        included.append(chrom)
        yield "{chrom}\t{startL}\t{endL}\n{chrom}\t{startR}\t{endR}\n".format(chrom=chrom,
                                                                        startL=startL,
                                                                        endL=endL,
                                                                        startR=startR,
                                                                        endR=endR,
                                                                       )

def cummulative_length(fasta, keys):
    """Return the sum of the lengths represented by `keys`."""
    keys = set(keys)

    lengths = []
    for name, seq in fasta.records.items():
        if name in keys:
            lengths.append(len(seq))

    return sum(lengths)



#### Do the real work ####

fas = Fasta(fasta_path)

# Write the bed file
with open(bed_path, 'w') as bed:

    bed.write("#beginnings and ends\n")

    excluded = []
    included = []
    for rec in make_ends_bed_recs(fasta=fas, window=window, excluded=excluded, included=included):
        bed.write(rec)

# Report the proportion of the genome retained and other information

excluded_bp = cummulative_length(fasta=fas, keys=excluded)
included_bp = cummulative_length(fasta=fas, keys=included)
total_bp = sum([len(seq) for seq in fas.records.values()])

assert excluded_bp + included_bp == total_bp

logging.info(" {excluded} of {total} scaffolds were excluded from the bed file.".format(excluded=len(excluded),
                                                                                           total=len(fas.records)))

logging.info(" {frac} of {total} bp are represented by the scaffolds that remained.".format(frac=included_bp/total_bp,total=total_bp))

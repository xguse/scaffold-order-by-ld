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



def make_ends_bed_recs(fasta, window, excluded):
    """Use the scaffold fasta to produce a bed file representing the ends of each scaffold.

    Scaffolds are only included if they are at least as 2x as long as the supplied window.
    """
    for name, seq in list(dict(fas.records).items())[:10]:
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

        yield "{chrom}\t{startL}\t{endL}\n{chrom}\t{startR}\t{endR}\n".format(chrom=chrom,
                                                                        startL=startL,
                                                                        endL=endL,
                                                                        startR=startR,
                                                                        endR=endR,
                                                                       )


#### Do the real work ####

fas = Fasta(fasta_path)

# WHY IS THIS ONLY OUTPUTTING 3 to 6 CHROMS and ENDING?!

with open(bed_path, 'w') as bed:

    excluded = []
    for rec in make_ends_bed_recs(fasta=fasta_path, window=window, excluded=excluded):
        bed.write(rec)

    logging.info(" {excluded} of {total} scaffolds were excluded from the bed file.".format(excluded=len(excluded),
                                                                                           total=len(fas.records)))

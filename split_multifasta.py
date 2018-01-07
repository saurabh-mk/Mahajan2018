## Author: Saurabh Mahajan
## Created: Oct 16, 2016
## Description: This script breaks a fasta file into several chunks. It is called by the shell script "codon_count_calculations_script.sh". 
## Input arguments: genome id (label for cds file), data directory path, output directory path, and batch size i.e. how many cds to put in each chunk
## Depends on: cds file (fasta), biopython

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

from Bio import SeqIO
import sys

GCF_id = sys.argv[1]
data_dir = sys.argv[2]
output_dir = sys.argv[3]
batch_size = int(sys.argv[4])

fn_oi = data_dir+"/"+GCF_id+"_cds_from_genomic.fna"
record_iter = SeqIO.parse(open(fn_oi),"fasta")
for i, batch in enumerate(batch_iterator(record_iter, batch_size)):
    w_fn = output_dir+"/"+GCF_id+"_cds_from_genomic_%i.fna" % (i + 1)
    with open(w_fn, "w") as w_fh:
        count = SeqIO.write(batch, w_fh, "fasta")
    print("Wrote %i records to %s" % (count, w_fn))


from __future__ import absolute_import
from collections import OrderedDict
from Bio.Seq import Seq
from Bio import SeqIO
from io import open
import six


def read_from_fasta(path, invert = False):
    '''
    Read in a FASTA file and create a (id:seq) or (seq:id) ordered dict.
    '''
    return OrderedDict(
        (str(r.seq), r.id) if invert else (r.id, str(r.seq))
        for r in SeqIO.parse(path, "fasta")
    )

def translate(s):
    '''
    Assume we are in frame and translate DNA to amino acids.
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))])
    return str(coding_dna.translate())

def write_to_fasta(d, filename):
    '''
    Write a FASTA dict to file.
    '''
    with open(filename, 'w') as f:
        for k, v in six.iteritems(d):
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))

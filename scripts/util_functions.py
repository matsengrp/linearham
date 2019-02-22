from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def translate(s):
    '''
    Assume we are in frame and translate DNA to amino acids.
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], IUPAC.ambiguous_dna)
    return str(coding_dna.translate())

def write_to_fasta(d, filename):
    '''
    Write a FASTA dict to file.
    '''
    with open(filename, 'w') as f:
        for k, v in d.items():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))

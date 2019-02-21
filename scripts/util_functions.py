from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


def translate(s):
    '''
    Assume we are in frame and translate DNA to amino acids.
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], IUPAC.ambiguous_dna)
    return str(coding_dna.translate())

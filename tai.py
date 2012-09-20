from biochem import *
from utils import *

class TAI(object):

    def __init__(self,trnas):
        """Initialize a TAI object from a dictionary or list (sorted
        according to sorted_codons) of tRNA copy numbers"""
        if type(trnas) is list:
            print trnas
            self.trnas = trnas
            self.trna_dict = {codon:trnas[i] for i,codon in enumerate(sorted_codons)}
        elif isinstance(trnas,dict):
            self.trnas = [trnas[codon] for codon in sorted_codons]
            self.trna_dict = trnas
        else:
            raise Exception("trnas is of type:",type(trnas))
    def W(self,codon_i):
        return sum([(1-self.s(codon_i,codon_j)) * self.tGCN(codon_i,codon_j)
                    for codon_j in wobble_neighbors(codon_i)])

    def s(self,codon_i,codon_j):
        if codon_i == codon_j:
            return 0
        elif codon_i[:2] == codon_j[:2]:
            return 0.5
        else:
            return 1

    def tGCN(self,codon_i,codon_j):
        return self.trna_dict[codon_j]

    def tai(self,seq):
        """Given a list of codons, return the tAI index"""
        seq_codons = [codon.lower() for codon in group_codons(seq)]
        Ws = map(self.W,seq_codons)
        L = len(seq_codons)
        Wmax = max(Ws)
        ws = [W/Wmax for W in Ws if W > 0]
        return (product(ws)) ** (1/float(L))

print "loaded tai"

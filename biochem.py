"""This module contains code describing biochemical relationships
re: nucleotides and amino acids"""

import random
from utils import *
from Bio.SubsMat import MatrixInfo as matlist

BLOSUM62 = matlist.blosum62

def blosum(a,b):
    try:
        score = BLOSUM62[(a,b)]
    except KeyError:
        score = BLOSUM62[(b,a)]
    return score

delta = "acgt"
codons = [x+y+z for x in delta for y in delta for z in delta]

def neighboring_codons(codon):
    return [c for c in codons
            if sum(zipWith(operator.eq,codon,c)) >= 2]

def wobble_neighbors(codon):
    return [codon[:2] + c for c in delta]

def neighboring_aas(aa):
    self_codons = translation_table[aa]
    neighbors = list(set([neighbor for codon in self_codons
                          for neighbor in neighboring_codons(codon)]))
    return list(set([translate(neighbor) for neighbor in neighbors]))


def translate(codon):
    if len(codon) == 3:
        for aa in translation_table:
            if codon in translation_table[aa]:
                return aa
    else:
        return "".join(map(translate,group_codons(codon)))

translation_table = {"A":["gca","gcc","gcg","gct"],
                     "R": ["aga","agg","cga","cgc","cgg","cgt"],
                     "N": ["aat", "aac"],
                     "D": ["gat", "gac"],
                     "C": ["tgt", "tgc"],
                     "E": ["gaa", "gag"],
                     "Q": ["caa", "cag"],
                     "G": ["gga","ggc","ggg","ggt"],
                     "H": ["cat", "cac"],
                     "I": ["att", "atc", "ata"],
                     "L": ["tta", "ttg", "cta","ctc","ctg","ctt"],
                     "K": ["aaa", "aag"],
                     "M": ["atg"],
                     "F": ["ttt", "ttc"],
                     "P": ["cca","ccc","ccg","cct"],
                     "S": ["agt", "agc", "tca","tcc","tcg","tct"],
                     "T": ["aca","acc","acg","act"],
                     "W": ["tgg"],
                     "Y": ["tat", "tac"],
                     "V": ["gta","gtc","gtg","gtt"],
                     "X": ["taa", "tag","tga"]}

six_box_translation_table = {"A":["gca","gcc","gcg","gct"],
                             "R2": ["aga","agg"],
                             "R4": ["cga","cgc","cgg","cgt"],
                             "N": ["aat", "aac"],
                             "D": ["gat", "gac"],
                             "C": ["tgt", "tgc"],
                             "E": ["gaa", "gag"],
                             "Q": ["caa", "cag"],
                             "G": ["gga","ggc","ggg","ggt"],
                             "H": ["cat", "cac"],
                             "I": ["att", "atc", "ata"],
                             "L2": ["tta", "ttg"],
                             "L4": ["cta","ctc","ctg","ctt"],
                             "K": ["aaa", "aag"],
                             "M": ["atg"],
                             "F": ["ttt", "ttc"],
                             "P": ["cca","ccc","ccg","cct"],
                             "S2": ["agt", "agc"],
                             "S4": ["tca","tcc","tcg","tct"],
                             "T": ["aca","acc","acg","act"],
                             "W": ["tgg"],
                             "Y": ["tat", "tac"],
                             "V": ["gta","gtc","gtg","gtt"],
                             "X": ["taa", "tag","tga"]}

four_box_codons = [codon for codon in codons
                   if len(translation_table[(translate(codon))]) == 4]

two_box_codons = [codon for codon in codons
                   if len(translation_table[(translate(codon))]) == 2]

sorted_codons = [codon for aa in translation_table
                 for codon in translation_table[aa]]

six_box_sorted_codons = [codon for aa in translation_table
                         for codon in translation_table[aa]]

table_codons = [x + y + z for x in "tcag" for y in "tcag" for z in "tcag"]

def synonymous_codons(codon,break_down_six_boxes=False):
    tt = six_box_translation_table if break_down_six_boxes else translation_table
    print codon
    return head(translation_table.values(),lambda codons:codon in codons)

aas = [k for k in translation_table.keys() if not k == 'X']
def substitutability(aa,new_aa):
    aa = aa.upper()
    new_aa = new_aa.upper()
    if aa == "X" or new_aa == 'X':
        score = -1000
    else:
        score = blosum(new_aa,aa)
    return score

def substitution_run(codon):
    history = []
    aa = translate(codon)
    while not aa == 'X':
        codon = mutate_codon(codon)
        new_aa = translate(codon)
        score = substitutability(aa,new_aa)
        history.append(score)
        aa = new_aa
    return history[:len(history) - 1]

def codon_exp(codon,n=1000):
    aa = translate(codon)
    for i in range(n):
        new_codon = mutate_codon(codon)
        new_aa = translate(new_codon)
        print new_codon,new_aa
        score = substitutability(aa,new_aa)
        print "%s -> %s, %s -> %s, %s" %(codon,new_codon,aa,new_aa,score)
        codon = new_codon
        aa = new_aa

def mutate_codon(codon):
    i = random.randint(0,2)
    nuc = codon[i]
    return (codon[:i] +
            random.choice([d for d in delta if not d == nuc]) +
            codon[i+ 1:])

def codon2index(codon):
    c = lambda n: "ACGT".index(n)
    return sum(map(lambda ((i, n)):(4**(2 - i))*c(n), enumerate(codon)))
    
def codon_summary(s):
    seq = s.tostring()
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    summary = [0 for i in range(64)]
    for codon in codons:
        summary[codon2index(codon)] += 1
    return summary

def wc(xs):
    def pair(c):
        return {"A":"T","T":"A","G":"C","C":"G"}[c]
    return "".join((map(pair,xs))[::-1])

print("loaded biochem")

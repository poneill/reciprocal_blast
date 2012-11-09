"""This module provides an interface to the COG database as abstracted
by the file whog.txt"""

import re
from collections import defaultdict
def parse_whog():
    filename = "cog/whog.txt"
    with open(filename) as f:
 	lines = [line.strip() for line in f.readlines()]
    cogs = defaultdict(list)
    heading_pattern = r"\[(.)\] (COG[0-9]+) (.+)"
    locus_tag_pattern = r":  (.*)"
    for line in lines:
        possible_heading = re.match(heading_pattern,line)
        if possible_heading:
            category = possible_heading.groups()[0]
        else:
            possible_gene = re.search(locus_tag_pattern,line)
            if possible_gene:
                gene = possible_gene.groups()[0].split(" ")
                cogs[category].extend(gene)
    return cogs

# cogs = parse_whog()

def get_cog(locus_tag):
    for cat in cogs:
        if locus_tag in cogs[cat]:
            return cat
    return False

cog_definitions = {"A":	"RNA processing and modification",
                   "B":	"Chromatin Structure and dynamics",
                   "C":	"Energy production and conversion",
                   "D":	"Cell cycle control and mitosis",
                   "E":	"Amino Acid metabolis and transport",
                   "F":	"Nucleotide metabolism and transport",
                   "G":	"Carbohydrate metabolism and transport",
                   "H":	"Coenzyme metabolis",
                   "I":	"Lipid metabolism",
                   "J":	"Tranlsation",
                   "K":	"Transcription",
                   "L":	"Replication and repair",
                   "M":	"Cell wall/membrane/envelop biogenesis",
                   "N":	"Cell motility",
                   "O":	"Post-translational modification, protein turnover, chaperone functions",
                   "P":	"Inorganic ion transport and metabolism",
                   "Q":	"Secondary Structure",
                   "T":	"Signal Transduction",
                   "U":	"Intracellular trafficing and secretion",
                   "V":	"Defense mechanisms",
                   "Y":	"Nuclear structure",
                   "Z":	"Cytoskeleton",
                   "R":	"General Functional Prediction only",
                   "S":	"Function Unknown"}

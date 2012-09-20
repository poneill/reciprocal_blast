#!/usr/bin/env python

# Author: Sefa Kilic
# Modifications: PON
"""
This module contains functions for common tasks related to the
analysis of homologous coding sequences through reciprocal BLAST.

A typical workflow begins by constructing localized BLAST databases
for a list of organisms of interest.  Assume that this list is
represented as a standard Python list of strings called 'orgs'.  We
can construct the databases with the populate_dbs function:

populate_dbs(orgs)

which creates a database for each organism consisting of all amino
acid sequences contained in its FASTA file.  NB: this function assumes
the NCBI blast suite to be installed in the directory specified by
BLAST_BIN_PATH in the source.  The FASTA file is assumed to live in
data/$orgdir/, where $orgdir is the directory name in
ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ from which the genome was
downloaded.  $orgdir is assumed to contain both the genus and species
name, but is otherwise free to vary.  The resulting database is built
in /blast_db.

Next, we wish to BLAST each protein in each genome against every other
genome.  In practice, the massively parallel nature of this task makes
cluster computing a more natural and efficient setting for this task.
The primary resource for parallel BLASTing can be found on
tara.rs.umbc.edu in the folder erill_common/cub, where the Python
script setup_slurm.py will automatically setup and process BLAST jobs
on tara's compute nodes.  In the anecdotal experience of the author,
the task of computing pairwise reciprocal BLASTs for ~20 organisms can
be accelerated from ~25hrs to ~10min by parallelizing, although
compressing, scp-ing and uncompressing the results files adds a
roughly O(log(n)) constant to the wall-time.  Since a typical job runs
in approximately 2 min, a reasonable rule of thumb is to parallelize
if running more than 15 jobs.  Otherwise, one can run the jobs locally with:

reciprocal_blasts2(orgs)

In any case, assume that the desired BLAST results are present in the
directory blast_results.  The next step is to collate the one-way
blast hits into reciprocal results.  We do this by calling:

collate_reciprocal_blasts(orgs)

which generates, for each pair of organisms, a tab-separated file
whose lines consist of a pair of locus tags, one from each organism,
if the locus tags are reciprocal blast hits.

One often wishes to find proteins conserved between more than two
organisms.  We solve this problem efficiently by constructing a graph
whose vertices are all locus tags for all organisms.  The graph
contains an edge between two tags if their respective proteins have
been identified as reciprocal blast hits.  The problem of identifying
proteins conserved over multiple organisms is therefore reduced to the
problem of finding n-vertex cliques, where n = |orgs|.  Although the
general complexity of this problem is daunting, it is quite tractable
in practice.  We construct the graph as follows:

all_results = load_reciprocals(orgs)
g = make_graph(all_results)

and find n-cliques by calling the find_full_cliques function:

cliques = find_full_cliques(g,orgs)

Ultimately, we will wish to record the results on disk.  We construct
an annotation dictionary of the form d[org][locus_tag] = annotation by calling:

anno_dict = annotation_dict(orgs)

We then construct locus tag dictionary of the form
d[org][gi] = locus_tag by calling:

named_ltds = make_named_ltds(orgs)

Now we can generate annotations for a group of cliques by computing:

clique_dict = analyze_cliques(cliques,orgs,anno_dict,named_ltds)

of the form clique_dict[i] = [annotations], where i is the index of
the clique in cliques and the list of annotations contains the
annotation for each locus tag in the clique.

Finally, the results are stored with:

cliques2csv(cliques,orgs,named_ltds,clique_dict,filename,index)

which represents the cliques in a tab-separated file.  A typical line
of this file consists of:

org1_locus_tag org1_index ... orgn_locus_tag orgn_index annotation

and has a header row giving the organism names.
"""
from math import *
import os
import time
import sys
import re
import csv
import operator
import string
import scipy.stats
import networkx as nx
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2

from orgs import *
from utils import *
from biochem import *
import tai
import cog

from collections import Counter, defaultdict
from string import Template
from matplotlib import pyplot as plt
import pylab

pairwise2.MAX_ALIGNMENTS = 1
BLAST_BIN_PATH = "/home/poneill/ncbi-blast-2.2.26+/bin/"
DB_PATH = 'blastdb'
SITE_LENGTH = 82
ORG_PATH = "data"
epsilon = 10**-10

def org2nc_id(org):
    """Get NCid associated with org"""
    dir_name = org2dirname(org)
    file_name = head(os.listdir(os.path.join("data", dir_name)))
    (NCid, _) = os.path.splitext(file_name)
    return NCid

def org2dirname(org):
    """Get directory name associated with org"""
    dir_contents = os.listdir("data")
    dir_name = head([dc for dc in dir_contents if org_matches_dir(org, dc)])
    return dir_name
    
def makedb(sp):
    """Make database using fasta files of species"""
    os.system(BLAST_BIN_PATH + 'makeblastdb \
            -in %s -title %s -dbtype prot -out %s -logfile %s'
            % (sp['genome_file'], sp['name'], sp['db'], sp['name']+'.log'))


def blastn(seq, sp):
    """Run blastn on given seq and db"""
    qfile = 'query.tmp' #query file
    ofile = 'my_blast.xml' # output file
    with open(qfile, 'w') as f:
        f.write('%s' % seq)
    os.system(BLAST_BIN_PATH + 'blastn \
            -query %s -task blastn -db %s -out %s -outfmt 5'
            % (qfile, sp['db'], ofile))
    # parse output file
    result_handle = open(ofile)
    blast_record = NCBIXML.read(result_handle)
    # if no hit, return None
    if len(blast_record.alignments) == 0:
        return None
    # get first alignment
    alignment = blast_record.alignments[0]
    # make sure that hsps are sorted
    alignment.hsps.sort(key=lambda x: x.score, reverse=True)
    os.remove(qfile) # delete query file
    os.remove(ofile) # delete output file
    return alignment.hsps[0] # return best hit

def extend_site(genome, sbjct_start, sbjct_end, query_start, query_end):
    '''Extend site'''
    assert query_start < query_end
    if sbjct_start < sbjct_end: # strand:1
        if sbjct_start < 50 and sbjct_end > len(genome.seq)-50:
            return None
        ext_start = sbjct_start - (query_start-1) # extend start to left
        ext_end = sbjct_end + (SITE_LENGTH - query_end) # extend end to right
        retseq = genome.seq[ext_start-1 : ext_end]
    else: # strand:0
        if sbjct_start > len(genome.seq)-50 and sbjct_end < 50:
            return None
        ext_start = sbjct_start + (query_start-1)
        ext_end = sbjct_end - (SITE_LENGTH - query_end)
        retseq = genome.seq[ext_end-1 : ext_start].reverse_complement()
    return retseq
#


def populate_dbs(orgs):
    """given a list of organism names,  find directories,  grab files
    and construct databases"""
    org_dirs = os.listdir(ORG_PATH)
    file_ext = ".faa"
    for org in orgs:
        print "constructing db for ",  org
        org_dir = head([od for od in org_dirs if org_matches_dir(org, od)])
        full_org_dir = os.path.join(ORG_PATH, org_dir)
        fasta_file = head([f for f in os.listdir(full_org_dir)
                           if f.endswith(file_ext)])
        print "using fasta file", fasta_file
        sp = {'name':org, 
              'genome_file':os.path.join(full_org_dir, fasta_file), 
              'db': os.path.join(DB_PATH, org)}
        makedb(sp)
    
def reciprocal_blasts2(orgs):
    """Compute reciprocal blast hits for orgs"""
    org_dirs = os.listdir(ORG_PATH)
    results_contents = os.listdir('blast_results')
    file_ext = '.faa'
    for org1 in orgs:
        for org2 in orgs:
            print "starting on: ", org1, org2, "at", time.ctime()
            out_file = "results_%s_%s.txt" % (org1, org2)
            if org1 == org2:
                print "skipping", org1, org2
                continue 
            org_dir = head([od for od in org_dirs if org_matches_dir(org1, od)])
            full_org_dir = os.path.join(ORG_PATH, org_dir)
            fasta_file = head([f for f in os.listdir(full_org_dir)
                                            if f.endswith(file_ext)])
            full_fasta_file = os.path.join(full_org_dir, fasta_file)
            full_db_path = os.path.join(DB_PATH, org2)
            if out_file in results_contents:
                print "skipping", org1, org2
                continue
            os.system(os.path.join(BLAST_BIN_PATH, 
                      ('blastp -query %s -task blastp -db %s -out %s -outfmt 5'
                       % (full_fasta_file, full_db_path, 
                          os.path.join("blast_results", out_file)))))
            print "finished", org1, org2, "at", time.ctime()

def find_missing_results_files(orgs):
    """Print a list of expected files not found in blast_results"""
    dir_contents = os.listdir('blast_results')
    for org1 in (orgs):
        for org2 in orgs:
            if org1 == org2:
                continue
            filename1 = ("results_%s_%s.txt" % (org1, org2))
            filename2 = ("results_%s_%s.txt" % (org2, org1))
            for fn in [filename1, filename2]:
                if not fn in dir_contents:
                    print fn
                    

def collate_reciprocal_blasts(orgs):
    full_path = lambda f: os.path.join("blast_results", f)
    out_path = lambda f: os.path.join("reciprocal_results", f)
    locus_tag_dicts = {}
    dir_contents = os.listdir('reciprocal_results')
    n = len(orgs)
    total = n * (n - 1)/2
    for i, (org1,org2) in enumerate(choose2(orgs)): 
            out_filename = "reciprocals_%s_%s.csv" % (org1, org2)
            out_file = out_path(out_filename)
            if out_filename in dir_contents:
                print "found %s, moving on" % out_filename
                continue
            print "collating recips for %s,%s at %s" % (org1, org2, time.ctime())
            if not org1 in locus_tag_dicts:
                locus_tag_dicts[org1] = make_locus_tag_dict(get_genome_filename(org1,"gbk"))
            if not org2 in locus_tag_dicts:
                locus_tag_dicts[org2] = make_locus_tag_dict(get_genome_filename(org2,"gbk"))
            org1_ltd, org2_ltd = locus_tag_dicts[org1], locus_tag_dicts[org2]
            filename1 = full_path("results_%s_%s.txt" % (org1, org2))
            filename2 = full_path("results_%s_%s.txt" % (org2, org1))
            print "parsing results for %s, %s" % (org1, org2)
            print "using:", filename1, filename2
            hits1 = parse_results(filename1, org1_ltd, org2_ltd)
            print "parsing results for %s, %s" % (org2, org1)
            hits2 = parse_results(filename2, org2_ltd, org1_ltd)
            print "finding reciprocals"
            reciprocals = find_reciprocals(hits1, hits2)
            print "writing file %s out of %s" % (i,total)
            with open(out_file, 'w') as f:
                text = "\n".join(["\t".join(tup) for tup in reciprocals]) + "\n"
                f.write(text)

def make_locus_tag_dict(gbk_filename):
    print gbk_filename
    genome = SeqIO.read(gbk_filename, 'genbank')
    print "finding CDSs"
    CDSs = ([feature for feature in genome.features if feature.type == 'CDS'])
    print "findings gis"
    gis = [cds.qualifiers['db_xref'][0][3:] for cds in CDSs]
    print "finding locus tags"
    locus_tags = [cds.qualifiers['locus_tag'][0] for cds in CDSs]
    return {gi:locus_tag for (gi, locus_tag) in zip(gis, locus_tags)}

def make_named_ltds(orgs):
    return {org:make_locus_tag_dict(get_genome_filename(org,"gbk")) for org in orgs}

master_lt2o_dict = {}

def locus_tag2org(tag, named_ltds):
    if not tag in master_lt2o_dict:
        master_lt2o_dict[tag] = locus_tag2org_inner(tag, named_ltds)
    return master_lt2o_dict[tag]

def locus_tag2org_inner(tag, named_ltds):
    prefix = re.match(r"^[A-Za-z]+", tag).group(0)
    for org in named_ltds:
        first_tag = named_ltds[org].values()[0]
        if prefix == locus_tag_prefix(first_tag):
            return org
    #otherwise, search thoroughly
    for org in named_ltds:
        if any(value.startswith(tag) for value in named_ltds[org].values()):
            return org
    raise Exception("Couldn't find locus tag %s" % tag)

def get_genome_filename(org_name,ext):
    dir_contents = os.listdir('data')
    dir_name = os.path.join("data", 
                            head(filter(lambda d: org_matches_dir(org_name, d), 
                                        dir_contents)))
    fn = head(filter(lambda f: f.endswith('.' + ext), os.listdir(dir_name)))
    return os.path.join(dir_name, fn)
    
def parse_results(filename, query_locus_dict, target_locus_dict):
    """Accept a file containing the results of blasting query_org against
    target_org and return a dictionary of the form {protein in query_org: first
    blast hit in target_org}"""
    hits = {}
    query_def_pattern = re.compile(r"""<Iteration_query-def>
                                       gi\|(.*)\|ref.*
                                       </Iteration_query-def>""", re.X)
    hit_def_pattern = re.compile(r"<Hit_def>gi\|(.*?)\|ref.*?</Hit_def>")
    evalue_pattern =  re.compile(r"<Hsp_evalue>(.*?)</Hsp_evalue>")
    cutoff = 1e-10
    found_best_hit = False
    with open(filename) as f:
        for line in f:
            query_def_match = query_def_pattern.search(line)
            hit_def_match = hit_def_pattern.search(line)
            evalue_match = evalue_pattern.search(line)
            if query_def_match:
                query_name = query_def_match.group(1)
                found_best_hit = False
#                print "found query_name: %s" % query_name
            elif hit_def_match and not found_best_hit: 
                hit_name = hit_def_match.group(1)
#                print "found hit_name: %s" % hit_name
            elif evalue_match and not found_best_hit:
                evalue = float(evalue_match.group(1))
 #               print "found evalue: %s" % str(evalue)
                if evalue < cutoff:
                    query_locus = query_locus_dict[query_name]
                    target_locus = target_locus_dict[hit_name]
                    hits[query_locus] = (target_locus, cutoff)
                    found_best_hit = True 
    return hits

def find_reciprocals(d1, d2):
    """Given two dictionaries d1 and d1 collecting the matches
    for org1 against org2 and vice versa, return a list of tuples
    [(prot1, prot2)] such that prot1:prot2 is in hits1 and prot2:prot1
    is in hits2"""
    hits1 = {k:v[0] for k, v in d1.iteritems()}
    hits2 = {k:v[0] for k, v in d2.iteritems()}
    reciprocals = [(h1, h2) for h1 in hits1 for h2 in hits2
                   if h2 in hits1[h1] and h1 in hits2[h2]]
    return reciprocals

def load_reciprocals(orgs):
    """Read the reciprocal blast results and return a dictionary of form
    {results_filename:[(geneA, geneB)]}"""
#    dir_contents = os.listdir("reciprocal_results")
    def get(org1,org2):
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        lines = open(os.path.join("reciprocal_results", filename)).readlines()
        results = [re.search(r"(.*)\t(.*)\n", line).groups() for line in lines]
        return results
#    dir_contents = ["reciprocals_%s_%s.csv" % tup for tup in choose2(orgs)]
    all_results = {}
    for org1,org2 in choose2(orgs):
        #for filename in dir_contents:
        print org1,org2
        try:
            results = get(org1,org2)
        except IOError:
            reversed_results = get(org2,org1)
            results = map(reverse,reversed_results)
        filename = "reciprocals_%s_%s.csv" % (org1,org2)
        all_results[filename] = results
    return all_results

def make_graph(all_results):
    """Given the a dictionary of form
    {results_filename:[(geneA, geneB)]}, construct a graph whose
    vertices are genes and whose edges are reciprocal blast hits"""
    G = nx.Graph()
    vertices = set([tup[0] for result in all_results.values()
                    for tup in result])
    for v in vertices:
        G.add_node(v)
    for result in all_results.values():
        for tup in result:
            G.add_edge(*tup)
    return G

def get_genome(org):
    return SeqIO.read(get_genome_filename(org,"gbk"), 'genbank')

def gene_name2locus_tag(gene, genome):
    print gene
    cdss = ([feature for feature in genome.features
                 if feature.type == 'CDS'])
    for cds in cdss:
        if ('gene' in cds.qualifiers
            and gene in cds.qualifiers['gene']):
            return cds.qualifiers['locus_tag'][0]
    
def annotation_dict(orgs):
    """return a dictionary of form anno_dict[org][locus_tag] == annotation"""
    anno_dict = {}
    for org in orgs:
        print org
        anno_dict[org] = {}
        genome = get_genome(org)
        CDSs = ([feature for feature in genome.features
                 if feature.type == 'CDS'])
        for cds in CDSs:
            description = (head(cds.qualifiers['product'])
                           if 'product' in cds.qualifiers
                           else 'None')
            anno_dict[org][head(cds.qualifiers['locus_tag'])] = description
    return anno_dict

def analyze_cliques(cliques, orgs, anno_dict, named_ltds):
#    clique_dict = {i:[] for i, cl in enumerate(cliques)}
    clique_dict = {}
    for i, cl in enumerate(cliques):
        clique_dict[i] = [anno_dict[locus_tag2org(tag, named_ltds)][tag]
                          for tag in cl]
    return clique_dict


def get_cdss(genome):
    return ([feature for feature in genome.features if feature.type == 'CDS'])

def cdss_from_org(org):
        return get_cdss(get_genome(org))
    
        
def maj_annotation(annotations):
    """Given a list of annotations (i.e. clique_dict[i]), return the
    most common"""
    return Counter(annotations).most_common(1)[0][0]

def find_full_cliques(G, orgs):
    return [cl for cl in nx.find_cliques(G) if len(cl) == len(orgs)]

def find_specific_cliques(G, orgs,named_ltds):
    return [cl for cl in nx.find_cliques(G) if clique_in_orgs(cl)]

def locus_tag_prefix(tag):
    return re.match("[A-Za-z]+", tag).group()

def locus_tag_prefixes(org, named_ltds):
    return set(locus_tag_prefix(tag) for tag in named_ltds[org].itervalues())

def cliques2csv(cliques, orgs, named_ltds, clique_dict, filename,index):
    print "finding all_prefixes"
    all_prefixes = set([locus_tag_prefix(tag) for clique in cliques
                        for tag in clique])
    print "sorting prefixes"
    sorted_prefixes = sort_on(all_prefixes, orgs, 
                              lambda tag: locus_tag2org(tag, named_ltds))
#    row_header = map(locus_tag_prefix, sorted_tags)
    print "indices"
    indices = index_dicts(orgs,index)
    print "sorted cliques"
    sorted_cliques = [sort_on(clique, sorted_prefixes, locus_tag_prefix)
                      for clique in cliques]
    print "rows"
    print orgs
    rows = [zipWith(lambda tag, org:(tag, indices[org][tag]), sc, orgs)
            for sc in sorted_cliques]
    annotations = [maj_annotation(clique_dict[i]) for i in range(len(cliques))]
    formatted_indices = ["\t".join(map(str, sum(row, ()))) for row in rows]
    formatted_rows = "\n".join(zipWith(lambda row, anno: row + '\t' + anno, 
                                       formatted_indices, annotations))
    formatted_header = ("\t".join(sum(zip(orgs, 
                                        ["" for i in range(len(orgs))]), ())) +
                        "\tAnnotation\n")
    with open(filename, 'w') as f:
        f.write(formatted_header)
        f.write(formatted_rows)

def generate_pairwise_spreadsheets(orgs, named_ltds, anno_dict):
    for org1, org2 in choose2(orgs):
        print "starting", org1, org2
        filename = "pairwise_correlations_%s_%s.csv" % (org1, org2)
        full_dir = os.path.join("clique_csvs", "pairwise_correlations")
        full_path = os.path.join(full_dir, filename)
        if filename in os.listdir(full_dir):
            print "found", filename, "skipping"
            continue
        pair = [org1, org2]
        print "making graph"
        g = make_graph(load_reciprocals(pair))
        print "finding cliques"
        cliques = find_full_cliques(g, pair)
        print "making clique_dict"
        clique_dict = analyze_cliques(cliques, pair, anno_dict, named_ltds)
        print "writing csv"
        cliques2csv(cliques, pair, named_ltds, clique_dict, full_path,"nRCA")

def read_genome_info_txt(filename):
    lines = open(filename).readlines()
    return [line.split("\t") for line in lines]

def find_indices_over_cliques(genome_info, cliques):
    return [float(line[0]) for clique in cliques
            for line in genome_info if line[1] in clique]

def nc_template(index):
    return "%s-1.0_rcc_RCA" if index in ["nRCA","RCA"]  else "%s-1.0-tAI"
    
def correlate_index_over_cliques(cliques, orgs,index):
    """Return a list of lists containing, for each organism, the index value of
    that organism's contribution to each clique"""
    genome_info_filenames = [genome_info_filename(org,index) for org in orgs]
    return [find_indices_over_cliques(read_genome_info_txt(gif), cliques)
            for gif in verbose_gen(genome_info_filenames)]

def genome_info_filename(org,index):
    "return the filepath for genome_info.txt belonging to org"
    if index in ["nRCA","RCA,""CAI"]:
        d = "%s-1.0_rcc_%s" % (org2nc_id(org), index)
    else:
        d = "%s-%s" % (org2nc_id(org), index)
    return os.path.join("index_results", d, "genome_info.txt")
    
def index_dict(org,index):
    """return a dictionary of the form {locus_tag:index} over org"""
    print org
    genome_info = read_genome_info_txt(genome_info_filename(org,index))
    return {row[1]:float(row[0]) for row in genome_info}

def index_dicts(orgs,index):
    return {org:index_dict(org,index) for org in orgs}

def check_correlation_in_pairwise_cliques():
    dirname = os.path.join("clique_csvs", "pairwise_correlations")
    fs = [f for f in os.listdir(dirname)
          if f.startswith("pairwise_correlations") and f.endswith(".csv")]
    results = []
    for f in fs:
        filename = os.path.join(dirname, f)
        lines = [line for line in csv.reader(open(filename), delimiter='\t')]
        xs = [float(line[1]) for line in lines[1:]]
        ys = [float(line[3]) for line in lines[1:]]
        results.append((f, scipy.stats.pearsonr(xs, ys)))
    return sorted(results, key=lambda result: result[1][0], reverse=True)

conditional_log_transform=True
def exp_dict(filename):
#    print filename
    MAX_LOG_VALUE = 1000000
    lines = [line for line in csv.reader(open(filename), delimiter=',')]
    d = {}
    for line in lines:
        locus_tag = line[0]
        data = map(float, line[1:])
        if not locus_tag in d:
            d[locus_tag] = data
        else:
            d[locus_tag].extend(data)
    if conditional_log_transform:
        vals = sum(d.values(),[])
        _max = max(vals)
        _min = abs(min(vals))
        def safe_log(x):
#            print x, _min
            return log(x + _min + 1)
        if _max > MAX_LOG_VALUE:
            d = {k:map(safe_log,d[k]) for k in d}
    return d

def exp_dicts(orgs):
    dicts = {}
    for org in orgs:
        try:
            dicts[org] = exp_dict(exp_filename(org))
            print "exp data for:",org
        except:
            print "no exp data for:",org
            continue
    return dicts

def pairwise_correlations(org1, org2):
    dirname = os.path.join("clique_csvs", "pairwise_correlations")
    f = "pairwise_correlations_%s_%s.csv" % (org1, org2)
    try:
        file_handle = open(os.path.join(dirname, f))
        rev = False
    except IOError:
        f = "pairwise_correlations_%s_%s.csv" % (org2, org1)
        file_handle = open(os.path.join(dirname, f))
        rev = True
    lines = [line for line in csv.reader(file_handle, delimiter='\t')]
    if not rev:
        pairs =  [(line[0], float(line[1]), line[2], float(line[3]))
                  for line in lines[1:]]
    else:
        pairs =  [(line[2], float(line[3]),line[0], float(line[1]))
                  for line in lines[1:]]
    return pairs

def exp_filename(org):
    return os.path.join("exp_csvs", org+"_exp.csv")
    
def expression_report_locals(org1, org2, p=0, upper=True):
    correlations = pairwise_correlations(org1, org2)
    (raw_exp_dict1, raw_exp_dict2) = map(lambda org:exp_dict(exp_filename(org)), 
                                        [org1, org2])
    exp_dict1_cutoff = percentile(p, raw_exp_dict1.values(), upper=upper)
    exp_dict2_cutoff = percentile(p, raw_exp_dict2.values(), upper=upper)
    meets_cutoff = operator.ge if upper else operator.le
    exp_dict1 = {k:raw_exp_dict1[k] for k in raw_exp_dict1
                 if meets_cutoff(raw_exp_dict1[k], exp_dict1_cutoff)}
    exp_dict2 = {k:raw_exp_dict2[k] for k in raw_exp_dict2
                 if meets_cutoff(raw_exp_dict2[k], exp_dict2_cutoff)}
    exp_data_for = lambda line: (line[0] in exp_dict1
                                 and line[2] in exp_dict2)
    exp_pairs = [map(mean, (exp_dict1[line[0]], exp_dict2[line[2]]))
                 for line in correlations
                 if exp_data_for(line)]
    (org1_exp, org2_exp) = zip(*exp_pairs) if len(exp_pairs) else ([], [])
    num_cliques = len(correlations)
    index_pairs = [(line[1], line[3]) for line in correlations
                  if exp_data_for(line)]
    (org1_indices, org2_indices) = zip(*index_pairs)
    org1_index_vs_org2_index = pearson(org1_indices, org2_indices)
    conserved_in_org1 = {line[0]:line[1] for line in correlations}
    conserved_in_org2 = {line[2]:line[3] for line in correlations}
    org1_exp_vs_org2_exp = pearson(org1_exp, org2_exp)
    num_expressed = len(exp_pairs)
    org1_index_vs_org1_exp = pearson(org1_indices, org1_exp)
    org2_index_vs_org2_exp = pearson(org2_indices, org2_exp)
    org1_index_vs_org2_exp = pearson(org1_indices, org2_exp)
    org1_exp_vs_org2_index = pearson(org2_indices, org1_exp)
    return locals() #sketchy!

def expression_vs_percentile_single(org, ps, upper=True):
    rs = [pearson(*expression_vs_percentile_single_data(org, p, upper))
          for p in ps]
    plt.plot(ps, rs)
    plt.show()

def expression_vs_percentile_single_data(org, index,p, upper=True):
    indices = index_dict(org,index)
    exp_path = lambda f: os.path.join("exp_csvs", f+"_exp.csv")
    exps = exp_dict(exp_path(org))
    cutoff = percentile(p, exps.values(), upper=upper)
    indices_exps = [(indices[tag], mean(exps[tag])) for tag in indices
                  if tag in exps and exps[tag] >= cutoff]
    final_indices, final_exps = zip(*indices_exps)
    return final_indices, final_exps
    
def expression_vs_percentile(org1, org2, ps, upper=True):
    vals = [expression_report_locals(org1, org2, p, upper=upper) for p in ps]
    variables = ["org1_index_vs_org2_index", 
                 "org1_exp_vs_org2_exp", 
                 "num_expressed", 
                 "org1_index_vs_org1_exp", 
                 "org2_index_vs_org2_exp", 
                 "org1_index_vs_org2_exp", 
                 "org1_exp_vs_org2_index"]
    for i, var in enumerate(variables):
        exec("%s_list=[v['%s'] for v in vals]" % (var, var))
        exec("p%s, =plt.plot(ps, %s_list)" % (i, var))
    plt.xlabel("Percentile cutoff")
    plt.ylabel("Pearson correlation")
    plt.legend([eval("p%s" % i) for i in range(len(variables))], variables, loc=3)
    plt.show()
    
def expression_report(org1, org2, p=0):
    d = expression_report_locals(org1, org2, p)
    template = Template("""Comparing: $org1, $org2
Found $num_cliques conserved proteins
Using $num_expressed found in both expression datasets
org1 index vs org2 index:\t\t$org1_index_vs_org2_index
org1 exp vs org2 exp:\t\t$org1_exp_vs_org2_exp
org1 self-expression:\t\t$org1_index_vs_org1_exp
org2 self-expression:\t\t$org2_index_vs_org2_exp
org1 index vs. org2 expression\t$org1_index_vs_org2_exp
org2 index vs. org1 expression\t$org1_exp_vs_org2_index
NB: Pearson r-values, not r^2
""")
    return template.substitute(d)

def all_expression_reports(orgs):
    for org1 in orgs:
        for org2 in orgs:
            try:
                print expression_report(org1, org2)
            except:
                continue

def expression_report_self_cross_comparison(orgs):
    selfs = defaultdict(list)
    crosses = defaultdict(list)
    for org1 in orgs:
        for org2 in orgs:
            print org1, org2
            try:
                print "trying"
                d = expression_report_locals(org1, org2)
                print "got dict for ", org1, org2
                selfs[org1].append(d["org%s_index_vs_org%s_exp" % (1, 1)])
                selfs[org2].append(d["org%s_index_vs_org%s_exp" % (2, 2)])
                print "extended selves"
                crosses[(org1, org2)].append(d["org1_index_vs_org2_exp"])
                crosses[(org2, org1)].append(d["org1_exp_vs_org2_index"])
                print "extended crosses"
            except:
                continue
    return selfs, crosses

def expression_report_graph(orgs):
    g = nx.MultiGraph()
    for org1 in orgs:
        for org2 in orgs:
            try:
                d = expression_report_locals(org1, org2)
                locals().update(d)
                org1_index = org1+"_index"
                org1_exp = org1+"_exp"
                org2_index = org2+"_index"
                org2_exp = org2+"_exp"
                if not org1 in g:
                    print "adding node", org1_index
                    g.add_node(org1_index)
                    print "adding node", org1_exp
                    g.add_node(org1_exp)
                if not org1 in g:
                    print "adding node", org2_index
                    g.add_node(org2_index)
                    print "adding node", org2_exp
                    g.add_node(org2_exp)
                var_pairs = choose2(["".join(pair)
                                     for pair in
                                     cart_product(["org1_", "org2_"], 
                                                  ["index", "exp"])])
                for var1, var2 in var_pairs:
                    label = eval("%s_vs_%s" % (var1, var2))
                    node1 = var1.replace("org1", org1).replace("org2", org2)
                    node2 = var2.replace("org2", org2).replace("org1", org1)
                    print "adding edge", node1, node2, label
                    g.add_edge(node1, node2, weight=label)
            except:
                pass
    return g

def collate_cdss_with_expression(org):
    genome = get_genome(org)
    cdss = get_cdss(genome)
    exp_d = exp_dict(exp_filename(org))
    return {cds.extract(genome).seq:exp_d[head(cds.qualifiers['locus_tag'])]
            for cds in verbose_gen(cdss)
            if head(cds.qualifiers['locus_tag']) in exp_d}

def tag2aas(tag, cdss, genome):
    print tag
    cds = head([cds for cds in cdss if tag in cds.qualifiers['locus_tag']])
    #drop symbol for stop codon
    return cds.extract(genome).seq.translate().tostring()[:-1] 

def analyze_index_vs_conservation(org1, org2, n=None, length_normalize=False):
    print "reading genomes"
    genome1, genome2 = map(get_genome, [org1, org2])
    print "finding CDSs"
    cdss1, cdss2 = map(get_cdss, [genome1, genome2])
    print "reading index dicts"
    index_dict1, index_dict2 = map(lambda org:index_dict(org,index), [org1, org2])
    print "reading correlations"
    correlations = pairwise_correlations(org1, org2)[:n]
    tag2aas1 = lambda tag: tag2aas(tag, cdss1, genome1)
    tag2aas2 = lambda tag: tag2aas(tag, cdss2, genome2)
    def nw_score(seq1, seq2):
        denom = float(seq1 + seq2)/2 if length_normalize else 1
        return head(pairwise2.align.globaldx(seq1, seq2, BLOSUM62))[2]/denom
    print "computing alignments"
    scores = [nw_score(tag2aas1(line[0]), tag2aas2(line[2]))
              for line in verbose_gen(correlations)]
    return [(line[0], index_dict1[line[0]], line[2], index_dict2[line[2]], scores[i])
            for i, line in enumerate(correlations)]

def test_gc(k, n):
    aa_alphabet = 'ABCDEFGHIKLMNPQRSTVWXYZ'
    for i in range(n):
        print i
        seq1 = "".join([random.choice(aa_alphabet) for i in range(k)])
        seq2 = "".join([random.choice(aa_alphabet) for i in range(k)])
        pairwise2.align.globaldx(seq1, seq2, BLOSUM62)

def parse_iterations_verbose(org):
    """parse iterationsVerbose.txt for org and return three
    dictionaries, the first of the form {codon:genome_wide_codon
    freq}, the second of the form {codon:refset_codon_freq}, and
    the third of the form {codon:w-value}"""
    ncid = org2nc_id(org)
    filename = os.path.join("index_results", ncid+"-1.0_rcc_nRCA", 
                            "iterationVerbose.txt")
    with open(filename) as f:
        lines = [filter(iota, line) for line in csv.reader(f, delimiter="\t")]
    genomic_freq_lines = lines[3:11] #magic numbers by inspection
    refset_freq_lines = lines[len(lines)-36:len(lines)-28]
    w_lines = lines[len(lines)-9:]
    def lines2dict(lines):
        d = {}
        for line in lines:
            for field in line:
                codon, val_string = field.split(':')
                val = float(val_string)
                d[codon] = val
        return d
    results = map(lines2dict, [genomic_freq_lines, refset_freq_lines, w_lines])
    return results

def gc_content(seq):
    return len([n for n in seq if n in "gcGC"])/float(len(seq))

def codon_dict2gc(refset_freqs):
    return sum([gc_content(codon) * refset_freqs[codon] for codon in refset_freqs])
        
def org2refset_codons(org):
    refset = 1
    return parse_iterations_verbose(org)[refset]

def org2genomic_codons(org):
    genomic = 0
    return parse_iterations_verbose(org)[genomic]

def org2w_codons(org):
    w = 2
    return parse_iterations_verbose(org)[w]

def org2refset_enrichment(org):
    epsilon = 10**-10
    org_refset_codons = org2refset_codons(org)
    org_genomic_codons = org2genomic_codons(org)
    return {codon:org_refset_codons[codon]/(org_genomic_codons[codon] + epsilon)
            for codon in codons}
    
def extract_boxes(codon_dict,box_length,break_down_six_boxes = False):
    tt = translation_table if not break_down_six_boxes else six_box_translation_table
    n_box_aas = [aa for aa in tt if len(tt[aa]) == box_length]
    n_box_codons = concat(map(lambda codon: tt[codon],n_box_aas))
    n_box_dict = {codon:codon_dict[codon] for codon in codon_dict
                     if codon in n_box_codons}
    return n_box_dict


def separate_nbox_into_t_ending(n_box_dict,break_down_six_boxes,n):
    """Given a dictionary of n_box codons (see
    extract_boxes),return a dictionary giving frequency of t_ending
    codons, non_t_ending codons.  Format: {amino
    acid:(t_ending_freq,non_t_ending_freq)}"""
    tt = translation_table if not break_down_six_boxes else six_box_translation_table
    def freq_tuple(aa):
        codons = tt[aa]
        t_ending_freqs = [n_box_dict[codon] for codon in codons
                          if codon.endswith('t')]
        non_t_ending_freqs = [n_box_dict[codon] for codon in codons
                               if not codon.endswith('t')]
        return (t_ending_freqs,non_t_ending_freqs)
    return {aa:freq_tuple(aa) for aa in tt
            if len(tt[aa]) == n}

def n_box_dict_separated(orgs,break_down_six_boxes,n):
    """Given a list of orgs, return a dictionary of the form:
    {org:{aa:([t_ending_frequency],[non_t_ending_frequencies])}}"""
    n_box_dict = {org:extract_boxes(org2refset_enrichment(org),n,break_down_six_boxes) for org in orgs}
    fbd_separated = {org:separate_nbox_into_t_ending(n_box_dict[org],break_down_six_boxes,n) for org in orgs}
    return fbd_separated
    
def n_box_dict_separated2csv(n_box_dict_separated,filename):
    two_boxing = len(n_box_dict_separated.values()[0].values())
    aas = head(n_box_dict_separated.values()).keys()
    break_down_six_boxes = True if len(aas) in [8,11] else False
    tt = translation_table if not break_down_six_boxes else six_box_translation_table
    first_header = "," + ''.join(["%s,,,," % aa for aa in aas])
    def second_header_for_aa(aa):
        codons = tt[aa]
        t_ending_codon = [codon for codon in codons
                          if codon.endswith('t')]
        non_t_ending_codon = [codon for codon in codons
                              if not codon.endswith('t')]
        return ",".join(t_ending_codon + non_t_ending_codon)
    second_header = "," + ",".join([second_header_for_aa(aa) for aa in aas])
    rows = "\n".join([",".join([org] + map(str,concat([concat(n_box_dict_separated[org][aa])
                            for aa in n_box_dict_separated[org]])))
            for org in n_box_dict_separated])
    with open(filename,'w') as f:
        f.write(first_header + "\n")
        f.write(second_header + "\n")
        f.write(rows + "\n")
        

def generate_refset_gc_content(group_name):
    orgs = eval(group_name)
    genomic = 0
    refset = 1
    gc_contents = [codon_dict2gc(parse_iterations_verbose(org)[refset]) for org in orgs]
    gc_contents_pseudos = [codon_dict2gc(parse_iterations_verbose(org)[refset]) for org in pseudos]
    gc_contents_psychros = [codon_dict2gc(parse_iterations_verbose(org)[refset]) for org in psychros]
    group_mean = mean(gc_contents)
    group_sd = sd(gc_contents)
    group_mean_psychros = mean(gc_contents_psychros)
    group_sd_psychros = sd(gc_contents_psychros)
    group_mean_pseudos = mean(gc_contents_pseudos)
    group_sd_pseudos = sd(gc_contents_pseudos)
    data2csv([gc_contents,["group mean",group_mean],["group sd",group_sd],
              ["group mean (pseudos)",group_mean_pseudos],
              ["group sd (pseudos)",group_sd_pseudos],
              ["group mean (psychros)",group_mean_psychros],
              ["group sd (psychros)",group_sd_psychros]],
             filename="refset_gcs/refset_gc_content_%s.csv" % group_name,header=orgs)


def all_pairwise_indices(orgs):
    return sorted([(all_correlations(org1, org2)[0], org1, org2)
                   for (org1, org2) in verbose_gen(choose2(orgs))], 
                  key = lambda tup: tup[0])

def get_codon_dicts(orgs,genomic=True):
    """Return a dictionary of codon usage biases index by org"""
    f = genomic_codon_freqs if genomic else refset_codon_freqs
    return {org:f(org) for org in orgs}

#helper functions...
def genomic_codon_freqs(org):
    return parse_iterations_verbose(org)[0]

def refset_codon_freqs(org):
    return parse_iterations_verbose(org)[1]

# def genomic_codon_counts(org):
#     genome = get_genome(org)
#     genomic_cdss = cdss_from_org(org)
#     counts = defaultdict(int)
#     for cds in genomic_cdss:
#         for codon in group_codons

def trna_codon_freqs(org):
    trnas = trna_counts(org)
    total = float(sum(trnas.values()))
    return {codon:(trnas[codon]/total if codon in trnas else 0)
            for codon in codons}

def trna_counts(org):
    return parse_trnas(org2trna_filename(org))
            
def all_correlations(org1, org2, n=None,correlation_func=pearson):
    """Compare the correlation between index values to the correlation
    between RCA weights"""
#    genome1, genome2 = map(get_genome, [org1, org2])
#    cdss1, cdss2 = map(get_cdss, [genome1, genome2])
#    cdss1, cdss2 = map(cdss_from_org, [org1,org2])
    index_dict1, index_dict2 = map(lambda org:index_dict(org,"nRCA"), [org1, org2])
    correlations = pairwise_correlations(org1, org2)[:n]
    index_table1 = [index_dict1[line[0]] for line in correlations]
    index_table2 = [index_dict2[line[2]] for line in correlations]
    gen_freqs1, refset_freqs1, ws1 = parse_iterations_verbose(org1)
    gen_freqs2, refset_freqs2, ws2 = parse_iterations_verbose(org2)
    gen_codon_table1 = [gen_freqs1[codon] for codon in codons]
    gen_codon_table2 = [gen_freqs2[codon] for codon in codons]
    refset_codon_table1 = [refset_freqs1[codon] for codon in codons]
    refset_codon_table2 = [refset_freqs2[codon] for codon in codons]
    w_table1 = [ws1[codon] for codon in codons]
    w_table2 = [ws2[codon] for codon in codons]
    results = map(lambda (x,y): correlation_func(x,y),[(index_table1, index_table2),
                                    (gen_codon_table1, gen_codon_table2),
                                    (refset_codon_table1, refset_codon_table2),
                                    (w_table1, w_table2)])
    index_corr, gen_codon_corr, refset_codon_corr, w_corr = results
    return (org1,org_group(org1),org2,org_group(org2),index_corr, gen_codon_corr, refset_codon_corr, w_corr)

def all_correlations_all_orgs(orgs,cf=pearson):
    d = all_correlations_all_orgs_dict(orgs,cf)
    return dictmap(d,choose2(orgs))
    
def all_correlations_all_orgs_dict(orgs,cf=pearson):
    return {(org1,org2):all_correlations(org1,org2,correlation_func = cf)
            for (org1,org2) in verbose_gen(choose2(orgs))}


def parse_ribo_str_crit(filename):
    with open(filename) as f:
        lines = f.readlines()
    results = [float(line.split(": ")[1]) for line in lines]
    (ribosomal_crit, strength_crit, content_crit) = results
    return (ribosomal_crit, strength_crit, content_crit)

def init_analysis_summary(orgs, outfile="initial_analysis_summary.csv"):
    def get_file(org, index):
        return os.path.join("index_results", org2nc_id(org) + "-1.0_rcc_" + index, 
                            "ribo_strength_criterion.txt")
    header = ("Organism, " +
              "RCA Ribosomal Crit, RCA Strength Crit, RCA Content Crit, "
              + "CAI Ribosomal Crit, CAI Strength Crit, CAI Content Crit")
    lines = [", ".join(map(str, (org, ) +
              parse_ribo_str_crit(get_file(org, "nRCA")) +
              parse_ribo_str_crit(get_file(org, "CAI")))) for org in orgs]
    with open(outfile, 'w') as f:
        f.write(header)
        f.write("\n".join(lines))

def pairwise_summary(orgs, outfile, f=pearson):
    """Summarize the pairwise_correlation files by providing csv
    tables that report the Pearson and Spearman correlations, and
    number of genes, in each pairwise_correlation file. """
    header = ", "+", ".join(orgs)
    def get_lines(org1, org2):
        filename = os.path.join("clique_csvs", "pairwise_correlations", 
                                "pairwise_correlations_%s_%s.csv" % (org1, org2))
        with open(filename) as fn:
            #skim off the header
            lines = [line for line in csv.reader(fn, delimiter="\t")][1:] 
        return lines 
    def extract(org1, org2, content):
        if org1 == org2:
            return float(content == pearson or content == spearman)
        try:
            lines = get_lines(org1, org2)
        except IOError:
            try:
                lines = get_lines(org2, org1)
            except IOError:
                raise IOError(org1, org2)
        vals = [map(float, (line[1], line[3])) for line in lines]
        org1_vals, org2_vals = zip(*vals)
        return f(org1_vals, org2_vals)
    
    outlines = [[org1] + [extract(org1, org2, f)
                       for org2 in orgs] for org1 in orgs]
    with open(outfile, 'w') as outf:
        outf.write(header+"\n")
        outf.write("\n".join([", ".join(map(str, line)) for line in outlines]))

def escape_spaces(path):
    return path.replace(' ', '\ ')
def run_trna_scan(orgs):
    home_dir = "/home/poneill"
    trna_scan_path = os.path.join(home_dir, "bin")
    relative_home = os.getcwd()
    esc_relative_home = escape_spaces(relative_home)
    home_path = lambda fn: os.path.join(esc_relative_home, fn)
    trnas_path = lambda org: os.path.join(esc_relative_home, 
                                          "trnas", org + "_trnas")
    outfiles = os.listdir("trnas")
    for org in orgs:
        print org, time.ctime()
        print os.getcwd()
        org_dir = os.path.join("data", org2dirname(org))
        fn = head(os.listdir(org_dir), lambda f: f.endswith(".fna"))
        if org + "_trnas" in outfiles:
            print "skipping:", org
            continue
        fna_filename = home_path(os.path.join(org_dir, fn)) 
        outfile = trnas_path(org)
        os.chdir(trna_scan_path)
        command = "tRNAscan-SE -BH %s -f %s" % (fna_filename, outfile)
        print command
        os.system(command)
        os.chdir(relative_home)

def reverse_complement(seq):
    table = string.maketrans("ACGTacgt", 'TGCAtgca')
    return seq.translate(table)[::-1]

def org2trna_filename(org):
        return os.path.join("trnas", "%s_trnas" % org)

def trna_dict2trna_list(trna_dict,codon_order = "lexicographic"):
    lookup = {"lexicographic":codons,
              "sorted":sorted_codons,
              "table":table_codons}
    cs = lookup[codon_order]
    return [trna_dict[c] for c in cs]

def write_tai_data_for_R_script(org,overwrite=False):
    print org
    w_list = w_list_for_R_script(org)
    orf_freqs = orf_frequencies_for_R_script(org)
    data2csv([w_list],os.path.join("tai_script",org+"_ws.csv"),overwrite=overwrite)
    data2csv(orf_freqs,os.path.join("tai_script",org+"_orf_freqs.csv"),
             overwrite=overwrite)
                  
def w_list_for_R_script(org):
    return trna_dict2trna_list(genomic_codon_freqs(org),codon_order = "table")

def process_output_for_R_script(org,anno_dict,overwrite=False):
    print(org)
    with open(os.path.join("tai_script",org+"_tais.csv")) as f:
        lines = [line for line in csv.reader(f)][1:]
    tais = [float(head(line)) for line in lines]
    org_locus_tags = locus_tags(org)
    annos = [anno_dict[org][tag] for tag in org_locus_tags]
    genome_info = transpose([tais,org_locus_tags,annos])
    dirname = org2nc_id(org)+"-tAI"
    if not dirname in os.listdir('index_results'):
        os.mkdir(os.path.join("index_results",dirname))
    filename = os.path.join("index_results",dirname,"genome_info.txt")
    data2csv(genome_info,filename=filename,sep="\t",overwrite=overwrite)
    
def orf_frequencies_for_R_script(org):
    sequences = coding_sequences(org)
    m = []
    for sequence in verbose_gen(sequences):
        m.append(trna_dict2trna_list(Counter(group_codons(str(sequence.seq))),
                                     codon_order="table"))
    return m
    

def get_trna_dicts(orgs):
    """Return a dictionary of trna counts.
    Format: {org:{codon:[scores]}}"""
    return {org:parse_trnas(org2trna_filename(org))
            for org in orgs}

def parse_trnas(filename,cutoff=None):
#    type_pat = re.compile(r"Type: ([A-Za-z]{3})")
    anticodon_pat = re.compile(r"Anticodon: ([ATGC]{3})")
    score_pat = re.compile(r"Score: ([0-9.]+)")
    with open(filename) as f:
        lines = f.readlines()
    trnas = ["".join(lines[7*i:7*(i+1)]) for i in range(len(lines)/7)]
    codon_keys = defaultdict(list)
    for trna in trnas:
        if "Undet" in trna or "???" in trna or "pseudo" in trna:
            continue
        anticodon = re.search(anticodon_pat, trna).groups(0)[0]
        score = float(re.search(score_pat, trna).groups(0)[0])
        codon = reverse_complement(anticodon).lower()
        codon_keys[codon].append(score)
    return defaultdict(int,{codon:len(codon_keys[codon])
                             for codon in verbose_gen(codon_keys)})
#    return codon_keys

def my_extract(cds,genome):
    if cds.sub_features:
        return "".join([my_extract(sf,genome) for sf in cds.sub_features])
    else:
        seq = str(genome[cds.location.start:cds.location.end].seq)
        return (seq if cds.strand == 1
                else reverse_complement(seq))

def coding_sequences(org):
    try:
        fn = get_genome_filename(org,"ffn")
        cdss = SeqIO.parse(fn,'fasta')
        print "found ffn"
        return [cds.lower() for cds in cdss]
    except IndexError:
        print "found gbk"
        genome = get_genome(org)
        cdss = get_cdss(genome)
        return [str(cds.extract(genome).seq).lower() for cds in verbose_gen(cdss)]

def coding_sequences_dict(org):
    genome = get_genome(org)
    cdss = get_cdss(genome)
    return {head(cds.qualifiers['locus_tag']):str(cds.extract(genome).seq).lower()
                    for cds in verbose_gen(cdss)}

def locus_tags(org):
    "Return a list of coding sequences in the order encountered in gbk file"
    genome = get_genome(org)
    cdss = get_cdss(genome)
    return [head(cds.qualifiers['locus_tag']) for cds in cdss]

def codon_usage(genome):
    codons = defaultdict(int)
    cdss = get_cdss(genome)
    for cds in verbose_gen(cdss):
#        print "extracting"
#        seq = cds.extract(genome)
        seq = my_extract(cds,genome)
#       print "finished extracting"
        for codon in group_codons(seq):
            codons[codon] += 1
    return codons

def codon_trna_correlation(codon_dicts,trna_dicts):
    assert(codon_dicts.keys() == trna_dicts.keys())
    orgs = codon_dicts.keys()
    for codon in codons:
        ts = [trna_dicts[org][codon] for org in orgs]
        print "ts:",ts
        cs = [codon_dicts[org][codon] for org in all_orgs]
        print "cs:",cs
        print codon,pearson(ts,cs)

def trna_analysis(gen_codon_dicts,refset_codon_dicts,trna_dicts):
    for org in all_orgs:
 	for codon in codons:
            gen_freq = gen_codon_dicts[org][codon]
            refset_freq = refset_codon_dicts[org][codon]
            copy_num = trna_dicts[org][codon]
            print "%s,%s,%s" % (refset_freq, gen_freq, copy_num)

def synonymous_codon_analysis(codon_dict,normalized=True):
    """Given a codon frequency dict, return a normcd alized dictionary of
    synonymous codons"""
    cond_norm = lambda xs:map(truncate,normalize(xs)) if normalized else iota
    return {aa:cond_norm([codon_dict[codon] for codon in translation_table[aa]])
            for aa in aas}
    
def codon_behavior(gen_codon_dicts,refset_codon_dicts,trna_dicts):
    for codon in codons:
        gen_freqs = [gen_codon_dicts[org][codon] for org in all_orgs]
        refset_freqs = [refset_codon_dicts[org][codon] for org in all_orgs]
        avg_cn = mean([trna_dicts[org][codon] for org in all_orgs])
        avg_diff = mean(zipWith(operator.sub,refset_freqs,gen_freqs))
        print codon,avg_diff,scipy.stats.wilcoxon(refset_freqs,gen_freqs)[1]<.05,avg_cn


def generate_correlations(): 
    def write_correlations(org_name,method):
        filename = "nrca_freq_correlations/%s_%s_correlations.csv" % (org_name,method.func_name)
        print filename
        if filename in os.listdir("nrca_freq_correlations"):
            print "found ",filename
            return
        return data2csv(all_correlations_all_orgs(eval(org_name),method),
                        filename,
                        header=["org1","org1_group","org2","org2_group","nRCAs","genfreqs","refsetfreqs","ws"])
    org_names = ["last_orgs","all_orgs","actinos","firmicutes","gammas","pseudos",
                 "psychros","pvp_orgs","enteros"]
#    org_names = ["pseudos","psychros"]
    methods = [spearman] #[pearson,spearman,l2]
    for org_name,method in cart_product(org_names,methods):
        write_correlations(org_name,method)

def normalize_codon_dict(d,break_down_six_boxes):
    """Take a dictionary of the form {codon:numeric_val} and return a
    dictionary normalized by amino acid bias"""
    d_copy = d.copy()
    break_down_six_boxes = any("2" in codon or "4" in codon for codon in d)
    def normalization(codon):
        print codon
        denom = float(sum([d_copy[c]
                           for c in synonymous_codons(codon,break_down_six_boxes)
                           if c in d_copy]))
        return denom if denom else 1
    print "defined"
    return defaultdict(int,{codon:d[codon]/normalization(codon) for codon in d})
    
def write_codon_trna_spreadsheet(filename,six_box=False):
    tt = translation_table if not six_box else six_box_translation_table
    sorted_codons = [codon for aa in tt for codon in tt[aa]]
    print "trna_dicts"
    trna_dicts = map_over_vals(normalize_codon_dict,get_trna_dicts(pvp_orgs))
    print "genomic_dicts"
    genomic_dicts = map_over_vals(normalize_codon_dict,get_codon_dicts(pvp_orgs,genomic=True))
    print "refset_dicts"
    refset_dicts = map_over_vals(normalize_codon_dict,get_codon_dicts(pvp_orgs,genomic=False))
    num_fields = 3 #refset freqs, genomic freqs, trna copy number
    rearrange_comma = lambda xs: "," + xs[:len(xs) - 1]
    aa_header = rearrange_comma("".join([aa + "," * len(tt[aa]) * num_fields
                                         for aa in tt]))
    codon_header = rearrange_comma("".join(codon + "," * num_fields
                                           for codon in sorted_codons))
    field_header = rearrange_comma("".join(["refset,genomic,trna,"
                                            for codon in sorted_codons]))
    def line(org):
        print org
        return [org] + sum([[refset_dicts[org][codon],
                             genomic_dicts[org][codon],
                         trna_dicts[org][codon]] for codon in sorted_codons],[])
    pseudo_lines = map(line,pseudos)
    psychro_lines = map(line,psychros)
    def moment_lines(lines,name,f):
        return [name] + map(lambda xs: truncate(f(xs),9),transpose(lines)[1:])#drop
                                                                                 #name
                                                                                 #column
    pseudo_avg = moment_lines(pseudo_lines,"Pseudomonas average",mean)
    psychro_avg = moment_lines(psychro_lines,"Psychrobacter average",mean)
    pseudo_std = moment_lines(pseudo_lines,"Pseudomonas standard dev",sd)
    psychro_std = moment_lines(psychro_lines,"Psychrobacter standard dev",sd)
    with open(filename,'w') as f:
        f.write(aa_header+"\n")
        f.write(codon_header+"\n")
        f.write(field_header+"\n")
        f.write("\n".join([",".join(map(str,line)) for line in pseudo_lines]) + "\n")
        f.write(",".join(map(str,pseudo_avg)) + "\n")
        f.write(",".join(map(str,pseudo_std)) + "\n")
        f.write("\n".join([",".join(map(str,line)) for line in psychro_lines]) + "\n")
        f.write(",".join(map(str,psychro_avg)) + "\n")
        f.write(",".join(map(str,psychro_std)) + "\n")
print "loaded"

def codon_trna_barchart(filename):
    with open(filename) as f:
        lines = [line for line in csv.reader(f, delimiter=",")]
    pseudo_averages = map(float,lines[12][1:])
    psychro_averages = map(float,lines[16][1:])
    pseudo_refset = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 0]
    psychro_refset = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 0]
    pseudo_genomic = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 1]
    psychro_genomic = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 1]
    pseudo_trnas = [field for (i,field) in enumerate(pseudo_averages) if i % 3 == 2]
    psychro_trnas = [field for (i,field) in enumerate(psychro_averages) if i % 3 == 2]
    ys = list(sum(zip(pseudo_refset,psychro_refset,pseudo_genomic,
             psychro_genomic,pseudo_trnas,psychro_trnas),()))
    fig = pylab.figure()
    ax = fig.add_subplot(1,1,1)
    indices = range(len(ys))
    sorted_codons = [codon for aa in tt for codon in translation_table[aa]]
    color_labels = ["%s %s %s" % (codon,species,typ) for codon in sorted_codons
                    for typ in ["refset","genomic","trnas"]
                    for species in ["pseudos","psychros"]]
    labels = ["%s (%s)" % (codon,translate(codon)) for codon in sorted_codons
              for typ in ["refset","genomic","trnas"]
              for species in ["pseudos","psychros"]]
    def color(label):
        if "genomic" in label and "pseudos" in label:
            return "magenta"
        elif "genomic" in label and "psychros" in label:
            return "cyan"
        elif "refset" in label and "pseudos" in label:
            return "red"
        elif "refset" in label and "psychros" in label:
            return "blue"
        elif "trna" in label and "pseudos" in label:
            return "yellow"
        elif "trna" in label and "psychros" in label:
            return "green"
        else:
            print(label)
    for index,y,label,color_label in zip(indices,ys,labels,color_labels):
        ax.bar(index,y,color=color(color_label))
    ax.set_xticks(indices)
    ax.set_xticklabels(labels,rotation='vertical')
    pylab.show()

def run_tai(orgs):
    """Make tai folders in index_results"""
    print "trna dicts"
    trna_dicts = get_trna_dicts(orgs)
    print "sorting"
    sorted_trna_dicts = {org:[trna_dicts[org][codon]
                              for codon in sorted_codons]
                         for org in orgs}
    print "annotations"
    anno_dict = annotation_dict(orgs)
    for org in orgs:
        print "taing ",org
        folder = org2nc_id(org) + "-tAI"
        dirname = os.path.join("index_results",folder)
        if not folder in os.listdir('index_results'):
            print "making folder for ",org
            os.mkdir(dirname)
        filename = "genome_info.txt"
        if filename in os.listdir(dirname):
            print "found ",filename, "for ", org
            continue
        org_tai = tai.TAI(sorted_trna_dicts[org])
        genome = get_genome(org)
        cdss = get_cdss(genome)
        lines = []
        for cds in cdss:
            sequence = str(cds.extract(genome).seq)
            locus_tag = head(cds.qualifiers['locus_tag'])
            description = anno_dict[org][locus_tag]
            tai_value = org_tai.tai(sequence)
            lines.append([tai_value,locus_tag,description])
            print locus_tag
        data2csv(lines,os.path.join(dirname,filename),sep="\t")
        
def refset_info(org,index):
    nc = org2nc_id(org)
    dirname = nc + "-1.0_rcc_%s" % index
    full_path = os.path.join("index_results",dirname,"refset_info.txt")
    with open (full_path) as f:
        lines = [line for line in csv.reader(f,delimiter="\t")]
    return lines

def top_scoring_annotations(org,index,n):
    nc = org2nc_id(org)
    dirname = nc + "-1.0_rcc_%s" % index
    full_path = os.path.join("index_results",dirname,"genome_info.txt")
    with open (full_path) as f:
        lines = [line for line in csv.reader(f,delimiter="\t")]
    return [line[5] for line in lines][:n]

def refset_tags(org,index):
    return [line[0] for line in refset_info(org,index)]

def cliques_in_refsets(orgs,all_results,anno_dict,named_ltds):
    print "loading reciprocals"
    print "making graph"
    g = make_graph(all_results)
    print "finding cliques"
    cliques_of_length = lambda clique_size: [clique for clique in nx.find_cliques(g)
                                             if len(clique) == clique_size]
    "loading refset dicts"
    refset_dicts = {org:refset_tags(org,"nRCA") for org in orgs}
    def num_refsets(clique):
        return sum([tag in refset_dicts[locus_tag2org(tag,named_ltds)]
                    for tag in clique])
    locus_tag2anno = locus_tag2anno_factory(anno_dict,named_ltds)
    _locus_tag2org = lambda tag:locus_tag2org(tag,named_ltds)
    header_clique = head(nx.find_cliques(g),lambda clique:len(clique)==len(orgs))
    sorted_orgs = [locus_tag2org(tag,named_ltds) for tag in sorted(header_clique)]
    sort_clique = lambda clique : [head(clique, lambda tag: _locus_tag2org(tag) == org)
                                   for org in sorted_orgs]
    return ([sorted_orgs + ["Majority Annotation","Refset Count","Clique Size"]] +
            [(sort_clique(clique) +
              [plurality(map(locus_tag2anno,clique))] +
              [num_refsets(clique)] + [len(clique)])
            for i in range(2,len(header_clique) + 1)
             for clique in verbose_gen(cliques_of_length(i))])

def generate_cliques_in_refsets_csvs():
    for group_name in ["pseudos","psychros","pvp_orgs","gammas",
                  "firmicutes","actinos","enteros","last_orgs"]:
        print group_name
        group = eval(group_name)
        group_results = load_reciprocals(group)
        data2csv(cliques_in_refsets(group,group_results,anno_dict,named_ltds),
                 filename = ("clique_csvs/%s_cliques_of_all_lengths_in_refsets,csv"
                             % group_name),
                 sep="\t",overwrite=False)
    
def refset_annotations(org,index,anno_dict):
    def locus_tag(line):
        #This function exists to cope with trailing whitespace, which
        #arbitrarily litters refset_info files.  OOOORRRR!
        return line[0]
    return [anno_dict[org][locus_tag(line)] for line in refset_info(org,index)]

def get_refsets(orgs,index,anno_dict):
    return {org:refset_annotations(org,index,anno_dict) for org in orgs}

def get_refset_tags(orgs,index):
    return {org:refset_tags(org,index) for org in orgs}

def refset_words(orgs,index,anno_dict):
    refsets = get_refsets(orgs,index,anno_dict)
    return (sum([term.replace(","," ").split(" ")
                         for org in orgs
                         for term in refsets[org]],[]))

def genomic_words(orgs,anno_dict):
    return [word
            for org in orgs
            for annotation in (anno_dict[org].values())
            for word in annotation.replace(","," ").split(" ")]

def weighted_words(orgs,anno_dict,indices):
    return [word
            for org in orgs
            for tag in anno_dict[org]
            for word in anno_dict[org][tag].replace(","," ").split(" ")
            if indices[org][tag] > random.random()]

def print_enrichment(refset_words,genomic_words):
    refset_counter = Counter(refset_words)
    genomic_counter = Counter(genomic_words)
    refset_length = len(refset_words)
    genomic_length = len(genomic_words)
    enrichments = {word:((refset_counter[word]/float(refset_length)) /
                         (genomic_counter[word]/float(genomic_length)))
                   for word in refset_counter}
    return [(enrichments[word],refset_counter[word],genomic_counter[word],word)
            for word in sorted(enrichments,key=lambda w:enrichments[w],reverse=True)]

def write_enrichments(org_name):
    orgs = eval(org_name)
    data = print_enrichment(refset_words(orgs,"nRCA",anno_dict),
                            genomic_words(orgs,anno_dict))
    data2csv(data,"refset_annotations/%s_refset_annotations.csv" % org_name,
             header=["Enrichment","refset_count","genomic_count","term"])
    
def refset_summary(orgs,anno_dict,filename,overwrite=False):
    data = []
    for org in orgs:
        data.append([org])
        data.append([])
        for annotation in refset_annotations(org,"nRCA",anno_dict):
            data.append([annotation])
    data2csv(data,filename,overwrite=overwrite)

def top_scoring_summary(orgs,filename,n):
    data = []
    for org in orgs:
        data.append([org])
        data.append([])
        for annotation in top_scoring_annotations(org,"nRCA",n):
            data.append([annotation])
    data2csv(data,filename)

def setup(orgs):
    """Construct databases, blast and collate reciprocals for orgs"""
    populate_dbs(orgs)
    reciprocal_blasts2(orgs)
    collate_reciprocal_blasts(orgs)

def compute_cliques(orgs,anno_dict,named_ltds):
    print "loading reciprocals"
    all_results = load_reciprocals(orgs)
    print "making graph"
    g = make_graph(all_results)
    print "finding cliques"
    cliques = map(sorted,find_full_cliques(g,orgs))
    print "analyzing cliques"
    return cliques,analyze_cliques(cliques,orgs,anno_dict,named_ltds).values()

def generate_all_full_clique_csvs():
    for group_name in group_names:
        cliques,_ = compute_cliques(eval(group_name),anno_dict,named_ltds)
        data2csv(map(sorted,cliques),filename="clique_csvs_"+group_name + "_full_cliques.csv")

# def ribosomal_cliques(orgs,anno_dict,named_ltds):
#     cliques,clique_dict = compute_cliques(orgs,anno_dict,named_ltds)
#     header_clique = head(cliques,lambda clique:len(clique)==len(orgs))
#     sorted_orgs = [locus_tag2org(tag,named_ltds) for tag in sorted(header_clique)]
#     return ([sorted_orgs] +
#             [[indices[locus_tag2org(tag,named_ltds)][tag] for tag in clique]
#              for clique in [clique for i,clique in enumerate(cliques)
#                             if any(("ribosomal" in term
#                                     for term in clique_dict[i]))]])

def generate_ribosomal_clique_csvs():
    for group_name in group_names:
         filename = "clique_csvs/%s_ribosomals.csv" % group_name
         try:
             open(filename)
             print "found",filename
             continue
         except IOError:
             data2csv(phrase_cliques(eval(group_name),anno_dict,named_ltds,"ribosomal"),
                      filename,
                      sep="\t",overwrite=False)

def generate_full_clique_csvs():
    for group_name in group_names:
         filename = "clique_csvs/%s_full_cliques.csv" % group_name
         try:
             open(filename)
             print "found",filename
             continue
         except IOError:
             data2csv(phrase_cliques(eval(group_name),anno_dict,named_ltds,""),
                      filename,
                      sep="\t",overwrite=False)

def phrase_cliques(orgs,anno_dict,named_ltds,phrase):
    cliques,clique_dict = compute_cliques(orgs,anno_dict,named_ltds)
    header_clique = head(cliques,lambda clique:len(clique)==len(orgs))
    sorted_orgs = [locus_tag2org(tag,named_ltds) for tag in sorted(header_clique)]
    _locus_tag2org = lambda tag:locus_tag2org(tag,named_ltds)
    sort_clique = lambda clique : [head(clique, lambda tag: _locus_tag2org(tag) == org)
                                   for org in sorted_orgs]
    return ([sorted_orgs] +
            [[indices[locus_tag2org(tag,named_ltds)][tag]
              for tag in sort_clique(clique)]
             for clique in [clique
                            for i,clique in enumerate(cliques)
                            if any((phrase in term for term in clique_dict[i]))]])

def compute_clique_crossings(cliques):
    """Given a list of nRCA values for orthologs in a group of cliques
    (structured as [[clique1],[clique2]...]) compute the number of
    _crossings_, or occurences of the form nRCA(org1,clique_i) <
    nRCA(org1,clique_j) but nRCA(org2,clique_i) < nRCA(org2,clique_j)"""
    orgs = transpose(cliques) #[[cliques in org1],[cliques in org2]...]
    crossings = 0
    for org_i,org_j in choose2(orgs):
        for x,y in choose2(range(len(org_i))):
            if not (cmp(org_i[x],org_i[y]) == cmp(org_j[x],org_j[y])):
                crossings += 1
    return crossings

def random_clique_matrix(num_cliques,num_orgs):
    return [[random.random() for org in range(num_orgs)]
            for clique in range(num_cliques)]

def locus_tag2index_factory(indices,named_ltds):
    return lambda tag: indices[locus_tag2org(tag,named_ltds)][tag]

def locus_tag2exp_factory(exp_dicts,named_ltds):
    def locus_tag2exp(tag):
        try:
            return exp_dicts[locus_tag2org(tag,named_ltds)][tag]
        except:
            return None
    return locus_tag2exp

def locus_tag2anno_factory(anno_dict,named_ltds):
    def locus_tag2anno(tag):
        return anno_dict[locus_tag2org(tag,named_ltds)][tag]
    return locus_tag2anno

def rca_vs_tai_comparison():
    all_exp_dicts = exp_dicts(all_orgs)
    
    rca_indices = index_dicts(all_orgs,"nRCA")
    tai_indices = index_dicts(all_orgs,"tAI")
    for org in all_exp_dicts:
	print org
        # tags,rcas,tais,exps = zip(*[(tag,indices[org][tag],
        #                              tai_indices[org][tag],
        #                              mean(all_exp_dicts[org][tag]))
        #                             for tag in all_exp_dicts[org]
        #                             if (tag in indices[org] and
        #                                 tag in tai_indices[org])])
        
        lines = ([(tag,indices[org][tag],
                       tai_indices[org][tag],
                       mean(all_exp_dicts[org][tag]))
                      for tag in all_exp_dicts[org]
                      if (tag in indices[org] and
                          tag in tai_indices[org])])
        tags,rcas,tais,exps = transpose(lines)
        data2csv(lines,"%s_rca_vs_tai_vs_exp_comparison.csv" % org,header = ["tag","rca","tai","exp"])
	print "rca vs exp:",spearman(rcas,exps),"tai vs exp",spearman(tais,exps)

def exp_vs_exp_pvp_plot(pvp_cliques,exp_dicts,named_ltds):
    locus_tag2exp = locus_tag2exp_factory(exp_dicts,named_ltds)
    expressables = ["Pseudomonas_putida",
                    "Pseudomonas_fluorescens",
                    "Pseudomonas_aeruginosa",
                    "Psychrobacter_arcticus_273"]
    tags = {org:[tag for clique in pvp_cliques for tag in clique
                 if locus_tag2org(tag,named_ltds) == org] for org in expressables}
    def all_data_for_ith_tag(i):
        return all ([mean(locus_tag2exp(other_tag))
                     for other_tag in [tags[other_org][i]
                                       for other_org in expressables]])
    tags_with_no_missing_data = {org:[tag for (i,tag) in enumerate(tags[org])
                                      if all_data_for_ith_tag(i)]
                                 for org in expressables}
    
    exps = {org:map(lambda tag: mean(locus_tag2exp(tag)),
                    tags_with_no_missing_data[org]) for org in expressables}
    psych_exps = exps[expressables[3]]
    pseudo_exps = map(mean,zip(*[exps[pseudo] for pseudo in expressables[:3]]))
    
def cliques_over(orgs):
    all_results = load_reciprocals(orgs)
    g = make_graph(all_results)
    cliques = find_full_cliques(g,orgs)
    return cliques

def clique_annos(cliques,locus_tag2anno):
    return [plurality(map(locus_tag2anno,clique)) for clique in cliques]

def build_zero_crossings_graph():
    pvp_orgs_full_cliques = [line for line in csv.reader(open("clique_csvs/pvp_orgs_full_cliques.csv"),delimiter="\t")]
    pvp_orgs_full_cliques = [map(float,line) for line in pvp_orgs_full_cliques[1:]]
    g = nx.Graph()
    for (i,j) in choose2(range(361)):
        line1 = pvp_orgs_full_cliques[i]
        line2 = pvp_orgs_full_cliques[j]
        if(reduce(lambda x,y:x==y,zipWith(lambda u,v: u < v,line1,line2))):
            g.add_edge(i,j)
            print i,j
    return g

def read_cliques(org_name):
    """Read csv of clique nRCA values, a la in R"""
    with open("clique_csvs/%s_full_cliques.csv" % org_name) as f:
        cliques = mapmap(float,list(csv.reader(f,delimiter="\t"))[1:])
    return cliques

def distance_graph(cliques,inverse=False):
    """Given a matrix of cliques, return a weighted graph whose nodes
    are cliques and whose edge weights are the Euclidean distances
    between their rows.  If inverse if false, weight=distance.  If
    True, weight=1/distance"""
    def distance(xs,ys):
        return sqrt(sum(zipWith(lambda x,y: (x-y)**2,xs,ys)))
    def cond_inv(x):
        return 1/x if inverse else x
    g = nx.Graph()
    n = len(cliques)
    for (i,j) in choose2(range(n)):
        g.add_edge(i,j,weight=cond_inv(distance(cliques[i],cliques[j])))
    return g

def normalize_by_synonymous_codons(codon_dict):
    def normalization(codon):
        denom = float(sum(dictmap(codon_dict, synonymous_codons(codon,False))))
        return 1/denom if denom else 0
    return {codon:codon_dict[codon]*normalization(codon) for codon in codons}


def pvp_barchart_unnormalized_by_aa(orgs,codon_set=four_box_codons,testing=False):
    """Compute the ending position frequencies for every pvp in the
    refset and genome, compare using MWU.
    Wed Aug  8 14:08:53 EDT 2012 """
    mannwhitneyu = scipy.stats.mannwhitneyu
    pseudo_genomic = defaultdict(list)
    pseudo_refset = defaultdict(list)
    pseudo_trna = defaultdict(list)
    psychro_genomic = defaultdict(list)
    psychro_refset = defaultdict(list)
    psychro_trna = defaultdict(list)
    def ending_freqs_of(c,freqs):
        denom = sum([freqs[codon] for codon in codon_set])
        return ([freqs[codon]/denom
                 for codon in codon_set
                 if codon.endswith(c)])
    for org in orgs:
        genomic_freqs = genomic_codon_freqs(org)
        refset_freqs = refset_codon_freqs(org)
        trna_freqs = trna_codon_freqs(org)
        if "Pseudo" in org:
            for c in "acgt":
                pseudo_genomic[c].append(sum(ending_freqs_of(c,genomic_freqs)))
                pseudo_refset[c].append(sum(ending_freqs_of(c,refset_freqs)))
                pseudo_trna[c].append(sum(ending_freqs_of(c,trna_freqs)))
        elif "Psychro" in org or "Acineto" in org:
            for c in "acgt":
                psychro_genomic[c].append(sum(ending_freqs_of(c,genomic_freqs)))
                psychro_refset[c].append(sum(ending_freqs_of(c,refset_freqs)))
                psychro_trna[c].append(sum(ending_freqs_of(c,trna_freqs)))
    for c in "agct":
        if testing:
            print "pseudo refset vs. genomic %s ending:" %c, mannwhitneyu(pseudo_refset[c],
                                                                     pseudo_genomic[c])[1]
            print"psychro refset vs. genomic %s ending:" %c, mannwhitneyu(psychro_refset[c],
                                                                      psychro_genomic[c])[1]
    return (pseudo_genomic,pseudo_refset,pseudo_trna,
            psychro_genomic,psychro_refset,psychro_trna)

def n_ending_freqs(orgs,codon_set=four_box_codons):
    genomic = defaultdict(list)
    refset = defaultdict(list)
    trna = defaultdict(list)
    def ending_freqs_of(c,freqs):
        denom = sum([freqs[codon] for codon in codon_set])
        return ([freqs[codon]/denom
                 for codon in codon_set
                 if codon.endswith(c)])
    for org in orgs:
        genomic_freqs = genomic_codon_freqs(org)
        refset_freqs = refset_codon_freqs(org)
        trna_freqs = trna_codon_freqs(org)
        for c in "acgt":
            genomic[c].append(sum(ending_freqs_of(c,genomic_freqs)))
            refset[c].append(sum(ending_freqs_of(c,refset_freqs)))
            trna[c].append(sum(ending_freqs_of(c,trna_freqs)))
    return (genomic,refset,trna)

def summarize_pvp_barchart(barchart_results,filename,group1_name,group2_name):
    (pseudo_genomic,pseudo_refset,pseudo_trna,
            psychro_genomic,psychro_refset,psychro_trna) = barchart_results
    data = [["VB" + c.upper(),aspect, mean(eval("pseudo_" + aspect)[c]),
            se(eval("pseudo_" + aspect)[c]),mean(eval("psychro_" + aspect)[c]),
            se(eval("psychro_" + aspect)[c])]
            for c in "acgt"
            for aspect in ["genomic","refset","trna"]]
    data2csv(data,filename="boxes_analysis/" + filename,
             header=["","",group1_name + " mean",group1_name + " se",group2_name + " mean",group2_name + " se"])
    
def cog_analysis():
    pa_lts = named_ltds["Pseudomonas_aeruginosa"].values()
    pa_cogs = map(cog.get_cog,pa_lts)
    pvp_clique_cogs = [cog.get_cog(head(clique,lambda t:t.startswith("PA")))
                       for clique in pvp_cliques]
    #count pa_cogs
    for c in cog.cogs:
        print c, cog.cog_definitions[c], pa_cogs.count(c)

def succ_diffs(w_seq):
    ps = pairs(w_seq)
    return sum(map(lambda(x,y): abs(x-y),ps))

def succ_diff_perm_test(w_seq,n):
    """Compare succ_diff statistic of w_seq against permutations of
    w_seq, return proportion of control statistics lower than
    succ_diffs(w_seq)"""
    statistic = succ_diffs(w_seq)
    length = len(w_seq)
    controls = [succ_diffs(random.sample(w_seq,length)) for i in range(n)]
    return len(filter(lambda control: control < statistic,controls))/float(n)
    
def w_sequence(cds,w_codons):
    epsilon = 10**-10
    cds = cds.lower()
    codons = group_codons(cds)
    if codons[-1] in ["tag","tga","taa"]:
        codons = codons[:len(codons) - 1]
    lookups = dictmap(w_codons,codons)
    return lookups
        
def score_cds(cds,w_codons):
    lookups = w_sequence(cds,w_codons)
    log_lookups = map(lambda w:log(w + epsilon),lookups) 
    nrca = exp(sum(log_lookups)/float(len(codons)))
    return  nrca
    
def split_cds(cds,fraction):
    cds_codons = group_codons(cds)
    cut = int(floor(fraction * len(cds_codons)))
    return ("".join(cds_codons[:cut]),"".join(cds_codons[cut:]))


def sliding_windows(cds_codons,n):
    """return a list containing sliding windows of length n"""
    return ["".join(cds_codons[i:i+n]) for i in range(len(cds_codons)-n +1)]


def sliding_score(cds,w_codons,frac=.1):
    print "frac:",frac
    cds_codons = group_codons(cds)
    n = int(frac * len(cds_codons))
    print "n:",n
    windows = sliding_windows(cds_codons,n)
    print "num_windows:",len(windows)
    return map(lambda seq:score_cds(seq,w_codons),windows)

def plot_sliding_windows(cds,w_codons,frac=.1):
    scores = sliding_score(cds,w_codons,frac)
    plt.plot(*zip(*[(3*i,score) for (i,score) in enumerate(scores)]))

def control_sequence(genomic_freqs,length):
    return "".join([sample(genomic_freqs.keys(),genomic_freqs.values())
                    for i in range(length)])

def control_score(genomic_freqs,ws,length):
    return score_cds(control_sequence(genomic_freqs,length),ws)

def generate_nrcas_vs_exps(orgs):
    exps = exp_dicts(orgs)
    nrcas = index_dicts(orgs,"nRCA")
    for org in exps:
        print org
        data = [(tag,nrcas[org][tag],mean(exps[org][tag]))
                for tag in exps[org] if tag in nrcas[org]]
        data2csv(data,filename="nrca_vs_exp/" + "nrca_vs_exp_%s" % org,
                 header = ["Locus tag","nRCA","exp"])

org1 = "Pseudomonas_aeruginosa"
org2 = "Psychrobacter_arcticus_273"
def compare_ws_in_cliques(org1,org2):
    """Are w-values conserved in orthologues?"""
    genome1 = get_genome(org1)
    genome2 = get_genome(org2)
    ws1 = org2w_codons(org1)
    ws2 = org2w_codons(org2)
    cdss1 = get_cdss(genome1)
    cdss2 = get_cdss(genome2)
    orthologues = [(tup[0],tup[2]) for tup in pairwise_correlations(org1,org2)]
    data = []
    for i,orthologue in enumerate(orthologues[320:]):
        gene1,gene2 = orthologue
        cds1 = head(cdss1,lambda cds: gene1 in cds.qualifiers['locus_tag'])
        seq1 = str(cds1.extract(genome1).seq)
        w_seq1 = w_sequence(seq1,ws1)
        protein1 = translate(seq1.lower())[1:]
        cds2 = head(cdss2,lambda cds: gene2 in cds.qualifiers['locus_tag'])
        seq2 = str(cds2.extract(genome2).seq)
        w_seq2 = w_sequence(seq2,ws2)
        protein2 = translate(seq2.lower())[1:]#remove ATG
        if len(protein1) > 2000:
            print "protein1:",len(protein1),"skipping"
            continue
        alignment = pairwise2.align.globalds(protein1,protein2,BLOSUM62,-7,-1)[0]
        align1,align2 = alignment[:2]
        w_seq1_cp = w_seq1[:]#copy for destructive pops
        w_seq2_cp = w_seq2[:]
        altered_align1 = [w_seq1_cp.pop(0) if a in string.ascii_uppercase else None
                          for a in align1]
        altered_align2 = [w_seq2_cp.pop(0) if a in string.ascii_uppercase else None
                          for a in align2]
        tuples = (filter(lambda tup: all(tup),
                            zip(altered_align1,altered_align2)))
        aligned_spearman = spearman(*transpose(tuples))
        m = min(len(w_seq1),len(w_seq2))
        raw_spearman = spearman(w_seq1[:m],w_seq2[:m])
        # data.append(w_seq1)o
        # data.append(w_seq2)
        print i, aligned_spearman, raw_spearman, len(tuples),mean([len(w_seq1),len(w_seq2)])
        del(seq1)
        del(seq2)

def welch_analysis():
    """Welch et al. 2009 report that expression is better correlated
    with codons whose tRNAs are preferentially charged during
    starvation, within synthetic genes for a given protein.  """

    #preliminaries:
    exps = exp_dict(exp_path(ecoli))
    exps = {k:mean(exps[k]) for k in exps}
    cdss = cdss_from_org(ecoli)
    genome = get_genome(ecoli)
    coding_seqs = coding_sequences_dict(ecoli)
    for locus_tag in coding_seqs:
        if locus_tag in exps:
            codons = group_codons(coding_seqs[locus_tag])
            serine_codons = translation_table("S")
            num_serine = len([c for c in codons if translate(c) == 'S'])
            serine_codon_percents = {codon:len([c for c in codons
                                          if c == codon])/float(num_serine+epsilon)
                                     for codon in serine_codons}
            
            print locus_tag,serine_codon_percents["tct"],exps[locus_tag]
            
    #how does usage of codon GCA (Ala) vary with expression?
    
    
main = False
if main:
    anno_dict = annotation_dict(all_orgs)
    named_ltds = make_named_ltds(all_orgs)
    lt2anno = locus_tag2anno_factory(anno_dict,named_ltds)
    pvp_cliques = cliques_over(pvp_orgs)

def refset_conservation(org1,org2):
    """For each org, find the number of genes in that org's refset
    which have orthologues in the other refset,divided by the total
    number of refset genes with orthologues.  Return the max of the
    two."""
    pcs = pairwise_correlations(org1,org2)
    refset1 = refset_tags(org1,"nRCA")
    refset2 = refset_tags(org2,"nRCA")
    refset1_to_genome2 = [lt1 for lt1 in refset1 if any(lt1 in pc for pc in pcs)]
    refset1_to_refset2 = [lt1 for lt1 in refset1 if any((lt1 in pc and lt2 in pc)
                                                        for pc in pcs
                                                        for lt2 in refset2)]
    refset2_to_genome1 = [lt2 for lt2 in refset2 if any(lt2 in pc for pc in pcs)]
    refset2_to_refset1 = [lt2 for lt2 in refset2 if any((lt2 in pc and lt1 in pc)
                                                        for pc in pcs
                                                        for lt1 in refset1)]
    r1r2 = len(refset1_to_refset2)
    r1g2 = len(refset1_to_genome2)
    r2r1 = len(refset2_to_refset1)
    r2g1 = len(refset2_to_genome1)
    ratio1 = r1r2/(r1g2 + epsilon)
    ratio2 = r2r1/(r2g1 + epsilon)
    ratio = max(ratio1,ratio2)
    print org1,org2, r1r2,r1g2,ratio1,r2r1,r2g1,ratio2,ratio
    
def refset_conservation_exp(orgs):
    for org1,org2 in choose2(orgs):
        refset_conservation(org1,org2)

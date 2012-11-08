#!/usr/bin/python2.6

"""This script accepts a list species names and sets up slurm jobs to
perform pair-wise reciprocal BLAST searches between them.

Usage:

$ ./setup_slurms.py list_of_species_names.txt

Preconditions:

"""

import os,sys,time
from string import Template
from utils import *

ORG_PATH = "../data"
BLAST_PATH = "../../ncbi-blast-2.2.26+/bin/blastp"
BLAST_RESULTS_PATH = '../blast_results'
DB_PATH = '../blastdb'
PARTITION = "batch"

def populate_dbs(orgs):
    """given a list of organism names, find directories, grab files
    and construct databases"""
    org_dirs = os.listdir(ORG_PATH)
    file_ext = ".faa"
    for org in orgs:
        print "constructing db for ", org
        org_dir = head([od for od in org_dirs if org_matches_dir(org,od)])
        full_org_dir = os.path.join(ORG_PATH,org_dir)
        fasta_file = head([f for f in os.listdir(full_org_dir)
                           if f.endswith(file_ext)])
        print "using fasta file", fasta_file
        sp = {'name':org,
              'genome_file':os.path.join(full_org_dir,fasta_file),
              'db': os.path.join(DB_PATH,org)}
        makedb(sp)

def reciprocal_blasts2(orgs,new_orgs=None):
    print orgs
    print new_orgs
    all_orgs = orgs + new_orgs if new_orgs else orgs
    org_dirs = os.listdir(ORG_PATH)
    results_contents = os.listdir(BLAST_RESULTS_PATH)
    file_ext = '.faa'
    for org1 in all_orgs:
        for org2 in all_orgs:
            print "starting on: ",org1, org2, "at", time.ctime()
            out_file = "results_%s_%s.txt" % (org1,org2)
            same_org = org1 == org2
            already_done = (new_orgs and not (org1 in new_orgs or
                                              org2 in new_orgs))
            if same_org or already_done:
                print "skipping", org1, org2
                continue 
#            org_dir = head([od for od in org_dirs if org1 in od])
            print "looking for data directory for ", org1
            org_dir = head([od for od in org_dirs if org_matches_dir(org1,od)])
            full_org_dir = os.path.join(ORG_PATH,org_dir)
            fasta_file = head([f for f in os.listdir(full_org_dir)
                                            if f.endswith(file_ext)])
            full_fasta_file = os.path.join(full_org_dir,fasta_file)
            full_db_path = os.path.join(DB_PATH,org2)
            full_out_file = os.path.join(BLAST_RESULTS_PATH,out_file)
            if out_file in results_contents:
                print "found results, skipping", org1, org2
                continue
            slurm_file = "%s_%s.slurm" % (org1,org2)
            with open(slurm_file,'w') as f:
                f.write(slurm_template.substitute(job_name="%s_%s" % (org1,org2),
                                                  partition=PARTITION,
                                                  query=full_fasta_file,
                                                  db=full_db_path,
                                                  outfile=full_out_file,
                                                  blast_path=BLAST_PATH))
            os.system("sbatch %s" % slurm_file)
            print "finished",org1, org2,"at", time.ctime()

def org_matches_dir(org,org_dir):
    return all(word.lower() in org_dir.lower() for word in org.split('_'))

slurm_template = Template("""#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err
#SBATCH --partition=${partition}
#SBATCH --qos=short

${blast_path} -query ${query} -db ${db} -out ${outfile} -outfmt 5 -task blastp""")


if __name__ == '__main__':
    orgs_file = sys.argv[1]
    orgs = [org.strip() for org in open(orgs_file).readlines()]
    try:
        new_orgs_file = sys.argv[2]
        new_orgs = [new_org.strip() for new_org in open(new_orgs_file).readlines()]
    except IndexError:
        new_orgs = None
    print ""
    reciprocal_blasts2(orgs,new_orgs)

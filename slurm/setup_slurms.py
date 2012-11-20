#!/usr/bin/python2.6

"""This script accepts a list species names and sets up slurm jobs to
perform pair-wise reciprocal BLAST searches between them.

Usage:

$ ./setup_slurms.py list_of_species_names.txt

Preconditions:

Postconditions:

"""

import os,sys,time
import subprocess
from string import Template
sys.path.append("..")
from utils import *
from orgs import *

ORG_PATH = "../data"
BLAST_RESULTS_PATH = '../blast_results'
DB_PATH = '../blastdb'
username = "pon2"
job_limit = 8

def populate_dbs(orgs):
    """given a list of organism names, find directories, grab files
    and construct databases"""
    org_dirs = os.listdir(ORG_PATH)
    file_ext = ".faa"
    for org in orgs:
        print "constructing db for ", org
        org_dir = head([od for od in org_dirs if org_matches_dir(org,od)])
        print ORG_PATH,org_dir
        full_org_dir = os.path.join(ORG_PATH,org_dir)
        fasta_file = head([f for f in os.listdir(full_org_dir)
                           if f.endswith(file_ext)])
        print "using fasta file", fasta_file
        sp = {'name':org,
              'genome_file':os.path.join(full_org_dir,fasta_file),
              'db': os.path.join(DB_PATH,org)}
        makedb(sp)

def watch_squeue_until_free():
    while True:
        result = subprocess.Popen("squeue",shell=True,stdout=subprocess.PIPE)
        job_count = result.stdout.read().count(username)
        if job_count < job_limit:
            return
        else:
            print "job queue full, sleeping"
            time.sleep(10)
    
def reciprocal_blasts2(orgs,new_orgs=[],program="blastp",intra_new_orgs=True,one_way=False):
    evalue = 1e-10
    BLAST_PATH = "../../ncbi-blast-2.2.26+/bin/" + program
    print orgs
    print new_orgs
    all_orgs = orgs + new_orgs
    org_dirs = os.listdir(ORG_PATH)
    results_contents = os.listdir(BLAST_RESULTS_PATH)
    if program == "blastp":
        file_ext = lambda(filename):filename.endswith('.faa')
    else:
        file_ext = lambda(filename):(filename.endswith(".fna") or filename.endswith("ffn"))
    for org1 in all_orgs:
        for org2 in all_orgs:
            print "starting on: ",org1, org2, "at", time.ctime()
            out_file = "results_%s_%s.txt" % (org1,org2)
            same_org = (org1 == org2)
            skip_since_two_old_orgs = new_orgs and (org1 in orgs and org2 in orgs)
            skip_since_two_new_orgs = ((org1 in new_orgs) and (org2 in new_orgs)
                                       and not intra_new_orgs)
            skip_since_one_way = one_way and not (org1 in orgs and org2 in new_orgs)
            if (same_org or skip_since_two_old_orgs or skip_since_two_new_orgs
                or skip_since_one_way):
                print "skipping", org1, org2
                continue 
#            org_dir = head([od for od in org_dirs if org1 in od])
            print "looking for data directory for ", org1
            org_dir = head([od for od in org_dirs if org_matches_dir(org1,od)])
            print ORG_PATH,org_dir
            full_org_dir = os.path.join(ORG_PATH,org_dir)
            fasta_file = head(os.listdir(full_org_dir), file_ext)
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
                                                  blast_path=BLAST_PATH,
                                                  evalue=evalue))
            watch_squeue_until_free()
            os.system("sbatch %s" % slurm_file)
            print "finished",org1, org2,"at", time.ctime()

def org_matches_dir(org,org_dir):
    return all(word.lower() in org_dir.lower() for word in org.split('_'))

slurm_template = Template("""#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err
#SBATCH --partition=${partition}
#SBATCH --qos=medium

${blast_path} -query ${query} -db ${db} -out ${outfile} -outfmt 5 -evalue ${evalue}""")


if __name__ == '__main__':
    orgs = eval(sys.argv[1])
    num_args = len(sys.argv)
    new_orgs = []
    program = "blastp"
    PARTITION = "batch"
    one_way = False
    intra_new_orgs = True
    #usage: ./setup_slurms orgs new_orgs program intra_new_orgs one_way partition
    if num_args >= 3: #including self
        new_orgs= eval(sys.argv[2])
    if num_args >= 4:
        program = sys.argv[3]
    if num_args >= 5:
        intra_new_orgs = bool(sys.argv[4])
    if num_args >= 6:
        one_way = bool(sys.argv[5])
    if num_args == 7:
        PARTITION = sys.argv[6]
    reciprocal_blasts2(orgs,new_orgs,program,intra_new_orgs)

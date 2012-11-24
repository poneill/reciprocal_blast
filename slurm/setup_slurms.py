#!/usr/bin/python2.6

"""This script accepts a list species names and sets up slurm jobs to
perform pair-wise reciprocal BLAST searches between them.

Usage:

$ ./setup_slurms.py list_of_species_names.txt

Preconditions:

Postconditions:

"""

import os,sys,time
import subprocess,getopt
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

def watch_until_free(slurm):
    if slurm:
        command = "squeue"
        search_term = username
    else:
        command = "ps aux"
        search_term = "blast"
    while True:
        result = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
        job_count = result.stdout.read().count(search_term)
        print "currently ", job_count, "jobs in queue" 
        if job_count < job_limit:
            return
        else:
            print "job queue full, sleeping"
            time.sleep(10)

     
    
def reciprocal_blasts2(orgs,new_orgs=[],program="blastp",intra_new_orgs=True,one_way=False,slurm=False):
    evalue = 1e-10
    BLAST_PATH = program
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
            blast_string = Template("""${blast_path} -query ${query} -db ${db}\
 -out ${outfile} -outfmt 5 -evalue ${evalue}""").substitute(blast_path=BLAST_PATH,
                                                            query=full_fasta_file,
                                                            db=full_db_path,
                                                            outfile=full_out_file,
                                                            evalue=evalue)
            job_file = write_blast_job(blast_string,org1,org2,slurm)
            watch_until_free(slurm)
            run_job(job_file,slurm)

def write_blast_job(blast_string,org1,org2,slurm):
    print "slurm:",slurm
    slurm_header_template = Template("""#SBATCH --job-name=${job_name}
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err
#SBATCH --partition=${partition}
#SBATCH --qos=medium""")
    job_template = Template("""#!/bin/bash
${slurm_headers}
${blast_string}
""")
    job_name = ("%s_%s" % (org1,org2))
    job_file = job_name + (".slurm" if slurm else ".sh")
    print "job_file:",job_file
    if slurm:
        slurm_header_string = slurm_header_template.substitute(job_name=job_name,
                                                               partition=PARTITION)
        job_string = job_template.substitute(blast_string=blast_string,
                                             slurm_headers=slurm_header_string)
    else:
        job_string = job_template.substitute(blast_string=blast_string,
                                             slurm_headers="")
    with open(job_file,'w') as f:
                f.write(job_string)
    return job_file    

def run_job(job_file,slurm):
    print "running job: ",job_file
    print "slurm:",slurm
    command = "sbatch" if slurm else "bash"
    print "command:",command
    # subprocess.Popen([command,job_file],stdin=None,stdout=None,stderr=None,
    #                 close_fds=None)
    os.system(command + " " + job_file)

    
def org_matches_dir(org,org_dir):
    return all(word.lower() in org_dir.lower() for word in org.split('_'))

if __name__ == '__main__':
    program = "../../ncbi-blast-2.2.26+/bin/blastp"
    PARTITION = "batch"
    one_way = False
    intra_new_orgs = True
    slurm = True
    #usage: ./setup_slurms --orgs=orgs --new_orgs=new_orgs --program=program --no_intra_new_orgs --one_way --develop
    args = ["orgs=","new_orgs=","program=","no_intra_new_orgs",
            "one_way","develop","no_slurm"]
    options,remainder = getopt.getopt(sys.argv[1:],"",args)
    for opt,arg in options:
        print "assigning opt:",opt,"arg:",arg 
        if opt == "--orgs":
            orgs = eval(arg)
        elif opt == "--new_orgs":
            new_orgs = eval(arg)
        elif opt == "--program":
            program = arg
        elif opt == "--no_intra_new_orgs":
            intra_new_orgs = False
        elif opt == "--one_way":
            one_way = True
        elif opt == "--develop":
            PARTITION = "develop"
        elif opt == "--no_slurm":
            slurm = False
    for arg in ["orgs","new_orgs","program","intra_new_orgs",
                "one_way","slurm","PARTITION"]:
        print arg,eval(arg)
    # if num_args >= 3: #including self
    #     new_orgs= eval(sys.argv[2])
    # if num_args >= 4:
    #     program = sys.argv[3]
    # if num_args >= 5:
    #     intra_new_orgs = bool(sys.argv[4])
    # if num_args >= 6:
    #     one_way = bool(sys.argv[5])
    # if num_args >= 7:
    #     PARTITION = sys.argv[6]
    # if num_args == 8:
    #     slurm = sys.argv[7]
    # print "one_way:",one_way
    reciprocal_blasts2(orgs,new_orgs,program,intra_new_orgs=intra_new_orgs,
                       one_way=one_way,slurm=slurm)

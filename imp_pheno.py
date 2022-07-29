"""
SYNOPSIS
    python3 imp_pheno.py <filename>
"""

import csv
import os
import sys

inds = {}
cl = {}

base_name = input('Insert name without extension: ') if len(sys.argv)<2 else sys.argv[1]

fam_file = base_name + '.fam'
if not os.path.isfile(fam_file):
    print("This file doesn't exist!")
    quit()

csv_file = 'Research_Group_List.csv'
if not os.path.isfile(csv_file):
    print("Please store in this folder a csv named 'Research_Group_List.csv' containing phenotypes")
    quit()

out_file = base_name + '_subj_list.txt'
status_phenotype = {
    "CN": 1,
    "EMCI": 2,
    "MCI": 2,
    "LMCI": 2,
    "AD": 3    
}

with open(fam_file, 'r') as in_fam:
    for r in in_fam.readlines():
        fid, pid = r.split(' ')[:2]
        inds[pid] = fid
        
file_reader = csv.reader(open(csv_file), delimiter = ',')

# skip the header...
next(file_reader)
for row in file_reader:
    pid = row[0].replace(' ', '_')
    study = row[1].replace(' ', "_")
    sex = row[2].replace(' ', '_')
    pheno = status_phenotype[row[3].replace(' ', '_')]
    
    cl[pid] = [study, sex, pheno]

# NOTE: using a dictionary will keep only the latest read entry for each patient
with open(out_file, 'w') as f:
    for p in cl.keys():
        if p in inds.keys():
            f.write(f"{inds[p]}\t{p}\t{cl[p][0]}\t{cl[p][1]}\t{cl[p][2]}\n")
    f.close()

print(f"*** File stored as {out_file} ***")
print(f'*** Changing phenotype in {base_name} through plink ***')
os.system(f"plink --bfile {base_name} --pheno {out_file} --mpheno 3 --make-bed --out imp_{base_name}")
fam_file = 'imp_' + fam_file
print("*** Generating lists ***")
for cat in status_phenotype.keys():
    # NOTE in cl.txt the phenotype is in the 5th column, in the .fam in the 6th
    os.system("awk '{ if ($6 == \"" + str(status_phenotype[cat]) + "\") print $1, $2, \"" + cat + "\"}' " + fam_file + " > list_" + cat + ".txt")
    print(f"# {cat} subjects: ")
    os.system(f"cat list_{cat}.txt | wc -l")
    
try:
    prune_pheno = input("Do you want to prune missing phenotypes [Y/n]?")
except Exception:
    prune_pheno = 'Y'

if prune_pheno != 'n' and prune_pheno != 'N':
    os.system(f"plink --bfile imp_{base_name} --prune --make-bed --out prun_{base_name}")

print('*** Done ***')
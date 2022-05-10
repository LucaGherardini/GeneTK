import os
import sys 
from subprocess import Popen, PIPE, STDOUT

try:
    stem = input('Insert stem: ')
except Exception as e:
    print(e)
    quit()
    
output = Popen(f"find . -name \'{stem}*.bed\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).replace('b\'', '').replace('\'', '').replace('.bed', '').split('\\n')
print(files)
for f in files:
    os.system(f"plink --bfile {f} --snps-only --make-bed --out snp_{f}")
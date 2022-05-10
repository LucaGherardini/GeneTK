import os
import sys 
from subprocess import Popen, PIPE, STDOUT

stem = input('Insert stem: ') if len(sys.argv) < 2 else sys.argv[1]
    
output = Popen(f"find . -name \'{stem}*.bed\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).replace('b\'', '').replace("\'", '').replace('./', '').replace('.bed', '').removesuffix('\\n').split('\\n')
print(files)
for f in files:
    os.system(f"plink --bfile {f} --snps-only --make-bed --out snp_{f}")
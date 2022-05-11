import os
import sys 
from subprocess import Popen, PIPE, STDOUT

stem = input('Insert stem: ') if len(sys.argv) < 2 else sys.argv[1]
threads = int(input('Insert number of threads: ')) if len(sys.argv) < 3 else sys.argv[2]

output = Popen(f"find . -name \'{stem}*.bed\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).replace('b\'', '').replace("\'", '').replace('./', '').replace('.bed', '').removesuffix('\\n').split('\\n')
print(files)
for f in files:
    os.system(f"plink --threads {threads} --bfile {f} --snps-only --make-bed --out ext_{f}")
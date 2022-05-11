import os
import sys 
from subprocess import Popen, PIPE, STDOUT
    
part_number = 1
threads = int(input('Insert number of threads: ')) if len(sys.argv) < 3 else sys.argv[2]

output = Popen(f"find . -name \'*.vcf\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).removeprefix('b\'').removesuffix('\'').removesuffix('\\n').split('\\n')
print(files)
for f in files:
    os.system(f"plink --threads {threads} --vcf {f} --make-bed --out p_{part_number}")
    # NOTE: you can't overwrite directly the input file
    part_number += 1
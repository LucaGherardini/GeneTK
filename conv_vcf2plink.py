import os
import sys 
from subprocess import Popen, PIPE, STDOUT
    
part_number = 1

output = Popen(f"find . -name \'*.vcf\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).removeprefix('b\'').removesuffix('\'').removesuffix('\\n').split('\\n')
print(files)
for f in files:
    os.system(f"plink --vcf {f} --make-bed --out p_{part_number}")
    # NOTE: you can't overwrite directly the input file
    part_number += 1
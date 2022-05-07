import os
import sys 
from subprocess import Popen, PIPE, STDOUT
    
part_number = 0

output = Popen(f"find . -name \'*.vcf\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).removeprefix('b\'').removesuffix('\'').removesuffix('\\n').split('\\n')
print(files)
for f in files:
    os.system(f"plink --vcf {f} --make-bed --out p{part_number}")
    # NOTE: you can't overwrite directly the input file
    os.system(f"sed -r \'s/rs\S+/./\t\' p{part_number}.bim > tmp.bim && mv tmp.bim p{part_number}.bim")
    os.system(f"plink --bfile p{part_number} --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out p_{part_number}")
    os.system(f"rm p{part_number}.*")
    part_number += 1
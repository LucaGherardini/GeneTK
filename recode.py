import os
import sys 
from subprocess import Popen, PIPE, STDOUT

stem = input('Insert stem: ') if len(sys.argv) < 2 else sys.argv[1]
    
output = Popen(f"find . -name \'{stem}*.bed\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).replace('b\'', '').replace("\'", '').replace('./', '').replace('.bed', '').removesuffix('\\n').split('\\n')

print(files)
for f in files:
    os.system(f"sed -r \'s/rs\S+/./\t\' {f}.bim > tmp.bim && mv tmp.bim {f}.bim")
    os.system(f"plink --bfile {f} --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out recoded_{f}")
import os
import sys
from subprocess import Popen, PIPE, STDOUT


stem = input('Insert stem: ') if len(sys.argv) < 2 else sys.argv[1]
threads = int(input('Insert number of threads: ')) if len(sys.argv) < 3 else sys.argv[2]

output = Popen(f"find . -name \'{stem}*.bed\'", shell=True, stdout=PIPE)
files = str(output.stdout.read()).replace('b\'', '').replace("\'", '').replace('./', '').replace('.bed', '').removesuffix('\\n').split('\\n')
if len(files)< 2:
    print("Less than 2 files found, nothing to merge")
    quit()
    
merge_list = open('merge_list.txt', 'w')
for k in range(0, len(files)):
    merge_list.write(files[k] + '\n')
merge_list.close()    
    
os.system(f"plink --threads {threads} --merge-list merge_list.txt --make-bed --out merge")
print("Done! File 'merge' wrote")
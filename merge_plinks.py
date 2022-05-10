import os
from subprocess import Popen, PIPE, STDOUT
    
stem = 'Insert stem: '
output = Popen(f"find . -name \'{stem}*.bed\' | wc -l", shell=True, stdout=PIPE)
max_number = int(str(output.stdout.read()).removeprefix('b\'').removesuffix('\'').removesuffix('\\n').split('\\n')[0])
print(max_number)
base = stem+"1"
merge_number = 1
for k in range(1, max_number-1):
    os.system(f"plink --bfile {base} --bmerge {stem}{k} --make-bed --out m{merge_number}")
    if os.path.isfile('m_'+str(merge_number)+'-merge.missnp'):
        os.system(f"plink --bfile {base} --flip m_{merge_number}-merge.missnp --make-bed --out {base}_f")
        os.system(f"plink --bfile {base}_f --bmerge p_{k} --make-bed --out m_{merge_number}")

    #os.system(f"rm p_{k+1} tmp* {base}*")
    base = "m_"+str(merge_number)
    merge_number += 1
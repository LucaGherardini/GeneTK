from datetime import datetime
import os 
import sys
import multiprocessing
from subprocess import Popen, PIPE, STDOUT
import logging


def command(command):
    # if a command doesn't redirect output, redirect on log to keep terminal clean
    print(command)
    os.system(command)
    return

def remove(file):
    command(f"rm {file}.bim")
    command(f"rm {file}.bed")
    command(f"rm {file}.fam")
    return

def print_info(phase, file):
    logging.info(f'*****[{phase}]*****')
    output = Popen(f"cat {file}.fam | wc -l", shell=True, stdout=PIPE)
    output = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
    logging.info(f"Individuals: {output}")
    output = Popen(f"cat {file}.bim | wc -l", shell=True, stdout=PIPE)
    output = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
    logging.info(f"SNPs: {output}")
    
def mkdir(phase):    
    print(f"Phase {phase}")
    if not os.path.isdir(f"QC_{phase}"):
        os.system(f"mkdir QC_{phase}")

def gather_inputs():
    if len(sys.argv)<2:
        try:
            f = input("Insert baseline: ")
        except Exception as e:
            print(e)
            quit()
    else:
        f = sys.argv[1]
        
    if not os.path.isfile(f+'.bim'):
        print(f"Error: {f} doesn't exist")
        quit()

    try:
        if len(sys.argv) < 3:
            threads = int(input("Number of threads [default all available]: "))
        else:
            threads = int(sys.argv[2])
        
        if threads <= 0: raise Exception
    except Exception as e:
        threads = multiprocessing.cpu_count()
        print(f"{threads} threads will be spawn")
            
    try:
        if len(sys.argv) < 4:
            pop_f = input('Insert the auxiliary file for population stratification (PLINK format)[optional]: ')
        else:
            pop_f = sys.argv[3]
    except Exception as e:
        pop_f = ''
        
    if not os.path.isfile(pop_f+'.bim'):
        print("Population file not found, skipping...")
        pop_f = ''
        
    try:
        if len(sys.argv) < 5:
            panel_f = input('Insert the auxiliary file for population information related to population file (PLINK format)[optional]: ')
        else:
            panel_f = sys.argv[4]
    except Exception as e:
        panel_f = ''

    if not os.path.isfile(panel_f):
        print("Population information file not found, skipping...")
        panel_f = ''

    if not os.path.isfile('inversion.txt'): 
        print('Put "inversion.txt" in the parent directory before proceeding')    
        quit()
        
    rscripts = True
    if not os.path.isdir('../Rscripts/'):
        print('Set a "Rscript" folder containing R scripts used for intermediate results (optional)')
        rscripts = False
        
    try:
        rm_tmp = input("Do you want to remove intermediate PLINK file during QC(.bim/.bed/.fam)[Y/n]?")
        if len(rm_tmp) < 1 or rm_tmp.upper() != 'N': raise Exception
    except Exception as e:
        rm_tmp = 'Y'
    rm_tmp = rm_tmp.upper()     
    
    return f, threads, pop_f, panel_f, rscripts, rm_tmp

def missingness(step, rm_tmp):   
    # Analysis of Missingness
    if rscripts:
        command(f"plink --threads {threads} --bfile {f} --missing")
        command("Rscript ../Rscripts/hist_miss.R")
        command(f"mv plink.log plink.hh plink.imiss plink.lmiss histlmiss.pdf histimiss.pdf QC_{phase}")

    # Delete SNPs with missingness >0.2
    print("--- STEP " +str(step))
    print("Delete SNPs with missingness > 0.2")
    command(f"plink --threads {threads} --bfile {f} --geno 0.2 --make-bed --out {f}_{step}")
    command(f"mv {f}_{step}.hh {f}_{step}.log QC_{phase}/")

    # Delete Individuals with missingness >0.2 
    print("--- STEP " +str(step))
    print("Delete Individuals with missingness > 0.2")
    command(f"plink --threads {threads} --bfile {f}_{step} --mind 0.2 --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.log QC_{phase}/")

    # Delete SNPs with missingness >0.02
    print("--- STEP " +str(step))
    print("Delete SNPs with missingness > 0.02")
    command(f"plink --threads {threads} --bfile {f}_{step} --geno 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.log QC_{phase}/")

    # Delete individuals with missingness >0.02
    print("--- STEP " +str(step))
    print("Delete Individuals with missingness > 0.02")
    command(f"plink --threads {threads} --bfile {f}_{step} --mind 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.irem {f}_{step}.log QC_{phase}/")
    
    return step

def sex_discrepancy(step, rm_tmp):
    try:
        choice = int(input("Do you want to try to impute missing sex it [0, default]? "))
    except Exception as e:
        choice = 0
    
    if choice == 0:
        print("--- STEP " +str(step))
        print("Imputing sex discrepancies")
        command(f"plink --threads {threads} --bfile {f}_{step} --impute-sex --make-bed --out {f}_{step+1}")   
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1 
        command(f"mv {f}_{step}.log {f}_{step}.hh QC_{phase}/")
    
    # check for still ambiguous sex
    command(f"plink --threads {threads} --bfile {f}_{step} --check-sex")
    if rscripts:
        command("Rscript ../Rscripts/gender_check.R")
        command(f"mv Gender_check.pdf Men_check.pdf Women_check.pdf  QC_{phase}/")

    command("grep \"PROBLEM\" plink.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt")
    print("--- STEP " +str(step))
    print("Remotion of sex discrepancies")
    command(f"plink --threads {threads} --bfile {f}_{step} --remove sex_discrepancy.txt --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log {f}_{step}.hh sex_discrepancy.txt plink.sexcheck QC_{phase}/")
    
    return step

def maf(step, rm_tmp):
    print("--- STEP " +str(step))
    command("awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' " + f+'_'+str(step)+'.bim' + " > snp_1_22.txt")
    print("Extraction of SNPs chr 1~22")
    command(f"plink --threads {threads} --bfile {f}_{step} --extract snp_1_22.txt --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log snp_1_22.txt QC_{phase}/")
    
    try:
        maf_threshold = float(input("Insert the MAF threshold you want to use [0, default]: "))
    except Exception as e:
        maf_threshold = 0.0
        
    if rscripts:
        command(f"plink --threads {threads} --bfile {f}_{step} --freq --out MAF_check")
        command("Rscript ../Rscripts/MAF_check.R")

    if maf_threshold > 0.0:
        print("--- STEP " +str(step))
        print(f"MAF filtering {maf_threshold}")
        command(f"plink --threads {threads} --bfile {f}_{step} --maf {maf_threshold} --make-bed --out {f}_{step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
    command(f"mv {f}_{step}.log MAF_* QC_{phase}/")
    
    return step

def hwe(step, rm_tmp):
    if rscripts:
        command(f"plink --threads {threads} --bfile {f}_{step} --hardy")
        command("awk '{ if ($9 < 0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe")
        command("Rscript ../Rscripts/hwe.R")
    print("--- STEP " +str(step))
    print("HWE filtering (p<=1e-10)")
    command(f"plink --threads {threads} --bfile {f}_{step} --hwe 1e-10 include-nonctrl --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log histhwe* plink.hwe plinkzoomhwe.hwe QC_{phase}/")
    
    return step

def heterozygosity(step, rm_tmp):
    command(f"plink --threads {threads} --bfile {f}_{step} --exclude range inversion.txt --indep-pairwise 50 5 0.2 --out indepSNP")
    command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --het --out R_check")
    if rscripts:
        command("Rscript ../Rscripts/check_heterozygosity_rate.R")
        command("Rscript ../Rscripts/heterozygosity_outliers_list.R")
    command("sed 's/\"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt")
    print("--- STEP " +str(step))
    print("Heterozygosity filtering (- 3 SD from mean)")
    command(f"plink --threads {threads} --bfile {f}_{step} --remove het_fail_ind.txt --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1 
    command(f"mv {f}_{step}.log fail-het-qc.txt heterozygosity* het_fail_ind.txt indepSNP.log indepSNP.prune.out R_check* QC_{phase}/")
    
    return step

def cryptic_relatedness(step, rm_tmp):
    command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2")
    command("awk '{ if ($8 > 0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome")
    if rscripts:
        command("Rscript ../Rscripts/Relatedness.R")
        
    print("--- STEP " +str(step))
    print("Cryptic relatedness filtration")
    command(f"plink --threads {threads} --bfile {f}_{step} --filter-founders --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders")
    command(f"plink --threads {threads} --bfile {f}_{step} --missing")
    command("awk '{ print $1, $2, $3, $4}' pihat_min0.2_in_founders.genome>pairs.txt")
    pairs = open("pairs.txt", "r")
    low_call = open("0.2_low_call_rate_pihat.txt", "w")
    for line in pairs.readlines():
        # check you're not reading the header
        if 'IID1' in line:
            continue
        
        sub1_fam, sub1_id, sub2_fam, sub2_id = line.split(' ')
        output = Popen(f"cat plink.imiss | grep {sub1_fam+sub1_id}", shell=True, stdout=PIPE)
        sub1_fmiss = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
        
        output = Popen(f"cat plink.imiss | grep {sub2_fam+sub2_id}", shell=True, stdout=PIPE)
        sub2_fmiss = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
        if sub1_fmiss > sub2_fmiss:
            low_call.write(sub1_fam + ' ' + sub1_id)
        else:
            low_call.write(sub2_fam + ' ' + sub2_id)
        low_call.write('\n')
    low_call.close()    

    command(f"mv pihat* relatedness.pdf zoom_pihat.genome {f}_{step}.log plink.* pairs.txt QC_{phase}/")

    print("--- STEP " +str(step))
    print("Remotion of related individuals with lowest call rates")
    command(f"plink --threads {threads} --bfile {f}_{step} --remove 0.2_low_call_rate_pihat.txt --make-bed --out {f}_{step+1}")
    if rm_tmp == 'Y': remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log  0.2_low_call_rate_pihat.txt QC_{phase}/")
    
    return step

def population_stratification(step, rm_tmp):
    pop_step = 1
    original_f = f+'_'+str(step)
    
    # NOTE: without R scripts, it's useless to merge different populations
    if rscripts:
        print("Preprocessing of population file")
        command(f"plink --threads {threads} --bfile {pop_f} --geno 0.2 --allow-no-sex --make-bed --out {pop_f}_{pop_step}")
        command(f"mv {pop_f}_{pop_step}.log QC_{phase}/")
        
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.2 --allow-no-sex --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log QC_{phase}/")
        
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --geno 0.02 --allow-no-sex --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log QC_{phase}/")
        
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.02 --allow-no-sex --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log QC_{phase}/")
        
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --maf 0.05 --allow-no-sex --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log QC_{phase}/")
        
        print("Uniforming SNPs between the two populations")
        command("awk '{print$2}' " + f + '_' + str(step) + ".bim > " + f + '_' + str(step) + "_SNPs.txt")
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --extract {f}_{step}_SNPs.txt --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {f}_{step}_SNPs.txt {pop_f}_{pop_step}.log QC_{phase}/")
        
        command("awk '{print$2}' " + pop_f + '_' + str(pop_step) + ".bim > " + pop_f + '_' + str(pop_step) + '_SNPs.txt')
        command(f"plink --threads {threads} --bfile {f}_{step} --extract {pop_f}_{pop_step}_SNPs.txt --recode --make-bed --out {f}_{step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log {pop_f}_{pop_step}_SNPs.txt QC_{phase}/")
        
        print("Updating build")
        command("awk '{print$2,$4}' " + f + '_' + str(step) + '.map > buildmap.txt')
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --update-map buildmap.txt --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log buildmap.txt QC_{phase}/")
        
        print("Setting reference genome")
        command("awk '{print$2,$5}' " + pop_f + '_' + str(pop_step) + '.bim > pop_ref_list.txt')
        command(f"plink --threads {threads} --bfile {f}_{step} --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")    
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log QC_{phase}/")
        
        print("Checking for potential strand issues")
        command("awk '{print$2,$5,$6}' " + pop_f + '_' + str(pop_step) + '.bim > ' + pop_f + '_tmp')
        command("awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
        command(f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > all_differences.txt")
        command(f"mv {f}_{step}_tmp QC_{phase}/")
        
        print("Flipping SNPs for resolving strand issues")
        command("awk '{print$1}' all_differences.txt | sort -u > flip_list.txt")
        command(f"plink --threads {threads} --bfile {f}_{step} --flip flip_list.txt --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        command(f"mv all_differences.txt pop_ref_list.txt flip_list.txt {f}_{step}.log QC_{phase}/")
        
        print("Checking still problematic SNPs")
        command("awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
        command(f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > uncorresponding_SNPs.txt")
        command(f"mv {f}_{step}_tmp QC_{phase}/")
        
        print("Removing problematic SNPs")
        command("awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt")
        command(f"plink --threads {threads} --bfile {f}_{step} --exclude SNPs_for_exclusion.txt --make-bed --out {f}_{step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        command(f"plink --threads {threads} --bfile {pop_f}_{pop_step} --exclude SNPs_for_exclusion.txt --make-bed --out {pop_f}_{pop_step+1}")
        command(f"mv {f}_{step}.log uncorresponding_SNPs.txt SNPs_for_exclusion.txt QC_{phase}/")
        if rm_tmp == 'Y': remove(f"{f}_{pop_step}")
        pop_step += 1
        
        print("Merge")
        command(f"plink --threads {threads} --bfile {f}_{step} --bmerge {pop_f}_{pop_step}.bed {pop_f}_{pop_step}.bim {pop_f}_{pop_step}.fam --allow-no-sex --make-bed --out {f}_{step+1}")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log QC_{phase}/")
        
        print("MDS")
        # NOTE: step is not incremented because it's just producing genome file
        command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
        command(f"plink --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
        command(f"mv {f}_{step}.log QC_{phase}/")
        
        print("MDS Plot")
        panel_step = 1
        command("awk '{print$1,$1,$2}' " + panel_f + ' > race_1kG.txt')
        command(f" sed 's/JPT/ASN/g' race_1kG.txt > race_1kG{panel_step}.txt")
        command(f" sed 's/ASW/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/CEU/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/CHB/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/CHD/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/YRI/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/LWK/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/TSI/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/MXL/AMR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/GBR/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/FIN/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/CHS/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        command(f" sed 's/PUR/AMR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
        panel_step += 1
        
        # Racefile for original population
        command("awk '{print$1,$2,\"OWN\"}' " + f + '_' + str(step) + '.fam > racefile_own.txt')

        # Concatenate racefiles
        command(f"cat race_1kG{panel_step}.txt racefile_own.txt | sed -e '1i\\FID IID race' > racefile.txt")
    
        command(f"Rscript ../Rscripts/MDS_merged.R {f}_{step}.mds pop.pdf")
        print('Population MDS displayed in "pop.pdf" file, please check and specify desired cut-offs')
        while True:
            try:
                low_lim = float(input(" Insert upper threshold on MDS component 1: "))
                high_lim = float(input("Insert lower threshold on MDS component 2: "))
            except Exception as e:
                print(e)
                print("Retry")
                continue
            break
        command("awk '{if($4 < " + str(low_lim) + " && $5 > " + str(high_lim) + ") print $1,$2}' " + f + '_' + str(step) + ".mds > filtered_subj")  
        command(f"mv pop.pdf {f}_{step}.mds {f}_{step}.cluster* {f}_{step}.genome {f}_{step}.nosex QC_{phase}/")  
        
        print("Filtering population from the original dataset")
        # Remove individuals from the pre-merged dataset (original_f)
        command(f"plink --threads {threads} --bfile {original_f} --keep filtered_subj --make-bed --out {f}_{step+1}")
        command(f"mv filtered_subj QC_{phase}/")
        if rm_tmp == 'Y': remove(f"{f}_{step}")
        step += 1
        
        # Showing updated (original) population
        command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
        command(f"plink --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
        command("awk '{print$1,$2,\"OWN\"}' " + f + '_' + str(step) + '.fam > racefile_own.txt')
        command(f"cat race_1kG{panel_step}.txt racefile_own.txt | sed -e '1i\\FID IID race' > racefile.txt")
        command(f"Rscript ../Rscripts/MDS_merged.R {f}_{step}.mds filtered.pdf")
        command(f"mv {f}_{step}.mds {f}_{step}.log {f}_{step}.cluster* filtered.pdf {f}_{step}.mds race*txt QC_{phase}/")
        print('Filtered original population displayed in "filtered.pdf" file')
    
    return step

def mds_covariates(step):
    # NOTE: producing genome the step counter is not incremented 
    command(f"plink --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
    command(f"plink --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
    command("awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' " + f + '_' + str(step) + '.mds > covar_mds.txt')
    command(f"mv {f}_{step}.cluster* {f}_{step}.genome {f}_{step}.log {f}_{step}.mds QC_{phase}/")
    
    return step

if __name__ == '__main__':
    start_time = datetime.today()
    logging.basicConfig(format='%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s', datefmt='%Y-%m-%d,%H:%M:%S', level=logging.INFO, filename = f"trace_{start_time.strftime('%Y-%m-%d-%H:%M:%S')}.log")
    f, threads, pop_f, panel_f, rscripts, rm_tmp = gather_inputs()

    step = 1

    ## PHASE 1 (Missingness)
    phase = 1
    mkdir(phase)
    step = missingness(step, rm_tmp)
    print_info(phase, f"{f}_{step}")

    ## PHASE 2 (Sex discrepancy)
    phase = 2
    mkdir(phase)
    step = sex_discrepancy(step, rm_tmp)
    print_info(phase, f"{f}_{step}")


    ## PHASE 3 (MAF filtration)
    phase = 3
    mkdir(phase)
    step = maf(step, rm_tmp)
    print_info(phase, f"{f}_{step}")

    ## PHASE 4 (Hardy-Weinberg Equilibrium)
    phase = 4
    mkdir(phase)
    step = hwe(step, rm_tmp)
    print_info(phase, f"{f}_{step}")

    ## PHASE 5 (Heterozygosity rate)
    phase = 5
    mkdir(phase)
    step = heterozygosity(step, rm_tmp)
    print_info(phase, f"{f}_{step}")

    ## PHASE 6 (Cryptic Relatedness)
    phase = 6
    mkdir(phase)
    step = cryptic_relatedness(step, rm_tmp)
    print_info(phase, f"{f}_{step}")

    if len(pop_f)>0:
        ## PHASE 7 (Population Stratification)
        phase = 7
        mkdir(phase)
        step = population_stratification(step, rm_tmp)
        print_info(phase, f"{f}_{step}")
            
    ## PHASE 8 (MDS Covariates)
    phase = 8
    mkdir(phase)
    step = mds_covariates(step, rm_tmp)  
    print_info(phase, f"{f}_{step}")
        
    command(f"mv {f}_{step}.fam {f}-QC.fam")
    command(f"mv {f}_{step}.bim {f}-QC.bim")
    command(f"mv {f}_{step}.bed {f}-QC.bed")
    print(f'QC completed, the final file is {f}-QC')
    
    if rm_tmp == 'Y' and len(pop_f)>0:
        # created during merging
        command(f"rm {f}_*.map")
        command(f"rm {f}_*.ped")
        command(f"rm {pop_f}_*")          
    
    print("DONE! You can proceed to imputation with 'imp_pheno.py'")
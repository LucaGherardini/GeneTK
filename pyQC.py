"""
SYNOPSYS:
python3 pyQC.py <file> <thrds> <pop_file> <pop_panel> <discard_sex> <MAF_thr> <rm_tmp>
"""

from datetime import datetime
import os
import sys
import multiprocessing
from subprocess import Popen, PIPE, STDOUT
import logging
import yaml

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
    if not os.path.isdir(f"Phase_{phase}"):
        os.system(f"mkdir Phase_{phase}")

def collect_inputs():
    if not os.path.isfile('inversion.txt'): 
        print('Put "inversion.txt" in the parent directory before proceeding')    
        quit()
   
    rscripts = True
    if not os.path.isdir('Rscripts'):
        print('Set a "Rscript" folder containing R scripts used for intermediate results (optional)')
        rscripts = False    
         
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    
    try:
        f = config['baseline']
    except Exception as e:
        print(e)
        print('Error during reading of baseline name')
        quit()
    if not os.path.isfile(f+'.bim'):
        print(f"Error: {f} doesn't exist")
        quit()    
    
    try:
        threads = config['cores']
    except Exception as e:
        print(e)
        print('Error during reading of cores number, using default (all available)')
        threads = multiprocessing.cpu_count()
            
    try:
        pop_f = config['population_file']
        panel_f = config['population_panel']
    except Exception as e:
        print(e)
        print('Error when reading population files')
        pop_f = ''
        panel_f = ''
        
    if not os.path.isfile(pop_f+'.bim') or not os.path.isfile(panel_f):
        print("Population files not found, skipping...")
        pop_f = ''
        panel_f = ''
        
    if (pop_f == '' or panel_f == '') and config['Population_filtering']:
        print('You requested population filtering but files have not been found, aborting...')
        quit()

    try:
        discard_sex = config['sex_discrepancy_discard']
    except Exception as e:
        print(e)
        print('Error when reading "sex_discrepancy_discard" option, using False (default)')
        discard_sex = False
   
    try:
        maf_threshold = float(config['MAF_threshold'])
    except Exception as e:
        print(e)
        print('Error when reading "MAF_threshold" option, using "0.01" (default)')
        maf_threshold = 0.01
       
    try:
        rm_tmp = bool(config['rm_tmp'])
    except Exception as e:
        rm_tmp = True
        
    try:
        sup_lim = float(config['superior_threshold_MDS_1'])
        inf_lim = float(config['inferior_threshold_MDS_2'])
    except Exception as e:
        print(e)
        print("MDS thresholds will be chosen at runtime")
        sup_lim = ''
        inf_lim = ''
    
    logging.info('****************')
    logging.info(f"Baseline: {f}")
    logging.info(f"Cores: {threads}")
    logging.info(f"Population file: {pop_f}")
    logging.info(f"Population panel: {panel_f}")
    logging.info(f"Discard sex: {discard_sex}")
    logging.info(f"MAF threshold: {maf_threshold}")
    logging.info(f"Remove temporary files: {rm_tmp}")
    logging.info(f"Superior and inferior thresholds: [{sup_lim}, {inf_lim}]")
    logging.info('****************')
    logging.info(f"Missingness: {config['Missingness']}")
    logging.info(f"Sex discrepancy: {config['Sex_discrepancy']}")
    logging.info(f"Autosomal SNPs and MAF filtering: {config['Autos_MAF']}")
    logging.info(f"Hardy-Weinberg Equilibrium filtering: {config['Hardy-Weinberg']}")
    logging.info(f"Heterozygosity rate: {config['Heterozygosity_rate']}")
    logging.info(f"Cryptic relatedness: {config['Cryptic_relatedness']}")
    logging.info(f"Population filtering: {config['Population_filtering']}")
    logging.info(f"MDS_covariates: {config['MDS_covariates']}")
    
    return f, threads, pop_f, panel_f, discard_sex, maf_threshold, rscripts, rm_tmp, sup_lim, inf_lim, config

def missingness(step, rm_tmp):   
    # Analysis of Missingness
    if rscripts:
        command(f"plink --threads {threads} --bfile {f}_{step} --missing")
        command("Rscript Rscripts/hist_miss.R")
        command(f"mv plink.log plink.hh plink.imiss plink.lmiss histlmiss.pdf histimiss.pdf Phase_{phase}")

    # Delete SNPs with missingness >0.2
    print("--- STEP " +str(step))
    print("Delete SNPs with missingness > 0.2")
    command(f"plink --threads {threads} --bfile {f}_{step} --geno 0.2 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.log Phase{phase}/")

    # Delete Individuals with missingness >0.2 
    print("--- STEP " +str(step))
    print("Delete Individuals with missingness > 0.2")
    command(f"plink --threads {threads} --bfile {f}_{step} --mind 0.2 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.log Phase_{phase}/")

    # Delete SNPs with missingness >0.02
    print("--- STEP " +str(step))
    print("Delete SNPs with missingness > 0.02")
    command(f"plink --threads {threads} --bfile {f}_{step} --geno 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.log Phase_{phase}/")

    # Delete individuals with missingness >0.02
    print("--- STEP " +str(step))
    print("Delete Individuals with missingness > 0.02")
    command(f"plink --threads {threads} --bfile {f}_{step} --mind 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.hh {f}_{step}.irem {f}_{step}.log Phase_{phase}/")
    
    return step

def sex_discrepancy(step, rm_tmp, discard_sex):
    # impute sex    
    print("--- STEP " +str(step))
    print("Imputing sex discrepancies")
    command(f"plink --threads {threads} --bfile {f}_{step} --impute-sex --make-bed --out {f}_{step+1}")   
    if rm_tmp: remove(f"{f}_{step}")
    step += 1 
    command(f"mv {f}_{step}.log {f}_{step}.hh Phase{phase}/")
    
    # Sex is checked after imputation or before discarding discrepancies    
    command(f"plink --threads {threads} --bfile {f}_{step} --check-sex")
    if rscripts:
        command("Rscript Rscripts/gender_check.R")
        command(f"mv Gender_check.pdf Men_check.pdf Women_check.pdf Phase_{phase}/")
    
    # check for ambiguous sex
    if discard_sex:
        command("grep \"PROBLEM\" plink.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt")
        print("--- STEP " +str(step))
        print("Remotion of sex discrepancies")
        command(f"plink --threads {threads} --bfile {f}_{step} --remove sex_discrepancy.txt --make-bed --out {f}_{step+1}")
        
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log {f}_{step}.hh sex_discrepancy.txt plink.sexcheck Phase_{phase}/")
    
    return step

def maf(step, rm_tmp, maf_threshold):
    print("--- STEP " +str(step))
    command("awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' " + f+'_'+str(step)+'.bim' + " > snp_1_22.txt")
    print("Extraction of SNPs chr 1~22")
    command(f"plink --allow-no-sex --threads {threads} --bfile {f}_{step} --extract snp_1_22.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log snp_1_22.txt Phase_{phase}/")
        
    if rscripts:
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --freq --out MAF_check")
        command("Rscript Rscripts/MAF_check.R")

    if maf_threshold > 0.0:
        print("--- STEP " +str(step))
        print(f"MAF filtering {maf_threshold}")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --maf {maf_threshold} --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
    command(f"mv {f}_{step}.log MAF_* Phase_{phase}/")
    
    return step

def hwe(step, rm_tmp):
    if rscripts:
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --hardy")
        command("awk '{ if ($9 < 0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe")
        command("Rscript Rscripts/hwe.R")
    print("--- STEP " +str(step))
    print("HWE filtering (p<=1e-10)")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --hwe 1e-10 include-nonctrl --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log histhwe* plink.hwe plinkzoomhwe.hwe Phase_{phase}/")
    
    return step

def heterozygosity(step, rm_tmp):
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --exclude range inversion.txt --indep-pairwise 50 5 0.2 --out indepSNP")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --het --out R_check")
    if rscripts:
        command("Rscript Rscripts/check_heterozygosity_rate.R")
        command("Rscript Rscripts/heterozygosity_outliers_list.R")
    command("sed 's/\"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt")
    print("--- STEP " +str(step))
    print("Heterozygosity filtering (- 3 SD from mean)")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --remove het_fail_ind.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1 
    command(f"mv {f}_{step}.log fail-het-qc.txt heterozygosity* het_fail_ind.txt indepSNP.log indepSNP.prune.out R_check* Phase_{phase}/")
    
    return step

def cryptic_relatedness(step, rm_tmp):
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2")
    command("awk '{ if ($8 > 0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome")
    if rscripts:
        command("Rscript Rscripts/Relatedness.R")
        
    print("--- STEP " +str(step))
    print("Cryptic relatedness filtering")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --filter-founders --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --missing")
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

    command(f"mv pihat* relatedness.pdf zoom_pihat.genome {f}_{step}.log plink.* pairs.txt Phase_{phase}/")

    print("--- STEP " +str(step))
    print("Remotion of related individuals with lowest call rates")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --remove 0.2_low_call_rate_pihat.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(f"{f}_{step}")
    step += 1
    command(f"mv {f}_{step}.log  0.2_low_call_rate_pihat.txt Phase_{phase}/")
    
    return step

def population_stratification(step, rm_tmp, sup_lim, inf_lim):
    pop_step = 1
    original_f = f+'_'+str(step)
    # a copy is made to maintain the original untouched
    os.system(f"cp {f}_{step}.bed {f}_{step+1}.bed")
    os.system(f"cp {f}_{step}.bim {f}_{step+1}.bim")
    os.system(f"cp {f}_{step}.fam {f}_{step+1}.fam")
    step += 1
    
    # NOTE: without R scripts, it's useless to merge different populations
    if rscripts:
        print("Preprocessing of population file")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f} --geno 0.2 --make-bed --out {pop_f}_{pop_step}")
        command(f"mv {pop_f}_{pop_step}.log Phase{phase}/")
        
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.2 --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log Phase_{phase}/")
        
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --geno 0.02 --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log Phase_{phase}/")
        
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.02 --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log Phase_{phase}/")
        
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --maf 0.05 --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log Phase_{phase}/")
        
        print("Uniforming SNPs between the two populations")
        command("awk '{print$2}' " + f + '_' + str(step) + ".bim > " + f + '_' + str(step) + "_SNPs.txt")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --extract {f}_{step}_SNPs.txt --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {f}_{step}_SNPs.txt {pop_f}_{pop_step}.log Phase_{phase}/")
        
        command("awk '{print$2}' " + pop_f + '_' + str(pop_step) + ".bim > " + pop_f + '_' + str(pop_step) + '_SNPs.txt')
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract {pop_f}_{pop_step}_SNPs.txt --recode --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log {pop_f}_{pop_step}_SNPs.txt Phase_{phase}/")
        
        print("Updating build")
        command("awk '{print$2,$4}' " + f + '_' + str(step) + '.map > buildmap.txt')
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --update-map buildmap.txt --make-bed --out {pop_f}_{pop_step+1}")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        command(f"mv {pop_f}_{pop_step}.log buildmap.txt Phase_{phase}/")
        
        print("Setting reference genome")
        command("awk '{print$2,$5}' " + pop_f + '_' + str(pop_step) + '.bim > pop_ref_list.txt')
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")    
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log Phase_{phase}/")
        
        print("Checking for potential strand issues")
        command("awk '{print$2,$5,$6}' " + pop_f + '_' + str(pop_step) + '.bim > ' + pop_f + '_tmp')
        command("awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
        command(f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > all_differences.txt")
        command(f"mv {f}_{step}_tmp Phase_{phase}/")
        
        print("Flipping SNPs for resolving strand issues")
        command("awk '{print$1}' all_differences.txt | sort -u > flip_list.txt")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --flip flip_list.txt --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"mv all_differences.txt pop_ref_list.txt flip_list.txt {f}_{step}.log Phase_{phase}/")
        
        print("Checking still problematic SNPs")
        command("awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
        command(f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > uncorresponding_SNPs.txt")
        command(f"mv {f}_{step}_tmp Phase_{phase}/")
        
        print("Removing problematic SNPs")
        command("awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --exclude SNPs_for_exclusion.txt --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"plink --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --exclude SNPs_for_exclusion.txt --make-bed --out {pop_f}_{pop_step+1}")
        command(f"mv {pop_f}_{pop_step}.log uncorresponding_SNPs.txt SNPs_for_exclusion.txt Phase_{phase}/")
        if rm_tmp: remove(f"{pop_f}_{pop_step}")
        pop_step += 1
        
        print("Merge")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --bmerge {pop_f}_{pop_step}.bed {pop_f}_{pop_step}.bim {pop_f}_{pop_step}.fam --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(f"{f}_{step}")
        step += 1
        command(f"mv {f}_{step}.log Phase_{phase}/")
        
        print("MDS")
        # NOTE: step is not incremented because it's just producing genome file
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
        command(f"mv {f}_{step}.log Phase_{phase}/")
        
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
    
        command(f"Rscript Rscripts/MDS_merged.R {f}_{step}.mds pop.pdf")
        print('Population MDS displayed in "pop.pdf" file')
        if sup_lim == '' or inf_lim == '':
            print('Please check and specify desired cut-offs')
            while True:
                try:
                    sup_lim = float(input(" Insert superior threshold on MDS component 1: "))
                    inf_lim = float(input("Insert inferior threshold on MDS component 2: "))
                except Exception as e:
                    print(e)
                    print("Retry")
                    continue
                break
            logging.info(f"Superior and inferior thresholds: [{sup_lim}, {inf_lim}]")
        
        command("awk '{if($4 < " + str(sup_lim) + " && $5 > " + str(inf_lim) + ") print $1,$2}' " + f + '_' + str(step) + ".mds > filtered_subj")  
        command(f"mv pop.pdf {f}_{step}.mds {f}_{step}.cluster* {f}_{step}.genome {f}_{step}.nosex Phase_{phase}/")  
        
        print("Filtering population from the original dataset")
        # Remove individuals from the pre-merged dataset (original_f)
        command(f"plink --allow-no-sex  --threads {threads} --bfile {original_f} --keep filtered_subj --make-bed --out {f}_{step+1}")
        command(f"mv filtered_subj Phase_{phase}/")
        if rm_tmp: 
            remove(f"{f}_{step}")
            remove(f"{original_f}")
        step += 1
        
        # Showing updated (original) population
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
        command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
        command("awk '{print$1,$2,\"OWN\"}' " + f + '_' + str(step) + '.fam > racefile_own.txt')
        command(f"cat race_1kG{panel_step}.txt racefile_own.txt | sed -e '1i\\FID IID race' > racefile.txt")
        command(f"Rscript Rscripts/MDS_merged.R {f}_{step}.mds filtered.pdf")
        command(f"mv {f}_{step}.mds {f}_{step}.log {f}_{step}.cluster* filtered.pdf {f}_{step}.mds race*txt Phase_{phase}/")
        print('Filtered original population displayed in "filtered.pdf" file')
    
    return step

def mds_covariates(step):
    # NOTE: producing genome the step counter is not incremented 
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
    command(f"plink --allow-no-sex  --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
    command("awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' " + f + '_' + str(step) + '.mds > covar_mds.txt')
    command(f"mv {f}_{step}.cluster* covar_mds.txt {f}_{step}.genome {f}_{step}.log {f}_{step}.mds Phase_{phase}/")
    
    return step

if __name__ == '__main__':
    start_time = datetime.today()
    logging.basicConfig(format='%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s', datefmt='%Y-%m-%d,%H:%M:%S', level=logging.INFO, filename = f"trace_{start_time.strftime('%Y-%m-%d-%H:%M:%S')}.log")
    f, threads, pop_f, panel_f, discard_sex, maf_threshold, rscripts, rm_tmp, sup_lim, inf_lim, config = collect_inputs()

    step = 1
    
    # this ensures a valid input file for whatever first step will be done
    os.system(f'cp {f}.bim {f}_{step}.bim')
    os.system(f'cp {f}.bed {f}_{step}.bed')
    os.system(f'cp {f}.fam {f}_{step}.fam')

    # PHASE 1 (Missingness)
    if config['Missingness'] == True:
        phase = 1
        mkdir(phase)
        step = missingness(step, rm_tmp)
        print_info(phase, f"{f}_{step}")

    # PHASE 2 (Sex discrepancy)
    if config['Sex_discrepancy'] == True:
        phase = 2
        mkdir(phase)
        step = sex_discrepancy(step, rm_tmp, discard_sex)
        print_info(phase, f"{f}_{step}")

    # PHASE 3 (MAF filtering)
    if config['Autos_MAF'] == True:
        phase = 3
        mkdir(phase)
        step = maf(step, rm_tmp, maf_threshold)
        print_info(phase, f"{f}_{step}")

    # PHASE 4 (Hardy-Weinberg Equilibrium)
    if config['Hardy-Weinberg'] == True:
        phase = 4
        mkdir(phase)
        step = hwe(step, rm_tmp)
        print_info(phase, f"{f}_{step}")

    # PHASE 5 (Heterozygosity rate)
    if config['Heterozygosity_rate'] == True:
        phase = 5
        mkdir(phase)
        step = heterozygosity(step, rm_tmp)
        print_info(phase, f"{f}_{step}")

    # PHASE 6 (Cryptic relatedness)
    if config['Cryptic_relatedness'] == True:
        phase = 6
        mkdir(phase)
        step = cryptic_relatedness(step, rm_tmp)
        print_info(phase, f"{f}_{step}")

    if len(pop_f)>0:
        # PHASE 7 (Population filtering)
        if config['Population_filtering'] == True:
            phase = 7
            mkdir(phase)
            step = population_stratification(step, rm_tmp, sup_lim, inf_lim)
            print_info(phase, f"{f}_{step}")
            
    # PHASE 8 (MDS Covariates)
    if config['MDS_covariates'] == True:
        phase = 8
        mkdir(phase)
        step = mds_covariates(step)  
        print_info(phase, f"{f}_{step}")
        
    command(f"mv {f}_{step}.fam QC_{f}.fam")
    command(f"mv {f}_{step}.bim QC_{f}.bim")
    command(f"mv {f}_{step}.bed QC_{f}.bed")
    print(f'QC completed, the final files are QC_{f}')
    command(f"rm {f}_*.sex")
    command(f"rm {f}_*.sexcheck")
    command(f"rm {f}*log")
    
    if rm_tmp and len(pop_f)>0:
        # created during merging
        command(f"rm {f}_*.map")
        command(f"rm {f}_*.ped")
        command(f"rm {pop_f}_*.*")  
        command(f"rm {pop_f}_tmp")           
    
    print("DONE! You can proceed to imputation with 'imp_pheno.py'")
"""
SYNOPSYS:
python3 pyQC.py <file> <thrds> <pop_file> <pop_panel> <discard_sex> <MAF_thr> <rm_tmp>
"""

from datetime import datetime
import os
import multiprocessing
import subprocess
from subprocess import Popen, PIPE, STDOUT
import logging
import yaml

def command(step, command, exit_on_error=True):
    # if a command doesn't redirect output, redirect on log to keep terminal clean
    print(command)
    logging.info(f"Step {step}: {command}")
    try:
        p = subprocess.run(command, shell=True, capture_output=True)
    except Exception as e:
        logging.error(f"During command {command} an error occurred: {e}")
        return
    
    if exit_on_error and p.returncode != 0:
        logging.error(f"Error: {p.stderr.decode('utf-8')}")
        print(f"Error: {p.stderr.decode('utf-8')}")
        quit()       
    
    return

def remove(step, file):
    command(step, f"rm {file}.bim", False)
    command(step, f"rm {file}.bed", False)
    command(step, f"rm {file}.fam", False)
    return

def print_info(phase, file):
    output = Popen(f"cat {file}.fam | wc -l", shell=True, stdout=PIPE)
    output = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
    logging.info(f"Individuals: {output}")
    print(f"Individuals: {output}")
    output = Popen(f"cat {file}.bim | wc -l", shell=True, stdout=PIPE)
    output = str(output.stdout.read()).lstrip('b\'').rstrip('\'').rstrip('\\n').split(' ')[-1]
    logging.info(f"SNPs: {output}")
    print(f"SNPs: {output}")
    logging.info(f'*****[Phase {phase} is over]*****')
    print(f'*****[Phase {phase} is over]*****')

def mkdir(phase):
    print(f"### Phase {phase}")
    if not os.path.isdir(f"Phase_{phase}"):
        os.system(f"mkdir Phase_{phase}")

def collect_inputs():
    if not os.path.isfile('inversion.txt'): 
        print('Put "inversion.txt" in the parent directory before proceeding')    
        quit()
   
    rscripts = True
    if not os.path.isdir('Rscripts'):
        print('You chose to use R scripts, but a folder named "Rscripts" has not been found in the current position. Set a "Rscript" folder containing R scripts used for intermediate results')
        quit()    
         
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
        threads = int(config['threads'])
    except Exception as e:
        print(e)
        print('Error during reading of threads number, using default (all cores available)')
        threads = multiprocessing.cpu_count()

    try:
        memory = int(config['memory'])
    except Exception as e:
        print(e)
        print('Error during reading of memory, using default (16384 MB = 16 GB)')
        memory = 16384
            
    try:
        pop_f = config['population_file']
        panel_f = config['population_panel']
    except Exception as e:
        print(e)
        print('Error when reading population files')
        pop_f = ''
        panel_f = ''
        
    if not os.path.isfile(pop_f+'.bim') or not os.path.isfile(panel_f):
        print("### Population files not found, skipping...")
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
        inf_lim_mds_1 = float(config['inferior_threshold_MDS_1'])
        sup_lim_mds_1 = float(config['superior_threshold_MDS_1'])
        inf_lim_mds_2 = float(config['inferior_threshold_MDS_2'])
        sup_lim_mds_2 = float(config['superior_threshold_MDS_2'])
    except Exception as e:
        print(e)
        print("### MDS thresholds will be chosen at runtime")
        inf_lim_mds_1 = ''
        sup_lim_mds_1 = ''
        inf_lim_mds_2 = ''
        sup_lim_mds_2 = ''
    
    logging.info('****************')
    logging.info(f"Baseline: {f}")
    logging.info(f"Cores: {threads}")
    logging.info(f"Population file: {pop_f}")
    logging.info(f"Population panel: {panel_f}")
    logging.info(f"Discard sex: {discard_sex}")
    logging.info(f"MAF threshold: {maf_threshold}")
    logging.info(f"Remove temporary files: {rm_tmp}")
    logging.info(f"Inferior and Superior thresholds on MDS component 1: [{inf_lim_mds_1}, {sup_lim_mds_1}]")
    logging.info(f"Inferior and Superior thresholds on MDS component 2: [{inf_lim_mds_2}, {sup_lim_mds_2}]")
    logging.info('****************')
    logging.info(f"Missingness: {config['Missingness']}")
    logging.info(f"Sex discrepancy: {config['Sex_discrepancy']}")
    logging.info(f"Autosomal SNPs and MAF filtering: {config['Autos_MAF']}")
    logging.info(f"Hardy-Weinberg Equilibrium filtering: {config['Hardy-Weinberg']}")
    logging.info(f"Heterozygosity rate: {config['Heterozygosity_rate']}")
    logging.info(f"Cryptic relatedness: {config['Cryptic_relatedness']}")
    logging.info(f"Population filtering: {config['Population_filtering']}")
    logging.info(f"MDS_covariates: {config['MDS_covariates']}")
    
    return f, threads, memory, pop_f, panel_f, discard_sex, maf_threshold, rscripts, rm_tmp, inf_lim_mds_1, sup_lim_mds_1, inf_lim_mds_2, sup_lim_mds_2, config

def missingness(step, rm_tmp):   
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
        
    # Analysis of Missingness
    if rscripts:
        command(step, f"plink --memory {memory} --threads {threads} --bfile {base} --missing")
        command(step, "Rscript Rscripts/hist_miss.R", False)
        command(step, f"mv plink.log plink.hh plink.imiss plink.lmiss histlmiss.pdf histimiss.pdf Phase_{phase}", False)    

    # Delete SNPs with missingness >0.2
    command(step, f"plink --memory {memory} --threads {threads} --bfile {base} --geno 0.2 --make-bed --out {f}_{step}")
    command(step, f"mv {f}_{step}.hh {f}_{step}.log Phase_{phase}/", False)
    if rm_tmp: remove(step, f"{base}")

    # Delete Individuals with missingness >0.2 
    command(step, f"plink --memory {memory} --threads {threads} --bfile {f}_{step} --mind 0.2 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.hh {f}_{step}.log Phase_{phase}/", False)

    # Delete SNPs with missingness >0.02
    command(step, f"plink --memory {memory} --threads {threads} --bfile {f}_{step} --geno 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.hh {f}_{step}.log Phase_{phase}/", False)

    # Delete individuals with missingness >0.02
    command(step, f"plink --memory {memory} --threads {threads} --bfile {f}_{step} --mind 0.02 --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.hh {f}_{step}.irem {f}_{step}.log Phase_{phase}/", False)
    
    return step

def sex_discrepancy(step, rm_tmp, discard_sex):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
        
    # impute sex    
    command(step, f"plink --memory {memory} --threads {threads} --bfile {base} --impute-sex --make-bed --out {f}_{step}") 
    if rm_tmp: remove(step, f"{base}")  
    command(step, f"mv {f}_{step}.log {f}_{step}.hh Phase_{phase}/", False)
    
    # Sex is checked after imputation or before discarding discrepancies 
    command(step, f"plink --memory {memory} --threads {threads} --bfile {f}_{step} --check-sex")
    if rscripts:
        command(step, "Rscript Rscripts/gender_check.R", False)
        command(step, f"mv Gender_check.pdf Men_check.pdf Women_check.pdf Phase_{phase}/", False)
    
    # check for ambiguous sex
    if discard_sex:
        command(step, "grep \"PROBLEM\" plink.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt")
        command(step, f"plink --memory {memory} --threads {threads} --bfile {f}_{step} --remove sex_discrepancy.txt --make-bed --out {f}_{step+1}")
        
        if rm_tmp: remove(step, f"{f}_{step}")
        step += 1
        command(step, f"mv {f}_{step}.log {f}_{step}.hh sex_discrepancy.txt plink.sexcheck Phase_{phase}/", False)
    
    return step

def maf(step, rm_tmp, maf_threshold):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
        
    command(step, "awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' " + base +'.bim' + " > snp_1_22.txt")
    command(step, f"plink --memory {memory} --allow-no-sex --threads {threads} --bfile {base} --extract snp_1_22.txt --make-bed --out {f}_{step}")
    if rm_tmp: remove(step, f"{base}")
    command(step, f"mv {f}_{step}.log snp_1_22.txt Phase_{phase}/", False)
        
    if rscripts:
        command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --freq --out MAF_check")
        command(step, "Rscript Rscripts/MAF_check.R", False)

    if maf_threshold > 0.0:
        command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --maf {maf_threshold} --make-bed --out {f}_{step+1}")
        if rm_tmp: remove(step, f"{f}_{step}")
        step += 1
    command(step, f"mv {f}_{step}.log MAF_* Phase_{phase}/", False)
    
    return step

def hwe(step, rm_tmp):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
    
    if rscripts:
        command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --hardy")
        command(step, "awk '{ if ($9 < 0.00001) print $0 }' plink.hwe > plinkzoomhwe.hwe")
        command(step, "Rscript Rscripts/hwe.R", False)
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --hwe 1e-10 include-nonctrl --make-bed --out {f}_{step}")
    if rm_tmp: remove(step, f"{base}")
    command(step, f"mv {f}_{step}.log histhwe* plink.hwe plinkzoomhwe.hwe Phase_{phase}/", False)
    
    return step

def heterozygosity(step, rm_tmp):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --exclude range inversion.txt --indep-pairwise 50 5 0.2 --out indepSNP")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --extract indepSNP.prune.in --het --out R_check")
    if rscripts:
        logging.info(f"Step {step}: check heterozygosity")
        command(step, "Rscript Rscripts/check_heterozygosity_rate.R", False)
        command(step, "Rscript Rscripts/heterozygosity_outliers_list.R", False)
    command(step, "sed 's/\"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --remove het_fail_ind.txt --make-bed --out {f}_{step}")
    if rm_tmp: remove(step, f"{base}")
    command(step, f"mv {f}_{step}.log fail-het-qc.txt heterozygosity* het_fail_ind.txt indepSNP.log indepSNP.prune.out R_check* Phase_{phase}/", False)
    
    return step

def cryptic_relatedness(step, rm_tmp):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
        
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2")
    command(step, "awk '{ if ($8 > 0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome")
    if rscripts:
        command(step, "Rscript Rscripts/Relatedness.R", False)
        
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --filter-founders --make-bed --out {f}_{step}")
    if rm_tmp: remove(step, f"{base}")
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --missing")
    command(step, "awk '{ print $1, $2, $3, $4}' pihat_min0.2_in_founders.genome>pairs.txt")
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

    command(step, f"mv pihat* relatedness.pdf zoom_pihat.genome {f}_{step}.log plink.* pairs.txt Phase_{phase}/", False)

    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --remove 0.2_low_call_rate_pihat.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.log  0.2_low_call_rate_pihat.txt Phase_{phase}/", False)
    
    return step

def population_stratification(step, rm_tmp, inf_lim_mds_1, sup_lim_mds_1, inf_lim_mds_2, sup_lim_mds_2):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        step += 1
        
    pop_step = 1
    
    command(step, f"sed -r \'s/rs\S+/./\t\' {pop_f}.bim > tmp.bim && mv tmp.bim {pop_f}.bim")
    command(step, f"sed -r \'s/,/_/\t\' {pop_f}.bim > tmp.bim && mv tmp.bim {pop_f}.bim")
    # NOTE: the character '_' is used as allele separator because of bimbam compatibility
    command(step, f"plink --threads {threads} --bfile {pop_f} --set-missing-var-ids @:#[b37]\$1_\$2 --make-bed --out {pop_f}_{step}")


    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{step} --geno 0.2 --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.2 --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --geno 0.02 --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --mind 0.02 --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --maf 0.05 --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, "awk '{print$2}' " + base + ".bim > " + f + '_' + str(step) + "_SNPs.txt")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --extract {f}_{step}_SNPs.txt --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {f}_{step}_SNPs.txt {pop_f}_{pop_step}.log Phase_{phase}/", False)
    
    command(step, "awk '{print$2}' " + pop_f + '_' + str(pop_step) + ".bim > " + pop_f + '_' + str(pop_step) + '_SNPs.txt')
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --extract {pop_f}_{pop_step}_SNPs.txt --recode --make-bed --out {f}_{step}")
    # NOTE: don't remore this file (it's what "original" is pointing to)
    #if rm_tmp: remove(step, f"{f}_{step}")
    command(step, f"mv {f}_{step}.log {pop_f}_{pop_step}_SNPs.txt Phase_{phase}/", False)
    
    command(step, "awk '{print$2,$4}' " + f + '_' + str(step) + '.map > buildmap.txt')
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --update-map buildmap.txt --make-bed --out {pop_f}_{pop_step+1}")
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    command(step, f"mv {pop_f}_{pop_step}.log buildmap.txt Phase_{phase}/", False)
    
    command(step, "awk '{print$2,$5}' " + pop_f + '_' + str(pop_step) + '.bim > pop_ref_list.txt')
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")    
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.log Phase_{phase}/", False)
    
    command(step, "awk '{print$2,$5,$6}' " + pop_f + '_' + str(pop_step) + '.bim > ' + pop_f + '_tmp')
    command(step, "awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
    command(step, f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > all_differences.txt")
    command(step, f"mv {f}_{step}_tmp Phase_{phase}/", False)
    
    command(step, "awk '{print$1}' all_differences.txt | sort -u > flip_list.txt")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --flip flip_list.txt --reference-allele pop_ref_list.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv all_differences.txt pop_ref_list.txt flip_list.txt {f}_{step}.log Phase_{phase}/", False)
    
    command(step, "awk '{print$2,$5,$6}' " + f + '_' + str(step) + '.bim > ' + f + '_' + str(step) + '_tmp')
    command(step, f"sort {pop_f}_tmp {f}_{step}_tmp | uniq -u > uncorresponding_SNPs.txt")
    command(step, f"mv {f}_{step}_tmp Phase_{phase}/", False)
    
    command(step, "awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exclusion.txt")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --exclude SNPs_for_exclusion.txt --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {pop_f}_{pop_step} --exclude SNPs_for_exclusion.txt --make-bed --out {pop_f}_{pop_step+1}")
    command(step, f"mv {pop_f}_{pop_step}.log uncorresponding_SNPs.txt SNPs_for_exclusion.txt Phase_{phase}/", False)
    if rm_tmp: remove(step, f"{pop_f}_{pop_step}")
    pop_step += 1
    
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --bmerge {pop_f}_{pop_step}.bed {pop_f}_{pop_step}.bim {pop_f}_{pop_step}.fam --make-bed --out {f}_{step+1}")
    if rm_tmp: remove(step, f"{f}_{step}")
    step += 1
    command(step, f"mv {f}_{step}.log Phase_{phase}/", False)
    
    # NOTE: step is not incremented because it's just producing genome file
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
    command(step, f"mv {f}_{step}.log Phase_{phase}/", False)
    
    panel_step = 1
    command(step, "awk '{print$1,$1,$2}' " + panel_f + ' > race_1kG.txt')
    command(step, f" sed 's/JPT/ASN/g' race_1kG.txt > race_1kG{panel_step}.txt")
    command(step, f" sed 's/ASW/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/CEU/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/CHB/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/CHD/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/YRI/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/LWK/AFR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/TSI/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/MXL/AMR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/GBR/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/FIN/EUR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/CHS/ASN/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    command(step, f" sed 's/PUR/AMR/g' race_1kG{panel_step}.txt > race_1kG{panel_step+1}.txt")
    panel_step += 1
    
    # Racefile for original population
    command(step, "awk '{print$1,$2,\"OWN\"}' " + f + '_' + str(step) + '.fam > racefile_own.txt')

    # Concatenate racefiles
    command(step, f"cat race_1kG{panel_step}.txt racefile_own.txt | sed -e '1i\\FID IID race' > racefile.txt")

    logging.info(f"Step {step}: population analysis")
    command(step, f"Rscript Rscripts/MDS_merged.R {f}_{step}.mds pop.pdf")
    print('Population MDS displayed in "pop.pdf" file')
    if inf_lim_mds_1 == '' or sup_lim_mds_1 == '' or inf_lim_mds_2 == '' or sup_lim_mds_2 == '':
        print('Please check and specify desired cut-offs')
        while True:
            try:
                inf_lim_mds_1 = float(input(" Insert inferior threshold on MDS component 1: "))
                sup_lim_mds_1 = float(input(" Insert superior threshold on MDS component 1: "))
                inf_lim_mds_2 = float(input(" Insert inferior threshold on MDS component 2: ")) 
                sup_lim_mds_2 = float(input(" Insert superior threshold on MDS component 2: "))
            except Exception as e:
                print(e)
                print("Retry")
                continue
            break
        logging.info(f"Inferior and Superior thresholds on MDS component 1: [{inf_lim_mds_1}, {sup_lim_mds_1}]")
        logging.info(f"Inferior and Superior thresholds on MDS component 2: [{inf_lim_mds_2}, {sup_lim_mds_2}]")
    
    # Create header
    command(step, 'echo "FID       IID    SOL           C1           C2           C3           C4           C5           C6           C7           C8           C9          C10" > filtered_subj')
    command(step, "awk '{if( " + str(inf_lim_mds_1) + " <= $4 && $4 <= " + str(sup_lim_mds_1) + " && " + str(inf_lim_mds_2) + " <= $5 && $5 <= " + str(sup_lim_mds_2) + " ) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' " + f + '_' + str(step) + ".mds >> filtered_subj")  
    command(step, f"mv pop.pdf {f}_{step}.mds {f}_{step}.cluster* {f}_{step}.genome {f}_{step}.nosex Phase_{phase}/", False)  
    
    # Remove individuals from the pre-merged dataset (original_f)
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --keep filtered_subj --make-bed --out {f}_{step+1}")
    
    command(step, f"Rscript Rscripts/MDS_merged.R filtered_subj filtered.pdf")
    command(step, f"mv filtered_subj race*.txt Phase_{phase}/", False)
    if rm_tmp: 
        remove(step, f"{f}_{step}")
        remove(step, f"{base}")
    step += 1
    
    # Showing updated (original) population
    '''
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --extract indepSNP.prune.in --genome --out {f}_{step}")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {f}_{step} --read-genome {f}_{step}.genome --cluster --mds-plot 10 --out {f}_{step}")
    command(step, "awk '{print$1,$2,\"OWN\"}' " + f + '_' + str(step) + '.fam > racefile_own.txt')
    command(step, f"cat race_1kG{panel_step}.txt racefile_own.txt | sed -e '1i\\FID IID race' > racefile.txt")
    command(step, f"Rscript Rscripts/MDS_merged.R {f}_{step}.mds filtered.pdf")
    '''
    command(step, f"mv {f}_{step}.mds {f}_{step}.log {f}_{step}.cluster* filtered.pdf {f}_{step}.mds race*txt Phase_{phase}/", False)
    print('# Filtered population displayed in "filtered.pdf" file in Phase_{phase}/ folder')

    return step

def mds_covariates(step):
    if step == 0:
        base = f
    else:
        base = f"{f}_{step}"
        
    # NOTE: producing genome the step counter is not incremented 
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --extract indepSNP.prune.in --genome --out {base}")
    command(step, f"plink --memory {memory} --allow-no-sex  --threads {threads} --bfile {base} --read-genome {base}.genome --cluster --mds-plot 10 --out {base}")
    command(step, "awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' " + base + '.mds > covar_mds.txt')
    command(step, f"mv {base}.cluster* covar_mds.txt {base}.genome {base}.log {base}.mds Phase_{phase}/", False)
    
    return step

if __name__ == '__main__':
    start_time = datetime.today()
    logging.basicConfig(format='%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s', datefmt='%Y-%m-%d,%H:%M:%S', level=logging.INFO, filename = f"trace_{start_time.strftime('%Y-%m-%d-%H:%M:%S')}.log")
    f, threads, memory, pop_f, panel_f, discard_sex, maf_threshold, rscripts, rm_tmp, inf_lim_mds_1, sup_lim_mds_1, inf_lim_mds_2, sup_lim_mds_2, config = collect_inputs()

    # NOTE: every step will check if step is equal to 0 to know if it exists a file "f_step" or if "f" must be taken as first input
    step = 0

    # PHASE 1 (Missingness)    
    if config['Missingness'] == True:
        phase = 1
        mkdir(phase)
        try:
            step = missingness(step, rm_tmp)
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
        print_info(phase, f"{f}_{step}")

    # PHASE 2 (Sex discrepancy)
    if config['Sex_discrepancy'] == True:
        phase = 2
        mkdir(phase)
        try:
            step = sex_discrepancy(step, rm_tmp, discard_sex)
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
        print_info(phase, f"{f}_{step}")

    # PHASE 3 (MAF filtering)
    if config['Autos_MAF'] == True:
        phase = 3
        mkdir(phase)
        try:
            step = maf(step, rm_tmp, maf_threshold)
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
        print_info(phase, f"{f}_{step}")

    # PHASE 4 (Hardy-Weinberg Equilibrium)
    if config['Hardy-Weinberg'] == True:
        phase = 4
        mkdir(phase)
        try:
            step = hwe(step, rm_tmp)
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
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
        try:
            step = cryptic_relatedness(step, rm_tmp)
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
        print_info(phase, f"{f}_{step}")

    if len(pop_f)>0:
        # PHASE 7 (Population filtering)
        if config['Population_filtering'] == True:
            phase = 7
            mkdir(phase)
            try:
                step = population_stratification(step, rm_tmp, inf_lim_mds_1, sup_lim_mds_1, inf_lim_mds_2, sup_lim_mds_2)
            except Exception as e:
                logging.error(f"Error during phase {phase}. Traceback:")
                logging.error(e)
            print_info(phase, f"{f}_{step}")
            
    # PHASE 8 (MDS Covariates)
    if config['MDS_covariates'] == True:
        phase = 8
        mkdir(phase)
        try:
            step = mds_covariates(step)  
        except Exception as e:
            logging.error(f"Error during phase {phase}. Traceback:")
            logging.error(e)
        print_info(phase, f"{f}_{step}")
        
    command(step, f"mv {f}_{step}.fam QC_{f}.fam")
    command(step, f"mv {f}_{step}.bim QC_{f}.bim")
    command(step, f"mv {f}_{step}.bed QC_{f}.bed")
    print(f'### QC completed, the final files are QC_{f}')
    command(step, f"rm {f}_*.sex", False)
    command(step, f"rm {f}_*.sexcheck", False)
    command(step, f"rm {f}*log", False)
    
    if rm_tmp and len(pop_f)>0:
        # created during merging
        command(step, f"rm {f}_*.map", False)
        command(step, f"rm {f}_*.ped", False)
        command(step, f"rm {pop_f}_*.*", False)  
        command(step, f"rm {pop_f}_tmp", False)           
    
    print("### DONE! You can proceed to imputation with 'imp_pheno.py'")

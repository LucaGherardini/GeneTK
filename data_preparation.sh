python3 conv_vcf2plink.py 
echo "### Conversion done"
python3 extract_SNPs.py p_
echo "### SNPs extraction done"
python3 recode.py snp_
echo "### ID recoding done"
python3 merge_plinks.py recoded_
echo "### Merging done"
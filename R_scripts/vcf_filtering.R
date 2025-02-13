#VCF filtering

##NOtes:
#- more functions in dartR, better if vcf was created using freebayes

gl2vcf(data_gl_filtered, outfile = "final_filtered", plink.bin.path = 'C:/Users/gerar/Desktop/plink_win64_20231018')
#50405 variants and 217 people pass filters and QC.


#Run below commands in Ubuntu

#quality command removed as not present
# vcftools --vcf final_filtered.vcf --max-missing 0.5 --mac 3 --recode --recode-INFO-all --out raw.g5mac3.recode.vcf
# 
# vcftools --vcf raw.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out raw.g5mac3dp3 
# 
# vcftools --vcf raw.g5mac3dp3.recode.vcf --missing-indv
# 
# cat out.imiss
# 
# mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
# 
# vcftools --vcf raw.g5mac3dp3.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out raw.g5mac3dplm
# 
# vcftools --vcf raw.g5mac3dplm.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out DP3g95maf05 --min-meanDP 20


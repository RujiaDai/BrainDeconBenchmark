#prepare for genotype file
for i in {1..22}
do
#selcet sample
bcftools view -S ../RNA_WGS2.txt --force-samples ROSMAP/WGS/"DEJ_11898_B01_GRM_WGS_2017-05-15_"$i.recalibrated_variants.vcf.gz -Oz > rna.chr$i.vcf.gz
gunzip rna.chr$i.vcf.gz
vcftools --vcf rna.chr$i.vcf --plink --out rna.chr$i
rm rna.chr$i.vcf
plink --file rna.chr$i --make-bed --out rna.chr$i
plink --bfile rna.chr$i --geno 0.05  --maf 0.05 --hwe 1e-6 --recode vcf-iid --out qc.chr$i
rm rna.chr$i.*
bgzip qc.chr$i.vcf
tabix -p vcf qc.chr$i.vcf.gz
#annotate rsid for variants
java -jar snpEff/SnpSift.jar annotate /zs32/home/tychu/tools/snpEff/data/All_20170710.vcf qc.chr$i.vcf.gz > geno/snp.chr$i.vcf
bgzip geno/snp.chr$i.vcf
rm qc.chr$i.*
done


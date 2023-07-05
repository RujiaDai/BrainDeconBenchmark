#eqtl
#nominal
for i in Ast End Mic Ex In Oli Opc Per
do
echo $i
mkdir -p nominal/$i
parallel -j 6 QTLtools cis --vcf geno/snp.{}.vcf.gz \
--bed pheno/$i/$i.chr{}.bed.gz \
--cov cov/$i.xls \
--out nominal/$i/$i.chr{} --nominal 1 --std-err  ::: {1..22}
done

#permutation
for i in Ast End Mic Ex In Oli Opc Per
do
echo $i
mkdir -p permutation/$i
parallel -j 8 QTLtools cis --vcf geno/snp.{}.vcf.gz \
--bed pheno/$i/$i.chr{}.bed.gz \
--cov cov/$i.xls \
--out permutation/$i/$i.chr{} --permute 1000  ::: {1..22} 
done

for i in Ast End Mic Ex In Oli Opc Per
do
for j  in {1..22}
do
cat permutation/$i/$i.chr$j >> $i.all
done
Rscript qtltools/scripts/qtltools_runFDR_cis.R $i.all 0.05 $i
done

#conditional
for i in Ast End Mic Ex In Oli Opc Per
do
echo $i
parallel -j 12 QTLtools cis --vcf geno/snp.{}.vcf.gz \
--bed pheno/$i/$i.chr{}.bed.gz \
--cov cov/$i.xls \
--mapping permutation/$i.thresholds.txt --normal --std-err --out condi/$i/$i.chr{} ::: {1..22} 
cat condi/$i/$i.chr* > condi/$i/$i.all
rm condi/$i/$i.chr*
done
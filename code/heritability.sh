#sldsc
#prepare annotation file
for i in Ast End Ex In Mic Oli Opc Per union sc_Ast sc_End sc_Ex sc_In sc_Mic sc_Oli sc_Opc sc_Per sc_union Bulk
do
ldfile=/zs32/data-analysis/liucy_group/jiangyi/database/ld.haploreg/LD_EUR.gt0.2.gz
#snp 的chr,pos path:/vg_sklmgdata_hw_01/data/tychu/eqtl/sldsc_eqtl/new
eqtl=/vg_sklmgdata_hw_01/data/tychu/eqtl/sldsc_eqtl/new/$i".eqtl"
outdir=/vg_sklmgdata_hw_01/data/tychu/eqtl/sldsc_multi/$i
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1"_"$2]++}NR>FNR&&($1"_"$2 in a){split($4,s,":");print $1,s[1],s[2]}' $eqtl <(zcat $ldfile) |awk '!a[$0]++' > $outdir/a.bed
cat $outdir/a.bed | perl -F"\t" -lane '$start=($F[2-1]-0)>0?($F[2-1]-0):0;print $F[1-1],"\t",$start,"\t",$F[3-1]+0' | sort -k1,1V -k2,2n | grep -P "^\d" > $outdir/a.bed2
bedtools intersect -a /zs32/data-analysis/liucy_group/jiangyi/database/1000G/baseline.annot.snps -b $outdir/a.bed2 -wa -c | perl -lane 'if($F[5]>1){$F[5]=1};print join "\t",@F' > $outdir/a.annot
awk '{print $0 > "'$outdir'/a.annot."$1}' $outdir/a.annot
parallel -j 11 cat $outdir/a.annot.{} \| cut -f 1,2,4,5,6 \| sed \''1i\CHR\tBP\tSNP\tCM\tinput'\' \| gzip \> $outdir/a.annot.{}.annot.gz ::: $(seq 1 22) 
parallel -j 10 ~/tools/ldsc-master/ldsc.py --l2 --bfile /zs32/data-analysis/liucy_group/jiangyi/database/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.{} \
--ld-wind-cm 1 --annot $outdir/a.annot.{}.annot.gz --out $outdir/a.annot.{} \
--print-snp  /zs32/home/yjiang/softwares/ldsc/hapmap3_snps/hm.{}.snp ::: $(seq 1 22)
rm $outdir/a.annot
done


#joint model for multiple annotation
for g in 25056061.2014Nature.sczVScontrol 29483656.2018NG.sczVScontrol 29906448.2018Cell.sczVScontrol; 
do  
nohup ~/tools/ldsc-master/ldsc.py --h2 /zs32/data-analysis/liucy_group/jiangyi/psychENCODE/heritability/gwas.sumstats/$g.sumstats.gz \
--ref-ld-chr /zs32/data-analysis/liucy_group/jiangyi/database/1000G/baseline_v1.1/baseline.,Ast/a.annot.,End/a.annot.,Ex/a.annot.,\
In/a.annot.,Mic/a.annot.,Oli/a.annot.,Opc/a.annot.,Per/a.annot.,union/a.annot.,sc_Ast/a.annot.,sc_End/a.annot.,\
sc_Ex/a.annot.,sc_In/a.annot.,sc_Mic/a.annot.,sc_Oli/a.annot.,sc_Opc/a.annot.,sc_Per/a.annot.,sc_union/a.annot.,Bulk/a.annot. \
--frqfile-chr /zs32/data-analysis/liucy_group/jiangyi/database/1000G/1000G_Phase3_frq/1000G.EUR.QC.  \
--w-ld-chr /zs32/data-analysis/liucy_group/jiangyi/database/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-coefficient \
--out /vg_sklmgdata_hw_01/data/tychu/eqtl/sldsc_multi/a.h2.$g &
done


awk 'FNR==1{print "GWAS""\t"$0}' a.h2.25056061.2014Nature.sczVScontrol.results > h2enrichment.txt
grep ^L2 a.h2.*.results |sed  -e's/\/a.h2./\t/' -e's/.results//' | tr ":" "\t" >> h2enrichment.txt
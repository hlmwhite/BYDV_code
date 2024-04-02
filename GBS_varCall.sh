#!/bin/bash

REF_GENOME=$(readlink -f $1 )
ORG_OUT=$2
SAMPLES_LIST=$(readlink -f $3)
PLOIDY=$4
THIN_VALUE=$5



mkdir "$ORG_OUT"_vcf
sleep 1s

echo 'running mpileup...'
sleep 1s

if [ -f "$ORG_OUT"_vcf/"$ORG_OUT"_mpileup.vcf ]; then
        echo 'found mpileup, skipping....'
else

bcftools mpileup --output-type z --skip-indels --annotate AD,DP --fasta-ref $REF_GENOME --min-MQ 20 --min-BQ 20  --no-version -b $SAMPLES_LIST -o "$ORG_OUT"_vcf/"$ORG_OUT"_mpileup.vcf

fi


echo 'running call...'
sleep 1s
if [ $PLOIDY == "1" ]; then
        bcftools call --multiallelic-caller --variants-only --ploidy $PLOIDY --no-version "$ORG_OUT"_vcf/"$ORG_OUT"_mpileup.vc* > "$ORG_OUT"_vcf/"$ORG_OUT"_raw.vcf 
elif [ $PLOIDY == "2" ]; then
        bcftools call --multiallelic-caller --variants-only --no-version "$ORG_OUT"_vcf/"$ORG_OUT"_mpileup.vc* > "$ORG_OUT"_vcf/"$ORG_OUT"_raw.vcf
else
        echo 'ploidys greater than 2 are prohibited, exiting...'
        exit
fi

echo 'running filter...'
sleep 1s
cat "$ORG_OUT"_vcf/"$ORG_OUT"_raw.vcf | java -jar SnpSift.jar filter "( DP > 2)" > "$ORG_OUT"_vcf/mod_"$ORG_OUT"_minDepth_2.vcf

#vcftools --vcf mod_BuchRp_minDepth_2.vcf --maf 0.05 --out maf0.05_PASSED_SNP.Anno.vcf --recode --recode-INFO-all --stdout > maf_0.05_PASSED_SNP.Anno.vcf
vcftools --vcf "$ORG_OUT"_vcf/mod_"$ORG_OUT"_minDepth_2.vcf --maf 0.05 --recode --recode-INFO-all --stdout > "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.vcf

echo 'gathering shared variant loci across all samples...'
get_genotypes.sh "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.vcf $PLOIDY > "$ORG_OUT"_vcf/genotypes.tbl

cat "$ORG_OUT"_vcf/genotypes.tbl | cut -f 2- | nl -b a | tail -n +2 > "$ORG_OUT"_vcf/genotypes.tbl_ln

awk 'NR == FNR{a[$0]; next};FNR in a' <(grep -v '\.' "$ORG_OUT"_vcf/genotypes.tbl_ln | cut -f 1 | sed 's/\s//g') "$ORG_OUT"_vcf/genotypes.tbl > "$ORG_OUT"_vcf/shared_genotypes.list

bcftools view -I "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.vcf -O z -o "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.zipped.vcf.bgz
bcftools index "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.zipped.vcf.bgz

bcftools view -R <(cat "$ORG_OUT"_vcf/shared_genotypes.list | cut -f 1,2) "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.zipped.vcf.bgz > "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.vcf

echo 'thinning vcf...'
sleep 1s
vcftools --vcf "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.vcf --stdout --thin $THIN_VALUE --recode > "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.thin5000.vcf

num_samples=$(cat $SAMPLES_LIST | wc -l)

echo 'generating PCA...'
plink --vcf "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.thin5000.vcf --double-id --allow-extra-chr --make-bed --pca $num_samples --out "$ORG_OUT"_vcf/plink

Rscript PCA.Rscript "$ORG_OUT"_vcf

echo 'generating MDS...'
plink --vcf "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.thin5000.vcf --double-id --allow-extra-chr --make-bed --cluster --mds-plot $num_samples --out "$ORG_OUT"_vcf/plink_mds

cat "$ORG_OUT"_vcf/plink_mds.mds | tr '/' '\t' | awk '{print $1}' | sed 's/FID/sample_label/' | paste - "$ORG_OUT"_vcf/plink_mds.mds > "$ORG_OUT"_vcf/mod.plink.mds

Rscript MDS.Rscript "$ORG_OUT"_vcf


echo 'generating phylogenetic tree...'
vk phylo tree nj "$ORG_OUT"_vcf/"$ORG_OUT"_minDepth_2.shared.thin5000.vcf | tr '\n' ' ' | sed -e 's/ //'g -e 's/$/\n/' > "$ORG_OUT"_vcf/"$ORG_OUT".nwk

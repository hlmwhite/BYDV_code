#!/bin/bash

####
# Author: Mark Whitehead, hlmwhite@liverpool.ac.uk
####


if [ "$#" == 0 ]; then

        echo ''
        echo 'usage: ./get_genotypes.sh <vcf> > <outfile> <ploidy>'
        echo ''
        echo '0/0 - the sample is homozygous reference'
echo '0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles'
echo '1/1 - the sample is homozygous alternate'

echo ''

echo 'more info at: https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it'

echo ''
exit

fi

PLOIDY=$2

SUBSTR_NUM=$(expr $PLOIDY \* 2 + 1)
cat $1 | grep -v '#' | cut -f 1,2,4,5,10- | awk -v ploid=$SUBSTR_NUM '{for(i=5;i<=NF;i++) $i=substr($i,1,ploid)}1' | sed 's/:.//g' | cat <(grep '#CHROM' $1 | cut -f 1,2,4,5,10- ) - | tr ' ' '\t' 

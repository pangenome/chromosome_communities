#!/bin/bash

# inputs described here https://github.com/human-pangenomics/hpp_pangenome_resources#pggb

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/vcfs/hprc-v1.0-pggb.chm13.1-22%2BX.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/vcfs/hprc-v1.0-pggb.grch38.1-22%2BX.vcf.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/vcfs/hprc-v1.0-pggb.grch38.Y.vcf.gz

zcat pggb.wgg.88.chm13.1-22+X.vcf.gz | cut -f -8 | awk 'length($4) == 1 && length($5) == 1 || /^#/' | sed 's/^chm13#//' | vcf2tsv | cut -f 1,2,8,13 | awk 'NR == 1 { print $0, "callset" } NR > 1 { print $0, "pggb.wgg.88.chm13" }' | tr ' ' '\t' | gzip >pggb.wgg.88.chm13.1-22+X.biallelic_snps.tsv.gz

zcat pggb.wgg.88.grch38.1-22+X.vcf.gz | cut -f -8 | awk 'length($4) == 1 && length($5) == 1 || /^#/' | sed 's/^grch38#//' | vcf2tsv | cut -f 1,2,8,13 | awk 'NR == 1 { print $0, "callset" } NR > 1 { print $0, "pggb.wgg.88.grch38" }' | tr ' ' '\t' | gzip >pggb.wgg.88.grch38.1-22+X.biallelic_snps.tsv.gz

zcat pggb.wgg.88.grch38.Y.vcf.gz | cut -f -8 | awk 'length($4) == 1 && length($5) == 1 || /^#/' | sed 's/^chm13#//' | vcf2tsv | cut -f 1,2,8,13 | awk 'NR == 1 { print $0, "callset" } NR > 1 { print $0, "pggb.wgg.88.chm13" }' | tr ' ' '\t' | gzip >pggb.wgg.88.grch38.Y.biallelic_snps.tsv.gz

( echo chrom pos ac lv graph; zcat pggb.wgg.88.chm13.1-22+X.biallelic_snps.tsv.gz | tail -n+2 ; zcat pggb.wgg.88.grch38.1-22+X.biallelic_snps.tsv.gz | tail -n+2 ; zcat pggb.wgg.88.grch38.Y.biallelic_snps.tsv.gz
     | tail -n+2 ) | tr ' ' '\t' | pigz >bisnps.tsv.gz

Rscript plot_bisnp_dens.R

pdfunite $(ls *.lv_by.pdf | sort -V | grep -v M) pggb.wgg.88.bisnips.chm13_vs_grch38.pdf

### Fst

```shell
cd ~/tools
git clone --recursive https://github.com/ksamuk/pixy.git

#guix install python-numcodecs
#guix install python-multiprocess
#guix install python-scipy
```

```shell
THREADS=48

#run_pixy="python3 ~/tools/pixy/pixy/__main__.py"
run_pixy="python3 /home/guarracino/git/pixy/pixy/__main__.py"
```

Choose the VCF file:

```shell
PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n84/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.vcf.gz
#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n84/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.grch38.vcf.gz
#PATH_REF_FASTA=xxxxxx

#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.vcf.gz
#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.grch38.vcf.gz

#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.chm13.vcf.gz
#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.grch38.vcf.gz
#PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.hg002.vcf.gz


OUTPUT_DIR=/lizardfs/guarracino/chromosome_communities/fst
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

VCF_NAME=$(basename $PATH_VCF_GZ .vcf.gz)
PATH_SNV_VCF_GZ=$OUTPUT_DIR/$VCF_NAME.snps.vcf.gz

# Get reference names
zgrep '^#' $PATH_VCF_GZ -v | cut -f 1 | uniq | sort | uniq > ref_names.txt
```

Take SNPs:

```shell

# Take SNPs and normalize
vcftools --gzvcf $PATH_VCF_GZ --remove-indels --recode --recode-INFO-all --stdout |\
  #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
  bgzip -@ $THREADS > $PATH_SNV_VCF_GZ
tabix $PATH_SNV_VCF_GZ
rm out.log
```

```shell
# Prepare populations
touch pop.txt
(seq 13 15; seq 21 22) | while read f; do
  PATH_FASTA_FAI=/lizardfs/erikg/HPRC/year1v2genbank/parts/chr$f.pan.fa.fai
  
  # Intersect all samples with samples really present in the VCF (removing the references)
  comm -12 \
    <(cut -f 1 $PATH_FASTA_FAI | sort) \
    <(zgrep '^#CHROM' $PATH_SNV_VCF_GZ -m 1 | cut -f 10- | tr '\t' '\n' | grep chr -v | sort) |\
     awk -v OFS="\t" -v chr=chr$f '{print ($0,chr)}' >> pop.txt
done

#FST estimator to use, one of either: 
#- 'wc' (Weir and Cockerham 1984) or
#- 'hudson' (Hudson 1992, Bhatia et al. 2013) 
for FST_TYPE in wc hudson; do
  for WIN_SIZE in 2000 5000 10000 20000; do
    $run_pixy --n_cores $THREADS \
      --vcf $PATH_SNV_VCF_GZ \
      --stats 'fst' --fst_type $FST_TYPE \
      --pop pop.txt \
      --output_prefix $FST_TYPE.$WIN_SIZE \
      --bypass_invariant_check 'yes' \
      --window_size $WIN_SIZE
  done
done
```


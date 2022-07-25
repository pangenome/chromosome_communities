# Linkage Disequilibrium

#### 0.Download files

```
wget http://hypervolu.me/~erik/chrcommunity/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.no_HG002-hifi.snv.vcf.gz

wget http://hypervolu.me/~guarracino/chrcommunity/chrARCO_25-Jul-22_PPRRs.bed
```

#### 1. Remove CONFLICT regions from vcf

```
zcat chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.no_HG002-hifi.snv.vcf.gz | grep -v "CONFLICT" > chrACRO+refs.vcf

bgzip chrACRO+refs.vcf

tabix -p vcf chrACRO+refs.vcf.gz
```

#### 2. Compute linkage disequilibrium

p arm
```
cat data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 > chrNames.txt
cat data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 | sed 's/chm13#//g' > acrocentrics.txt
mkdir ld

paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do plink --vcf chrACRO+refs.vcf.gz --allow-extra-chr  --chr-set -5 --chr $f1 --r2  --extract 'range' data/annotation/chm13.p_arms.approximate.acros.bed --ld-window-kb 10  --ld-window-r2 0 --out ld/$f2.pArm
```
q arm
```
paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do plink --vcf chrACRO+refs.vcf.gz --allow-extra-chr  --chr-set -5 --chr $f1 --r2  --extract 'range' data/annotation/chm13.q_arms.approximate.acros.bed --ld-window-kb 10  --ld-window-r2 0 --out ld/$f2.qArm
```

putative recombinant
```
paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do plink --vcf chrACRO+refs.vcf.gz --allow-extra-chr  --chr-set -5 --chr $f1 --r2  --extract 'range' chrARCO_25-Jul-22_PPRRs.bed --ld-window-kb 10  --ld-window-r2 0 --out ld/$f2.recombinant
```

#### 3. Plot

```
Rscript plot_ld.R 'ld/*' acrocentrics_ld.pdf

```

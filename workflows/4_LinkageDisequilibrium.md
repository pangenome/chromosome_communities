# Linkage Disequilibrium

## Data preparation

Download the files:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/linkage_disequilibrium
cd /lizardfs/guarracino/chromosome_communities/linkage_disequilibrium

wget http://hypervolu.me/~erik/chrcommunity/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.no_HG002-hifi.snv.vcf.gz
wget http://hypervolu.me/~guarracino/chrcommunity/chrARCO_25-Jul-22_PPRRs.bed
```

Remove CONFLICT regions from the VCF file:

```shell
zgrep -v "CONFLICT" chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.no_HG002-hifi.snv.vcf.gz | \
    bgzip -c -@ 48 > chrACRO+refs.vcf.gz
tabix -p vcf chrACRO+refs.vcf.gz
```

## Compute linkage disequilibrium (LD)

Prepare files:

```shell
cat ../data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 > chrNames.txt
cat ../data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 | sed 's/chm13#//g' > acrocentrics.txt
sed -i 's/^/chm13#/g' chrARCO_25-Jul-22_PPRRs.bed
mkdir ld
```

Compute LD:

```shell
# p arm
paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do
  plink \
    --vcf chrACRO+refs.vcf.gz \
    --allow-extra-chr \
    --chr-set -5 \
    --chr $f1 \
    --r2 \
    --extract 'range' ../data/annotation/chm13.p_arms.approximate.acros.bed \
    --ld-window-kb 10 --ld-window-r2 0 \
    --out ld/$f2.pArm;
done

# q arm
paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do
  plink \
    --vcf chrACRO+refs.vcf.gz \
    --allow-extra-chr \
    --chr-set -5 \
    --chr $f1 \
    --r2 \
    --extract 'range' ../data/annotation/chm13.q_arms.approximate.acros.bed \
    --ld-window-kb 10 --ld-window-r2 0\
    --out ld/$f2.qArm;
done

# putative recombinant
paste -d'\n' chrNames.txt acrocentrics.txt | while read f1 && read f2 ; do\
  plink \
    --vcf chrACRO+refs.vcf.gz \
    --allow-extra-chr \
    --chr-set -5 \
    --chr $f1 \
    --r2 \
    --extract 'range' chrARCO_25-Jul-22_PPRRs.bed \
    --ld-window-kb 10 --ld-window-r2 0 \
    --out ld/$f2.recombinant;
done
```

#### Visualization

```
Rscript ../scripts/plot_ld.R 'ld/*.ld' acrocentrics_ld.pdf
```

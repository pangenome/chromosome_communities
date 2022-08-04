# Linkage Disequilibrium

## Data preparation

Download the files:

```shell
mkdir -p ../linkage_disequilibrium
cd ../linkage_disequilibrium

wget http://hypervolu.me/~guarracino/chrcommunity/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.snv.norm.no_HG002-hifi.vcf.gz
wget http://hypervolu.me/~guarracino/chrcommunity/chrACRO_29-Jul-22_PHRs.bed
```

Remove CONFLICT regions from the VCF file:

```shell
zgrep -v "CONFLICT" chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.snv.norm.no_HG002-hifi.vcf.gz | \
    bgzip -c -@ 48 > chrACRO+refs.vcf.gz
tabix -p vcf chrACRO+refs.vcf.gz
```
## Data exploration 

```shell
mkdir stats 
```

#### Determine fraction onf missing sites per locus 

```shell
plink \
--vcf chrACRO+refs.vcf.gz \
--missing  \
--chr-set -5 \
--allow-extra-chr \
--extract 'range' ../data/annotation/chm13.p_arms.approximate.acros.bed \
--out stats/chrACRO+refs.pArm \

plink \
--vcf chrACRO+refs.vcf.gz \
--missing  \
--chr-set -5 \
--allow-extra-chr \
--extract 'range' ../data/annotation/chm13.q_arms.approximate.acros.bed \
--out stats/chrACRO+refs.qpArm \

plink \
--vcf chrACRO+refs.vcf.gz \
    --chr-set -5 \
--missing  \
--allow-extra-chr \
--extract 'range' chrACRO_29-Jul-22_PHRs.bed \
--out stats/chrACRO+refs.recombinant \
```

#### Calculate allele frequencies 
```shell
plink \
--vcf chrACRO+refs.vcf.gz \
--freq  \
--chr-set -5 \
--allow-extra-chr \
--out stats/chrACRO+refs
```

#### Make a file of chr positions 
```shell
zcat chrACRO+refs.vcf.gz | grep -v '##' | cut -f1,2 > stats/chrpos.txt
```

#### visualization - supplementary figure 

Rscript ../scripts/plot_stats.R '../stats/*.lmiss' '../stats/chrACRO+refs.frq' '../chrACRO_29-Jul-22_PHRs.bed' '../stats/chrpos.txt'  statsSupp.pdf


## Compute linkage disequilibrium (LD)

Prepare files:

```shell
cat ../data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 > chrNames.txt
cat ../data/annotation/chm13.p_arms.approximate.acros.bed | cut -f 1 | sed 's/chm13#//g' > acrocentrics.txt
sed -i 's/^/chm13#/g' chrACRO_29-Jul-22_PHRs.bed
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
    --extract 'range' chrACRO_29-Jul-22_PHRs.bed \
    --ld-window-kb 10 --ld-window-r2 0 \
    --out ld/$f2.recombinant;
done
```

#### Visualization

```shell
Rscript ../scripts/plot_ld.R 'ld/*.ld' acrocentrics_ld.pdf
```

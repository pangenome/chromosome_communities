# Sex graphs

## Tools

```shell
mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
git pull
git checkout d5a8de4de4d5bd16f683f19c5708f845f51b2f9a
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-d5a8de4de4d5bd16f683f19c5708f845f51b2f9a
cd ..
```


Collect assemblies from CHM13 and HG002 (GRCh38's chrY is particularly incomplete) and put everything together.

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

# Get HG002 chromosome X and Y
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
gunzip HG002.mat.cur.20211005.fasta.gz
samtools faidx HG002.mat.cur.20211005.fasta
samtools faidx HG002.mat.cur.20211005.fasta SX | sed 's/SX/HG002#MAT#chrX/g' |\
  bgzip -@ 48 -c > hg002.chrX.fa.gz

~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002#PAT#' /lizardfs/erikg/T2T/liftover/split_chrY/chm13v2.0_chrY.fasta |\
  bgzip -@ 48 -c > hg002.chrY.fa.gz
  
# Get CHM13 chromosome X too
samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa chm13#chrX | bgzip -@ 48 -c > chm13.chrX.fa.gz

# Put all together
zcat chm13.chrX.fa.gz hg002.chrX.fa.gz hg002.chrY.fa.gz | bgzip -@ 48 -c > CHM13+HG002.chrXY.fa.gz
samtools faidx CHM13+HG002.chrXY.fa.gz
```

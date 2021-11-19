# Traces of recombination

Create the main folder.

```
mkdir /lizardfs/chromosome_communities/
cd /lizardfs/chromosome_communities/
```

### Build graphs

Apply `pggb` on the chromosome-partitioned HPRC dataset. We make two graphs considering the following chromosome
jointly:

- acrocentric chromosomes (chr13, chr14, chr15, chr21, and chr22);
- sex chromosomes (chrX and chrY).

```
( echo A S | tr ' ' '\n') | while read i; do sbatch -p workers -c 48 --wrap 'hostname; cd /scratch && pggb -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -o chr'$i'.pan -t 48 -p 98 -s 100000 -n 90 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# -Z ; mv /scratch/chr'$i'.pan /lizardfs/chromosome_communities/'; done >>pggb.jobids
```

### Untangling
ToDo

### ...
ToDo

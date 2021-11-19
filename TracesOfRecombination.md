# Traces of recombination

```
(echo pggb wfmash seqwish smoothxg odgi gfaffix | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '

/gnu/store/2rwrch6gc5r5pazikhfc25j2am2rh22a-pggb-0.2.0+531f85f-1/bin/pggb
/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
/gnu/store/1v50lsaj5wjblznm3r8f83ccsw084wkr-seqwish-0.7.1+2ab95d7-6/bin/seqwish
/gnu/store/hyww54rrsga53xfr9gaixk08jq387zhj-smoothxg-0.6.0+d39280b-9/bin/smoothxg
/gnu/store/aw7x195pl3v83jk2ncni8xv1h82gxmwx-odgi-0.6.2+b04c7a9-11/bin/odgi
/gnu/store/31b6l8wxb6qgm4i8447233jv0z8bra5h-gfaffix-0.1.2.2/bin/gfaffix
```

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

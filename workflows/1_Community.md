# Community detection

We evaluate the chromosome separation among the HPRCy1v2 assemblies.

We first apply all-to-all approximate mapping:

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9
PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank+refs.fa.gz

mkdir -p /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

# Mappings without references included
for s in 50k 20k 10k; do
  for p in 98 95 90; do
    # There are 94 haplotypes, so `-n` is equal to `94 - 1`
    sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n 93 -m > HPRCy1v2genbank.self.s'$s'.p'$p'.paf && mv HPRCy1v2genbank.self.s'$s'.p'$p'.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
  done
done

# Mappings with references included
for s in 50k 20k 10k; do
  for p in 98 95 90; do
    # There are 96 haplotypes, so `-n` is equal to `96 - 1`
    sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n 95 -m > HPRCy1v2genbank+refs.self.s'$s'.p'$p'.paf && mv HPRCy1v2genbank+refs.self.s'$s'.p'$p'.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
  done
done

# More stringent mapping for visualization
for s in 50k; do
  for p in 95 90; do
    sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n 3 -m > HPRCy1v2genbank+refs.self.s'$s'.p'$p'.n3.paf && mv HPRCy1v2genbank+refs.self.s'$s'.p'$p'.n3.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
  done
done
```

Then, we collect mappings between sequences > `l` kbp long:

```shell
cd /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

for s in 50k ; do
  for p in 98 95; do
    for l in 0 500000 1000000; do
      for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
        echo $s $p $l $prefix
        awk -v len=$l '$2 > len && $7 > len' $prefix.self.s$s.p$p.paf > $prefix.self.s$s.p$p.l$l.paf
      done
    done
  done
done
```

To evaluate chromosome communities, we build an "alignment graph" from our mappings using `paf2net` scripts delivered in
`pggb` repository. In this graph, nodes are contigs and edges are mappings between them.

```shell
for s in 50k ; do
  for p in 98 95; do
    for l in 0 500000 1000000; do
      for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
        echo $s $p $l $prefix
        python3 /home/guarracino/tools/pggb/scripts/paf2net.py -p $prefix.self.s$s.p$p.l$l.paf
      done
    done
  done
done
```

Then we obtain the ``Leiden`` communities:

```shell
# Prepare chromosome info
( seq 1 22; echo X; echo Y; echo M ) | while read i; do
  cut -f 1 /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$i.pan.fa.fai | grep chr -v | grep MT -v | awk -v OFS='\t' -v chr=chr$i '{print($0,chr)}' >> contig2chr.tsv
done

for s in 50k ; do
  for p in 98 95; do
    for l in 0 500000 1000000; do
      for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
        echo $s $p $l $prefix
        
        PAF=$prefix.self.s$s.p$p.l$l.paf
        
        ID2NAME=$PAF.vertices.id2name.txt
        if [ $prefix == "HPRCy1v2genbank" ]; then
          # Prepare labels by contig
          join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
            <(sort -k 2,2 $PAF.vertices.id2name.txt) \
            <(sort -k 1,1 /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_with_pq/*.partitioning_with_pq.tsv) | cut -f 1,3 -d ' ' > $PAF.vertices.id2name.chr.txt
#            <(sort -k 1,1 contig2chr.tsv) | cut -f 1,3 -d ' ' > $PAF.vertices.id2name.chr.txt
            
          ID2NAME=$PAF.vertices.id2name.chr.txt
        fi
        
        # guix install python-igraph
        python3 /home/guarracino/tools/pggb/scripts/net2communities.py \
          -e $PAF.edges.list.txt \
          -w $PAF.edges.weights.txt \
          -n $ID2NAME
      done
    done
  done
done

ls HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done
ls HPRCy1v2genbank+refs.self.s50k.p95.l0.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; grep chr $f; done
```

Finally, we look for groupings of reference chromosomes in the same community.

```shell
join <(tr ' ' '\t' <$o.louvain | sort -k 1b,1 ) \
     <(tr ' ' '\t' <$o.names | sort -k 1b,1) \
    | grep 'chm13\|grch38' | awk '{ print $2, $3 }' | grep -v _ \
    | sort -n | awk '{ c=$1 ; if (c == p) { x=x" "$2 } else { print x; x=$2 }; p = c; }'
```


Check how chromosomes where partitioned:

```shell


for s in 50k ; do
  for p in  95; do
    for l in 1000000; do
      echo $s $p $l
      
      PAF=HPRCy1v2genbank.self.s$s.p$p.l$l.paf
      
      ls $PAF.edges.weights.txt.community.*.txt | while read f; do
        echo $f
        cat $f | sort | uniq -c
#        join -a 1 -e 'unmapped' -o '1.1 2.2' \
#          <(sort -k 1,1 $f) \
#          <(sort -k 1,1 contig2chr.tsv) | cut -f 2 -d ' '| sort | uniq -c
      done
    done
  done
done
```

```shell
( echo "Id,Label,Chromosome"; \
  join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
    <(sort -k 2,2 HPRCy1v2genbank.self.s50k.p95.l1000000.paf.vertices.id2name.txt) \
    <(sort -k 1,1 contig2chr.tsv) | tr ' ' ',' \
) > HPRCy1v2genbank.self.s50k.p95.l1000000.nodes.csv
  
( echo "Source,Target,Weight"; \
  paste -d ' ' HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.list.txt HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.weights.txt | tr ' ' ','
) > HPRCy1v2genbank.self.s50k.p95.l1000000.edges.csv
```


HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.list.txt
HPRCy1v2genbank.self.s50k.p95.l1000000.paf.vertices.id2name.txt

join <(tr ' ' '\t' < uff | sort -k 1b,1 ) <(tr ' ' '\t' < HPRCy1v2genbank.self.s50k.p95.l1000000.paf.vertices.id2name.txt | sort -k 1b,1) |\
awk '{ print $2, $3 }' | grep -v _ | sort -n | awk '{ c=$1 ; if (c == p) { x=x" "$2 } else { print x; x=$2 }; p = c; }'

join <(tr ' ' '\t' < uff | sort -k 1b,1 ) <(tr ' ' '\t' < HPRCy1v2genbank.self.s50k.p95.l1000000.paf.vertices.id2name.txt | sort -k 1b,1) \
| grep 'chm13\|grch38' | awk '{ print $2, $3 }' | grep -v _ \
| sort -n | awk '{ c=$1 ; if (c == p) { x=x" "$2 } else { print x; x=$2 }; p = c; }'





For the given example, these are as follows:

```shell
chm13#chr1 grch38#chr1
chm13#chr2 grch38#chr2
chm13#chr3 grch38#chr3
chm13#chr4 grch38#chr4
chm13#chr5 grch38#chr5
chm13#chr6 grch38#chr6
chm13#chr7 grch38#chr7
chm13#chr8 grch38#chr8
chm13#chr10 grch38#chr10
chm13#chr11 grch38#chr11
chm13#chr12 grch38#chr12
chm13#chr13 chm13#chr21 grch38#chr13 grch38#chr21
chm13#chr14 chm13#chr22 grch38#chr14 grch38#chr22
chm13#chr15 grch38#chr15
chm13#chr16 grch38#chr16
chm13#chr17 grch38#chr17
chm13#chr18 grch38#chr18
chm13#chr19 grch38#chr19
chm13#chr20 grch38#chr20
chm13#chrX grch38#chrX grch38#chrY
```

A script to reproduce this is in `do_louvain.sh` in this repo.

# Community detection

We evaluate the chromosome separation among the HPRCy1v2 assemblies.

We first apply all-to-all approximate mapping:

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

mkdir -p /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

# Mappings without references included
PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

for s in 50k 10k; do
  for p in 95; do
    # There are 94 haplotypes, so `-n` is equal to `94 - 1`
    sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n 93 -H 0.001 -m > HPRCy1v2genbank.self.s'$s'.p'$p'.n93.h0001.paf && mv HPRCy1v2genbank.self.s'$s'.p'$p'.n93.h0001.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
  done
done

# More stringent mapping for visualization
for s in 50k; do
  for p in 95; do
    for n in 5 7; do
      sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n '$n' -H 0.001 -m > HPRCy1v2genbank.self.s'$s'.p'$p'.n'$n'.h0001.paf && mv HPRCy1v2genbank.self.s'$s'.p'$p'.n'$n'.h0001.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
    done
  done
done

## Mappings with references included
#PATH_HPRCY1_REFS_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank+refs.fa.gz
#for s in 50k 10k; do
#  for p in 98 95; do
#    # There are 96 haplotypes, so `-n` is equal to `96 - 1`
#    sbatch -p workers -c 48 --job-name all-vs-all --wrap 'hostname; cd /scratch; \time -v '$RUN_WFMASH' '$PATH_HPRCY1_REFS_FA_GZ' -t 48 -Y "#" -p '$p' -s '$s' -n 95 -H 0.001 -m > HPRCy1v2genbank+refs.self.s'$s'.p'$p'.n95.h0001.paf && mv HPRCy1v2genbank+refs.self.s'$s'.p'$p'.n95.h0001.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/'
#  done
#done
```

Then, we collect mappings between sequences > `l` kbp long:

```shell
cd /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

for s in 50k 10k; do
  for p in 95; do
    for l in 0 500000 1000000; do
      #for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
      for prefix in HPRCy1v2genbank; do
        echo $s $p $l $prefix
        awk -v len=$l '$2 > len && $7 > len' $prefix.self.s$s.p$p.n93.h0001.paf > $prefix.self.s$s.p$p.n93.h0001.l$l.paf
      done
    done
  done
done
```

To evaluate chromosome communities, we build an "alignment graph" from our mappings using `paf2net` scripts delivered in
`pggb` repository. In this graph, nodes are contigs and edges are mappings between them.

```shell
for s in 50k 10k; do
  for p in 95; do
    for l in 0 500000 1000000; do
      #for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
      for prefix in HPRCy1v2genbank; do
        echo $s $p $l $prefix
        python3 /home/guarracino/tools/pggb/scripts/paf2net.py -p $prefix.self.s$s.p$p.n93.h0001.l$l.paf
      done
    done
  done
done
```

Then we obtain the ``Leiden`` communities:

```shell
# Prepare chromosome info
#( seq 1 22; echo X; echo Y; echo M ) | while read i; do
#  cut -f 1 /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$i.pan.fa.fai | grep chr -v | grep MT -v | awk -v OFS='\t' -v chr=chr$i '{print($0,chr)}' >> contig2chr.tsv
#done

for s in 50k 10k; do
  for p in 95; do
    for l in 0 500000 1000000; do
      #for prefix in HPRCy1v2genbank+refs HPRCy1v2genbank; do
      for prefix in HPRCy1v2genbank; do
        echo $s $p $l $prefix
        
        PAF=$prefix.self.s$s.p$p.n93.h0001.l$l.paf
        
        ID2NAME=$PAF.vertices.id2name.txt
        if [ $prefix == "HPRCy1v2genbank" ]; then
          # Prepare labels by contig
          join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
            <(sort -k 2,2 $PAF.vertices.id2name.txt) \
            <(sort -k 1,1 /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/pq_info/*.partitioning_with_pq.tsv) | cut -f 1,3 -d ' ' > $PAF.vertices.id2name.chr.txt
#            <(sort -k 1,1 contig2chr.tsv) | cut -f 1,3 -d ' ' > $PAF.vertices.id2name.chr.txt
            
          ID2NAME=$PAF.vertices.id2name.chr.txt
        fi
        
        # guix install python-igraph
        python3 /home/guarracino/tools/pggb/scripts/net2communities.py \
          -e $PAF.edges.list.txt \
          -w $PAF.edges.weights.txt \
          -n $ID2NAME --accurate-detection
      done
    done
  done
done
```

Finally, we look for groupings of reference chromosomes in the same community, or we check how chromosomes where partitioned:

```shell
ls HPRCy1v2genbank.self.s50k.p95.n93.h0001.l0.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done
ls HPRCy1v2genbank.self.s50k.p95.n93.h0001.l500000.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done
ls HPRCy1v2genbank.self.s50k.p95.n93.h0001.l1000000.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 5' | grep unmapped -v; done

ls HPRCy1v2genbank.self.s10k.p95.n93.h0001.l0.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done
ls HPRCy1v2genbank.self.s10k.p95.n93.h0001.l500000.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done
ls HPRCy1v2genbank.self.s10k.p95.n93.h0001.l1000000.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | sort | uniq -c | awk '$1 > 0' | grep unmapped -v; done


ls HPRCy1v2genbank+refs.self.s50k.p95.l0.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; grep chr $f; done
```


```shell
( echo "Id,Label,Chromosome"; \
  join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
    <(sort -k 2,2 HPRCy1v2genbank.self.s50k.p95.l1000000.paf.vertices.id2name.txt) \
    <(sort -k 1,1 /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_with_pq/*.partitioning_with_pq.tsv) | tr ' ' ',' \
) > HPRCy1v2genbank.self.s50k.p95.l1000000.nodes.csv
  
( echo "Source,Target,Weight"; \
  paste -d ' ' HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.list.txt HPRCy1v2genbank.self.s50k.p95.l1000000.paf.edges.weights.txt | tr ' ' ','
) > HPRCy1v2genbank.self.s50k.p95.l1000000.edges.csv
```

# Community detection

We evaluate the chromosome separation among the HPRCy1v2 assemblies.

We first apply all-to-all approximate mapping:

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

mkdir -p /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank

# Mappings without references included
PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

p=95
s=50k
l=250k
n=93 # There are 94 haplotypes, so the max `-n` should be equal to `94 - 1`
sbatch -p workers -c 48 --job-name all-vs-all --wrap "hostname; cd /scratch; \time -v $RUN_WFMASH $PATH_HPRCY1_FA_GZ -t 48 -Y '#' -p $p -s $s -l $l -n $n -H 0.001 -m > HPRCy1v2genbank.self.s$s.l$l.p$p.n$n.h0001.paf && mv HPRCy1v2genbank.self.s$s.l$l.p$p.n$n.h0001.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/"


# Mappings with chm13v2 included
PATH_HPRCY1_REFS_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank+chm13v2.fa.gz

w=800
for s in 50k; do
  for ln in 5; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' $ln | bc)
    l=${l_no_k}k
    
    for p in 98; do    
      if [[ $s == "100k" && $p == "98" ]]; then
        w=20000
      elif [[ $s == "100k" && $p == "95" ]]; then
        w=10000
      elif [[ $s == "50k" && $p == "98" ]]; then
        w=10000
      elif [[ $s == "50k" && $p == "95" ]]; then
        w=5000
      elif [[ $s == "20k" && $p == "98" ]]; then
        w=4000
      elif [[ $s == "20k" && $p == "95" ]]; then
        w=2000
      fi
    
      echo $s $p $l $w
      
      # Runs with stringent n values are for visualization
      for n in 5 3; do
        sbatch -p workers -c 48 --job-name all-vs-all --wrap "hostname; cd /scratch; \time -v $RUN_WFMASH $PATH_HPRCY1_REFS_FA_GZ -t 48 -Y '#' -p $p -s $s -l $l -n $n -H 0.001 -w $w -m > HPRCy1v2genbank+chm13v2.self.s$s.l$l.p$p.n$n.h0001.big_w.paf && mv HPRCy1v2genbank+chm13v2.self.s$s.l$l.p$p.n$n.h0001.big_w.paf /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/"
      done
    done
  done
done
```

Then, we collect mappings between sequences at least `l` kbp long:

```shell
cd /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

for s in 50k; do
  for ln in 5; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' $ln | bc)
    l=${l_no_k}k
    
    for p in 95; do  
      for n in 93 5 3; do
        for L in 1000000; do
          for prefix in HPRCy1v2genbank+chm13v2 HPRCy1v2genbank; do
            for big_y in Y N; do
                PREFIX=$prefix.self.s$s.l$l.p$p.n$n.h0001
                if [ $big_y == Y ]; then
                    PREFIX=$prefix.self.s$s.l$l.p$p.n$n.h0001.big_w
                fi
                if [[ -s $PREFIX.paf && ! -s $PREFIX.l$L.paf ]]; then
                  echo $s $p $n $l $prefix $big_y
                  awk -v len=$L '$2 >= len && $7 >= len' $PREFIX.paf > $PREFIX.l$L.paf
                fi
            done
          done
        done
      done
    done
  done
done
```

To evaluate chromosome communities, we build an "alignment graph" from our mappings using `paf2net` scripts delivered in
`pggb` repository. In this graph, nodes are contigs and edges are mappings between them.

```shell
for s in 50k; do
  for ln in 5; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' $ln | bc)
    l=${l_no_k}k
    
    for p in 95; do  
      for n in 93 5 3; do
        for L in 1000000; do
          for prefix in HPRCy1v2genbank+chm13v2 HPRCy1v2genbank; do
            for big_y in Y N; do
              PAF=$prefix.self.s$s.l$l.p$p.n$n.h0001.l$L.paf
              if [ $big_y == Y ]; then
                PAF=$prefix.self.s$s.l$l.p$p.n$n.h0001.big_w.l$L.paf
              fi

              if [[ -s $PAF && ! -s $PAF.edges.list.txt ]]; then
                echo $s $l $p $n $L $prefix $big_y
                python3 /home/guarracino/tools/pggb/scripts/paf2net.py -p $PAF
                
                ID2NAME=$PAF.vertices.id2name.txt
                if [ $prefix == "HPRCy1v2genbank" ]; then
                  # Prepare labels by contig
                  join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
                    <(sort -k 2,2 $PAF.vertices.id2name.txt) \
                    <(sort -k 1,1 /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/pq_info/*.partitioning_with_pq.tsv) | cut -f 1,2,3 -d ' ' | \
                    awk '{print($1,$2"-"$3)}' > $PAF.vertices.id2name.chr.txt
                    
                  ID2NAME=$PAF.vertices.id2name.chr.txt
                fi
              fi
            done
          done
        done
      done
    done
  done
done
```

Then we obtain the ``Leiden`` communities:

```shell
for s in 50k; do
  for ln in 5; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' $ln | bc)
    l=${l_no_k}k
    
    for p in 95; do  
      for n in 93; do
        for L in 1000000; do
          for prefix in HPRCy1v2genbank; do
            for big_y in Y N; do
                PAF=$prefix.self.s$s.l$l.p$p.n$n.h0001.l$L.paf
                if [ $big_y == Y ]; then
                    PAF=$prefix.self.s$s.l$l.p$p.n$n.h0001.big_w.l$L.paf
                fi

                if [[ -s $PAF && ! -s $PAF.edges.weights.txt.community.0.txt ]]; then
                  echo $s $l $p $n $L $prefix $big_y
                                 
                  ID2NAME=$PAF.vertices.id2name.txt
                  if [ $prefix == "HPRCy1v2genbank" ]; then                    
                    ID2NAME=$PAF.vertices.id2name.chr.txt
                  fi
                                 
                  # guix install python-igraph
                  python3 /home/guarracino/tools/pggb/scripts/net2communities.py \
                    -e $PAF.edges.list.txt \
                    -w $PAF.edges.weights.txt \
                    -n $ID2NAME --accurate-detection
                fi
            done          
          done
        done
      done
    done
  done
done
```

[//]: # (For the mappings with the reference&#40;s&#41;, we apply the `Lovain` algorithm too:)
[//]: # ()
[//]: # (```shell)
[//]: # (for s in 20k; do)
[//]: # (  for ln in 5; do)
[//]: # (    s_no_k=${s::-1})
[//]: # (    l_no_k=$&#40;echo $s_no_k '*' $ln | bc&#41;)
[//]: # (    l=${l_no_k}k)
[//]: # (    )
[//]: # (    for p in 98; do  )
[//]: # (      for n in 94 93 5 3; do)
[//]: # (        for L in 1000000; do)
[//]: # (          for prefix in HPRCy1v2genbank+chm13v2; do)
[//]: # (            PAF=$prefix.self.s$s.l$l.p$p.n$n.h0001.l$L.paf)
[//]: # (          )
[//]: # (            /home/guarracino/tools/pafnet/target/release/pafnet $PAF -e > $PAF.edgelist)
[//]: # (            /home/guarracino/tools/pafnet/target/release/pafnet $PAF -l > $PAF.names)
[//]: # (            /home/guarracino/tools/louvain/louvain $PAF.edgelist $PAF.louvain)
[//]: # (          done)
[//]: # (        done)
[//]: # (      done)
[//]: # (    done)
[//]: # (  done)
[//]: # (done)
[//]: # (```)

Finally, we check how chromosomes where partitioned:

```shell
cd /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

s=50k
l=250k
p=95
n=93
L=1000000
PAF=HPRCy1v2genbank.self.s$s.l$l.p$p.n$n.h0001.l$L.paf

python3 /lizardfs/guarracino/chromosome_communities/scripts/get_table_about_communities.py $PAF.edges.weights.txt.community.*.txt | cut -f 1-25,27 > $PAF.community.leiden.tsv
cat $PAF.community.leiden.tsv | column -t

ls $PAF.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | cut -f 2 -d '-' | cut -f 2 -d '#' | sort | uniq -c | sort -k 1nr | awk '$1 > 0' ; done

# File for mapping contigs to communities
ls $PAF.edges.weights.txt.community.*.txt | while read f; do
  n=$(echo $f | rev | cut -f 2 -d '.' | rev)
  cut -f 1 $f -d '-' | awk -v OFS='\t' -v n=$n '{print($1,n)}' >> $PAF.contig2community.tsv
done


# Total number of contigs
cat /lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz.fai | wc -l

# Total length of the contigs
cat /lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz.fai | \
  awk -F'\t' 'BEGIN{LEN=0}{ LEN+=$2 }END{print LEN}'

# Mapped contigs
cut -f 1,6 $PAF | tr '\t' '\n' | sort | uniq | wc -l

# Filtered contigs (1 Mbps)
cut -f 1,6 HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf | tr '\t' '\n' | sort | uniq | wc -l
16118

cat \
  <( cut -f 1,2 $PAF | sort | uniq ) \
  <( cut -f 6,7 $PAF | sort | uniq ) | sort | uniq | \
  awk -F'\t' 'BEGIN{LEN=0}{ LEN+=$2 }END{print LEN}'


# Contigs of the communities containing only not partitioned contigs
grep -f <(cat \
  HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf.edges.weights.txt.community.28.txt \
  HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf.edges.weights.txt.community.29.txt \
  HPRCy1v2genbank.self.s50k.l250k.p95.n93.h0001.l1000000.paf.edges.weights.txt.community.30.txt | cut -f 1 -d '-') \
  /lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz.fai 
```

[//]: # (```shell)
[//]: # (# OTHER STUFF)
[//]: # (s=50k)
[//]: # (l=250k)
[//]: # (p=95)
[//]: # (n=5)
[//]: # (L=1000000)
[//]: # (ls HPRCy1v2genbank+chm13v2.self.s$s.l$l.p$p.n$n.h0001.big_w.l$L.paf.edges.weights.txt.community.*.txt | while read f; do echo $f; cat $f | grep chr; done)
[//]: # ()
[//]: # (s=20k)
[//]: # (l=100k)
[//]: # (p=98)
[//]: # (n=94)
[//]: # (L=1000000)
[//]: # (PAF=HPRCy1v2genbank+chm13v2.self.s$s.l$l.p$p.n$n.h0001.l$L.paf)
[//]: # ()
[//]: # (# Info communities)
[//]: # (&#40;echo -e "id.community\tnum.contigs\tchm13.chromosome"; \)
[//]: # (join -a 1 -e 'EMPTY' \)
[//]: # (  <&#40; cut -f 2 -d ' ' $PAF.louvain | sort -n | uniq -c | awk -v OFS='\t' '{print $2,$1}' | sort &#41; \)
[//]: # (  <&#40; join <&#40;tr ' ' '\t' < $PAF.louvain | sort -k 1b,1&#41; <&#40;tr ' ' '\t' < $PAF.names | sort -k 1b,1&#41; | grep 'chm13\|grch38' | awk -v OFS='\t' '{ print $2, $3 }' | sort &#41; | \)
[//]: # (  sort -k 1,1n&#41; | tr ' ' '\t' > $PAF.louvain.tsv)
[//]: # ()
[//]: # (# OLD)
[//]: # (join <&#40;tr ' ' '\t' < $PAF.louvain | sort -k 1b,1 &#41; \)
[//]: # (     <&#40;tr ' ' '\t' < $PAF.names | sort -k 1b,1&#41; \)
[//]: # (    | grep 'chm13\|grch38' | awk '{ print $2, $3 }' | grep -v _ \)
[//]: # (    | sort -n | awk '{ c=$1 ; if &#40;c == p&#41; { x=x" "$2 } else { print x; x=$2 }; p = c; }')
[//]: # (```)

Obtain files for `gephi`:

```shell
cd /lizardfs/guarracino/chromosome_communities/mappings/HPRCy1v2genbank/

s=50k
l=250k
p=95
n=3
L=1000000
PAF=HPRCy1v2genbank+chm13v2.self.s$s.l$l.p$p.n$n.h0001.big_w.l$L.paf

# Add labels only to acrocentric chromosomes
( echo "Id,Label,Chromosome"; \
  join -1 2 -2 1 -a 1 -e 'unmapped' -o '1.1 1.2 2.2' \
    <(sort -k 2,2 $PAF.vertices.id2name.txt) \
    <(cat \
      <(sort -k 1,1 /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/pq_info/*.partitioning_with_pq.tsv |\
          sed 's/chm13#//' | sed 's/grch38#//') \
      <(echo -e chm13#chr1"\t"chr1_pq) \
      <(echo -e chm13#chr2"\t"chr2_pq) \
      <(echo -e chm13#chr3"\t"chr3_pq) \
      <(echo -e chm13#chr4"\t"chr4_pq) \
      <(echo -e chm13#chr5"\t"chr5_pq) \
      <(echo -e chm13#chr6"\t"chr6_pq) \
      <(echo -e chm13#chr7"\t"chr7_pq) \
      <(echo -e chm13#chr8"\t"chr8_pq) \
      <(echo -e chm13#chr9"\t"chr9_pq) \
      <(echo -e chm13#chr10"\t"chr10_pq) \
      <(echo -e chm13#chr11"\t"chr11_pq) \
      <(echo -e chm13#chr12"\t"chr12_pq) \
      <(echo -e chm13#chr13"\t"chr13_pq) \
      <(echo -e chm13#chr14"\t"chr14_pq) \
      <(echo -e chm13#chr15"\t"chr15_pq) \
      <(echo -e chm13#chr16"\t"chr16_pq) \
      <(echo -e chm13#chr17"\t"chr17_pq) \
      <(echo -e chm13#chr18"\t"chr18_pq) \
      <(echo -e chm13#chr19"\t"chr19_pq) \
      <(echo -e chm13#chr20"\t"chr20_pq) \
      <(echo -e chm13#chr21"\t"chr21_pq) \
      <(echo -e chm13#chr22"\t"chr22_pq) \
      <(echo -e chm13#chrX"\t"chrX_pq) \
      <(echo -e chm13#chrY"\t"chrY_pq) \
      <(echo -e grch38#chr1"\t"chr1_pq) \
      <(echo -e grch38#chr2"\t"chr2_pq) \
      <(echo -e grch38#chr3"\t"chr3_pq) \
      <(echo -e grch38#chr4"\t"chr4_pq) \
      <(echo -e grch38#chr5"\t"chr5_pq) \
      <(echo -e grch38#chr6"\t"chr6_pq) \
      <(echo -e grch38#chr7"\t"chr7_pq) \
      <(echo -e grch38#chr8"\t"chr8_pq) \
      <(echo -e grch38#chr9"\t"chr9_pq) \
      <(echo -e grch38#chr10"\t"chr10_pq) \
      <(echo -e grch38#chr11"\t"chr11_pq) \
      <(echo -e grch38#chr12"\t"chr12_pq) \
      <(echo -e grch38#chr13"\t"chr13_pq) \
      <(echo -e grch38#chr14"\t"chr14_pq) \
      <(echo -e grch38#chr15"\t"chr15_pq) \
      <(echo -e grch38#chr16"\t"chr16_pq) \
      <(echo -e grch38#chr17"\t"chr17_pq) \
      <(echo -e grch38#chr18"\t"chr18_pq) \
      <(echo -e grch38#chr19"\t"chr19_pq) \
      <(echo -e grch38#chr20"\t"chr20_pq) \
      <(echo -e grch38#chr21"\t"chr21_pq) \
      <(echo -e grch38#chr22"\t"chr22_pq) \
      <(echo -e grch38#chrX"\t"chrX_pq)  \
      <(echo -e grch38#chrY"\t"chrY_pq) | sort)  | \
        awk '{print($1,$2,$3)}' | tr ' ' ',' \
  ) > $PAF.nodes.tmp.csv

# Add community information
join -1 1 -2 2 \
  <(cut -f 1,2 HPRCy1v2genbank.self.s$s.l$l.p$p.n93.h0001.l$L.paf.community.leiden.tsv | sort -k 1) \
  <(sort -k 2 HPRCy1v2genbank.self.s$s.l$l.p$p.n93.h0001.l$L.paf.contig2community.tsv) |  awk -v OFS='\t' '{print $3,$2}' > $PAF.contig2namedCommunity.tsv

(echo "Id,Label,Chromosome,Community"; \
  join -1 2 -2 1 \
    <(sed '1d' $PAF.nodes.tmp.csv | tr ',' ' ' | sort -k 2) \
    <(sort -k 1 $PAF.contig2namedCommunity.tsv) | \
      awk -v OFS=',' '{print $2,$1,$3,$4}') > $PAF.nodes.tmp2.csv

# Add color column (for the "givecolortonodes" plugin [that doesn't work!])
(echo "Id,Label,Chromosome,Community,ColorOfNode"; \
  join -1 3 -2 1 \
    <(sed '1d' $PAF.nodes.tmp2.csv | tr ',' ' ' | sort -k 3) \
    <(sed '1d' /lizardfs/guarracino/chromosome_communities/data/chromosome.colors.csv | tr ',' ' ' | sort -k 1) | \
      awk -v OFS=',' '{print $2,$3,$1,$4,$5}') \
     > $PAF.nodes.csv

rm $PAF.contig2namedCommunity.tsv $PAF.nodes.tmp*.csv


( echo "Source,Target,Weight"; \
  paste -d ' ' $PAF.edges.list.txt $PAF.edges.weights.txt | tr ' ' ','
) > $PAF.edges.csv
```

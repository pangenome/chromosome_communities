# Chromosome communities in the human pangenome

This repo contains scripts to evaluate the chromosome separation among the HPRCy1v2 assemblies.

We first apply all-to-all approximate mapping:

```
wfmash -m -t 48 -Y "#" -p 98 -s 100000 -l 300000 -n 90 hprcy1v2.fa hprcy1v2.fa >hprcy1v2.self.paf
```

Then, we collect mappings between sequences > 1Mbp long:

```
awk '$2 > 1000000 && $7 > 1000000' hprc.y1v2.self.paf >hprc.y1v2.self.1M.paf
```

To evaluate chromosome communities, we build an "alignment graph" from our mappings using [pafnet](https://github.com/ekg/pafnet).
In this graph, nodes are contigs and edges are mappings between them.

```
o=hprc.y1v2.self.1M.paf
pafnet $f -e >$o.edgelist
pafnet $f -l >$o.names
```

Then we obtain the Louvain community assignment using [an implementation of the method](https://github.com/jlguillaume/louvain).

```
louvain $o.edgelist $o.louvain
```

Finally, we look for groupings of reference chromosomes in the same community.

```
join <(tr ' ' '\t' <$o.louvain | sort -k 1b,1 ) \
     <(tr ' ' '\t' <$o.names | sort -k 1b,1) \
    | grep 'chm13\|grch38' | awk '{ print $2, $3 }' | grep -v _ \
    | sort -n | awk '{ c=$1 ; if (c == p) { x=x" "$2 } else { print x; x=$2 }; p = c; }'
```

For the given example, these are as follows:

```
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

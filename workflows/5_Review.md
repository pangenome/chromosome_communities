# Repetitive region centered on the SST1 array (chr13, chr14, chr21)

Flip the acrocentric graph:

```shell
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-0b21b3525fc2ed9305c7df2386475a008a9337bd 

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

sbatch -p workers -c 48 --job-name acro-flip --wrap "hostname; cd /scratch && $RUN_ODGI flip -i $path_input_og -o ${prefix}.flip.og -t 48 -P; mv ${prefix}.flip.og /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/"
```

Save the path names with the information of the flipped path, and rename them by removing the '_inv' suffix:

```shell
path_flipped_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.og
prefix=$(basename $path_flipped_og .og)

$RUN_ODGI paths -i $path_flipped_og -L > $prefix.path_names.txt

$RUN_ODGI view -i $path_flipped_og -g | sed '/^P/ s/_inv//' > $prefix.gfa
rm $path_flipped_og
$RUN_ODGI build -g $prefix.gfa -o $path_flipped_og -t 48 -P
```

Plots:

```shell
#chm13#chr13  12301367  12440010  SST1#222222
path_flip_grounded_pq_touching_reliable_tsv_gz=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#chr13.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_flip_grounded_pq_touching_reliable_tsv_gz \
  11301367 13440010 \
  90 0.9 \
  1.0 \
  3 1 \
  13 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_402b1_506680_chr13_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX6.chr13.SST1.1Mbps.n1.nref1.pdf

#chm13#chr14  6960008 6988409 SST_Composite#222222
path_flip_grounded_pq_touching_reliable_tsv_gz=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#chr14.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_flip_grounded_pq_touching_reliable_tsv_gz \
  5960008 7988409 \
  90 0.9 \
  1.0 \
  3 1 \
  14 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_92a4_50bdd0_chr14_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX7.chr14.SST1.1Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_flip_grounded_pq_touching_reliable_tsv_gz=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_flip_grounded_pq_touching_reliable_tsv_gz \
  8375567 10453313 \
  90 0.9 \
  1.0 \
  3 1 \
  21 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_bf20_50ca50_chr21_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX8.chr21.SST1.1Mbps.n1.nref1.pdf
```

# Plot FLAGGER's annotations (for pq-contig quality control)

Take only reliable blocks [flagged with "Hh" or "Hc"](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components).

```shell
path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt

flipped=false
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

rm xyz.tsv z.tsv # Cleaning

for e in 50000; do
  for m in 1000 ; do
    cat $path_targets_txt | while read ref; do 
      echo "-e $e -m $m $ref"

      # Plot additional tracks for the annotations
      path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz

      zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
        SAMPLE=$( echo $CONTIG | cut -f 1 -d '#');
        
        path_unreliable_bed=/lizardfs/guarracino/HPRC/annotations/unreliable/$SAMPLE.hifi.flagger_final.simplified.unreliable_only.bed
        if [[ -s $path_unreliable_bed ]]; then
          echo $CONTIG "--->" $SAMPLE
           
          zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{ if($10 <= 1 && $15 <= 1) {print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16)} }' > x.$CONTIG.bed
          
          # To get only Unreliable intervals
          grep $CONTIG $path_unreliable_bed | cut -f 1,2,3,4 > y.$CONTIG.bed
          # To get both Reliable and Unreliable intervals
          #cat \
          #  <( grep $CONTIG $path_unreliable_bed | cut -f 1,2,3,4) \
          #  <( grep $CONTIG $path_unreliable_bed | bedtools sort | \
          #    bedtools complement -i - -g <(grep $CONTIG /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1,2) | awk -v OFS='\t' '{print($1,$2,$3,"Reliable")}' ) \
          #  > y.$CONTIG.bed
            
          # wao: write the original A and B entries plus the number of base pairs of overlap between the two features.
          bedtools intersect -a x.$CONTIG.bed -b y.$CONTIG.bed -wo >> z.tsv
          
          rm x.$CONTIG.bed y.$CONTIG.bed
        fi
      done
        
      cat \
        <( zcat $path_grounded_pq_touching_tsv_gz | sed '1d' ) \
        <( python3 /lizardfs/guarracino/chromosome_communities/scripts/get_annotation_track.py z.tsv ) | tr ' ' '\t' >> xyz.tsv
      rm z.tsv
    done
    
    ref=$( head -n 1 $path_targets_txt ) # The first (all the same for getting the header line)
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz
    cat \
      <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
       xyz.tsv | pigz -c > e$e.m$m.annot.tsv.gz
    rm xyz.tsv
  done
done

# Supplementary Figures
for e in 50000; do
  for m in 1000 ; do
    (seq 13 15; seq 21 22) | while read i; do
      echo "-e $e -m $m chr$i"
      
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_and_flagger.R \
        e$e.m$m.annot.tsv.gz \
        0 25000000 \
        90 \
        0.9 1.0 \
        1 1 \
        $i \
        0.9 \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        SupplementaryFigureX.e$e.m$m.annot.chr$i.pdf
    done
    
    # Merge chromosomes's PDF files
    #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
    #  e$e.m$m.annot.chr*.pdf e$e.m$m.annot.chrACRO.pdf
  done
done
```

```shell
# OLD CODE
e=50000
m=1000
j=0.8
j_str=$(echo $j | sed 's/\.//g')
n=1
ref=chm13#chr13

# Brutal: remove untangled regions if they overlap with unreliable regions
touch xyz.tsv
cat $path_targets_txt | while read ref; do
    echo $ref
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
    
    touch z.tsv
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
      SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')
    
      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed
      if [[ -s $path_unreliable_bed ]]; then
        echo $CONTIG "--->" $SAMPLE
        
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13)}' > x.bed
        grep $CONTIG $path_unreliable_bed > y.bed
        # -A: remove entire feature if any overlap
        bedtools subtract -a x.bed -b y.bed -A |\
          awk -v OFS='\t' '{split($4, a, "_"); print($1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10])}' >> xyz.tsv
      else
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n >> xyz.tsv
      fi
    done
    
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' >> xyz.tsv
done

path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
cat \
  <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
   xyz.tsv | pigz -c > xyz.tsv.gz
rm xyz.tsv

Rscript scripts/plot_untangle_all.R xyz.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 360 200

# Plot additional tracks for the annotations
touch xyz.tsv
cat $path_targets_txt | while read ref; do
    echo $ref
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
    
    touch z.tsv
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
      SAMPLE=$( echo $CONTIG | cut -f 1 -d '#');
    
      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed
      if [[ -s $path_unreliable_bed ]]; then
        echo $CONTIG "--->" $SAMPLE
        
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13)}' > x.bed
        grep $CONTIG $path_unreliable_bed > y.bed
        
        # wao: write the original A and B entries plus the number of base pairs of overlap between the two features.
        bedtools intersect -a x.bed -b y.bed -wo >> z.tsv
      fi
    done
    
    cat \
      <( zcat $path_grounded_pq_touching_tsv_gz | sed '1d' ) \
      <( python3 scripts/get_annotation_track.py z.tsv ) | tr ' ' '\t' >> xyz.tsv

    rm z.tsv
done

path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
cat \
  <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
   xyz.tsv | pigz -c > xyz.annot.tsv.gz
rm xyz.tsv

Rscript scripts/plot_untangle_all.R xyz.annot.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 120 200
```


# Inversion state of chr13, chr14, chr15 around the SST1 array in the untangle space

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
  90 0.9 12 \
  1.0 \
  3 1 \
  13 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr13_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX6.chr13.SST1.1Mbps.n1.nref1.pdf

#chm13#chr14  6960008 6988409 SST_Composite#222222
path_flip_grounded_pq_touching_reliable_tsv_gz=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#chr14.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_flip_grounded_pq_touching_reliable_tsv_gz \
  5960008 7988409 \
  90 0.9 12 \
  1.0 \
  3 1 \
  14 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr14_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX7.chr14.SST1.1Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_flip_grounded_pq_touching_reliable_tsv_gz=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_flip_grounded_pq_touching_reliable_tsv_gz \
  8375567 10453313 \
  90 0.9 12 \
  1.0 \
  3 1 \
  21 \
  0.9 \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_flip_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr21_SST1_1Mbps.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/SupplementaryFigureX8.chr21.SST1.1Mbps.n1.nref1.pdf
```


# Dot plots (to hypothesize recombination mechanisms for the graphical abstract)

## With Gepard

Extract regions centered in the SST1 array:

```shell
cd /lizardfs/guarracino/chromosome_communities

samtools faidx assemblies/chm13.chr13.fa.gz chm13#chr13:11301367-13440010 | sed 's/>chm13#/>/g' > chr13.SST1.1Mbp.fa
samtools faidx assemblies/chm13.chr13.fa.gz chm13#chr14:5960008-7988409 | sed 's/>chm13#/>/g' > chr14.SST1.1Mbp.fa
samtools faidx assemblies/chm13.chr13.fa.gz chm13#chr21:8375567-10453313 | sed 's/>chm13#/>/g' > chr21.SST1.1Mbp.fa
samtools faidx chr_13_14_21.p_arms.fa.gz chr14:5960008-7988409 >> chr_13_14_21.SST1.1Mbps.fa

java -jar Gepard-2.1.jar
```

Open `gepard`, load pairs of sequences, and use several `word length` values to change the level of detail:

```shell
java -jar Gepard-2.1.jar
```

## With odgi untangle's output

Dot plots from the untangle output:

```shell
ref=chr13
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_13_27_SST1" \
  12301367 12440010 \
  1000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.1Mbp

ref=chr14
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_14_39_SST_Composite" \
  6960008 6988409 \
  1000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.1Mbp

ref=chr21
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_21_45_SST1_Composite" \
  9175567 9653313 \
  1000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.1Mbp


ref=chr13
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_13_27_SST1" \
  12301367 12440010 \
  6000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.6Mbps

ref=chr14
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_14_39_SST_Composite" \
  6960008 6988409 \
  6000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.6Mbps

ref=chr21
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_21_45_SST1_Composite" \
  9175567 9653313 \
  6000000 \
  5 \
  <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai) \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}.6Mbps
```


# Evolutionary strata

```shell
e=50000
m=1000
eid=0.900
eid_str=$(echo $eid | sed 's/\.//g')
n=2
path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/chm13.SEX.target_paths.txt
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
prefix=$(basename $path_input_og .og)

path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}
        
#chrX#PAR1   0  2394410 su CHM13
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_BED_annotation.R \
  $path_entropy_match_order_tsv \
  0 2894410 \
  40 \
  'X' \
  /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chrX.PAR1.pdf
# The entropy (the metric we use to find PHRs) is greater than 0 up to 2435289 => 2435289−2394410 = 40879 bp

#chrX#PAR2	153925834	154259566
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_BED_annotation.R \
  $path_entropy_match_order_tsv \
  153425834 154259566 \
  40 \
  'X' \
  /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chrX.PAR2.pdf
Su PAR2 non sforiamo

#chrX#XTR	87642550	91570785
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_BED_annotation.R \
  $path_entropy_match_order_tsv \
  87392550 91820785 \
  40 \
  'X' \
  /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chrX.XTR.pdf
Su XTR sforiamo un poco di 2636 bp sulla destra.

path_entropy_match_order_tsv=/home/guarracino/Downloads/Pangenomics/chromosome_communities/Review1/evolutionary_strata/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n2.tsv
cat $path_entropy_match_order_tsv | sed '1d' | awk '$4 > 0 && $5 > 0' | \
  bedtools merge -i - -d 10000 -c 4,5 -o mean | \
  awk '$3 - $2 > 0' | sed 's/chm13#//' > chrSEX_29-Nov-22_PHRs.bed # Supplementary Table
```


# Robertsonian translocation breakpoints

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/robertsonian_translocation
cd /lizardfs/guarracino/chromosome_communities/robertsonian_translocation

# BAC clones from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4257996/
# Our results from PCR and FISH analysis showed that only the clones CR382285, CR382287,
# and a small fragment of CR382332 are retained in the examined ROBs. [...].
# Given our results, we propose localization of the breakpoints in or nearby to clone CR382332. 
# https://www.ncbi.nlm.nih.gov/nuccore/CR382285
# https://www.ncbi.nlm.nih.gov/nuccore/CR382287
# https://www.ncbi.nlm.nih.gov/nuccore/CR381572
# https://www.ncbi.nlm.nih.gov/nuccore/CR381535
# https://www.ncbi.nlm.nih.gov/nuccore/CR381653
# https://www.ncbi.nlm.nih.gov/nuccore/CR382332
# https://www.ncbi.nlm.nih.gov/nuccore/CR381670
# https://www.ncbi.nlm.nih.gov/nuccore/CR392039
# Download the BAC clones:
# - Put all IDs in a text file, one per line
# - Go to https://www.ncbi.nlm.nih.gov/sites/batchentrez
# - Upload the IDs list
# - Select all files
# - Download in FASTA format
# - mv sequence.fasta PMC4257996.clones.fasta
# - bgzip PMC4257996.clones.fasta
# - Upload: scp PMC4257996.clones.fasta.gz guarracino@octopus02:/lizardfs/guarracino/chromosome_communities/robertsonian_translocation

samtools faidx PMC4257996.clones.fasta.gz

RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9
$RUN_WFMASH /lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz PMC4257996.clones.fasta.gz -t 48 \
  -p 90 -s 1k -n 10 -N -m > PMC4257996.clones.vs.CHM13.p90.s1k.n10.N.paf

# Check which clone(s) are not mapped
comm -23 \
  <(cut -f 1 PMC4257996.clones.fasta.gz.fai | sort) \
  <(cut -f 1 PMC4257996.clones.vs.CHM13.p90.s1k.n10.N.paf | sort | uniq) > clones.unaligned.txt
#if [[ $(wc -l clones.unaligned.txt | cut -f 1 -d ' ' ) != 0 ]];
#    samtools faidx PMC4257996.clones.fasta.gz $(tr '\n' ' ' < clones.unaligned.txt) > clones.unaligned.fa
#    samtools faidx clones.unaligned.fa
#
#    # Allow query split
#    $RUN_WFMASH /lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz clones.unaligned.fa -t 48 \
#  -p 90 -s 1k -n 10 -m > clones.unaligned.vs.CHM13.p90.s1k.n10.N.paf#
#fi

# How much PHRs do the mappings cover?

min_id=99 #90 (all) and 99

A_and_B=$(bedtools intersect \
  -a <(cut -f 1,2,3 /lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed | grep '^chr13\|^chr14\|^chr15\|^chr21\|^chr22' | awk -v OFS='\t' '{print("chm13#"$0,"100.0","PHRs")}' | awk -v min_id=$min_id '$4 >= min_id' | bedtools sort -i -) \
  -b <(sed 's/id:f://' PMC4257996.clones.vs.CHM13.p90.s1k.n10.N.paf | awk -v OFS='\t' '{print($6,$8,$9,$13,$1)}' | grep '^chm13#chr13\|^chm13#chr14\|^chm13#chr15\|^chm13#chr21\|^chm13#chr22' | awk -v min_id=$min_id '$4 >= min_id' | bedtools sort -i -)  | awk '{sum+=$3-$2}END{print(sum)}')

A=$(cut -f 1,2,3 /lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed | grep '^chr13\|^chr14\|^chr15\|^chr21\|^chr22' | awk -v OFS='\t' '{print("chm13#"$0,"100.0","PHRs")}' | awk -v min_id=$min_id '$4 >= min_id' | awk '{sum+=$3-$2}END{print(sum)}')
B=$(sed 's/id:f://' PMC4257996.clones.vs.CHM13.p90.s1k.n10.N.paf | awk -v OFS='\t' '{print($6,$8,$9,$13,$1)}' | grep '^chm13#chr13\|^chm13#chr14\|^chm13#chr15\|^chm13#chr21\|^chm13#chr22' | awk -v min_id=$min_id '$4 >= min_id' | awk '{sum+=$3-$2}END{print(sum)}')

perc_PHRs_covered=$(echo "$A_and_B / $A" | bc -l)
perc_MAPs_covered=$(echo "$A_and_B / $B" | bc -l)
echo "A = $A"
echo "B = $B"
echo "A and B = $A_and_B"
echo "%PHRs covered = $perc_PHRs_covered"
echo "%MAPPING covered$ = $perc_MAPs_covered"

# min_id 90
A = 18329090
B = 6063106
A and B = 4714652
%PHRs covered = .25722237165074752756
%MAPPING covered$ = .77759682908397115273


# min_id 99
A = 18329090
B = 2675655
A and B = 1964707
%PHRs covered = .10719064612591241572
%MAPPING covered$ = .73429010840336291487
```

Plots:

```shell
# Prepare the BED with PHRs, mappings and annotations on the SST1\rDNA
cat \
  <(cut -f 1,2,3 /lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed | awk -v OFS='\t' '{print("chm13#"$0,"100.0","PHRs")}') \
  <(sed 's/id:f://' PMC4257996.clones.vs.CHM13.p90.s1k.n10.N.paf | awk -v OFS='\t' '{print($6,$8,$9,$13,$1)}') \
  <(grep SST /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.rDNA.acros.SST1.200kbps.approximate.white_grch38_Ns.bed | grep 2222 | cut -f 1,2 -d '#' | cut -f 1,2,3 | awk -v OFS='\t' '{print($0,"100.0","SST1")}') \
  <(grep rDNA /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.rDNA.acros.SST1.200kbps.approximate.white_grch38_Ns.bed | cut -f 1,2 -d '#' | cut -f 1,2,3 | awk -v OFS='\t' '{print($0,"100.0","rDNA")}') \
  > chrACRO_7-Dec-22_PHRs+PMC4257996+SST1+rDNA.bed

# Use the plot_PHRs_and_clone_CR382332_mappings.R script
```


# Recombination hotspots (PRDM9)

Prepare the tools:

```shell
cd /home/guarracino/tools/

# Download and install the MEME Suite
wget -c https://meme-suite.org/meme/meme-software/5.5.0/meme-5.5.0.tar.gz
tar -xf meme-5.5.0.tar.gz 
cd ~/tools/meme-5.5.0
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make

# Download and compile TideHunter
wget https://github.com/yangao07/TideHunter/releases/download/v1.5.4/TideHunter-v1.5.4.tar.gz
tar -zxvf TideHunter-v1.5.4.tar.gz
cd TideHunter-v1.5.4
make
```

Prepare PRDM9 motifs (downloaded from [here](https://pubmed.ncbi.nlm.nih.gov/29072575/):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination_hotspots
cd /lizardfs/guarracino/chromosome_communities/recombination_hotspots

# Download the PWMs for all PRDM9 motifs
wget -c https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5705219/bin/elife-28383-fig1-data2.txt
# All 17 motif logos returned by our motif-finding algorithm are listed, 
# along with histograms indicating their positions within the central 300 bp
# of our human PRDM9 peaks, as a measure of how centrally enriched they are
# (and therefore likely to represent true binding targets). 
# Only the seven motifs for which greater than 85% of occurrences within peaks
# are within 100 bp of the peak center were retained for downstream analyses.
# The remaining, less centrally enriched, motifs are either degenerate (as seen
# in mice containing the human allele: (Davies et al., 2016) or may arise as a 
# consequence of PRDM9 binding to promoter regions 

# Remove non-human motifs
n=$(grep 'MOTIF Chimp1' elife-28383-fig1-data2.txt -n | cut -f 1 -d ':')
n=$((n-2))
head -n $n elife-28383-fig1-data2.txt > PRDM9_motifs.human.txt
```

Search the PRDM9 motifs in the whole CHM13 v2.0:

```shell
samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz $(grep chm /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz.fai | cut -f 1) > chm13v2.fa

RUN_FIMO=/home/guarracino/tools/meme-5.5.0/src/fimo
sbatch -p workers -c 48 --job-name meme-PRDM9 --wrap "hostname; cd /scratch && $RUN_FIMO --oc /lizardfs/guarracino/chromosome_communities/recombination_hotspots/ --verbosity 1 --thresh 1.0E-4 /lizardfs/guarracino/chromosome_communities/recombination_hotspots/PRDM9_motifs.human.txt /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.fa"
```

Convert the output in BED format:

```shell
# Rename all files
mv fimo.gff chm13v2.PRDM9.gff
mv fimo.html chm13v2.PRDM9.html
mv fimo.tsv chm13v2.PRDM9.tsv
mv fimo.xml chm13v2.PRDM9.xml
mv cisml.xml cisml.chm13v2.PRDM9.xml

# https://meme-suite.org/meme/doc/fimo-output-format.html#tsv_results
#	The start position of the motif occurrence; 1-based sequence coordinates.
#	The end position of the motif occurrence; 1-based sequence coordinates.

# Remove the last lines, remove the header line, remove the last empty line, prepare the columns, and sort the BED file
grep '^#' chm13v2.PRDM9.tsv -v | sed '1d' | sed '/^$/d' | awk -v OFS='\t' '{print($2,$3-1,$4,$1,$6,$5,$7,$8,$9)}' | bedtools sort > chm13v2.PRDM9.bed
```

Counts the number of hits in windows:

```shell
#TODO IF NEEDED: separate in Human1 hits, Human2 hits, ...

samtools faidx /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.fa

max_qvalue=1
window_size=20000

rm chm13v2.PRDM9.chrACRO.w${window_size}.bed
(seq 13 15; seq 21 22) | while read i; do
  echo $i

  bedtools intersect \
    -a <(bedtools makewindows -g <(cat /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.fa.fai | grep "chm13#chr$i" | cut -f 1,2) -w $window_size) \
    -b <(grep chm13#chr$i chm13v2.PRDM9.bed | grep -P 'Human[1-7]*[0-9]\t' | awk -v max_qvalue=$max_qvalue '$8 <= max_qvalue') -c \
    >> chm13v2.PRDM9.chrACRO.w${window_size}.bed
done

rm chm13v2.PRDM9.chrSEX.w${window_size}.bed
(echo X; echo Y) | while read i; do
  echo $i

  bedtools intersect \
    -a <(bedtools makewindows -g <(cat /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.fa.fai | grep "chm13#chr$i" | cut -f 1,2) -w $window_size) \
    -b <(grep chm13#chr$i chm13v2.PRDM9.bed | grep -P 'Human[1-7]*[0-9]\t' | awk -v max_qvalue=$max_qvalue '$8 <= max_qvalue') -c \
    >> chm13v2.PRDM9.chrSEX.w${window_size}.bed
done
```

Plot the number of hits in each window across the whole chromosomes:

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_without_annotation.all_chromosomes.R \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.chrACRO.w${window_size}.bed \
  1000000 \
  'Position (Mbp)' \
  '' \
  35 \
  'Chromosome' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.chrACRO.w${window_size}.pdf

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_without_annotation.all_chromosomes.R \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.chrSEX.w${window_size}.bed \
  1000000 \
  'Position (Mbp)' \
  '' \
  35 \
  'Chromosome' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.chrSEX.w${window_size}.pdf
```

Plot the number of hits in each window across a chromosome region, with annotation on the top:

```shell
(seq 13 15; seq 21 22) | while read i; do
  Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_with_annotation.R \
    /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.chrACRO.w${window_size}.bed \
    0 25000000 \
    $i \
    35 \
    /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
    /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.chr$i.PRDM9.w${window_size}.pdf
done

# --delta for white space between the pieces
pdfjam --delta '0 7' --no-landscape --nup 1x5 \
  chm13v2.chr13.PRDM9.w${window_size}.pdf \
  chm13v2.chr14.PRDM9.w${window_size}.pdf \
  chm13v2.chr15.PRDM9.w${window_size}.pdf \
  chm13v2.chr21.PRDM9.w${window_size}.pdf \
  chm13v2.chr22.PRDM9.w${window_size}.pdf \
  --outfile output.pdf

pdfcrop --margins "1 1 1 1" output.pdf SupplementaryFigureX.chm13v2.PRDM9.chrACRO.with_annotation.w${window_size}.pdf
rm output.pdf chm13v2.chr1*.PRDM9.w${window_size}.pdf chm13v2.chr2*.PRDM9.w${window_size}.pdf
```

## rDNA/SST1 array repetitive unit

### rDNA

Obtain the repetitive unit of the array:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA
cd /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v1.1.rdna_units.bed

# Get all possible rDNA repetitie units with different lengths (> 40kbps)

cat chm13v1.1.rdna_units.bed | awk -v OFS='\t' '{print($0,($3-$2)/1000)}' | awk '!a[$1"-"$5]++' | grep \
  -f <(cat chm13v1.1.rdna_units.bed | awk -v OFS='\t' '{print($0,($3-$2)/1000))}' | cut -f 5 | awk '$1 > 0' | sort | uniq) \
  > chm13v1.1.rdna_units.unique.bed

bedtools getfasta -fi /lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa -bed <(cat chm13v1.1.rdna_units.unique.bed | cut -f 1,2,3 | sed 's/chr/chm13#chr/g') > chm13v1.1.rdna_units.unique.fa


# IT DOES NOT WORK
#samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr13:5770549-9348041 chm13#chr14:2099538-2817811 chm13#chr15:2506443-4707485 chm13#chr21:3108299-5612715 chm13#chr22:4793795-5720650 > chm13.rDNA.fa
#
# https://github.com/yangao07/TideHunter#tabular-format
#$RUN_TIDEHUNTER -f 2 chm13.rDNA.fa -t 48 -k 13 > chm13.rDNA.TideHunter.tsv
#awk '{print(">"$1"_"$7"\n"$11)}' < chm13.rDNA.TideHunter.tsv > chm13.rDNA.TideHunter.fa
```

From the PRDM9 motifs found in the whole chromosomes, take those fully covering the repetitive unit:

 ```shell
rm chm13v2.PRDM9.rdna_units.bed
bedtools sort -i /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v1.1.rdna_units.unique.bed | sed 's/chr/chm13#chr/g' | while read f; do
  chr=$(echo $f | cut -f 1 -d ' ')
  start=$(echo $f | cut -f 2 -d ' ')
  end=$(echo $f | cut -f 3 -d ' ')
  echo "$chr:$start-$end"

  bedtools intersect \
    -a <(grep '^chm13#chr13\|^chm13#chr14\|^chm13#chr15\|^chm13#chr21\|^chm13#chr22' /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.bed | bedtools sort) \
    -b <(echo $f | tr ' ' '\t') -f 1.0 | \
    awk -v chr=$chr -v start=$start -v end=$end -v OFS='\t' '{print(chr":"start"-"end,$2-start,$3-start,$4,$5,$6,$7,$8,$9)}' >> chm13v2.PRDM9.rdna_units.bed
done

#ALTERNATIVE
#bedtools intersect \
#  -a <(grep '^chm13#chr13\|^chm13#chr14\|^chm13#chr15\|^chm13#chr21\|^chm13#chr22' /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.bed | bedtools sort) \
#  -b <(bedtools sort -i /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v1.1.rdna_units.unique.bed | sed 's/chr/chm13#chr/g') -f 1.0 | \
#  > chm13v2.PRDM9.rdna_units.bed
```

-----------------------------------------------------------------------------------------------------------
ALTERNATIVE APPROACH: it gives much more results, I think because of the p-value correction.
Having less pvalues (hits) to correct, more values get an adjusted p-value lower than the threshold (1.0E-4).

```shell
# Search the PRDM9 motifs (stronger adjusted p-value threshold because, as we are not looking at motifs genome-wide, we are going to correct less p-values)
RUN_FIMO=/home/guarracino/tools/meme-5.5.0/src/fimo
$RUN_FIMO --oc /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/ --verbosity 1 --thresh 1.0E-9 /lizardfs/guarracino/chromosome_communities/recombination_hotspots/PRDM9_motifs.human.txt /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v1.1.rdna_units.unique.fa

# Convert the output in BED format
grep '^#' fimo.tsv -v | sed '1d' | sed '/^$/d' | awk -v OFS='\t' '{print($2,$3-1,$4,$1,$6,$5,$7,$8,$9)}' | bedtools sort > chm13v1.1.rdna_units.unique.PRDM9.bed
```
-----------------------------------------------------------------------------------------------------------


Counts the number of hits in windows:

```shell
#TODO IF NEEDED: separate in Human1 hits, Human2 hits, ...
samtools faidx chm13v1.1.rdna_units.unique.fa

max_qvalue=1
window_size=1

rm chm13v2.PRDM9.rdna_units.w${window_size}.bed
bedtools sort -i /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v1.1.rdna_units.unique.bed | sed 's/chr/chm13#chr/g' | while read f; do
  chr=$(echo $f | cut -f 1 -d ' ')
  start=$(echo $f | cut -f 2 -d ' ')
  end=$(echo $f | cut -f 3 -d ' ')
  echo "$chr:$start-$end"

  bedtools intersect \
    -a <(bedtools makewindows -g <(cat chm13v1.1.rdna_units.unique.fa.fai | grep "$chr:$start-$end" | cut -f 1,2) -w $window_size) \
    -b <(grep $chr /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v2.PRDM9.rdna_units.bed | grep -P 'Human[1-7]*[0-9]\t' | awk -v max_qvalue=$max_qvalue '$8 <= max_qvalue') -c \
    >> chm13v2.PRDM9.rdna_units.w${window_size}.bed
done
```

Plot the number of hits in each window across on the units:

```shell
window_size=1

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_without_annotation.all_chromosomes.R \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v2.PRDM9.rdna_units.w${window_size}.bed \
  1 \
  'Position (bp)' \
  '' \
  45 \
  'Chromosome' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/rDNA/chm13v2.PRDM9.rdna_units.w${window_size}.pdf
```

## Adam Phillippy's sequence

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518
cd /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518

# Download https://www.ncbi.nlm.nih.gov/nuccore/KY962518
# mv /home/guarracino/Downloads/sequence.fasta KY962518.fasta
# scp /home/guarracino/Downloads/KY962518.fasta guarracino@octopus02:/lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518

# Search the PRDM9 motifs (stronger adjusted p-value threshold because, as we are not looking at motifs genome-wide, we are going to correct less p-values)
RUN_FIMO=/home/guarracino/tools/meme-5.5.0/src/fimo
$RUN_FIMO --oc /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/ --verbosity 1 --thresh 1.0E-8 /lizardfs/guarracino/chromosome_communities/recombination_hotspots/PRDM9_motifs.human.txt /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/KY962518.fasta

# Convert the output in BED format
grep '^#' fimo.tsv -v | sed '1d' | sed '/^$/d' | awk -v OFS='\t' '{print($2,$3-1,$4,$1,$6,$5,$7,$8,$9)}' | bedtools sort > KY962518.PRDM9.bed



samtools faidx KY962518.fasta
max_qvalue=1
window_size=1

rm KY962518.PRDM9.w${window_size}.bed
(seq 13 15; seq 21 22) | while read i; do
  echo $i

  bedtools intersect \
    -a <(bedtools makewindows -g <(cat KY962518.fasta.fai | cut -f 1,2) -w $window_size) \
    -b <(grep -P 'Human[1-7]*[0-9]\t' KY962518.PRDM9.bed | awk -v max_qvalue=$max_qvalue '$8 <= max_qvalue') -c \
    >> KY962518.PRDM9.w${window_size}.bed
done

window_size=1
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_without_annotation.all_chromosomes.R \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/KY962518.PRDM9.w${window_size}.bed \
  1 \
  'Position (bp)' \
  '' \
  45 \
  'Sequence' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/KY962518.PRDM9.w${window_size}.png



Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_BED.all_chromosomes.R \
  <(cat /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/KY962518.PRDM9.bed | awk -v OFS='\t' '{print($1,$2,$3,$4,$6,$8)}' ) \
  45 \
  'Sequence' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/KY962518/KY962518.PRDM9.png
```


### SST1

Obtain the repetitive unit of the array:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1
cd /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1

samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr13:12301367-12440010 > chm13.SST1.fa
# Put the SST1 region of chr14 in reverse complement.
# This help in getting with TideHunter all the repeat unit starting from the same point
samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr14:6960008-6988409 -i >> chm13.SST1.fa
samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr21:9375567-9453313 >> chm13.SST1.fa

RUN_TIDEHUNTER=/home/guarracino/tools/TideHunter-v1.5.4/bin/TideHunter

# https://github.com/yangao07/TideHunter#tabular-format
$RUN_TIDEHUNTER -f 2 chm13.SST1.fa -t 48 -k 13 > chm13.SST1.TideHunter.tsv


# For making dotplots with Gepard (to check that all the repeat units have the same start/end)
awk '{print(">"$1"_"$7"\n"$11)}' < chm13.SST1.TideHunter.tsv > chm13.SST1.TideHunter.fa
samtools faidx chm13.SST1.TideHunter.fa chm13#chr13:12301367-12440010_1409 > SST1.chr13.fa
samtools faidx chm13.SST1.TideHunter.fa chm13#chr14:6960008-6988409/rc_1407 > SST1.chr14rc.fa
samtools faidx chm13.SST1.TideHunter.fa chm13#chr21:9375567-9453313_1406 > SST1.chr21.fa
```

From the PRDM9 motifs found in the whole chromosomes, take those fully covering the repetitive unit:

 ```shell
cut -f 1,10 chm13.SST1.TideHunter.tsv | cut -f 1,2 -d ','
#chm13#chr13:12301367-12440010	967,2371  => chm13#chr13 12301367+967-1  12301367+2370-1
#chm13#chr14:6960008-6988409/rc	758,2163  => chm13#chr14 6988409-2163+1  6988409-758+1
#chm13#chr21:9375567-9453313	770,2174    => chm13#chr21 9375567+770-1   9375567+2174-1

# To check that the outputs are identical
diff \
  <(samtools faidx chm13.SST1.fa chm13#chr13:12301367-12440010:967-2370) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr13:12302333-12303736)
diff \
  <(samtools faidx chm13.SST1.fa chm13#chr14:6960008-6988409/rc:758-2163) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr14:6986247-6987652 -i)
diff \
  <(samtools faidx chm13.SST1.fa chm13#chr21:9375567-9453313:770-2174) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz chm13#chr21:9376336-9377740)
 

#chm13#chr13:12301367-12440010	967,2371  => chm13#chr13 12301367+967-1  12301367+2370-1
#chm13#chr14:6960008-6988409/rc	758,2163  => chm13#chr14 6988409-2163+1  6988409-758+1
#chm13#chr21:9375567-9453313	770,2174    => chm13#chr21 9375567+770-1   9375567+2174-1

rm chm13v2.PRDM9.SST1.bed
#(echo chm13#chr13:12302333-12303736; echo chm13#chr14:6986247-6987652; echo chm13#chr21:9376336-9377740) | while read f; do
(echo chm13#chr13:12301367-12440010; echo chm13#chr14:6960008-6988409; echo chm13#chr21:9375567-9453313) | while read f; do
  chr=$(echo $f | cut -f 1 -d ':')
  start=$(echo $f | cut -f 2 -d ':' | cut -f 1 -d '-')
  end=$(echo $f | cut -f 2 -d ':' | cut -f 2 -d '-')
  echo "$chr:$start-$end"

  start1=$((start-50))
  end1=$(($end+50))

  bedtools intersect \
    -a <(grep '^chm13#chr13\|^chm13#chr14\|^chm13#chr15\|^chm13#chr21\|^chm13#chr22' /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.bed | bedtools sort) \
    -b <(echo $chr $start $end | tr ' ' '\t') -f 1.0 | \
    awk -v chr=$chr -v start=$start -v end=$end -v OFS='\t' '{print(chr":"start"-"end,$2-start,$3-start,$4,$5,$6,$7,$8,$9)}' >> chm13v2.PRDM9.SST1.bed
done
```

-----------------------------------------------------------------------------------------------------------
ALTERNATIVE APPROACH: it gives much more results, I think because of the p-value correction.
Having less pvalues (hits) to correct, more values get an adjusted p-value lower than the threshold (1.0E-4).

```shell
#Search the PRDM9 motifs
RUN_FIMO=/home/guarracino/tools/meme-5.5.0/src/fimo
$RUN_FIMO --oc /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/ --verbosity 1 --thresh 1.0E-8 /lizardfs/guarracino/chromosome_communities/recombination_hotspots/PRDM9_motifs.human.txt /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13.SST1.TideHunter.fa

#Convert the output in BED format
samtools faidx chm13.SST1.TideHunter.fa
cat chm13.SST1.TideHunter.fa.fai | awk -v OFS='\t' '{print($1,"0",$2,"SST1","+","0")}' > chm13.SST1.TideHunter.PRDM9.bed
grep '^Human' /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/fimo.tsv | awk -v OFS='\t' '{print($2,$3-1,$4,$1,$5,$8)}' >> chm13.SST1.TideHunter.PRDM9.bed
```
-----------------------------------------------------------------------------------------------------------

Counts the number of hits in each base pair:

```shell
samtools faidx chm13.SST1.fa

max_qvalue=1
window_size=1

rm chm13v2.PRDM9.SST1.w${window_size}.bed
cat chm13v2.PRDM9.SST1.bed | cut -f 1 | sort | uniq | while read f; do
  chr=$(echo $f | cut -f 1 -d ':')
  start=$(echo $f | cut -f 2 -d ':' | cut -f 1 -d '-')
  end=$(echo $f | cut -f 2 -d ':' | cut -f 2 -d '-')
  echo "$chr:$start-$end"

  bedtools intersect \
    -a <(bedtools makewindows -g <(sed 's/\/rc//g' chm13.SST1.fa.fai | grep "$chr:$start-$end" | cut -f 1,2) -w $window_size) \
    -b <(grep $chr /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13v2.PRDM9.SST1.bed | grep -P 'Human[1-7]*[0-9]\t' | awk -v max_qvalue=$max_qvalue '$8 <= max_qvalue') -c \
    >> chm13v2.PRDM9.SST1.w${window_size}.bed
done
```

Plot the number of hits in each window across on the units:

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_without_annotation.all_chromosomes.R \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13v2.PRDM9.SST1.w${window_size}.bed \
  1 \
  'Position (bp)' \
  '' \
  35 \
  'Chromosome' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13v2.PRDM9.SST1.w${window_size}.pdf
```

Show where the PRDM9 hits are on the SST1 units:

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_PRDM9_hits_BED.all_chromosomes.R \
  <(cat /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13v2.PRDM9.SST1.bed | awk -v OFS='\t' '{print($1,$2,$3,$4,$6,$8)}' ) \
  35 \
  'Chromosome' \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/repeat_unit/SST1/chm13v2.PRDM9.SST1.pdf
```


NOT USED
Alternative plots of the hits across in the untangle space:

```shell
# List of contigs with at least one match
###cut -f 1 fimo.bed | sort | uniq > fimo.contigs.txt

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt
flipped=false
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

PATH_FIMO_BED=/lizardfs/guarracino/chromosome_communities/recombination_hotspots/tmp/fimo.bed

rm xyz.tsv z.tsv # Cleaning

for e in 50000; do
  for m in 1000 ; do
    cat $path_targets_txt | while read ref; do 
      echo "-e $e -m $m $ref"

      # Plot additional tracks for the annotations
      path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz

      zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
        SAMPLE=$( echo $CONTIG | cut -f 1 -d '#');
        
        N=$(grep $CONTIG $PATH_FIMO_BED | wc -l)
        if (( $N > 0 )); then
          echo $CONTIG "--->" $SAMPLE
           
          zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{ if($10 <= 1 && $15 <= 1) {print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16)} }' > x.$CONTIG.bed
          
          # To get only Unreliable intervals
          grep $CONTIG $PATH_FIMO_BED | cut -f 1,2,3,4 > y.$CONTIG.bed
            
          # wao: write the original A and B entries plus the number of base pairs of overlap between the two features.
          bedtools intersect -a x.$CONTIG.bed -b y.$CONTIG.bed -wo >> z.tsv
          
          rm x.$CONTIG.bed y.$CONTIG.bed
        fi
      done
        
      cat \
        <( zcat $path_grounded_pq_touching_tsv_gz | sed '1d' ) \
        <( python3 /lizardfs/guarracino/chromosome_communities/scripts/get_annotation_track.py z.tsv ) | tr ' ' '\t' >> xyz.tsv
      rm z.tsv
    done
    
    ref=$( head -n 1 $path_targets_txt ) # The first (all the same for getting the header line)
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz
    cat \
      <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
       xyz.tsv | pigz -c > e$e.m$m.annot.tsv.gz
    rm xyz.tsv
  done
done


zcat e$e.m$m.annot.tsv.gz | sed 's/Human[0-9]\{1,\}/Err/' | pigz -c > e$e.m$m.annot.sed.tsv.gz

for e in 50000; do
  for m in 1000 ; do
    (seq 13 15; seq 21 22) | while read i; do
      echo "-e $e -m $m chr$i"
      
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_and_flagger.R \
        e$e.m$m.annot.sed.tsv.gz \
        0 25000000 \
        90 \
        0.9 1.0 \
        1 1 \
        $i \
        0.9 \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        e$e.m$m.annot.chr$i.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      e$e.m$m.annot.chr*.pdf e$e.m$m.annot.chrACRO.pdf
  done
done
```


# NOT USED: new figure 5 centered in the SST1 region

```shell
#chm13#chr13	12301367	12440010	SST1#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr13.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  11301367 13440010 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  13 \
  0.9 \
  <(echo 'chm13#chr13' 'grch38#chr13' 'HG002#MAT#chr13.prox' 'HG002#PAT#chr13.prox' 'HG01361#2#JAGYYW010000010.1' 'HG01978#1#JAGYVS010000056.1' 'HG02486#1#JAGYVM010000043.1' 'HG03540#2#JAGYVX010000153.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr13_SST1_1Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr13.SST1.1Mbps.n1.nref1.pdf
  
#chm13#chr14  6960008 6988409 SST_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr14.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  5960008 7988409 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  14 \
  0.9 \
  <(echo 'chm13#chr14' 'grch38#chr14' 'HG002#MAT#chr14.prox' 'HG002#PAT#chr14.prox' 'HG00735#1#JAHBCH010000039.1' 'HG00741#2#JAHALX010000038.1' 'HG01978#1#JAGYVS010000055.1' 'HG03453#2#JAGYVV010000008.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr14_SST1_1Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr14.SST1.1Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  8375567 10453313 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  21 \
  0.9 \
  <(echo 'chm13#chr21' 'grch38#chr21' 'HG002#MAT#chr21.prox' 'HG002#PAT#chr21.prox' 'HG00735#2#JAHBCG010000066.1' 'HG02886#1#JAHAOU010000106.1' 'NA18906#1#JAHEOO010000072.1' 'NA19240#2#JAHEOL010000065.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr21_SST1_1Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr21.SST1.1Mbps.n1.nref1.pdf


#chm13#chr13	12301367	12440010	SST1#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr13.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  9301367 15440010 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  13 \
  0.9 \
  <(echo 'chm13#chr13' 'grch38#chr13' 'HG002#MAT#chr13.prox' 'HG002#PAT#chr13.prox' 'HG01361#2#JAGYYW010000010.1' 'HG01978#1#JAGYVS010000056.1' 'HG02486#1#JAGYVM010000043.1' 'HG03540#2#JAGYVX010000153.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr13_SST1_3Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr13.SST1.3Mbps.n1.nref1.pdf
  
#chm13#chr14  6960008 6988409 SST_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr14.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  3960008 9988409 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  14 \
  0.9 \
  <(echo 'chm13#chr14' 'grch38#chr14' 'HG002#MAT#chr14.prox' 'HG002#PAT#chr14.prox' 'HG00735#1#JAHBCH010000039.1' 'HG00741#2#JAHALX010000038.1' 'HG01978#1#JAGYVS010000055.1' 'HG03453#2#JAGYVV010000008.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr14_SST1_3Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr14.SST1.3Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  6375567 12453313 \
  90 0.9 7.5 \
  1.0 \
  5 1 \
  21 \
  0.9 \
  <(echo 'chm13#chr21' 'grch38#chr21' 'HG002#MAT#chr21.prox' 'HG002#PAT#chr21.prox' 'HG00735#2#JAHBCG010000066.1' 'HG02886#1#JAHAOU010000106.1' 'NA18906#1#JAHEOO010000072.1' 'NA19240#2#JAHEOL010000065.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr21_SST1_3Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr21.SST1.3Mbps.n1.nref1.pdf
```

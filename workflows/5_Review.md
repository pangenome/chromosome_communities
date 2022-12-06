#### Plot FLAGGER's annotations (for pq-contig quality control)

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
       xyz.tsv | pigz -c > e$m.m$m.annot.tsv.gz
    rm xyz.tsv
  done
done

# Supplementary Figures
for e in 50000; do
  for m in 1000 ; do
    (seq 13 15; seq 21 22) | while read i; do
      echo "-e $e -m $m chr$i"
      
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_and_flagger.R \
        e$m.m$m.annot.tsv.gz \
        0 25000000 \
        90 \
        0.9 1.0 \
        1 1 \
        $i \
        0.9 \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        e$m.m$m.annot.chr$i.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      e$m.m$m.annot.chr*.pdf e$m.m$m.annot.chrACRO.pdf
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
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}

ref=chr14
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_14_39_SST_Composite" \
  6960008 6988409 \
  1000000 \
  5 \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}

ref=chr21
path_untangle_single_chr_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#${ref}.e50000.m1000.j0.n100.bed.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_SST1_dotplot.R \
  $path_untangle_single_chr_tsv_gz \
  "censat_21_45_SST1_Composite" \
  9175567 9653313 \
  1000000 \
  5 \
  <((
    cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.vs.refs.partitions.tsv | grep 'chr13$\|chr14$\|chr21$'; \
    zgrep 'HG002#.AT' pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | cut -f 1 | grep 'chr13\|chr14\|chr21'; \
    echo "chm13#chr13"; echo "chm13#chr14"; echo "chm13#chr21")) \
  /lizardfs/guarracino/chromosome_communities/untangle/sst1_region_dotplots/query_vs_${ref}
```


# New figure 5 centered in the SST1 region

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
  3 1 \
  14 \
  0.9 \
  <(echo 'chm13#chr14' 'grch38#chr14' 'HG002#MAT#chr14.prox' 'HG002#PAT#chr14.prox' 'HG00735#1#JAHBCH010000039.1' 'HG00741#2#JAHALX010000038.1' 'HG01978#1#JAGYVS010000055.1' 'HG02630#1#JAHAOQ010000067.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr14_SST1_1Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr14.SST1.1Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  8375567 10453313 \
  90 0.9 7.5 \
  1.0 \
  3 1 \
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
  
#chm13#chr14  3960008 6988409 SST_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr14.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  5960008 9988409 \
  90 0.9 7.5 \
  1.0 \
  3 1 \
  14 \
  0.9 \
  <(echo 'chm13#chr14' 'grch38#chr14' 'HG002#MAT#chr14.prox' 'HG002#PAT#chr14.prox' 'HG00735#1#JAHBCH010000039.1' 'HG00741#2#JAHALX010000038.1' 'HG01978#1#JAGYVS010000055.1' 'HG02630#1#JAHAOQ010000067.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr14_SST1_3bps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr14.SST1.3Mbps.n1.nref1.pdf

#chm13#chr21  9375567 9453313   SST1_Composite#222222
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chr21.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_SST1.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  6375567 12453313 \
  90 0.9 7.5 \
  1.0 \
  3 1 \
  21 \
  0.9 \
  <(echo 'chm13#chr21' 'grch38#chr21' 'HG002#MAT#chr21.prox' 'HG002#PAT#chr21.prox' 'HG00735#2#JAHBCG010000066.1' 'HG02886#1#JAHAOU010000106.1' 'NA18906#1#JAHEOO010000072.1' 'NA19240#2#JAHEOL010000065.1' | tr ' ' '\n') \
  /lizardfs/guarracino/chromosome_communities/data/annotation/genome_browser_chr21_SST1_3Mbps_CenSatAnnDense.png \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/Figure5.chr21.SST1.3Mbps.n1.nref1.pdf
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


# Robertsonian translocation

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

mkdir -p /lizardfs/guarracino/chromosome_communities/robertsonian_translocation
cd /lizardfs/guarracino/chromosome_communities/robertsonian_translocation

# Download the sequence from https://www.ncbi.nlm.nih.gov/nuccore/CR382332

$RUN_WFMASH /lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz CR382332.fasta.gz -t 48 \
  -p 98 -s 1k -n 10 -N -m > CR382332.vs.CHM13.p98.s1k.n10.N.paf
  

cat \
  <(cut -f 1,2,3 SupplementaryTable4.chrACRO_29-Jul-22_PHRs.bed | awk -v OFS='\t' '{print($0,"100.0","PHRs")}') \
  <(cut -f 6,8,9,13 CR382332.vs.CHM13.p97.s1k.n10.N.paf | sed 's/id:f://' | awk -v OFS='\t' '{print($0,"CR382332")}') \
  > PHRs+CR382332.bed
```


# Recombination hotspots

Prepare the tool:

```shell
# Download and install the MEME Suite
wget -c https://meme-suite.org/meme/meme-software/5.5.0/meme-5.5.0.tar.gz
tar -xf meme-5.5.0.tar.gz 
cd ~/tools/meme-5.5.0
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
```

Search PRDM9 motifs (downloaded from [here](https://pubmed.ncbi.nlm.nih.gov/29072575/)):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination_hotspots
cd /lizardfs/guarracino/chromosome_communities/recombination_hotspots

# Download the PWMs for all PRDM9 motifs
wget -c https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5705219/bin/elife-28383-fig1-data2.txt

# Remove non-human motifs
n=$(grep 'MOTIF Chimp1' elife-28383-fig1-data2.txt -n | cut -f 1 -d ':')
n=$((n-2))
head -n $n elife-28383-fig1-data2.txt > PRDM9_motifs.human.txt

zcat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz \
  > /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa

RUN_FIMO=/home/guarracino/tools/meme-5.5.0/src/fimo
sbatch -p workers -c 48 --job-name meme-PRDM9 --wrap "hostname; cd /scratch && $RUN_FIMO --oc /lizardfs/guarracino/chromosome_communities/recombination_hotspots/ --verbosity 1 --thresh 1.0E-4 /lizardfs/guarracino/chromosome_communities/recombination_hotspots/PRDM9_motifs.human.txt /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa"


```

Convert the output in BED format:

```shell
# Remove the last lines, remove the header line, remove the last empty line, prepare the columns, and sort the BED file
grep '^#' fimo.tsv -v | sed '1d' | sed '/^$/d' | awk -v OFS='\t' '{print($2,$3,$4,$1,$6,$5,$7,$8,$9)}' | bedtools sort > fimo.bedE
```

Plots:

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
       xyz.tsv | pigz -c > e$m.m$m.annot.tsv.gz
    rm xyz.tsv
  done
done


# Supplementary Figures
zcat e$m.m$m.annot.tsv.gz | sed 's/Human[0-9]\{1,\}/Err/' | pigz -c > e$m.m$m.annot.sed.tsv.gz

for e in 50000; do
  for m in 1000 ; do
    (seq 13 15; seq 21 22) | while read i; do
      echo "-e $e -m $m chr$i"
      
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation_and_flagger.R \
        e$m.m$m.annot.sed.tsv.gz \
        0 25000000 \
        90 \
        0.9 1.0 \
        1 1 \
        $i \
        0.9 \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        e$m.m$m.annot.chr$i.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      e$m.m$m.annot.chr*.pdf e$m.m$m.annot.chrACRO.pdf
  done
done


bedtools intersect -a <(bedtools makewindows -g <(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.fai | grep 'HG002#MAT#chr13' | cut -f 1,2) -w 20000) -b <(grep HG002#MAT#chr13 *bed) -c | column -t | less -S
```




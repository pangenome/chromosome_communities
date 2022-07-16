# Untangling

## Sex chromosome

Use HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

# Get HG002 chromosome X and Y

PATH_CHM13_FA_GZ=/lizardfs/erikg/human/chm13v2.0.fa.gz
PATH_GRCH38_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz

PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.partitions.tsv | awk -v OFS='\t' '($2 == "chrX" || $2 == "chrY") {print $1}' | sort > chrXY_contigs.txt

# HiFi-based assemblies + hg01978 verkko + hg002-bakeoff + hg002 verkko
cat \
  <( samtools faidx $PATH_CHM13_FA_GZ $(grep 'chrX\|chrY' $PATH_CHM13_FA_GZ.fai | cut -f 1) | sed 's/^>/>chm13#/g' ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chrX\|chrY' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <(samtools faidx $PATH_HPRCY1_FA_GZ $(comm -12 <(cut -f 1 "$PATH_HPRCY1_FA_GZ".fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa.fai | sort) <(cat chrXY_contigs.txt))) | \
  bgzip -@ 48 -c > chrSEX+refs.fa.gz
samtools faidx chrSEX+refs.fa.gz
```

### Pangenome building

Apply `pggb` on the chromosome-partitioned HPRC dataset:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

RUN_PGGB=/home/guarracino/tools/pggb/pggb-a4a6668d9ece42c80ce69dc354f0cb59a849286f

num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz.fai | sort | uniq | wc -l)
sbatch -p highmem -c 48 --job-name sexpggb --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz -o chrSEX+refs.s50k.l250k.p98.n${num_of_haplotypes} -t 48 -s 50k -l 250k -p 98 -n ${num_of_haplotypes} -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrSEX+refs.s50k.l250k.p98.n${num_of_haplotypes} /lizardfs/guarracino/chromosome_communities/graphs";
```

### Untangling

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/chm13.SEX.target_paths.txt
path_fasta_fai=/lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz.fai
grep chm13 $path_fasta_fai | cut -f 1 > $path_targets_txt


path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
prefix=$(basename $path_input_og .og)

RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-e2de6cbca0169b0720dca0c668743399305e92bd
```

```shell
# All references and emit cut points
for e in 50000; do
  for m in 1000 2000 5000; do
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.bed.gz
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.cut_points.txt
    
    if [[ ! -s ${path_cut_points_txt} ]]; then
      echo "-e $e -m $m"
      sbatch -p workers -c 16 --job-name sexuntangle --wrap "\time -v $RUN_ODGI untangle -t 16 -P -i $path_input_og -R $path_targets_txt -e $e -m $m --cut-points-output $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_bed_gz"
    fi;
  done
done
  
# Single reference by using the same cut points
for e in 50000; do
  for m in 1000 2000 5000; do
    echo "-e $e -m $m"
      
    cat $path_targets_txt | while read ref; do
      echo $ref
      
      path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
      if [[ ! -s ${path_ref_bed_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.cut_points.txt
        
        sbatch -p workers -c 8 --job-name sexuntangle --wrap "\time -v $RUN_ODGI untangle -t 8 -P -i $path_input_og -r $ref -e $e -m $m --cut-points-input $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_ref_bed_gz"
      fi
    done
  done
done


# Grounding
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded

for e in 50000; do
  for m in 1000; do
    cat $path_targets_txt | while read ref; do            
      path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.tsv.gz
            
      if [[ ! -s ${path_grounded_tsv_gz} ]]; then
        echo "-e $e -m $m $ref"

        # Grounding
        ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage nth.best ref ref.begin ref.end ref.jaccard ref.nth.best | tr ' ' '\t'
          join \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.bed.gz | awk -v j=0 -v n=5 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz | awk -v j=0 -v n=10 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
          tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -10,14-17,20 | sort -k 1,3 -k 7,7nr -k 10,10n -k 14,14nr -k 15,15n ) | tr ' ' '\t' | pigz -c > x.tsv.gz
                    
        # Contigs overlapping (or close at least 100kbps to) a PAR
        ref_chr=$(echo $ref | sed 's/chm13#//')
        bedtools intersect \
          -a <(zcat x.tsv.gz | awk -v OFS="\t" '{print $11,$12,$13,$1, "", "+"}' | sed '1d' | bedtools sort) \
          -b <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed |\
            bedtools sort |\
            bedtools slop -b 100000 -g /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.sizes |\
            awk -v OFS='\t' -v ref=$ref '{print(ref,$2,$3)}') | \
          #awk '$3-$2+1>=100000' | \
          cut -f 4 | \
          #Remove references to avoid grepping everything later (with zgrep -f)
          grep -v chr |\
          sort | uniq > $ref.tmp.txt
           
        # Add grounded.target column, re-add the references, and add annotation
        cat \
          <(zcat x.tsv.gz | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
          <(zgrep x.tsv.gz -f $ref.tmp.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
          <(zcat x.tsv.gz | awk -v OFS='\t' -v ref=$ref '$1 ~ /chr/ { print $0, ref}' ) \
          <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed | \
            awk -v OFS='\t' -v ref=$ref '{print "chm13#"$1,".",".",ref,".",".","1","+","1","1",ref,$2,$3,"1","1",ref}') |\
          pigz -c -9 > $path_grounded_tsv_gz
        rm x.tsv.gz
      fi;
    done
  done
done
```


Plot:

```shell
#PAR1/2/3
# https://link.springer.com/article/10.1007/s10142-013-0323-6/figures/1

for e in 50000  ; do
  for m in 1000  ; do
    for refn in 1 ; do
      (echo X; echo Y) | while read i; do      
        echo "-e $e -m $m -refn $refn chr$i"
    
        path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr$i.e$e.m$m.grounded.tsv.gz
        PREFIX=$(basename $path_grounded_tsv_gz .tsv.gz);
        
        if [[ $i == "X" ]]; then
            # chrX#PAR1:0-2396333
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              0 2446333 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.PAR1.pdf
              
            # chrX#PAR2:154012988-154349815
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              153962988 154399815 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.PAR2.pdf
              
            # chrX#PAR2:87723198-91647350
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              87673198 91697350 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.XTR.pdf
        
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              0 154259566 \
              200 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.pdf
        else
            # chrY#PAR1:0-2458348
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              0 2508348 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.PAR1.pdf
              
            # chrY#PAR2:62122794-62460029
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              62072794 62510029 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.PAR2.pdf
              
            # chrY#XTR1:2727073-59145619
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              2677073 59195619 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.XTR1.pdf
              
            # chrY#XTR2:6200825-6400875
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              6150825 6450875 \
              90 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.XTR2.pdf
              
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_tsv_gz \
              0 62460029 \
              80 0.01 \
              0 \
              1 $refn \
              $i \
              <(zcat $path_grounded_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n1.nref${refn}.pdf
        fi
      done
      
#      # Merge chromosomes's PDF files
#      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
#        /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.n1.nref${refn}.pdf \
#        /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.pq_touching.reliable.n1.nref${refn}.merged.pdf
#      rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.n1.nref${refn}.pdf
    done
  done
done
```


[//]: # (## Variant calling)
[//]: # ()
[//]: # (Call variants in a haploid setting:)
[//]: # ()
[//]: # (```shell)
[//]: # (path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.gfa)
[//]: # (path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.chm13.vcf.gz)
[//]: # (#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.grch38.vcf.gz)
[//]: # (path_hg002_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.hg002.vcf.gz)
[//]: # (sbatch -p workers -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz)
[//]: # (#sbatch -p workers -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz)
[//]: # (sbatch -p workers -c 48 --job-name vghg002 --wrap '\time -v vg deconstruct -P HG002 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_hg002_vcf_gz' && tabix '$path_hg002_vcf_gz)
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # (# In the VCF there are variants with all samples having missing genotype!)
[//]: # (#num_samples=`bcftools query -l $PATH_VCF_GZ | wc -l`)
[//]: # (#num_miss_gen=$&#40;echo $num_samples - 1 | bc&#41;)
[//]: # (#--max-missing-count $num_miss_gen \)
[//]: # (```)

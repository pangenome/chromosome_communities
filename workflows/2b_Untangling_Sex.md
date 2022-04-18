
## Sex chromosome

Use HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

# Get HG002 chromosome X and Y
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
gunzip HG002.mat.cur.20211005.fasta.gz
samtools faidx HG002.mat.cur.20211005.fasta
samtools faidx HG002.mat.cur.20211005.fasta SX | sed 's/SX/HG002#MAT#chrX/g' |\
  bgzip -@ 48 -c > hg002.chrX.fa.gz

~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002#PAT#' /lizardfs/erikg/T2T/liftover/split_chrY/chm13v2.0_chrY.fasta |\
  bgzip -@ 48 -c > hg002.chrY.fa.gz
```

Put HG002's chrX and chrY with the partitioned chrXs and chrYs.

```shell
# Prepare sequence order, with all references on the top
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai | cut -f 1 > sequence_order.txt
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai -v | cut -f 1 >> sequence_order.txt

cat \
  <(zcat hg002.chrX.fa.gz) \
  <(zcat hg002.chrY.fa.gz) \
  <(samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa $(cat sequence_order.txt)) |\
   bgzip -@ 48 > chrS.pan+HG002chrXY.fa.gz
samtools faidx chrS.pan+HG002chrXY.fa.gz

rm sequence_order.txt
```

### Pangenome building

Apply `pggb` on the chromosome-partitioned HPRC dataset:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz.fai | sort | uniq | wc -l)
sbatch -p highmem -c 48 --job-name sexpggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz -o chrS.pan+HG002chrXY.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrS.pan+HG002chrXY.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';
```

### Untangling

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.og
prefix=$(basename $path_input_og .og)

run_odgi=/home/guarracino/tools/odgi/bin/odgi-694948ccf31e7b565449cc056668e9dcc8cc0a3e
path_fasta_fai=/lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz.fai

# All references and emit cut points
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  grep $refpattern $path_fasta_fai | cut -f 1 > $path_targets_txt
    
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      echo "-e $e -m $m"
    
      path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.bed.gz
      path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.cut_points.txt
    
      if [[ ! -s ${path_cut_points_txt} ]]; then
        sbatch -p workers -c 24 --job-name sexuntangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -R '$path_targets_txt' -e '$e' -m '$m' --cut-points-output '$path_cut_points_txt' -j 0 -n 1 | pigz -c > '$path_bed_gz';'
      fi;
    done
  done
done

# Single reference by using the same cut points
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      echo "-e $e -m $m"
        
      cat $path_targets_txt | while read ref; do
        echo $ref
            
        path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
        if [[ ! -s ${path_ref_bed_gz} ]]; then
          path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.cut_points.txt
                    
          sbatch -p workers -c 24 --job-name sexuntangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -r '$ref' -e '$e' -m '$m' --cut-points-input '$path_cut_points_txt' -j 0 -n 100 | pigz -c > '$path_ref_bed_gz';'
        fi;
      done
    done
  done
done

# Grounding
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt

  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
          echo "-e $e -m $m -j $j -n $n"
                
          cat $path_targets_txt | while read ref; do
            echo -e "\t"$ref

            path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
            if [[ ! -s ${path_grounded_tsv_gz} ]]; then
              path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.bed.gz
              path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
              
              # Grounding
              ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage ref ref.begin ref.end | tr ' ' '\t'
                join \
                  <(zcat $path_bed_gz | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) \
                  <(zcat $path_ref_bed_gz | awk -v j=$j -v n=$n '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
                tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -9,14-16 ) | tr ' ' '\t' | pigz -c > x.tsv.gz
              
              # Contigs overlapping (or close at least 100kbps to) a PAR
              # Note that chrX PAR is from chm13, not HG002
              ref_chr=$(echo $ref | rev | cut -f 1 -d '#' | rev)
              bedtools intersect \
                -a <(zcat x.tsv.gz | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
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
                <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.approximate.bed | \
                  awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",ref,".",".",".",".",".",".",$2,$3, ref}') |\
                pigz -c > $path_grounded_tsv_gz

              rm x.tsv.gz
            fi;
          done
        done
      done
    done
  done
done
```

Plot:

```shell
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt

  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
          echo "-e $e -m $m -j $j -n $n"
    
          path_grounded_all_references_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
          if [[ ! -s ${path_grounded_all_references_tsv_gz} ]]; then
            # Merge single reference results
            cat \
              <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern*.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz | head -n 1) \
              <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern*.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz | grep query -v) |\
              pigz -c > x.tsv.gz
            # Rename after to avoid getting itself with the previous '*' expansion
            mv x.tsv.gz $path_grounded_all_references_tsv_gz
          fi;

          Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R $path_grounded_all_references_tsv_gz "-e $e -m $m -j $j -n $n" 0 155000000 120 800
        done
      done
    done
  done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

# Merge
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        echo "-e $e -m $m -j $j"
            
        j_str=$(echo $j | sed 's/\.//g')
        PDFs=$((seq 1 5; seq 10 10 50) | while read n; do \
          echo /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz.pdf
        done | tr '\n' ' ')
        #echo $PDFs

        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite $PDFs /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.merged.grounded.tsv.gz.pdf
      done
    done
  done
done

#PAR1/2/3
# https://link.springer.com/article/10.1007/s10142-013-0323-6/figures/1
```


## Variant calling

Call variants in a haploid setting:

```shell
path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.gfa
path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.chm13.vcf.gz
#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.grch38.vcf.gz
path_hg002_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.hg002.vcf.gz
sbatch -p workers -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz
#sbatch -p workers -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz
sbatch -p workers -c 48 --job-name vghg002 --wrap '\time -v vg deconstruct -P HG002 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_hg002_vcf_gz' && tabix '$path_hg002_vcf_gz





# In the VCF there are variants with all samples having missing genotype!
#num_samples=`bcftools query -l $PATH_VCF_GZ | wc -l`
#num_miss_gen=$(echo $num_samples - 1 | bc)
#--max-missing-count $num_miss_gen \
```

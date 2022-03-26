# History

## Preparation

### Variables

```shell
THREADS=48
RUN_MINIMAP2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2

# HACKED bamCaller.py to emit haploid VCF files --> if altAllele != "." and re.match("^[01]", genotypes):
RUN_HACKED_BAM_CALLER=/home/guarracino/tools/msmc-tools/bamCaller.py

# Using a hacked version, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#phased = True
#geno = "%s|%s" % (geno[0], geno[0])
RUN_HACKED_GENERATE_MULTIHETSEP=/home/guarracino/tools/msmc-tools/generate_multihetsep.py
```

### Steps

Split `chrS` by haplotype:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/haplotypes
cd /lizardfs/guarracino/chromosome_communities/haplotypes

cat /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai | cut -f 1,2 -d '#' | grep chr -v | sort | uniq > haplotypes.txt

( seq 13 15; seq 21 22 ) | while read i; do
  mkdir -p chr$i
  
  PATH_CHR_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/parts/chr$i.pan.fa
  cat haplotypes.txt | while read HAPLO; do
    echo "chr$i" $HAPLO
  
    samtools faidx $PATH_CHR_FASTA $(grep "^$HAPLO" $PATH_CHR_FASTA.fai | cut -f 1) |
      bgzip -@ THREADS -c > chr$i/chr$i.$HAPLO.fa.gz
    samtools faidx chr$i/chr$i.$HAPLO.fa.gz
  done
done
```

Map contigs against acrocentric chromosomes:

```shell
( seq 13 15; seq 21 22 ) | while read i; do
  cat haplotypes.txt | while read HAPLO; do
    echo "chr$i" $HAPLO
    
    ( seq 13 15; seq 21 22 ) | while read j; do
      PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz
      
      #$run_mp2 $PATH_REF_CHR_FASTA chr$i/chr$i.$HAPLO.fa.gz -t $THREADS -x asm20 -R "@RG\tID:$HAPLO\tSM:$HAPLO" -ca | samtools sort > chr$i/chr$i.$HAPLO.on.chr$j.bam
      sbatch -p workers -c 6 --job-name "$i-on-$j" --wrap 'hostname; \time -v '$RUN_MINIMAP2' '$PATH_REF_CHR_FASTA' /lizardfs/guarracino/chromosome_communities/haplotypes/chr'$i'/chr'$i'.'$HAPLO'.fa.gz -t 6 -x asm20 -R "@RG\tID:'$HAPLO'\tSM:'$HAPLO'" -ca | samtools sort > /lizardfs/guarracino/chromosome_communities/haplotypes/chr'$i'/chr'$i'.'$HAPLO'.on.chr'$j'.bam'
    done
  done
done
```

Generating VCF and Mask files:

```shell
( seq 13 15; seq 21 22 ) | while read i; do  
  cat haplotypes.txt | grep 'HG00438\|HG00621\|HG00673\|HG00733' | while read HAPLO; do
    echo "chr$i" $HAPLO
  
    ( seq 13 15; seq 21 22 ) | while read j; do
      PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz
      
      bcftools mpileup -B -q 0 -Q 0 -C 50 -f $PATH_REF_CHR_FASTA chr$i/chr$i.$HAPLO.on.chr$j.bam |\
        bcftools call -c --skip-variants indels --ploidy 1 --threads $THREADS |\
        $RUN_HACKED_BAM_CALLER 1 chr$i/chr$i.$HAPLO.on.chr$j.mask.bed.gz --minMapQ 1 |\
        bgzip -@ $THREADS -c > chr$i/chr$i.$HAPLO.on.chr$j.vcf.gz
    done
  done
done
```

Combining into one file multiple individuals from all possible pairs of populations:

```shell
( seq 13 15; seq 21 22 ) | while read r; do 
  ( seq 13 15; seq 21 22 ) | while read a; do  
    ( seq 13 15; seq 21 22 ) | while read b; do
      if (( $a < $b )); then
        echo "chr$a vs chr$b on chr$r"
        
        MASKsA=$( ls chr$a/chr$a.*.on.chr$r.mask.bed.gz | sed -z 's/\n/ --mask /g' | sed 's/ --mask $//g')
        VCFsA=$( ls chr$a/chr$a.*.on.chr$r.vcf.gz )
        
        MASKsB=$( ls chr$b/chr$b.*.on.chr$r.mask.bed.gz | sed -z 's/\n/ --mask /g' | sed 's/ --mask $//g')
        VCFsB=$( ls chr$b/chr$b.*.on.chr$r.vcf.gz )
        
        python3 ~/tools/msmc-tools/generate_multihetsep.py --mask $MASKsA --mask $MASKsB $VCFsA $VCFsB \
          > chr$a.vs.chr$b.on.chr$r.phased.hetsep
      fi
    done
  done
done


```
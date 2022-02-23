# python3 sample_contigs.py /lizardfs/guarracino/chromosome_communities/pq_contigs/chr13.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai /lizardfs/guarracino/chromosome_communities/pq_contigs/chr14.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai /lizardfs/guarracino/chromosome_communities/pq_contigs/chr15.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai /lizardfs/guarracino/chromosome_communities/pq_contigs/chr21.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai /lizardfs/guarracino/chromosome_communities/pq_contigs/chr22.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai 10

import sys

path_contigs_chr13_txt = sys.argv[1]
path_contigs_chr14_txt = sys.argv[2]
path_contigs_chr15_txt = sys.argv[3]
path_contigs_chr21_txt = sys.argv[4]
path_contigs_chr22_txt = sys.argv[5]

num_iteration=int(sys.argv[6])


contigs_list = [[], [], [], [], []]

for i, path_contigs_chrX_txt in enumerate([path_contigs_chr13_txt, path_contigs_chr14_txt, path_contigs_chr15_txt, path_contigs_chr21_txt, path_contigs_chr22_txt]):
	with open(path_contigs_chrX_txt) as f:
		for line in f:
			contig = line.split('\t')[0]
			if 'chr' not in contig:
				contigs_list[i].append(contig)



import random

for i in range(num_iteration):
	random_contig_list = []
	for x in contigs_list:
		random_contig_list.append(random.sample(x, 1)[0])
	print(','.join(random_contig_list))

#!/usr/bin/env python

import gzip
import sys


path_support_tsv_gz=sys.argv[1]


# Read untangling information
ground_2_group_2_query_2_pieces_dict = {}

last_start = None
last_end = None
last_line = ''

current_start = None
current_end = None
current_line = ''

print('\t'.join([
    'ground.target', 'start', 'end',
    'num.contigs.supporting.chr13',
    'num.contigs.supporting.chr14',
    'num.contigs.supporting.chr15',
    'num.contigs.supporting.chr21',
    'num.contigs.supporting.chr22'
]))


#ground.target, start, end, num_contigs_supporting_chr13, num_contigs_supporting_chr14, num_contigs_supporting_chr15, num_contigs_supporting_chr21, num_contigs_supporting_chr22
with gzip.open(path_support_tsv_gz, "rt") as f:
    f.readline()

    for line in f:
        ground_target, start, end, num_contigs_supporting_chr13, num_contigs_supporting_chr14, num_contigs_supporting_chr15, num_contigs_supporting_chr21, num_contigs_supporting_chr22 = line.strip().split('\t')

        current_start = start
        current_end = end
        current_line = '_'.join([ground_target, num_contigs_supporting_chr13, num_contigs_supporting_chr14, num_contigs_supporting_chr15, num_contigs_supporting_chr21, num_contigs_supporting_chr22])

        if last_line == '':
            last_start = current_start
            last_end = current_end
            last_line = current_line
        elif current_line != last_line:
            ground_target_x, num_contigs_supporting_chr13_x, num_contigs_supporting_chr14_x, num_contigs_supporting_chr15_x, num_contigs_supporting_chr21_x, num_contigs_supporting_chr22_x = last_line.split('_')
            print(
                '\t'.join([ground_target_x, last_start, last_end, num_contigs_supporting_chr13_x, num_contigs_supporting_chr14_x, num_contigs_supporting_chr15_x, num_contigs_supporting_chr21_x, num_contigs_supporting_chr22_x])
            )

            last_start = start
            last_end = end
            last_line = current_line
        else:
            last_end = end

    ground_target_x, num_contigs_supporting_chr13_x, num_contigs_supporting_chr14_x, num_contigs_supporting_chr15_x, num_contigs_supporting_chr21_x, num_contigs_supporting_chr22_x = last_line.split('_')
    print(
        '\t'.join([ground_target_x, last_start, last_end, num_contigs_supporting_chr13_x, num_contigs_supporting_chr14_x, num_contigs_supporting_chr15_x, num_contigs_supporting_chr21_x, num_contigs_supporting_chr22_x])
    )

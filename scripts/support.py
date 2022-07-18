#!/usr/bin/env python

import gzip
import sys

path_grounded_tsv_gz = sys.argv[1]
path_target_length_txt = sys.argv[2]
n = int(sys.argv[3])
refn = int(sys.argv[4])
path_query_to_consider_txt = sys.argv[5]

# Read chromosome lengths
ground_2_len_dict = {}

with open(path_target_length_txt) as f:
    for line in f:
        ground, target_len = line.strip().split('\t')
        ground_2_len_dict[ground] = int(target_len)

# Read queries to consider
query_to_consider_set = set()

with open(path_query_to_consider_txt) as f:
    for line in f:
        query = line.strip().split('\t')[0]
        query_to_consider_set.add(query)

# Read untangling information
ground_2_group_2_query_2_pieces_dict = {}

query_set = set()

# query, query.begin, query.end, target, target.begin, target.end, jaccard, strand, self.coverage, nth.best, ref, ref.begin, ref.end, ref.jaccard, ref.nth.best, grounded.target
with gzip.open(path_grounded_tsv_gz, "rt") as f:
    f.readline()

    for line in f:
        query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, nth_best, ref, ref_begin, ref_end, ref_jaccard, ref_nth_best, grounded_target = line.strip().split('\t')

        if query not in query_to_consider_set:
            continue

        nth_best = int(nth_best)
        ref_nth_best = int(ref_nth_best)
        if nth_best > n or ref_nth_best > refn:
            continue

        query_set.add(query)

        ref_begin = int(ref_begin)
        ref_end = int(ref_end)
        target_int = int(target.split('chm13#chr')[-1])

        # PATERNAL == 1, MATERNAL == 2
        group = 'PAT' if query.split('#')[1] in ['1', 'PAT'] else 'MAT'

        if grounded_target not in ground_2_group_2_query_2_pieces_dict:
            ground_2_group_2_query_2_pieces_dict[grounded_target] = {}
        if group not in ground_2_group_2_query_2_pieces_dict[grounded_target]:
            ground_2_group_2_query_2_pieces_dict[grounded_target][group] = {}
        if query not in ground_2_group_2_query_2_pieces_dict[grounded_target][group]:
            # `set` to avoid counting multiple times segment with self.coverage > 1
            ground_2_group_2_query_2_pieces_dict[grounded_target][group][query] = set()
        ground_2_group_2_query_2_pieces_dict[grounded_target][group][query].add((ref_begin, ref_end, target_int))

# Debugging print
# for ground, group_2_query_2_pieces_dict in ground_2_group_2_query_2_pieces_dict.items():
#     for group, query_2_pieces_dict in group_2_query_2_pieces_dict.items():
#         print(ground, group, query_2_pieces_dict)


acros2index_dict = {
    13: 0,
    14: 1,
    15: 2,
    21: 3,
    22: 4
}

num_query_computed = 0
print('\t'.join([
    'ground.target', 'start', 'end',
    'num.contigs.supporting.chr13',
    'num.contigs.supporting.chr14',
    'num.contigs.supporting.chr15',
    'num.contigs.supporting.chr21',
    'num.contigs.supporting.chr22'
]))
for ground_target, group_2_query_2_pieces_dict in ground_2_group_2_query_2_pieces_dict.items():
    print(ground_target, file=sys.stderr)
    ground_target_len = ground_2_len_dict[ground_target]

    # Prepare the counts for each base across the ground_target
    ground_target_list = [[0, 0, 0, 0, 0] for _ in range(ground_target_len)]

    for group, query_2_pieces_dict in group_2_query_2_pieces_dict.items():
        # print(ground_target, group, query_2_pieces_dict.keys())

        for query, pieces_list in query_2_pieces_dict.items():
            # print(ground_target, group, query, file=sys.stderr)

            for start, end, target_int in pieces_list:
                for pos in range(start, end):
                    # print(pos, ground_target_list[pos])
                    ground_target_list[pos][acros2index_dict[target_int]] += 1

            num_query_computed += 1
            print('Progress: {:3}%'.format(float(num_query_computed)/len(query_set)*100), file=sys.stderr)

    for pos, counter_list in enumerate(ground_target_list):
        print('\t'.join([str(x) for x in [ground_target, pos, pos + 1] + counter_list]))

    del ground_target_list

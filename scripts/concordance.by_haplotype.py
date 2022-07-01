#!/usr/bin/env python

import gzip
import sys

path_grounded_tsv_gz = sys.argv[1]
path_target_length_txt = sys.argv[2]
n = int(sys.argv[3])
refn = int(sys.argv[4])

# Read chromosome lengths
ground_2_len_dict = {}

with open(path_target_length_txt) as f:
    for line in f:
        ground, target_len = line.strip().split('\t')
        ground_2_len_dict[ground] = int(target_len)

# Read untangling information
ground_2_group_2_query_2_pieces_dict = {}

# query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, nth_best, ref, ref_begin, ref_end, ref_jaccard, ref_nth_best, grounded_target
with gzip.open(path_grounded_tsv_gz, "rt") as f:
    f.readline()

    for line in f:
        query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, nth_best, ref, ref_begin, ref_end, ref_jaccard, ref_nth_best, grounded_target = line.strip().split('\t')

        if query.startswith('HG002#'):
            nth_best = int(nth_best)
            ref_nth_best = int(ref_nth_best)
            if nth_best > n or ref_nth_best > refn:
                continue

            ref_begin = int(ref_begin)
            ref_end = int(ref_end)
            target_int = int(target.split('chm13#chr')[-1])

            # PATERNAL == 1, MATERNAL == 2
            group = 'PATERNAL' if query.split('#')[1] in ['1', 'PAT'] else 'MATERNAL'

            if grounded_target not in ground_2_group_2_query_2_pieces_dict:
                ground_2_group_2_query_2_pieces_dict[grounded_target] = {}
            if group not in ground_2_group_2_query_2_pieces_dict[grounded_target]:
                ground_2_group_2_query_2_pieces_dict[grounded_target][group] = {}
            if query not in ground_2_group_2_query_2_pieces_dict[grounded_target][group]:
                ground_2_group_2_query_2_pieces_dict[grounded_target][group][query] = set()
            ground_2_group_2_query_2_pieces_dict[grounded_target][group][query].add((ref_begin, ref_end, target_int))


# Scan each group for concordance
print('\t'.join([
    'grounded.target',
    'start.pos', 'end.pos', 'haplotype',
    'num.different.targets',
    'verkko.target',
    'all.targets', 'all.queries'
]))
for grounded_target, group_2_query_2_pieces_dict in ground_2_group_2_query_2_pieces_dict.items():
    grounded_target_len = ground_2_len_dict[grounded_target]

    for group, query_2_pieces_dict in group_2_query_2_pieces_dict.items():
        # print(grounded_target, group, query_2_pieces_dict.keys())

        query_2_filled_dict = {}
        for query, pieces_list in query_2_pieces_dict.items():
            # print(f'Fill {query} on {grounded_target}')

            # Fill (slowly!) targets over the ground target
            query_2_filled_dict[query] = [0] * grounded_target_len
            for start, end, target_int in pieces_list:
                for pos in range(start, end):
                    query_2_filled_dict[query][pos] = target_int

        last_pos = -1
        last_list = None
        last_query_list = None
        last_verkko = 0
        for pos in range(grounded_target_len):
            current_list = []
            current_query_list = []
            current_verkko = 0
            for query, filled_list in query_2_filled_dict.items():
                target = filled_list[pos]
                if target != 0:
                    current_list.append(target)
                    current_query_list.append(query)
                    if query.split('#')[1] in ['PAT', 'MAT']:
                        current_verkko = target  # The information is present in verkko's assembly

            if last_pos < 0:
                # First time
                last_pos = pos
                last_list = current_list.copy()
                last_query_list = current_query_list.copy()
                last_verkko = current_verkko
            elif last_list != current_list or last_query_list != current_query_list or last_verkko != current_verkko:
                if len(last_list) > 1 and last_verkko != 0:
                    print('\t'.join([
                        grounded_target,
                        str(last_pos), str(pos - 1), group,
                        str(len(set(last_list))),
                        str(last_verkko),
                        ','.join([str(x) for x in last_list]), ','.join(last_query_list)
                    ]))
                last_pos = pos
                last_list = current_list.copy()
                last_query_list = current_query_list.copy()
                last_verkko = current_verkko

        if len(last_list) > 1 and last_verkko != 0:
            print('\t'.join([
                grounded_target,
                str(last_pos), str(grounded_target_len - 1), group,
                str(len(set(last_list))),
                str(last_verkko),
                ','.join([str(x) for x in last_list]), ','.join(last_query_list)
            ]))

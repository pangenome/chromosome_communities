#!/usr/bin/env python

import gzip
import sys
import math

path_grounded_tsv_gz = sys.argv[1]
path_target_length_txt = sys.argv[2]
n = int(sys.argv[3])
refn = int(sys.argv[4])
estimated_identity_threshold = float(sys.argv[5])

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

            jaccard = float(jaccard)
            estimated_identity = math.exp((1.0 + math.log(2.0 * jaccard/(1.0+jaccard)))-1.0)
            if estimated_identity < estimated_identity_threshold:
                continue

            # PATERNAL == 1, MATERNAL == 2
            group = 'PAT' if query.split('#')[1] in ['1', 'PAT'] else 'MAT'

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
    'start.pos', 'end.pos',
    'contig',
    'contig.target',
    'verkko.target'
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
            for start, stop, target_int in pieces_list:
                for pos in range(start, stop):
                    query_2_filled_dict[query][pos] = target_int

        verkko_contig_name = ''
        for query in query_2_pieces_dict.keys():
            if query.split('#')[1] in ['PAT', 'MAT']:
                verkko_contig_name = query

        for query, filled_list in query_2_filled_dict.items():
            if query != verkko_contig_name:
                last_pos = -1
                last_query_target = ''
                last_verkko_target = ''
                pos = -1
                for pos in range(grounded_target_len):
                    current_query_target = filled_list[pos]
                    current_verkko_target = query_2_filled_dict[verkko_contig_name][pos]

                    if last_pos < 0:
                        # First time
                        last_pos = pos
                        last_query_target = current_query_target
                        last_verkko_target = current_verkko_target
                    elif last_query_target != current_query_target or last_verkko_target != current_verkko_target:
                        if last_query_target != 0 and last_verkko_target != 0:
                            print('\t'.join([
                                grounded_target,
                                str(last_pos), str(pos - 1),
                                query,
                                str(last_query_target),
                                str(last_verkko_target)
                            ]))
                        last_pos = pos
                        last_query_target = current_query_target
                        last_verkko_target = current_verkko_target

                # Last row
                if last_query_target != 0 and last_verkko_target != 0:
                    print('\t'.join([
                        grounded_target,
                        str(last_pos), str(pos - 1),
                        query,
                        str(last_query_target),
                        str(last_verkko_target)
                    ]))

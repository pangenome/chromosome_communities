#!/usr/bin/env python

import gzip
import sys
import math
import numpy as np

path_grounded_tsv_gz = sys.argv[1]
path_target_length_txt = sys.argv[2]
n = int(sys.argv[3])
estimated_identity_threshold = float(sys.argv[4])
path_query_to_consider_txt = sys.argv[5]

refn = 1  # In this way, the grounded interval is always the same for all `n` target hits of each segment

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
ground_2_group_2_query_2_segment_2_hits_dict = {}

ground_2_query_2_index_dict = dict()

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

        jaccard = float(jaccard)
        estimated_identity = math.exp((1.0 + math.log(2.0 * jaccard/(1.0+jaccard)))-1.0)
        if estimated_identity < estimated_identity_threshold:
            continue

        if grounded_target not in ground_2_query_2_index_dict:
            ground_2_query_2_index_dict[grounded_target] = {}
        if query not in ground_2_query_2_index_dict[grounded_target]:
            ground_2_query_2_index_dict[grounded_target][query] = len(ground_2_query_2_index_dict[grounded_target])

        nth_best = int(nth_best)
        ref_begin = int(ref_begin)
        ref_end = int(ref_end)
        target_int = int(target.split('chm13#chr')[-1])

        # PATERNAL == 1, MATERNAL == 2
        group = 'PAT' if query.split('#')[1] in ['1', 'PAT'] else 'MAT'

        if grounded_target not in ground_2_group_2_query_2_segment_2_hits_dict:
            ground_2_group_2_query_2_segment_2_hits_dict[grounded_target] = {}
        if group not in ground_2_group_2_query_2_segment_2_hits_dict[grounded_target]:
            ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group] = {}
        if query not in ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group]:
            # `set` to avoid counting multiple times segment with self.coverage > 1
            ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query] = {}
        if (query_begin, query_end) not in ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query]:
            ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query][(query_begin, query_end)] = [(ref_begin, ref_end), list()]
        ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query][(query_begin, query_end)][1].append((nth_best, target_int))

# Debugging print
# for ground_target, group_2_query_2_segment_2_hits_dict in ground_2_group_2_query_2_segment_2_hits_dict.items():
#     for group, query_2_segment_2_hits_dict in group_2_query_2_segment_2_hits_dict.items():
#         for query, segment_2_hits_dict in query_2_segment_2_hits_dict.items():
#             print(ground_target, group, query)
#             for (query_begin, query_end), range_and_hits_list in sorted(segment_2_hits_dict.items(), key=lambda item: item[1][0]):
#                 (ref_begin, ref_end), hit_list = range_and_hits_list
#                 hit_list = sorted(hit_list, key=lambda x: x[0])
#                 print('\t', (ref_begin, ref_end), hit_list)


for ground_target, group_2_query_2_segment_2_hits_dict in ground_2_group_2_query_2_segment_2_hits_dict.items():
    print(ground_target, file=sys.stderr)
    ground_target_len = ground_2_len_dict[ground_target]
    num_queries_on_ground = len(ground_2_query_2_index_dict[grounded_target])

    match_orders_np = np.zeros((num_queries_on_ground, 20000000, n), dtype=np.uint8) #todo to use ground_target_len as len, not 20000000

    for group, query_2_segment_2_hits_dict in group_2_query_2_segment_2_hits_dict.items():
        for query, segment_2_hits_dict in query_2_segment_2_hits_dict.items():
            query_index = ground_2_query_2_index_dict[grounded_target][query]
            print(ground, group, query, query_index, file=sys.stderr)

            for (query_begin, query_end), range_and_hits_list in sorted(segment_2_hits_dict.items(), key=lambda item: item[1][0]):
                (ref_begin, ref_end), hit_list = range_and_hits_list
                if ref_end <= 20000000: #todo to remove
                    hit_sorted_list = sorted(hit_list, key=lambda x: x[0])

                    # To get numpy array of the same size (n)
                    hit_sorted_filled_list = hit_sorted_list + [(0, 0)] * (n - len(hit_sorted_list))
                    target_sorted_np = np.array([target for n, target in hit_sorted_filled_list], dtype=np.uint8)
                    #print('\t', (ref_begin, ref_end), hit_sorted_list, target_sorted_np)

                    for pos in range(ref_begin, ref_end):
                        match_orders_np[query_index][pos] = target_sorted_np

    # Todo: for each position, compute the entropy across the queries and emit when something changes (compute and deduplicate at the same time?)
    for pos in range(ground_target_len):
        print(pos)
        #print(match_orders_np[0][11000000])
        #print(match_orders_np[1][11000000])

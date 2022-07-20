#!/usr/bin/env python
import collections
import gzip
import sys
import math
import numpy as np

def sdi(data):
    """ Given a hash { 'species': count } , returns the SDI

    >>> sdi({'a': 10, 'b': 20, 'c': 30,})
    1.0114042647073518"""

    from math import log as ln

    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return (float(n) / N) * ln(float(n) / N)

    N = sum(data.values())

    return -sum(p(n, N) for n in data.values() if n != 0)


path_grounded_tsv_gz = sys.argv[1]
path_target_length_txt = sys.argv[2]
n = int(sys.argv[3])  # TODO: if n == 0, detect automatically the max n and use it
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
ground_2_query_2_segment_2_hits_dict = {}

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

        if grounded_target not in ground_2_query_2_segment_2_hits_dict:
            ground_2_query_2_segment_2_hits_dict[grounded_target] = {}
        if query not in ground_2_query_2_segment_2_hits_dict[grounded_target]:
            # `set` to avoid counting multiple times segment with self.coverage > 1
            ground_2_query_2_segment_2_hits_dict[grounded_target][query] = {}
        if (query_begin, query_end) not in ground_2_query_2_segment_2_hits_dict[grounded_target][query]:
            ground_2_query_2_segment_2_hits_dict[grounded_target][query][(query_begin, query_end)] = [(ref_begin, ref_end), list()]
        ground_2_query_2_segment_2_hits_dict[grounded_target][query][(query_begin, query_end)][1].append((nth_best, target_int))


print('\t'.join(['ground.target', 'start', 'end', 'shannon_div_index', 'num.queries']))

tot_num_queries = sum([len(index_dict) for index_dict in [query_2_index_dict for query_2_index_dict in ground_2_query_2_index_dict.values()]])

num_query_computed = 0
for ground_target, query_2_segment_2_hits_dict in ground_2_query_2_segment_2_hits_dict.items():
    ground_target_len = ground_2_len_dict[ground_target]
    num_queries_on_ground = len(ground_2_query_2_index_dict[ground_target])

    match_orders_np = np.zeros((num_queries_on_ground, ground_target_len, n), dtype=np.uint8)

    for query, segment_2_hits_dict in query_2_segment_2_hits_dict.items():
        query_index = ground_2_query_2_index_dict[ground_target][query]

        for (query_begin, query_end), range_and_hits_list in sorted(segment_2_hits_dict.items(), key=lambda item: item[1][0]):
            (ref_begin, ref_end), hit_list = range_and_hits_list
            hit_sorted_list = sorted(hit_list, key=lambda x: x[0])

            # To get numpy array of the same size (n)
            hit_sorted_filled_list = hit_sorted_list + [(0, 0)] * (n - len(hit_sorted_list))
            target_sorted_np = np.array([target for n, target in hit_sorted_filled_list], dtype=np.uint8)

            for pos in range(ref_begin, ref_end):
                match_orders_np[query_index][pos] = target_sorted_np

        num_query_computed += 1
        print('Query preparation: {:.2f}%'.format(float(num_query_computed)/tot_num_queries*100), file=sys.stderr)

    num_queries = match_orders_np.shape[0]

    last_start = None
    last_end = None
    last_sdi = None
    last_count = None

    for pos in range(ground_target_len):
        if pos % 1000000 == 0:
            print('Deduplication {}: {:.2f}%'.format(ground_target, float(pos)/ground_target_len*100), file=sys.stderr)

        match_order_list = ['_'.join([str(x) for x in match_orders_np[i][pos]]) for i in range(num_queries) if sum(match_orders_np[i][pos]) > 0]
        if len(match_order_list) > 0:
            current_sdi = sdi(collections.Counter(match_order_list))
        else:
            current_sdi = -1

        current_count = len(match_order_list)

        if last_start is None:
            # It is the first info
            last_start = pos
            last_end = pos + 1
        elif last_sdi != current_sdi or last_count != current_count:
            # Something changed, so print the last line
            print('\t'.join([str(x) for x in [ground_target, last_start, last_end, last_sdi, last_count]]))

            last_start = pos
            last_end = pos + 1
        else:
            # Nothing changed, so update the end
            last_end = pos + 1

        last_sdi = current_sdi
        last_count = current_count

    print('\t'.join([str(x) for x in [ground_target, last_start, last_end, last_sdi, last_count]]))
    print('Deduplication {}: {:.2f}%'.format(ground_target, 100), file=sys.stderr)

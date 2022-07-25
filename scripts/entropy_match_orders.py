#!/usr/bin/env python

import collections
import gzip
import sys
import math


# Shannon Diversity Index
# http://en.wikipedia.org/wiki/Shannon_index

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
ground_2_segment_2_query_2_hits_dict = {}

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

        nth_best = int(nth_best)
        ref_begin = int(ref_begin)
        ref_end = int(ref_end)
        target = target.split('chm13#chr')[-1]

        # PATERNAL == 1, MATERNAL == 2
        group = 'PAT' if query.split('#')[1] in ['1', 'PAT'] else 'MAT'

        if grounded_target not in ground_2_segment_2_query_2_hits_dict:
            ground_2_segment_2_query_2_hits_dict[grounded_target] = {}
        if (ref_begin, ref_end) not in ground_2_segment_2_query_2_hits_dict[grounded_target]:
            ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)] = {}
        if query not in ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)]:
            ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)][query] = (ref_jaccard, query_begin, query_end, list())
        else:
            # The query is already present, check if the query-segment is the same
            r_jaccard, q_begin, q_end, _ = ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)][query]

            if query_begin == q_begin and query_end == q_end:
                pass  # Same query-segment: continue to fill the information
            else:
                # Different query-segment, grounded to the same target-segment.
                if ref_jaccard > r_jaccard:
                    # The current query-segment has a better grounding jaccard, so replace the old one
                    ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)][query] = (ref_jaccard, query_begin, query_end, list())
                else:
                    continue  # The current query-segment has not a better grounding

        ground_2_segment_2_query_2_hits_dict[grounded_target][(ref_begin, ref_end)][query][3].append((nth_best, target))

# for ground_target, segment_2_query_2_hits_dict in ground_2_segment_2_query_2_hits_dict.items():
#     for segment, query_2_hits_dict in sorted(segment_2_query_2_hits_dict.items(), key=lambda key: key):
#         print(ground_target, segment, len(query_2_hits_dict), query_2_hits_dict)
#         # for query, hits in query_2_hits_dict.items():
#         #     print(query, hits)


print('\t'.join(['ground.target', 'start', 'end', 'shannon_div_index', 'num.queries']))

for ground_target, segment_2_query_2_hits_dict in ground_2_segment_2_query_2_hits_dict.items():
    last_start = None
    last_end = None
    last_sdi = None
    last_count = None
    last_end_emitted = 0

    for (ref_begin, ref_end), query_2_hits_dict in sorted(segment_2_query_2_hits_dict.items(), key=lambda key: key):
        # print(ground_target, (ref_begin, ref_end), len(query_2_hits_dict))

        match_order_list = []
        for query, (_, _, _, hit_list) in query_2_hits_dict.items():
            sorted_hit_list = sorted(hit_list, key=lambda x: x[0])
            # print(query, sorted_hit_list)
            match_order_list.append('_'.join([target for nth_best, target in sorted_hit_list]))

        current_sdi = sdi(collections.Counter(match_order_list))
        current_count = len(match_order_list)

        # Deduplicate and emit
        if last_start is None:
            # It is the first info
            last_start = ref_begin
            last_end = ref_end
        elif last_end != ref_begin or last_sdi != current_sdi or last_count != current_count:
            # There is a hole, or something changed, so print the last line

            if last_end_emitted != last_start:
                # Fill ranges with missing information
                print('\t'.join([str(x) for x in [ground_target, last_end_emitted, last_start, -1, 0]]))

            print('\t'.join([str(x) for x in [ground_target, last_start, last_end, last_sdi, last_count]]))

            last_end_emitted = last_end
            last_start = ref_begin
            last_end = ref_end
        else:
            # Nothing changed, so update the end
            last_end = ref_end

        last_sdi = current_sdi
        last_count = current_count

    if last_end_emitted != last_start:
        # Fill ranges with missing information
        print('\t'.join([str(x) for x in [ground_target, last_end_emitted, last_start, -1, 0]]))

    print('\t'.join([str(x) for x in [ground_target, last_start, last_end, last_sdi, last_count]]))

    last_end_emitted = last_end

    if last_end != ground_2_len_dict[ground_target]:
        # Fill ranges with missing information
        print('\t'.join([str(x) for x in [ground_target, last_end_emitted, ground_2_len_dict[ground_target], -1, 0]]))

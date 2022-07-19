#!/usr/bin/env python

import gzip
import sys

path_grounded_tsv_gz = sys.argv[1]
jaccard_treshold = float(sys.argv[2])
self_cov_threshold = float(sys.argv[3])
path_query_to_consider_txt = sys.argv[4]

refn = 1  # In this way, the grounded interval is always the same for all `n` target hits of each segment

# Read queries to consider
query_to_consider_set = set()

with open(path_query_to_consider_txt) as f:
    for line in f:
        query = line.strip().split('\t')[0]
        query_to_consider_set.add(query)

# Read untangling information
ground_2_group_2_query_2_segment_2_hits_dict = {}

# query, query.begin, query.end, target, target.begin, target.end, jaccard, strand, self.coverage, nth.best, ref, ref.begin, ref.end, ref.jaccard, ref.nth.best, grounded.target
with gzip.open(path_grounded_tsv_gz, "rt") as f:
    f.readline()

    for line in f:
        query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, nth_best, ref, ref_begin, ref_end, ref_jaccard, ref_nth_best, grounded_target = line.strip().split('\t')

        if query not in query_to_consider_set:
            continue

        nth_best = int(nth_best)
        ref_nth_best = int(ref_nth_best)
        if ref_nth_best > refn:
            continue

        jaccard = float(jaccard)
        if jaccard < jaccard_treshold:
            continue

        self_coverage = float(self_coverage)
        if 0 < self_cov_threshold < self_coverage:
            continue

        query_begin = int(query_begin)
        query_end = int(query_end)
        ref_begin = int(ref_begin)
        ref_end = int(ref_end)
        target = target.split('chm13#chr')[-1]

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
            ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query][(query_begin, query_end)] = [(ref_begin, ref_end), set()]
        ground_2_group_2_query_2_segment_2_hits_dict[grounded_target][group][query][(query_begin, query_end)][1].add((jaccard, target))


print('\t'.join(['query', 'query.start', 'query.end', 'ground', 'ground.start', 'ground.end', 'different.targets']))
for ground, group_2_query_2_segment_2_hits_dict in ground_2_group_2_query_2_segment_2_hits_dict.items():
    for group, query_2_segment_2_hits_dict in group_2_query_2_segment_2_hits_dict.items():
        for query, segment_2_hits_dict in query_2_segment_2_hits_dict.items():

            for segment in sorted(segment_2_hits_dict.keys()):
                (ref_begin, ref_end), hit_set = segment_2_hits_dict[segment]
                target_set = set([target for (jaccard, target) in hit_set])

                if len(target_set) > 1:
                    query_begin, query_end = segment
                    print('\t'.join([query, str(query_begin), str(query_end), ground, str(ref_begin), str(ref_end), ','.join([str(x) for x in target_set])]))

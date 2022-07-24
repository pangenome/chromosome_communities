#!/usr/bin/env python

# Shannon Diversity Index
# http://en.wikipedia.org/wiki/Shannon_index

import collections
import gzip
import sys
import math

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
window_size = int(sys.argv[3])
n = int(sys.argv[4])
refn = int(sys.argv[5])
estimated_identity_threshold = float(sys.argv[6])
path_query_to_consider_txt = sys.argv[7]

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

# query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, nth_best, ref, ref_begin, ref_end, ref_jaccard, ref_nth_best, grounded_target
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

        query_set.add(query)

        ref_begin = int(ref_begin)
        ref_end = int(ref_end)
        target = target.split('chm13#chr')[-1]

        # PATERNAL == 1, MATERNAL == 2
        group = 'PATERNAL' if query.split('#')[1] in ['1', 'PAT'] else 'MATERNAL'

        if grounded_target not in ground_2_group_2_query_2_pieces_dict:
            ground_2_group_2_query_2_pieces_dict[grounded_target] = {}
        if group not in ground_2_group_2_query_2_pieces_dict[grounded_target]:
            ground_2_group_2_query_2_pieces_dict[grounded_target][group] = {}
        if query not in ground_2_group_2_query_2_pieces_dict[grounded_target][group]:
            ground_2_group_2_query_2_pieces_dict[grounded_target][group][query] = []
        ground_2_group_2_query_2_pieces_dict[grounded_target][group][query].append((ref_begin, ref_end, target))


num_query_computed = 0
print('\t'.join(['ground.target', 'start.pos', 'end.pos', 'contig', 'shannon_div_index']))
for ground_target, group_2_query_2_pieces_dict in ground_2_group_2_query_2_pieces_dict.items():
    print(ground_target, file=sys.stderr)
    ground_target_len = ground_2_len_dict[ground_target]

    for group, query_2_pieces_dict in group_2_query_2_pieces_dict.items():
        # print(ground_target, group, query_2_pieces_dict.keys())

        for query, pieces_list in query_2_pieces_dict.items():
            # print(ground_target, group, query, file=sys.stderr)

            # Fill (slowly!) targets over the ground target
            query_on_ref_list = [0] * ground_target_len
            for start, stop, target in pieces_list:
                for pos in range(start, stop):
                    query_on_ref_list[pos] = target

            # Compute entropy for each window
            for start in range(0, len(query_on_ref_list), window_size):
                end = start + window_size - 1
                if end > ground_2_len_dict[ground_target]:
                    end = ground_2_len_dict[ground_target]
                targets_list = [x for x in query_on_ref_list[start:end] if x != 0]
                if len(targets_list) > 0:
                    shannon_div_index = sdi(collections.Counter(targets_list))
                else:
                    shannon_div_index = -1

                print('\t'.join([ground_target, str(start), str(end + 1), query, str(shannon_div_index)]))

            num_query_computed += 1
            print('Progress: {:.2f}%'.format(float(num_query_computed)/len(query_set)*100), file=sys.stderr)

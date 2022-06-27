#!/usr/bin/env python

import gzip
import sys

path_untangle_tsv_gz = sys.argv[1]
path_contig_to_chr_tsv = sys.argv[2]

# Partitioning info
contig_to_chr = {}
with open(path_contig_to_chr_tsv) as f:
    for line in f:
        contig, chr = line.strip().split('\t')

        contig_to_chr[contig] = chr


# Read untangling information
key_2_info_dict = {}

# query.name, query.start, query.end, ref.name, ref.start, ref.end, score, inv, self.cov, nth.best
with gzip.open(path_untangle_tsv_gz, "rt") as f:
    print(f.readline().strip())

    for line in f:
        line = line.strip()
        query_name, query_start, query_end, ref_name, ref_start, ref_end, score, inv, self_cov, nth_best = line.split('\t')

        key = (query_name, query_start, query_end)
        if key not in key_2_info_dict:
            key_2_info_dict[key] = []
        key_2_info_dict[key].append([int(nth_best), float(score), ref_name, line])

for key, info_list in key_2_info_dict.items():
    query_name = key[0]

    info_list.sort(key=lambda x: x[0], reverse=False)  # sort by nth.best (rank)
    #print(key, [x[0:3] for x in info_list])

    if len(info_list) > 1:
        # Check if we have to fix the targets
        max_jaccard = info_list[0][1]

        info_with_max_jaccard_list = []
        targets_with_max_jaccard_set = set()
        info_without_max_jaccard_list = []
        for nth_best, jaccard, ref_name, line in info_list:
            if jaccard == max_jaccard:
                info_with_max_jaccard_list.append([nth_best, jaccard, ref_name, line])
                targets_with_max_jaccard_set.add(ref_name)
            else:
                info_without_max_jaccard_list.append([nth_best, jaccard, ref_name, line])

        # If there are multiple best-hits and there is the grounded target in one of them, and it is not the first one, then fix this
        if len(info_with_max_jaccard_list) > 1 and contig_to_chr[query_name] in targets_with_max_jaccard_set and contig_to_chr[query_name] != info_with_max_jaccard_list[0][2]:
            # Find the first case with the target we want (contig_to_chr[query_name])
            i = 0
            for nth_best, jaccard, ref_name, line in info_with_max_jaccard_list:
                if ref_name == contig_to_chr[query_name]:
                    break
                i += 1

            best_hit = info_with_max_jaccard_list.pop(i)

            # Put the best hit on the head
            info_with_max_jaccard_list = [best_hit] + info_with_max_jaccard_list

            del info_list
            info_list = []
            for i, (nth_best, jaccard, ref_name, line) in enumerate(info_with_max_jaccard_list):
                # Fix rank
                query_name, query_start, query_end, ref_name, ref_start, ref_end, score, inv, self_cov, nth_best = line.split('\t')
                nth_best = i+1
                line = '\t'.join([query_name, query_start, query_end, ref_name, ref_start, ref_end, score, inv, self_cov, str(nth_best)])
                info_list.append([nth_best, jaccard, ref_name, line])

            # Put the other rows, that have already correct ranks
            info_list.extend(info_without_max_jaccard_list)
            info_list.sort(key=lambda x: x[0], reverse=False)  # sort by nth.best (rank)

            # print the lines
            for x in info_list:
                print(x[3])

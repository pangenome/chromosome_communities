#!/usr/bin/env python

# Shannon Diversity Index
# http://en.wikipedia.org/wiki/Shannon_index

import collections
import sys

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
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n != 0)


import gzip
import sys


path_grounded_tsv_gz=sys.argv[1]
path_target_length_txt=sys.argv[2]


ground_2_len_dict = {}

with open(path_target_length_txt) as f:
    for line in f:
        ground, target_len = line.strip().split('\t')
        ground_2_len_dict[ground] = int(target_len)


ground_2_query_2_pieces_dict = {}

#query   query.begin query.end   target  target.begin    target.end  jaccard strand  self.coverage   ref ref.begin   ref.end grounded.target
with gzip.open(path_grounded_tsv_gz,"rt") as f:
    f.readline()

    for line in f:
        query, query_begin, query_end, target, target_begin, target_end, jaccard, strand, self_coverage, ref, ref_begin, ref_end, grounded_target = line.strip().split('\t')
        
        if not query.startswith('chm') and not query.startswith('grch'):
            ref_begin = int(ref_begin)
            ref_end = int(ref_end)

            if grounded_target not in ground_2_query_2_pieces_dict:
                ground_2_query_2_pieces_dict[grounded_target] = {}
            if query not in ground_2_query_2_pieces_dict[grounded_target]:
                ground_2_query_2_pieces_dict[grounded_target][query] = []
            ground_2_query_2_pieces_dict[grounded_target][query].append((ref_begin, ref_end, int(target.split('chm13#chr')[-1])))


#ground_2_query_2_entropy_dict = {}


window_size = 50000

print('\t'.join(['query', 'ground.target', 'start', 'end', 'shannon_div_index']))
for ground_target, query_2_pieces_dict in ground_2_query_2_pieces_dict.items():
    #ground_2_query_2_entropy_dict[ground_target] = {}

    for query, pieces_list in query_2_pieces_dict.items():
        #ground_2_query_2_entropy_dict[ground_target][query] = []
        #print(ground_target, query, sorted(pieces_list))
        
        # Fill (slowly!) targets over the ground target
        query_on_ref_list = [0] * ground_2_len_dict[ground_target]
        for start, stop, target in pieces_list:
            for pos in range(start, stop):
                query_on_ref_list[pos] = target
        
        # Compute entropy for each window
        for start in range(0, len(query_on_ref_list), window_size):
            end = start + window_size
            #if end > ground_2_len_dict[ground_target]:
            #    end = ground_2_len_dict[ground_target]
            shannon_div_index = sdi(collections.Counter(query_on_ref_list[start:end]))

            #ground_2_query_2_entropy_dict[ground_target][query].append(
            #    shannon_div_index
            #)

            print('\t'.join([query, ground_target, str(start), str(end), str(shannon_div_index)]))

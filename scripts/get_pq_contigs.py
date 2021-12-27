import os
import sys

# Notes:
# - acro-contigs (from chr13,14,15,21,22) longer than 10Mbps vs CHM13 chr13,14,15,21,22 (wfmash -s 10k/20k/50k/100k -p 90/95/98/99/99.5 -l 3*s/0)
# - for each run, count contigs that have trans / cys mappings
# - using CHM13 centromere annotations, counts contigs trans that go from p to q
# erikg
# - Try to measure the average size of the trans mappings for each mapping run
# - How long are the parts that don't map cys?
# - Also, the length of the longest (or average of longest)
# - To give an impression of scale

path_mappings_paf = sys.argv[1]
shift_coordinates = int(sys.argv[2])

chromosome2centromere_dict = {
	'chm13#chr13': {
		'p' : (16500000, 17700000),
		'q' : (17700000, 18767509),
	},
	'chm13#chr14': {
		'p' : (16100000, 17200000),
		'q' : (17200000, 18200000),
	},
	'chm13#chr15': {
		'p' : (17500000, 19000000),
		'q' : (19000000, 20500000),
	},
	'chm13#chr21': {
		'p' : (10900000, 12000000),
		'q' : (12000000, 12081303),
	},
	'chm13#chr22': {
		'p' : (13700000, 15000000),
		'q' : (15000000, 17400000),
	}
}

seq2len_dict = {}
query2target2info_dict = {}

with open(path_mappings_paf) as f:
    for line in f:
        query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len, _, est_id = line.strip().split('\t')

        # Skip reference sequences
        if 'chm13' in query or 'grch38' in query:
            continue

        query_len = int(query_len)
        query_start = int(query_start)
        query_end = int(query_end)
        target_len = int(target_len)
        target_start = int(target_start)
        target_end = int(target_end)

        seq2len_dict[query] = query_len
        seq2len_dict[target] = target_len

        #print(query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len, est_id)

        if query not in query2target2info_dict:
            query2target2info_dict[query] = {}
        if target not in query2target2info_dict[query]:
            query2target2info_dict[query][target] = list()
        query2target2info_dict[query][target].append(
            ((query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id)
        )

query2arm2target2info_dict = {}

for query, target2info_dict in query2target2info_dict.items():
    query2arm2target2info_dict[query] = {'p' : {}, 'q': {}, 'pq': {}, 'not_good': {}}

    #print('\t\t', query)
    for target, info_list in target2info_dict.items():
        #print('\t\t\t', target, sorted(info_list, key=lambda x: x[2][0]))

        for (query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id in info_list:
            #print('\t\t\t\t',(query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id)

            if target_end <= chromosome2centromere_dict[target]['p'][0] - shift_coordinates:
                arm = 'p'
            elif target_start <= chromosome2centromere_dict[target]['p'][0] - shift_coordinates and target_end >= chromosome2centromere_dict[target]['q'][1] + shift_coordinates:
                arm = 'pq'
            elif target_start >= chromosome2centromere_dict[target]['q'][1] + shift_coordinates:
                arm = 'q'
            else:
                arm = 'not_good'

            #print(arm, '\t\t\t\t',query, (query_start, query_end), strand, target, (target_start, target_end), num_matches, alignment_len, est_id)

            if target not in query2arm2target2info_dict[query][arm]:
                query2arm2target2info_dict[query][arm][target] = list()
            query2arm2target2info_dict[query][arm][target].append(
                ((query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id)
            )

for query, arm2target2info_dict in query2arm2target2info_dict.items():
    #print(query, arm2target2info_dict)
    if len(arm2target2info_dict['pq']) > 0 or (len(arm2target2info_dict['p']) > 0 and len(arm2target2info_dict['q']) > 0):
        print(query)

path_output = path_mappings_paf + '.pq_contigs.paf'
with open(path_output, 'w') as fw:
    for query, arm2target2info_dict in query2arm2target2info_dict.items():
        if len(arm2target2info_dict['pq']) > 0 or (len(arm2target2info_dict['p']) > 0 and len(arm2target2info_dict['q']) > 0):
            #print(query)
            for arm, target2info_dict in arm2target2info_dict.items():
                for target, info_list in target2info_dict.items():
                    for (query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id in info_list:
                        fw.write('\t'.join([str(x) for x in [query, seq2len_dict[query], query_start, query_end, strand, target, seq2len_dict[target], target_start, target_end, num_matches, alignment_len, 255, est_id]]) + '\n')

# num_contigs_trans = 0
# num_contigs_cys = 0
# for query, target2info_dict in query2target2info_dict.items():
#     if len(query2target2info_dict[query]) == 1:
#         num_contigs_cys += 1
#     else:
#         num_contigs_trans += 1
#
#         path_output =path_xxx_paf + f'.split.{query}.paf'
#         with open(path_output, 'w') as fw:
#             for target, info_list in target2info_dict.items():
#                 for (query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, est_id in info_list:
#                     fw.write('\t'.join([str(x) for x in [query, seq2len_dict[query], query_start, query_end, strand, target, seq2len_dict[target], target_start, target_end, num_matches, alignment_len, 255, est_id]]) + '\n')
#
#         print(f'\t\t {query} --> targets: {set(query2target2info_dict[query].keys())}')
#
# print('\t\t num. contigs trans|cys {:3}|{:3} --> ratio {:.3f}'.format(num_contigs_trans, num_contigs_cys, num_contigs_trans/num_contigs_cys))

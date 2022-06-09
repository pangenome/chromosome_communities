import sys

PATH_PAF=sys.argv[1]

query_2_target_2_wmap_dict = {}

with open(PATH_PAF) as f:
    for line in f:
        query, _, query_start, query_end, _, target, _, target_start, target_end, _, _, _, estimated_identity = line.strip().split('\t')                                                     

        num_mapped_bases = max(int(query_end) - int(query_start), int(target_end) - int(target_start))                                                                                       
        estimated_identity = float(estimated_identity.split('id:f:')[1]) / 100.0
        weight = float(num_mapped_bases) * estimated_identity

        #print(line)
        if query not in query_2_target_2_wmap_dict:
            query_2_target_2_wmap_dict[query] = {}
        if target not in query_2_target_2_wmap_dict[query]:
            query_2_target_2_wmap_dict[query][target] = 0.0
        query_2_target_2_wmap_dict[query][target] += weight



target_2_contig_dict = {}

for query, target_2_wmap_dict in query_2_target_2_wmap_dict.items():
    #print(query, dict(sorted(target_2_wmap_dict.items(), key=lambda item: item[1], reverse=True)))

    best_target = max(target_2_wmap_dict, key=target_2_wmap_dict.get)

    if best_target not in target_2_contig_dict:
        target_2_contig_dict[best_target] = []
    target_2_contig_dict[best_target].append(query)

for target, query_list in target_2_contig_dict.items():
    for query in query_list:
        print('\t'.join([query, target]))

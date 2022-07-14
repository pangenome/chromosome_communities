import sys


def set_chr_nr(chr):
    """ Sort by chromosome """
    if chr:
        new = chr[3:]
        if new == 'X':
            chr = 23
        elif new == 'Y':
            chr = 24
        elif new == 'M':
            chr = 25
        elif not chr.startswith('chr'):
            chr = 26
        else:
            chr = int(new)
    else:
        chr = 0
    return chr


labels_set = set()
id_2_labels_dict = {}

for i, filename in enumerate(sys.argv[1:]):
    id_2_labels_dict[i] = []
    with open(filename) as f:
        for line in f:
            label = line.strip().split('-')[-1].split('#')[-1].split('_')[0]
            label = label.replace('unmapped', 'not.partitioned')
            labels_set.add(label)
            id_2_labels_dict[i].append(label)

community_2_info_dict = {}

import collections

print('\t'.join(['community.of'] + [f'chr{x}' for x in range(1, 23)] + ['chrX', 'chrY', 'chrM', 'not.partitioned']))
for id, label_list in id_2_labels_dict.items():
    counter_list = [0] * 26  # Max number of chromosomes (including the `not.partitioned` set)

    chr_2_count_dict = collections.Counter(label_list)
    chr_with_max_count = max(chr_2_count_dict, key=chr_2_count_dict.get)
    for label, count in chr_2_count_dict.items():
        counter_list[set_chr_nr(label) - 1] = count

    counter_list = [str(x) for x in counter_list]
    if chr_with_max_count not in community_2_info_dict:
        community_2_info_dict[chr_with_max_count] = []
    community_2_info_dict[chr_with_max_count].append(counter_list)

    #print('\t'.join([chr_with_max_count] + counter_list))

chromosomes = community_2_info_dict.keys()
for chr in sorted(chromosomes, key=lambda x: set_chr_nr(x)):
    for i, info in enumerate(community_2_info_dict[chr]):
        suffix = ''
        if len(community_2_info_dict[chr]) > 1:
            suffix = f'_{i+1}'
        print('\t'.join([chr + suffix] + info))

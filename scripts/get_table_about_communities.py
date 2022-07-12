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
community_to_labels_dict = {}

for i, filename in enumerate(sys.argv[1:]):
    community_to_labels_dict[i] = []
    with open(filename) as f:
        for line in f:
            label = line.strip().split('#')[-1].split('_')[0]
            labels_set.add(label)
            community_to_labels_dict[i].append(label)

import collections

print('\t'.join(['community.of'] + [f'chr{x}' for x in range(1, 23)] + ['chrX', 'chrY', 'chrM', 'not.partitioned']))
for community, label_list in community_to_labels_dict.items():
    counter_list = [0] * 26  # Max number of chromosomes (including the `not.partitioned` set)

    chr_2_count_dict = collections.Counter(label_list)
    chr_with_max_count = max(chr_2_count_dict, key=chr_2_count_dict.get)
    for label, count in chr_2_count_dict.items():
        counter_list[set_chr_nr(label) - 1] = count

    counter_list = [str(x) for x in counter_list]
    print('\t'.join([chr_with_max_count] + counter_list))

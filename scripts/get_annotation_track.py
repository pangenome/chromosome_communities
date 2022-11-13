import sys

path_x_tsv = sys.argv[1]

with open(path_x_tsv) as f:
    for line in f:
        query, query_start, query_end, info, _, ann_start, ann_end, label, _, _, _, _, _, ann_len = line.strip().split('\t')

        info_list = info.split('_')
        ground_ref_start, ground_ref_end = info_list[8:10]

        if int(ann_start) <= int(query_start):
            new_start = int(ground_ref_start)
        else:
            new_start = int(ground_ref_start) + (int(ann_start) - int(query_start))

        new_end = new_start + int(ann_len)

        info_list[0] = label
        info_list[1] = '.'
        info_list[2] = '.'
        info_list[3] = '.'
        info_list[4] = '.'
        info_list[5] = '1'
        info_list[7] = str(new_start)
        info_list[8] = str(new_end)


        print(
            '\t'.join([
                query + '_',
                query_start if int(ann_start) <= int(query_start) else ann_start,
                ann_end if int(ann_end) <= int(query_end) else query_end,
            ] + info_list)
        )

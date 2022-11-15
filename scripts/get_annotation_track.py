import sys

path_x_tsv = sys.argv[1]

with open(path_x_tsv) as f:
    for line in f:
        # HG002#1#JAHKSE010000070.1
        # 6188
        # 10236
        # chm13#chr13_113563195_113564966_0.413751_-_1_1_chm13#chr13_113563195_113564966_0.413751_1_chm13#chr13
        # HG002#1#JAHKSE010000070.1	0
        # 18254
        # Err
        # 4048
        query, query_start, query_end, info, _, ann_start, ann_end, label, ann_len = line.strip().split('\t')

        info_list = info.split('_')
        ground_ref_start, ground_ref_end = info_list[8:10]

        if int(ann_start) <= int(query_start):
            new_start = int(ground_ref_start)
        else:
            new_start = int(ground_ref_start) + (int(ann_start) - int(query_start))

        new_end = new_start + int(ann_len)

        info_list[0] = label           # target
        info_list[1] = '.'             # target.begin
        info_list[2] = '.'             # target.end
                                       # jaccard
                                       # strand
                                       # self.coverage
                                       # nth.best
                                       # ref
        info_list[8] = str(new_start)  # ref.begin
        info_list[9] = str(new_end)    # ref.end
                                       # ref.jaccard
                                       # ref.nth.best
                                       # grounded.target

        # query, query.begin, query.end,
        # target, target.begin, target.end,
        # jaccard, strand, self.coverage, nth.best,
        # ref, ref.begin, ref.end, ref.jaccard, ref.nth.best,
        # grounded.target
        print(
            '\t'.join([
                query + '.',
                query_start if int(ann_start) <= int(query_start) else ann_start,
                ann_end if int(ann_end) <= int(query_end) else query_end,
            ] + info_list)
        )

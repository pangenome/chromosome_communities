import sys


path_concordance_tsv = sys.argv[1]

# info = [num.bp.ok, num.bo.not.ok, num.bp.ok.parm, num.bp.not.ok.parm, num.bp.ok.qarm, num.bp.ok.qarm]
ground_2_haplotype_2_info_dict = {}

# grounded.target, start.pos, end.pos, haplotype, num.different.targets, verkko.target, all.targets, all.queries
with open(path_concordance_tsv) as f:
    f.readline()

    for line in f:
        grounded_target, start_pos, end_pos, haplotype, num_different_targets = line.strip().split('\t')[0:5]
        start_pos = int(start_pos)
        end_pos = int(end_pos)

        if grounded_target not in ground_2_haplotype_2_info_dict:
            ground_2_haplotype_2_info_dict[grounded_target] = {}
        if haplotype not in ground_2_haplotype_2_info_dict[grounded_target]:
            ground_2_haplotype_2_info_dict[grounded_target][haplotype] = [0, 0]

        diff_end_start = end_pos - start_pos

        if int(num_different_targets) == 1:
            ground_2_haplotype_2_info_dict[grounded_target][haplotype][0] += diff_end_start
        else:
            ground_2_haplotype_2_info_dict[grounded_target][haplotype][1] += diff_end_start

print('\t'.join(['grounded.target', 'haplotype', 'num.bp.ok', 'num.bp.not.ok', 'percentage.bp.ok']))
for grounded_target, haplotype_2_info_dict in ground_2_haplotype_2_info_dict.items():
    for haplotype, info_list in haplotype_2_info_dict.items():
        print('\t'.join([grounded_target, haplotype] + [str(x) for x in info_list] + [str(info_list[0]/sum(info_list))]))

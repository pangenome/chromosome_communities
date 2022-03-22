import sys

diploid = sys.argv[1] == 'diploid'

for line in sys.stdin:
    if line.startswith('#'):
        print(line.strip())
    else:
        line_list = line.strip().split('\t')

        for i, gt in enumerate(line_list[9:]):
            new_gt = ''
            if diploid and len(gt) != 3:
                if gt in ['01']:
                    new_gt == f'{gt}|{gt}'
                else:
                    new_gt == '0|0'
            for j, c in enumerate(gt):
                if j % 2 == 0 and c == '.':
                    c = '0'
                new_gt += c
            
            line_list[9 + i] = new_gt
        
        print('\t'.join(line_list))

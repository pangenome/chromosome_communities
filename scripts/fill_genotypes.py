import sys

diploid = True if len(sys.argv) > 1 and sys.argv[1] == 'diploid' else False

for line in sys.stdin:
    if line.startswith('#'):
        print(line.strip())
    else:
        line_list = line.strip().split('\t')

        for i, gt in enumerate(line_list[9:]):
            new_gt = ''
            for x in gt.split(','):
                if diploid and len(x) != 3:
                    if x[0] in '01':
                        x = f'{x[0]}|{x[0]}'
                    else:
                        x = '0|0'
    
                for j, c in enumerate(x):
                    if j % 2 == 0 and c == '.':
                        c = '0'
                    new_gt += c
                
                new_gt += ','
            
            line_list[9 + i] = new_gt.strip(',')
        
        print('\t'.join(line_list))

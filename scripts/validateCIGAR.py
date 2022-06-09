import sys

path_seq = sys.argv[1]
path_out = sys.argv[2]

with open(path_seq) as f:
    query = f.readline().strip().lstrip('>')
    target = f.readline().strip().lstrip('<')

len_query = len(query)
len_target = len(target)

with open(path_out) as f:
    score, cigar = f.readline().strip().split('\t')

#print(len_query, query)
#print(len_target, target)
#print(score)
#print(cigar)

i = 0
j = 0

matches = 0
mismatches = 0
insertions = 0
insertions_bp = 0
deletions = 0
deletions_bp = 0

i_max = i + len_query
j_max = j + len_target

last_op_size = ""
ok = True
for c in cigar:
    if c in '0123456789':
        last_op_size += c
        continue
    #print(f'{last_op_size}{c}')
    if c == 'M':
        # check that we match
        for _ in range(int(last_op_size)):
            if query[i] != target[j]:
                print(f'ERROR: mismatch at {i}-{j} ({query[i]} != {target[j]})')
                ok = False
            if i >= i_max:
                print(f'ERROR: query out of bounds at {i}-{j}')
                ok = False
            if j >= j_max:
                print(f'ERROR: target out of bounds at {i}-{j}')
                ok = False
            i += 1
            j += 1

        matches += int(last_op_size)
    elif c == 'X':
        # check that we don't match
        for _ in range(int(last_op_size)):
            if query[i] == target[j]:
                print(f'ERROR: match at {i}-{j} ({query[i]} == {target[j]})')
                ok = False
            if i >= i_max:
                print(f'ERROR: query out of bounds at {i}-{j}')
                ok = False
            if j >= j_max:
                print(f'ERROR: target out of bounds at {i}-{j}')
                ok = False
            i += 1
            j += 1

        mismatches += int(last_op_size)
    elif c == 'I':
        j += int(last_op_size)

        insertions += 1
        insertions_bp += int(last_op_size)
    elif c == 'D'  :
        i += int(last_op_size)

        deletions += 1
        deletions_bp += int(last_op_size)
    else:
        print('ERROR: unrecognized CIGAR operator')
        ok = False
        break
    last_op_size = ""

gap_compressed_identity = 0.0
block_identity = 0.0

if not ok:
    matches = 0
    mismatches = 0
    insertions = 0
    insertions_bp = 0
    deletions = 0
    deletions_bp = 0
else:
    gap_compressed_identity = float(matches) / float(matches + mismatches + insertions + deletions)
    block_identity = float(matches) / float(matches + mismatches + insertions_bp + deletions_bp)

print('\t'.join(
    ['OK' if ok else 'INCORRECT'] + [
        str(x) for x in [
        gap_compressed_identity,
        block_identity,
        matches,
        mismatches,
        insertions,
        insertions_bp,
        deletions,
        deletions_bp]
    ]
))

#athaliana16/athaliana16.s50000.p99.sample100.paf.seq.pair_99.seq athaliana16/athaliana16.s50000.p99.sample100.paf.pair_99.med.out
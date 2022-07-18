import sys

prev_chrom = ''
prev_start = ''
prev_end = ''
prev_count = ''

# chrom, start_chrom, end_chrom, position, count
for line in sys.stdin:
	chrom, _, _, position, count = line.strip().split('\t')
	position = int(position) # It starts from 1
	count = int(count)

	if prev_chrom == '':
		prev_chrom = chrom
		prev_start = position
		prev_end = position
		prev_count = count
	else:
		if prev_chrom != chrom or prev_count != count:
			print(prev_chrom, prev_start-1, prev_end, prev_count)

			prev_start = position
			prev_chrom = chrom
			prev_count = count
		else:
			prev_end = position

print(prev_chrom, prev_start-1, prev_end, prev_count)

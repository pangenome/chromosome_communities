#! /bin/sh

seqName=$1
nn=$2
lpart=$3
pathLDhat=$4
jobId=$5
prefix=$6
output_dir=$7

#echo "seqName - $seqName"
#echo "nn - $nn"
#echo "lpart - $lpart"
#echo "pathLDhat - $pathLDhat"
#echo "jobId - $jobId"
#echo "out - $out"
#echo "output_dir - $output_dir"
#( pwd )

# sed -i "1i $nn $lpart 1" $seqName

awk -v n=1 -v s="$nn $lpart 1" 'NR == n {print s} {print}' "$seqName" > "${output_dir}"/"${prefix}".FastaTemp.${jobId}.fa
mv "${output_dir}"/"${prefix}".FastaTemp.${jobId}.fa "$seqName"

mkdir -p "${prefix}"${jobId} && cd "${prefix}"${jobId}
$pathLDhat/convert -seq ../"$seqName" | tail -n 7 > "${prefix}".temp.${jobId}.txt

no_output=$(grep 'No data to output' "${prefix}".temp.${jobId}.txt -c)
if [ $no_output != 0 ]; then
  printf "0\t0\t0\t0\t0\t0\t\n" > ../"$output_dir"/"${prefix}".sums_part_main.${jobId}.txt
else
  # Takes only the values (on the right of the '=' character), remove the last empty line,
  # trim spaces, and put in a tab-separated row, with a '\n' at the end
  (head "${prefix}".temp.${jobId}.txt -n 6 | tr -d ' ' | cut -f 2 -d '=' | tr '\n' '\t '; echo) \
    > ../"$output_dir"/"${prefix}".sums_part_main.${jobId}.txt
fi

cd .. && rm -rf "${prefix}"${jobId}

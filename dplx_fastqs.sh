dh=$1
#dh =/media

#dplx

for fh in $dh/dplxd/*qcd.fastq; do
	echo $i
	cutadapt -u -12 -o $fh.trimmed.fastq $fh
done

mv $dh/dplxd/*.trimmed.fastq $dh/trimmed
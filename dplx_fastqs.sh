dh=$1
#dh =/media

#dplx

for fh in $dh/dplxd/*qcd.fastq; do
	if [ ! -f $fh.trimmed.fastq ]; then
		echo $i
		cutadapt -u -12 -o $fh.trimmed.fastq $fh
	fi
done

mv $dh/dplxd/*.trimmed.fastq $dh/trimmed
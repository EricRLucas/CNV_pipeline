bamfilefolder=$1
bamfilepaths=$2
samplenum=$3 

samplename=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f1))
bampath=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f2))
bamfile=${bamfilefolder}/${samplename}.bam

echo Downloading bam file for $samplename from $bampath to $bamfile

# Download the required bam
curl ${bampath} > ${bamfile} 2> ${bamfilefolder}/${samplename}_download.log
curl ${bampath}.bai > ${bamfile}.bai 2>> ${bamfilefolder}/${samplename}_download.log


bamfilepaths=$1
joblimit=$2
study=$3

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch118/malaria/team112/personal/el10
outputfolder=$rootfolder/$study
bamfilefolder=$outputfolder/bamlinks
logfolder=$outputfolder/logfolders/coverage_and_diagnostic_reads
errorfolder=$outputfolder/errorfolders/coverage_and_diagnostic_reads

echo "Downloading bam files for samples from ${bamfilepaths} on `date`." >> $outputfolder/added_samples.log

# Create the output folders if necessary
mkdir -p $bamfilefolder
mkdir -p $logfolder
mkdir -p $errorfolder

# Get the number of bamfiles that need processing
numbams=($(wc -l $bamfilepaths))
echo "This sample set contains ${numbams} bam files." >> $outputfolder/added_samples.log

bsub -J "bamArray[1-$numbams]%$joblimit" \
     -R"select[mem>500] rusage[mem=500]" \
     -M500 \
     -o $logfolder/log_%J.%I.txt \
     -e $errorfolder/error_%J.%I.txt \
     ' '${scriptsfolder}'/get_windowed_coverage_and_diagnostic_reads.sh '${bamfilefolder}' '${bamfilepaths}' ${LSB_JOBINDEX} 


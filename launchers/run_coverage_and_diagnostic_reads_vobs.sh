bamfilepaths=$1
joblimit=$2
study=$3

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10
outputfolder=$rootfolder/$study
bamfilefolder=$outputfolder/bamlinks
coveragefolder=$outputfolder/coverage
SSFAfolder=$outputfolder/diagnostic_reads/SSFA
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
logfolder=$outputfolder/logfolders/coverage_and_diagnostic_reads
errorfolder=$outputfolder/errorfolders/coverage_and_diagnostic_reads
allchrom=(2L 2R 3L 3R X)

echo "Calculating coverage and identifying diagnostic reads for samples from ${bamfilepaths} on `date`." >> $outputfolder/added_samples.log

# Create the output folders if necessary
mkdir -p $bamfilefolder
mkdir -p $logfolder
mkdir -p $errorfolder
for chrom in ${allchrom[@]}
do
	mkdir -p $coveragefolder/$chrom/coveragelogs
done
mkdir -p $SSFAfolder/2R/Ace1_region/SSFAlogs
mkdir -p $SSFAfolder/2R/Cyp6_region/SSFAlogs
mkdir -p $SSFAfolder/3R/Cyp6zm_region/SSFAlogs
mkdir -p $SSFAfolder/3R/Gste_region/SSFAlogs
mkdir -p $SSFAfolder/X/Cyp9k1_region/SSFAlogs
mkdir -p $breakpointsfolder/2R/Ace1_region/breakpointlogs
mkdir -p $breakpointsfolder/2R/Cyp6_region/breakpointlogs
mkdir -p $breakpointsfolder/3R/Cyp6zm_region/breakpointlogs
mkdir -p $breakpointsfolder/3R/Gste_region/breakpointlogs
mkdir -p $breakpointsfolder/X/Cyp9k1_region/breakpointlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $bamfilepaths))
echo "This sample set contains ${numbams} bam files." >> $outputfolder/added_samples.log

bsub -J "bamArray[1-$numbams]%$joblimit" \
     -R"select[mem>500] rusage[mem=500]" \
     -M500 \
     -o $logfolder/log_%J.%I.txt \
     -e $errorfolder/error_%J.%I.txt \
     ' '${scriptsfolder}'/get_windowed_coverage_and_diagnostic_reads.sh '${bamfilefolder}' '${bamfilepaths}' ${LSB_JOBINDEX} '$outputfolder


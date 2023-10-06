bamfilepaths=$1
joblimit=$2
study=$3

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10
outputfolder=$rootfolder/$study
bamfilefolder=$outputfolder/bamlinks
SSFAfolder=$outputfolder/diagnostic_reads/SSFA
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
logfolder=$outputfolder/logfolders/coverage_and_diagnostic_reads
errorfolder=$outputfolder/errorfolders/coverage_and_diagnostic_reads

echo "Identifying diagnostic reads from Coeae clusters for samples from ${bamfilepaths} on `date`." >> $outputfolder/added_samples.log

# Create the output folders if necessary
mkdir -p $bamfilefolder
mkdir -p $logfolder
mkdir -p $errorfolder
mkdir -p $SSFAfolder/2L/Coeaexf_region/SSFAlogs
mkdir -p $SSFAfolder/2L/Coeaexg_region/SSFAlogs
mkdir -p $breakpointsfolder/2L/Coeaexf_region/breakpointlogs
mkdir -p $breakpointsfolder/2L/Coeaexg_region/breakpointlogs

# Get the number of bamfiles that need processing
numbams=($(wc -l $bamfilepaths))
echo "This sample set contains ${numbams} bam files." >> $outputfolder/added_samples.log

bsub -J "bamArray[1-$numbams]%$joblimit" \
     -R"select[mem>500] rusage[mem=500]" \
     -M500 \
     -o $logfolder/log_%J.%I.txt \
     -e $errorfolder/error_%J.%I.txt \
     ' '${scriptsfolder}'/get_diagnostic_reads_coeae_vobs.sh '${bamfilefolder}' '${bamfilepaths}' ${LSB_JOBINDEX} '$outputfolder


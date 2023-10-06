study=$1

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10
workingfolder=$rootfolder/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
ncores=4
logfolder=$diagnostic_reads_folder/logfolders/target_regions_analysis
errorfolder=$diagnostic_reads_folder/errorfolders/target_regions_analysis

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_analysis_output_%J.txt \
     -e $errorfolder/target_regions_analysis_error_%J.txt \
     -q long \
     -n $ncores \
     -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" \
     -M10000 \
     $scriptsfolder/target_regions_analysis_vobs_extras_1.sh $workingfolder

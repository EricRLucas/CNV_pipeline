samplelist=$1
species_id_file=$2
metadata_file=$3
study=$4

scriptsfolder=~/scripts/CNV_scripts/scripts
rootfolder=/lustre/scratch126/gsu/team112/personal/el10
workingfolder=$rootfolder/$study
diagnostic_reads_folder=$workingfolder/diagnostic_reads
coverage_variance_file=$workingfolder/coverage/coverage_variance_masked_09_05_all.csv
gene_coordinates_file=$rootfolder/phase3_data/tables/gene_regions.csv
ncores=7
logfolder=$diagnostic_reads_folder/logfolders/target_regions_analysis
errorfolder=$diagnostic_reads_folder/errorfolders/target_regions_analysis

mkdir -p $logfolder
mkdir -p $errorfolder

bsub -o $logfolder/target_regions_analysis_output_%J.txt \
     -e $errorfolder/target_regions_analysis_error_%J.txt \
     -q long \
     -n $ncores \
     -R"select[mem>50000] rusage[mem=50000] span[hosts=1]" \
     -M50000 \
     $scriptsfolder/target_regions_analysis_230911.sh $workingfolder \
                                                      $samplelist \
                                                      $species_id_file \
                                                      $coverage_variance_file \
                                                      $gene_coordinates_file \
                                                      $metadata_file \
                                                      $ncores

workingfolder=$1

scriptsfolder=~/scripts/CNV_scripts/scripts

cd $workingfolder

export R_LIBS_USER="~/R-modules:$R_LIBS_USER"

R-3.6.1 --version

R-3.6.1 --slave -f $scriptsfolder/target_regions_analysis_extras_1.r > target_regions_analysis/target_regions_analysis_extras_1.log 2>&1



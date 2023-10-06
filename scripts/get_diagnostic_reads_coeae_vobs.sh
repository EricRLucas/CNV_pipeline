bamfilefolder=$1
bamfilepaths=$2
samplenum=$3 
outputfolder=$4
scriptsfolder=~/scripts/CNV_scripts/scripts

samplename=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f1))
bampath=($(head -n$samplenum $bamfilepaths | tail -n1 | cut -f2))
bamfile=${bamfilefolder}/${samplename}.bam

echo Downloading bam file for $samplename from $bampath to $bamfile

# Download the required bam
curl ${bampath} > ${bamfile} 2> ${bamfilefolder}/${samplename}_download.log
curl ${bampath}.bai > ${bamfile}.bai 2>> ${bamfilefolder}/${samplename}_download.log

source activate cnv37 

echo Identifying discordant reads
# Get the discordant reads
SSFA_script=$scriptsfolder/SSFA.py
SSFAfolder=$outputfolder/diagnostic_reads/SSFA
python $SSFA_script $bamfile 2L 28490000:28600000 ${SSFAfolder}/2L/Coeaexf_region/${samplename}_Coeaexf_SSFA_output.csv 10 > ${SSFAfolder}/2L/Coeaexf_region/SSFAlogs/${samplename}_Coeaexf_SSFA_output.log 2>&1
python $SSFA_script $bamfile 2L 37250000:37350000 ${SSFAfolder}/2L/Coeaexg_region/${samplename}_Coeaexg_SSFA_output.csv 10 > ${SSFAfolder}/2L/Coeaexg_region/SSFAlogs/${samplename}_Coeaexg_SSFA_output.log 2>&1

# Get the soft clipped reads
breakpoints_script=$scriptsfolder/breakpoint_detector.py
breakpointsfolder=$outputfolder/diagnostic_reads/breakpoints
python $breakpoints_script $bamfile 2L 28490000:28600000 ${breakpointsfolder}/2L/Coeaexf_region/${samplename}_Coeaexf_breakpoints_output 10 > ${breakpointsfolder}/2L/Coeaexf_region/breakpointlogs/${samplename}_Coeaexf_breakpoints_output.log 2>&1
python $breakpoints_script $bamfile 2L 37250000:37350000 ${breakpointsfolder}/2L/Coeaexg_region/${samplename}_Coeaexg_breakpoints_output 10 > ${breakpointsfolder}/2L/Coeaexg_region/breakpointlogs/${samplename}_Coeaexg_breakpoints_output.log 2>&1

echo


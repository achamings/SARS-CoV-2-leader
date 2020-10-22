#! /bin/bash
#=======================================================================
#
#					sars_cov2_leader.sh
#
#					Version 1.0
#			Created by: Anthony Chamings
#			anthony.chamings@deakin.edu.au
#
#				Updated: 22-Oct-2020
#
#=======================================================================

#USER SPECIFIES INPUT BAM FILE, OUTPUT FILE, BED FILE, THREADS, Q VALUE
scriptSource="${BASH_SOURCE[0]%/*}"


function showUsage {
	echo "Measure the coverage from a BAM file over the SARS CoV2 genome" 
	echo "at a specific Q value or higher. Provide a BED file to restrict the coverage"
	echo "to the BED regions, otherwise the coverage across the entire genome will be calculated."
	echo "Will create a SAM and BAM file with reads which contain the partial SARS-CoV2 leader sequence."
	echo "Providing more threads will speed up samtools functions." 
	echo "If the output file exists, the data will be appended to the file."
	echo ""
	echo "sars_cov2_leader.sh -i <bam/sam file> -o <output csv file>"
	echo -e "\t -t <threads> -b <BED file> -q <Q value>" 
	echo -e "\t -r <reference name - default 2019-nCoV>"
	echo -e "\t -n <Name of this sample. Default: inputFile name>"
	exit 0
}

function cleanUp {
	rm -f "$tmpSAM"
	rm -f "$tmpBAM"
	
}

function cleanUpError {
	cleanUp
	rm -r -f "$tmpDIR"
}

trap cleanUpError EXIT








if [ "$#" == 1 ]; then
	#No arguments were supplied
	showUsage
	exit 0
fi



threads=""
outputFile=""
bedFile=""
tmpBAM=""
tmpSAM=""
qScore=0
referenceName="2019-nCoV"
sampleName=""
workingFolder="."
minLength=0
leaderSequence="GTAGATCTGTTCTCT"
while [ "$1" ]; do 

	case "$1" in
		-h|--help|\?) 
			showUsage
			;;
		-i|--input)
			shift
			inputFile="$1"
			;;
			
		-o|--output)
			shift
			outputFile="$1"
			;;
		-t|--threads)
			shift
			threads="$1"
			;;
		-b|--bedfile)
			shift
			bedFile="$1"
			;;
		-q|--quality)
			shift
			qScore="$1"
			;;
		-r|--reference)
			shift
			referenceName="$1"
			;;
		-n|--name)
			shift
			sampleName="$1"
			;;
		-w|--working-folder)
			shift
			workingFolder="$1"
			;;
		-l|--min-length) 
			shift
			minLength="$1"
			;;
		--leader)
			shift
			leaderSequence="$1"
			;;
	esac
	
	shift		
done

#CHECK USER INPUT
if [ ! -f "$inputFile" ]; then
	echo "Error: Cannot find input file"
	exit 1
fi

if [  $threads == "" ]; then
	threads=1
fi

if [ $outputFile == "" ]; then
	echo "Error: No output File nominated."
	exit 1
fi


#MAKE A TEMP DIRECTORY TO SAVE RESULTS IN


tmpDIR=`date +%D-%T-%N | md5sum | grep -o '^.*\b'`

while [ -d "$tmpDIR" ]; do 
	tmpDIR=`date +%D-%T-%N | md5sum | grep -o '^.*\b'`
done 

tmpDIR="$workingFolder/$tmpDIR"
mkdir -p "$tmpDIR"



#IF THE OUTPUT FILE DOES NOT EXIST CREATE A NEW ONE WITH A HEADER LINE
if [ ! -f "$outputFile" ]; then
	echo -n -e "Sample Name\tQ Score\tTotal Mapped Reads\tNo. Reads with Leader\tNon 5'UTR Reads with Leader\tPercentage non-5'UTR Reads with Leader\tAverage Coverage 5'UTR+ORF1ab\tAverage Log10 Coverage 5'UTR+ORF1ab\t" > "$outputFile"
	echo -n -e "Average coverage S Gene-3'UTR\tAverage Log10 coverage S Gene-3'UTR\tCoverage Ratio\tLog Difference\tConverted Log Difference\n" >> "$outputFile"
fi


#CHECK SAMPLE NAME
if [ "$sampleName" == "" ]; then
	# NO SAMPLE NAME GIVEN SO USE FILE NAME
	sampleName=${inputFile##*/}
	sampleName=${sampleName%.*}
fi
echo "Sample will be called: $sampleName"



#DETERMINE IF INPUT FILE IS BAM OR SAM FILE.

inputExt=${inputFile##*.}
inputExt=${inputExt,,}

echo "Input Extension: $inputExt"
 

if [ "$inputExt" == "bam" ]; then
	echo "BAM file supplied as input"
	tmpSAM=`date +%D-%T-%N | md5sum | grep -o '^.*\b'`
	tmpSAM="$tmpDIR"/"$tmpSAM.sam"
	samtools view -O SAM -@ "$threads" -h "$inputFile" > "$tmpSAM"
	workingBAM="$inputFile"
	workingSAM="$tmpSAM"
else
	#SAM FILE
	echo "SAM file supplied as input"
	tmpBAM=`date +%D-%T-%N | md5sum | grep -o '^.*\b'`
	tmpBAM="$tmpDIR/$tmpBAM.bam"
	samtools view -b -@ "$threads" "$inputFile" | samtools sort -@ "$threads" -o "$tmpBAM"
	samtools index "$tmpBAM"
	workingBAM="$tmpBAM"
	workingSAM="$inputFile"
fi 


echo ""
echo "Now calculating coverage over reference $referenceName in $workingBAM at Q $qScore..."
echo ""
if [ "$bedFile" == "" ]; then
	samtools depth -a -Q "$qScore" -m 0 -l "$minLength" -r "$referenceName" "$workingBAM" > "${tmpDIR}/${sampleName}_coverage_Q${qScore}.csv"
else
	#m FLAG IS THE MAXIMUM DEPTH - SET TO 0 OTHERWISE CAPS AT 8000
	
	samtools depth -a -Q "$qScore" -m 0 -l "$minLength" -b "$bedFile" -r "$referenceName" "$workingBAM" > "${tmpDIR}/${sampleName}_coverage_Q${qScore}.csv"
fi


#RUN AWK SCRIPT TO CALCULATE THE AVERAGE COVERAGE BEFORE AND AFTER
coverageString=$($scriptSource/SARS_COV2_coverage_ratio.awk "${tmpDIR}/${sampleName}_coverage_Q${qScore}.csv")

echo "done."


#FIND READS WITH THE LEADER SEQUENCE AT THE SPECIFIC Q VALUE

#PRODUCE A SAM AND BAM FILE WITH ONLY LEADER READS MAPPED AT Q 
echo ""
echo "Making Mapped BAM file with reads containing leader at Q $qScore..."
echo ""

if [ "$minLength" -lt 0 ]; then
	leaderReadFileName="${tmpDIR}/${sampleName}_leader_reads_Q${qScore}_len${minLength}"
else 
	leaderReadFileName="${tmpDIR}/${sampleName}_leader_reads_Q${qScore}"

	#leaderReadFileName="${tmpDIR}/${sampleName}_leader_reads"
fi


$scriptSource/identify_SARS_leader_reads_at_Q.awk "$qScore" "$referenceName" "$minLength" "$leaderSequence" "$workingSAM"  > "$leaderReadFileName.sam"

samtools view -b -@ threads "$leaderReadFileName.sam" | samtools sort -@ $threads -o "$leaderReadFileName.bam"
samtools index "$leaderReadFileName.bam"
echo -e "\ndone."
echo -e "\nCalculating number of reads with leader sequence...."

readsWithLeaderCount=$(awk '{ if ( $0 !~ /^@/ && $5 >= "$qScore" ) { print $0; } }' "$leaderReadFileName.sam" | wc -l) 
echo -e "\ndone."


echo -e "\nCalculating total number of mapped reads at Q $qScore..."

if [ "$bedFile" != "" ]; then
	#NO BED FILE SET
	
	mappedReadsAtQ=$(samtools view -O SAM -q $qScore -L $bedFile -@ $threads $workingSAM | wc -l )
else
	mappedReadsAtQ=$(samtools view -O SAM -q $qScore -@ $threads $workingSAM | wc -l )
fi
echo ""
echo "Found $readsWithLeaderCount reads with leader at Q $qScore." 
echo "Found $mappedReadsAtQ reads at Q of $qScore or higher mapped to bed File regions"

#FIND LEADER READS NOT BELONGING TO 5' UTR

non5UTRLeaderReads=$(awk '{ if ( $0 !~ /^@/) { if ( $4 > 100 && $5 >= "$qScore" ) { print $0; }  } }' "$leaderReadFileName.sam" | wc -l)

percentageLeaderReads=$(echo "scale=5;${non5UTRLeaderReads}/${mappedReadsAtQ}*100" | bc -l)



echo -e "$sampleName\t$qScore\t$mappedReadsAtQ\t$readsWithLeaderCount\t$non5UTRLeaderReads\t$percentageLeaderReads\t$coverageString" >> "$outputFile"



#CREATE AN OUTPUT DIRECTORY FOR THE DATA FILES AND COPY ALL DATA FILES ACROSS

outputDir="${sampleName}_leader_data"


cleanUp	#REMOVE THE BAM AND SAM TMP FILES
mv "$tmpDIR" "$outputDir"



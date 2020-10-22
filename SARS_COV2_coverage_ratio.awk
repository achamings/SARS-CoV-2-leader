#!/usr/bin/awk -f
#CALCULATES THE AVERAGE COVERAGE UP TO THE END OF THE ORF1AB GENE 
#THIS ENDS AT POSITION 21555 ON the 2019-nCoV


BEGIN {
	FS="\t";
	OFS="\t";
	#sampleName=ARGV[1];
	#delete ARGV[1];
}
{
	if ( $0 ~ /^\s*$/ )  {
	
		#SKIP BLANK LINES
		next;
	}	
	  
	#	$1 IS THE GENOME NAME 2019-nCoV
	#	$2 IS THE GENOME POSITION
	#	$3 IS THE NUCLEOTIDE COVERAGE
	if ( $2 <= 21555 ) {
		
		#print "ORF1ab",$3;
		orf1abCoverage[$2]=$3;
		orf1abLogCoverage[$2]=log($3+1)/log(10);
		#print orf1abCoverage[$2];
	} else {
		
		#print "Other",$3;
		remainingGenomeCoverage[$2]=$3;
		remainingGenomeLogCoverage[$2]=(log($3+1)/log(10));
	}
	#print $0
}
END {
	averageORF1ab=0;
	ORF1abCount=0;
	ORF1abTotal=0;
	ORF1abLogTotal=0;
	averageRemaining=0;
	remainingCount=0;
	remainingTotal=0;
	remainingLogTotal=0;
	
	
	for (cov1 in orf1abCoverage) {
		ORF1abCount++;
		#print ORF1abCount;
		#print orf1abCoverage[cov1];
		ORF1abTotal=ORF1abTotal+orf1abCoverage[cov1];

	}
	
	
	for (covLog1 in orf1abLogCoverage) {
		ORF1abLogTotal=ORF1abLogTotal+orf1abLogCoverage[covLog1];
	}
	
	
	
	
	for (cov2 in remainingGenomeCoverage) {
		remainingCount++;
		remainingTotal=remainingTotal+remainingGenomeCoverage[cov2];
	}
	
	
	for (covLog2 in remainingGenomeLogCoverage) {
		remainingLogTotal=remainingLogTotal+remainingGenomeLogCoverage[covLog2];
	}
	
	averageORF1ab=ORF1abTotal/ORF1abCount;
	averageRemaining=remainingTotal/remainingCount;
	
	NaturalRatio=averageRemaining/averageORF1ab;
	
	
	
	ORF1abLogAverage=ORF1abLogTotal/ORF1abCount;
	RemainingLogAverage=remainingLogTotal/remainingCount;
	
	logDifference=RemainingLogAverage-ORF1abLogAverage;
	
	convertedLogDifference=10^logDifference;
	
	
	#print sampleName,averageORF1ab,ORF1abLogAverage,averageRemaining,RemainingLogAverage,NaturalRatio,LogDifference,convertedLogDifference;
	print averageORF1ab,ORF1abLogAverage,averageRemaining,RemainingLogAverage,NaturalRatio,logDifference,convertedLogDifference;

}

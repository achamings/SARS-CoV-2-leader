#! /usr/bin/awk -f

#UPDATED 12 June 2020 TO ALLOW FOR ANY SEQUENCE PATTERN TO BE SEARCH FOR

#FIND READS CONTAINING THE LEADER SEQUENCE
#	 e.g.	GTAGATCTGTTCTCT for SARS-CoV2


BEGIN {
	qScore=ARGV[1];
	refName=ARGV[2];
	minLength=ARGV[3];
	leaderSeq=ARGV[4];
	delete ARGV[1];
	delete ARGV[2];
	delete ARGV[3];
	delete ARGV[4];
	
	totalReads=0;
	leaderReads=0;
	OFS="\t"
}
{
	
	if ( $0 !~ /^@/ ) {
		if ( $10 ~ leaderSeq && $5 >= qScore && $3 == refName && length($10) >= minLength ) {
		#if ( $10 ~ /GTAGATCTGTTCTCT/ && $5 >= qScore && $3 == refName && length($10) >= minLength ) {
			#if ( $10 ~ /GTAGATCTGTTCTCT/ && $3 == refName) {
				#LEADER SEQUENCE FOUND
				totalReads++;
				
				
						 
				
				#LABEL THE POSITION IN THE READ NAME
				leaderPos=match($10,leaderSeq)
				$1=$1 "_L:" leaderPos
				
				#CHECK TO SEE IF THE LEADER IS AT THE 5' END OF THE READ
				#JUST PRIOR TO THE START OF THE MAPPED REGION OR IS IN
				#THE START OF THE MAPPED REGION (SOME LEADERS ARE POORLY
				#MAPPED WITHIN A  READ	
					
				split($6,CIGARnums,/[A-Za-z]/);
				split($6,CIGARchars,/([0-9])+/);
				if (CIGARchars[2]=="S") {
					softClip=CIGARnums[1];
					
					if ( (softClip-(leaderPos+15)) < 20 ) {
					
						$1=$1 "S"
						print $0;
						leaderReads++;
						
					}			
							 
				}
				if (CIGARchars[2]=="M") {
					# READ IS MAPPED FROM THE START
					if (leaderPos < 40) {
						$1=$1 "M"
						print $0;
						leaderReads++;
						$1=$1 "M"
					} 
				}
				
				
		}
	} else {
		print $0;
		totalReads++;
	}
	
}

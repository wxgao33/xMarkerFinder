#!/bin/bash
helpdoc(){
	cat <<EOF
Description:

	This is a help document
	- Calculate Inter-Microbiota Correlation

Usage:
	
	$0 -w <workplace> -i <input microbiota file> -o <output file> -t <threads>

Option:
	-w workplace
	-i input microbiota file
	-o output file
	-t threads
EOF
}
#getopts

if [ "$1" == "-h" ]
then
	helpdoc
	exit 1
fi

if [ $# == 0 ]
then
	helpdoc
	exit 1
fi

while getopts ":W:i:o:t:" opt
do
    case $opt in
    	W)
		workplace=`echo $OPTARG`
		;;
        i)
        inputfile=`echo $OPTARG`
        ;;
        o)
		outputfile=`echo $OPTARG`
		;;
		t)
		threads=`echo $OPTARG`
		;;
        ?)
        echo "Undefined parameter"
        exit 1;;
    esac
done


fastspar --threads $threads --otu_table $workplace$inputfile --correlation ${workplace}median_correlation.tsv --covariance ${workplace}median_covariance.tsv

mkdir ${workplace}bootstrap_feat
fastspar_bootstrap --otu_table $workplace$inputfile --number 1000 --prefix ${workplace}bootstrap_feat/feat

mkdir ${workplace}$outputfile
parallel fastspar --otu_table $workplace$inputfile --correlation $workplace$outputfile/cor_{/} --covariance $workplace$outputfile/cov_{/} -i 5 ::: ${workplace}bootstrap_feat/*

fastspar_pvalues --otu_table $workplace$inputfile --correlation ${workplace}median_correlation.tsv --prefix $workplace$outputfile/cor_ --permutations 1000 --outfile ${workplace}pvalues.tsv

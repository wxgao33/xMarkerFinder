#!/bin/bash
helpdoc(){
	cat <<EOF
Description:

	This is a help document
	- Calculate Multi-omics Correlation

Usage:
	
	$0 -w <workplace> -i <input microbiota file1> -d <input microbiota file2> -o <output file>

Option:

	-w workplace
	-i input microbiota file 1
	-d input microbiota file 2
	-o output file prefix

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

while getopts ":w:i:d:o:" opt
do
    case $opt in
    	w)
		workplace=`echo $OPTARG`
		;;
        i)
        inputfile1=`echo $OPTARG`
        ;;
        d)
        inputfile2=`echo $OPTARG`
        ;;
        o)
		output=`echo $OPTARG`
		;;
        ?)
        echo "Undefined parameter"
        exit 1;;
    esac
done

halla -x $workplace$inputfile1 -y $workplace$inputfile2 -o $workplace$output -m spearman --fdr_alpha 0.05


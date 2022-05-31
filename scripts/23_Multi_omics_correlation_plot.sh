#!/bin/bash
helpdoc(){
    cat <<EOF
Description:

    This is a help document
    - Plot Multi-omics Correlation

Usage:
    
    $0 -i <input result file>
Option:

    -i input result file

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

while getopts ":W:i:" opt
do
    case $opt in
        W)
        workplace=`echo $OPTARG`
        ;;
        i)
        inputfile=`echo $OPTARG`
        ;;
        ?)
        echo "Undefined parameter"
        exit 1;;
    esac
done

hallagram \
    -i $workplace$inputfile/ \
    --x_dataset_label 'Signature1' \
    --y_dataset_label 'Signature2' \
    --block_border_width 2 \
    -c
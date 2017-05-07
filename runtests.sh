#!/bin/bash

usage() { 
echo "Usage: $0 -p <string> [-c <string>]" 1>&2; 
echo " -p <prog> parallel prog name" 1>&2;
echo " -c coprocessor name" 1>&2;
exit 1; 
}

mic=mic1

while getopts ":p:c:" o; do
    case "${o}" in
        p)
            progpar=${OPTARG}
            ;;
        c)
            mic=${OPTARG}
            ;;    

        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

# check if the setted params
if [ -z "$progpar" ]; then echo "-p is missing"; usage ; else make $progpar ; fi

file="$progpar.$mic.log"
if [ ! -e "$file" ] ; then
         # if not create the file
         touch $file
fi

echo "file copying "$mic
scp $progpar $mic:

echo "PAR tests con $progpar"

# for (( i=8; i <= 128; i*=2 ))
for iter in 10
do

#for w in 1 10 20 30 40 50 58 59 60 61 80 90 100 110 117 118 119 120 130 140 150 160 176 177 178 179 190 200 210 220 230 235 236 237 239
for w in {1..239..1} 
do 
# 60.mm 64.mm 128.mm 250.mm 256.mm 500.mesh 512.m 124.m 2k.mesh 4k.mesh 5k.mesh 10k.mesh 15k.mesh
for mesh in 512.m 124.m 2k.mesh 4k.mesh 5k.mesh 10k.mesh
do
	stringa=$(ssh $mic './'$progpar -w $w -i $iter -m $mesh)
	echo ${stringa/%(ms)} >>$file
	echo $stringa
done
done
done

ssh $mic rm $progpar  

echo "finished "$file
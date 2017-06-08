#!/bin/bash
# $1: GraphletAlign base folder
# $2: input folder
# $3: name
# $4: number of copies
# $5: min size of graphlet
# $6: max size of graphlet
# $7: file that stores # accumulated proteins in each species
# $8: output folder
echo "DIR=$1"
JAVA=java
Y=$5
Z=$6
params=$#              # Number of command-line parameters.
while [ "$Y" -le "$Z" ]
do
  echo "$JAVA -d64 -Xms1024M -Xmx6000M -cp \$DIR/class GraphletAlign $2/$3.hits $2/$3.adjl $2/$3.blastp.no \$DIR/data/$Y.adj \$DIR/data/$Y.auto $4 $7 $8"

  Y=$(($Y+1))
done

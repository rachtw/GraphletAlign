#!/bin/bash
# $1: config file

GRAPHLET_ALIGN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Graphlet alignment folder: $GRAPHLET_ALIGN_DIR"

DATA_DIR="$( cd "$( dirname "$1" )" && pwd )"
echo "Data folder: $DATA_DIR"

sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\r//g' $1 # remove CR at end of line
. $1
echo "Read configure file: $1"

COMBINED=${SPECIES[0]}
Y=1
while [ "$Y" -lt "${#SPECIES[@]}" ]
do
  COMBINED="$COMBINED"_"${SPECIES[$Y]}"
  Y=$(($Y+1))
done

cd $DATA_DIR
mv $CLUSTER_DIR "$CLUSTER_DIR"_bak # backup the output folder
rm -rf $CLUSTER_DIR
mkdir $CLUSTER_DIR
cd $TMP_DIR

if [ "$MULTIPLE_SPECIES" -eq 1 ]; then
  CMD="$GRAPHLET_ALIGN_DIR/bin/AlignmentClustering $COMBINED.name ../$CLUSTER_DIR $GRAPHLET_ALIGN_DIR/data false $CLUSTER_CUTOFF $CLUSTER_SIZE_PER_SPECIES_CUTOFF "
else
  CMD="$GRAPHLET_ALIGN_DIR/bin/AlignmentClustering $COMBINED.name ../$CLUSTER_DIR $GRAPHLET_ALIGN_DIR/data true $CLUSTER_CUTOFF $CLUSTER_SIZE_PER_SPECIES_CUTOFF "
fi

Y=$GRAPHLET_MIN_SIZE
while [ "$Y" -le "$GRAPHLET_MAX_SIZE" ]
do
  CMD=$CMD`ls $DATA_DIR/$ALIGN_DIR/*.$Y-*`" "
  Y=$(($Y+1))
done
echo $CMD
$CMD

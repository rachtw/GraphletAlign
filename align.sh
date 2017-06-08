#!/bin/bash
# $1: config file

GRAPHLET_ALIGN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Graphlet alignment folder: $GRAPHLET_ALIGN_DIR"

DATA_DIR="$( cd "$( dirname "$1" )" && pwd )"
echo "Data folder: $DATA_DIR"

sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\r//g' $1 # remove CR at end of line
. $1
echo "Read configure file: $1"

PWD=`pwd`
cd $DATA_DIR
BLASTP_FILE=`ls *blastp`
cd $PWD
echo "Blastp file: $BLASTP_FILE"

COMBINED=${SPECIES[0]}
echo "Net file: ${SPECIES[0]}.net"
Y=1
while [ "$Y" -lt "${#SPECIES[@]}" ]
do
  COMBINED="$COMBINED"_"${SPECIES[$Y]}"
  echo "Net file: ${SPECIES[$Y]}.net"
  Y=$(($Y+1))
done

cd $DATA_DIR
mv $ALIGN_DIR "$ALIGN_DIR"_bak # backup the output folder
rm -rf $ALIGN_DIR
mv $TMP_DIR "$TMP_DIR"_bak # backup the temparory folder
rm -rf $TMP_DIR
mkdir $TMP_DIR
mkdir $ALIGN_DIR
cd $TMP_DIR

Y=0
while [ "$Y" -lt "${#SPECIES[@]}" ]
do
  #echo "java -cp $GRAPHLET_ALIGN_DIR/class ParseNetwork ../${SPECIES[$Y]}.net ${SPECIES[$Y]}.adjl ${SPECIES[$Y]}.name"
  java -cp $GRAPHLET_ALIGN_DIR/class ParseNetwork ../${SPECIES[$Y]}.net ${SPECIES[$Y]}.adjl ${SPECIES[$Y]}.name
  Y=$(($Y+1))
done

if [ "$MULTIPLE_SPECIES" -eq 1 ]; then
  echo -n "" > $COMBINED.name
  Y=0
  while [ "$Y" -lt "${#SPECIES[@]}" ]
  do
    cat ${SPECIES[$Y]}.name >> $COMBINED.name
    Y=$(($Y+1))
  done
  #echo $COMBINED.name


  CMD="java -cp $GRAPHLET_ALIGN_DIR/class CombineAdjLists 0 "
  Y=0
  while [ "$Y" -lt "${#SPECIES[@]}" ]
  do
    CMD="$CMD${SPECIES[$Y]}.adjl "
    Y=$(($Y+1))
  done
  CMD="$CMD$COMBINED.adjl"
  #echo "$CMD"
  $CMD
fi
#echo "$GRAPHLET_ALIGN_DIR/bin/AdjList2 $COMBINED.adjl $COMBINED.adjl.2"
$GRAPHLET_ALIGN_DIR/bin/AdjList2 $COMBINED.adjl $COMBINED.adjl.2

#echo "sort ../$BLASTP_FILE > $COMBINED.blastp.sorted"
sort ../$BLASTP_FILE > $COMBINED.blastp.sorted
#echo "$GRAPHLET_ALIGN_DIR/bin/SortBlastp $COMBINED.blastp.sorted"
$GRAPHLET_ALIGN_DIR/bin/SortBlastp $COMBINED.blastp.sorted

CMD="$GRAPHLET_ALIGN_DIR/bin/BlastpTopHits "
Y=0
while [ "$Y" -lt "${#SPECIES[@]}" ]
do
  CMD="$CMD${SPECIES[$Y]}.name "
  Y=$(($Y+1))
done
CMD="$CMD$COMBINED.blastp.sorted.sorted $COMBINED.blastp.no $COMBINED.hits $BLASTP_EXP $TOP_HITS $BIDIRECTIONAL_HITS"
#echo "$CMD"
$CMD

echo -n "" > species_size
SUM_PROTEINS=0
Y=0
while [ "$Y" -lt "${#SPECIES[@]}" ]
do
  NUM_PROTEINS=`wc -l $DATA_DIR/$TMP_DIR/${SPECIES[$Y]}.name | cut -f 1 -d ' '`
  SUM_PROTEINS=$(($SUM_PROTEINS+$NUM_PROTEINS))
  echo $SUM_PROTEINS >> species_size
  Y=$(($Y+1))
done

if [ "$MULTIPLE_SPECIES" -eq 1 ]; then
  bash $GRAPHLET_ALIGN_DIR/script/GraphletAlignGen.sh $GRAPHLET_ALIGN_DIR . $COMBINED ${#SPECIES[@]} $GRAPHLET_MIN_SIZE $GRAPHLET_MAX_SIZE species_size $DATA_DIR/$ALIGN_DIR > GraphletAlign.sh
else
  bash $GRAPHLET_ALIGN_DIR/script/GraphletAlignGen.sh $GRAPHLET_ALIGN_DIR . $COMBINED 2 $GRAPHLET_MIN_SIZE $GRAPHLET_MAX_SIZE species_size $DATA_DIR/$ALIGN_DIR > GraphletAlign.sh
fi

#echo "bash GraphletAlign.sh"
bash GraphletAlign.sh

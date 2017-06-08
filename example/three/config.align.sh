#!/bin/bash
MULTIPLE_SPECIES=1 # 0: false; 1: true
TOP_HITS=35
BIDIRECTIONAL_HITS=1 # 0: false; 1: true
GRAPHLET_MAX_SIZE=4
SPECIES[0]=fly
SPECIES[1]=worm
SPECIES[2]=yeast

#------------------------default values
GRAPHLET_MIN_SIZE=3
BLASTP_EXP=1e-6

#------------------------system configuration
ALIGN_DIR=alignment
TMP_DIR=tmp
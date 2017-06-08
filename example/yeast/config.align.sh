#!/bin/bash
MULTIPLE_SPECIES=0 # 0: false; 1: true
TOP_HITS=10
BIDIRECTIONAL_HITS=1 # 0: false; 1: true
GRAPHLET_MAX_SIZE=4
SPECIES[0]=yeast

#------------------------default values
GRAPHLET_MIN_SIZE=3
BLASTP_EXP=1e-6

#------------------------system configuration
ALIGN_DIR=alignment
TMP_DIR=tmp
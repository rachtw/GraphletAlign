#!/bin/bash
MULTIPLE_SPECIES=1 # 0: false; 1: true
TOP_HITS=35
BIDIRECTIONAL_HITS=1 # 0: false; 1: true
GRAPHLET_MAX_SIZE=3
SPECIES[0]=Ecol
SPECIES[1]=Hpyl
SPECIES[2]=Styp
SPECIES[3]=Vcho

#------------------------default values
GRAPHLET_MIN_SIZE=3
BLASTP_EXP=1e-6

#------------------------system configuration
ALIGN_DIR=alignment
TMP_DIR=tmp
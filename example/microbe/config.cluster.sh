#!/bin/bash
MULTIPLE_SPECIES=1 # 0: false; 1: true
GRAPHLET_MAX_SIZE=3
SPECIES[0]=Ecol
SPECIES[1]=Hpyl
SPECIES[2]=Styp
SPECIES[3]=Vcho

#------------------------default values
GRAPHLET_MIN_SIZE=3
CLUSTER_CUTOFF=0.5
CLUSTER_SIZE_PER_SPECIES_CUTOFF=100

#------------------------system configuration
CLUSTER_DIR=cluster
ALIGN_DIR=alignment
TMP_DIR=tmp
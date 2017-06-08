#!/bin/bash
MULTIPLE_SPECIES=0 # 0: false; 1: true
GRAPHLET_MAX_SIZE=4
SPECIES[0]=yeast

#------------------------default values
GRAPHLET_MIN_SIZE=3
CLUSTER_CUTOFF=0.5
CLUSTER_SIZE_PER_SPECIES_CUTOFF=100

#------------------------system configuration
CLUSTER_DIR=cluster
ALIGN_DIR=alignment
TMP_DIR=tmp
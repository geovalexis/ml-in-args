#!/bin/bash

# Config variables
input_gff=$1
output_bed=$2

awk -F'\t' '$3 == "gene" {split($9, a, ";"); for (i in a) {if (a[i] ~ /^ID=/) {split(a[i], b, "="); print $1 "\t" $4 "\t" $5 "\t" b[2]}}}' \
$input_gff > $output_bed
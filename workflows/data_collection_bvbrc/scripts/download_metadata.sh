#!/bin/bash

# Config variables
query=$1
output=$2

# Download metadata from BV-BRC servers
curl -K -OJ 'https://www.bv-brc.org/api-for-website/genome_amr/?&http_accept=text/csv&http_download=true' \
  --data-urlencode "rql=$query" \
  --compressed -o "$output"
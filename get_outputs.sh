#!/bin/bash
output_file="output.txt"
> "$output_file"

for i in {1..30}; do
  row=$("$@" | tr -s '\t' '\t' | tr -d '\n')
  echo "$row" | tr '.' ',' >> "$output_file"
done
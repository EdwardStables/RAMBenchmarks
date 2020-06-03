#!/bin/bash
{
for filename in ~/Documents/julia_benchmarks/RAMBenchmarks/results/*.csv
do
    printf ${filename}
    printf "\n"
    csvlook ${filename} -I
    printf "\n"
done
} | less -S

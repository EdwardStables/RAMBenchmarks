#!/bin/bash
{
for filename in ./*.csv
do
    printf ${filename}
    printf "\n"
    #csvlook ${filename} -I
    column -s, -t ${filename}
    printf "\n"
done
} | less -S

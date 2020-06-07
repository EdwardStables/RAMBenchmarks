#!/bin/bash
initial_threading=$JULIA_NUM_THREADS

#exit if test file doesn't exist
[[ -e "$1" ]] || exit 1

#$1 : test file
#$2 : output file
#$3 : threaded 
#$4 : if threaded, how many iterations
#$5 : if threaded, max number of threads
#$6 : comp, run the comparison solvers

if [[ "$3" = "threaded" ]]; then
    export JULIA_NUM_THREADS=1
    julia1.4 --startup-file=no -O -- ./top_julia.jl "$4" "10" "$4" "threaded" "$6" $1 $2
    for ((t=5; t <= $5; t=t+5)); do
        export JULIA_NUM_THREADS=$t
        julia1.4 --startup-file=no -O -- ./top_julia.jl "$4" "10" "$4" "threaded" "$6" $1 $2
    done
else
    julia1.4 --startup-file=no -O -- ./top_julia.jl 10 10 50 "not_threaded" "$6" $1 $2
fi



#Reset number of threads to that at the start
export JULIA_NUM_THREADS=$initial_threading

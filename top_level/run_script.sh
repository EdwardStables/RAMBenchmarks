#!/bin/bash
initial_threading=$JULIA_NUM_THREADS

#exit if test file doesn't exist
[[ -e "$1" ]] || exit 1

#$1 : test file
#$2 : output file
#$3 : threaded 
#$4 : min iterations
#$5 : iteration step
#$6 : max iterations
#$7 : if threaded, max number of threads
#$8 : comp, run the comparison solvers

echo "Starting top level..."

if [[ "$3" = "threaded" ]]; then
    echo "Running threaded"
    export JULIA_NUM_THREADS=1
    julia1.4 --startup-file=no -O -- ./top_julia.jl "$4" "$5" "$6" "threaded" "$8" $1 $2
    for ((t=5; t <= $7; t=t+5)); do
        export JULIA_NUM_THREADS=$t
        julia1.4 --startup-file=no -O -- ./top_julia.jl "$4" "$5" "$6" "threaded" "$8" $1 $2
    done
else
    echo "Running not-threaded"
    julia1.4 --startup-file=no -O -- ./top_julia.jl "$4" "$5" "$6" "not_threaded" "$8" $1 $2
fi



#Reset number of threads to that at the start
export JULIA_NUM_THREADS=$initial_threading

using RAMBenchmarks

function main(args)
    rng = parse.(Int64, args[1:3])
    rng = rng[1]:rng[2]:rng[3]

    threads = args[4] == "threaded"

    isfile(args[end-1]) || error("file")
    in_file = args[end-1] 
    out_file = args[end] 

    #warm up the compiler
    run_ram("HS21")
    

    tests = []

    open(in_file, "r") do f
        for test in readlines(f)
            push!(tests, test)
        end
    end

    for test in tests
        run_range(test, [rng...], threads; out_file=out_file)
    end
end

main(ARGS)

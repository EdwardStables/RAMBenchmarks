using RAMBenchmarks

function main(args)
    rng = parse.(Int64, args[1:3])
    rng = rng[1]:rng[2]:rng[3]

    threads = args[4] == "threaded"
    comp = args[5] == "comp"

    isfile(args[end-1]) || error("file")
    in_file = args[end-1] 
    out_file = args[end] 

    #warm up the compiler
    bm = Benchmark("HS21", threads)
    run_ram_range(bm, [rng...], false)
    

    tests = []

    open(in_file, "r") do f
        for test in readlines(f)
            push!(tests, test)
        end
    end

    for test in tests
        bm = Benchmark(test, threads, out_file=out_file)
        if comp
            run_OSQP(bm, true)
            #run_COSMO(bm, true)
            run_Ipopt(bm, true)
        else
            run_ram_range(bm, [rng...], true)
        end
    end
end

main(ARGS)

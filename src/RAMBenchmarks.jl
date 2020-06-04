module RAMBenchmarks
export run_benchmarks, run_osqp, run_cosmo,run_ram, run_ipopt, RAMB, run_range

using RowActionMethods
using QPSReader
using SparseArrays
using COSMO
using JuMP
using OSQP
using Ipopt
using CSV
using Tables
using Dates

const RAMB = RAMBenchmarks

global valid_problems = ["AUG2D", "AUG2DC", "AUG2DCQP", "AUG2DQP", "AUG3D", "AUG3DC",
						 "AUG3DCQP", "AUG3DQP", "BOYD1", "CONT-050", "CONT-100", "CONT-101",
						 "CONT-200", "CONT-201", "CONT-300", "HS118", "HS21", "HUES-MOD",
						 "HUESTIS", "KSIP", "LISWET1", "LISWET10", "LISWET11", "LISWET12",
						 "LISWET2", "3LISWET3", "LISWET4", "LISWET5", "LISWET6", "LISWET7",
						 "LISWET8", "LISWET9", "POWELL20", "QPCBLEND", "QPCBOEI1", "QPCBOEI2",
						 "QPCSTAIR", "YAO"]
#known good: HS21, HS118, , LISWET1, LISWET2, AUG2DC, 
#poor performance: CONT-100
#running problems: AUG2D,

function run_range(f::String, iters::Vector{Int}, threads::Bool)
    #run_osqp(f)
    #run_ipopt(f)
    for i in iters
        run_ram(f;iterations=i,threads=threads, write=true)
    end
end

#TODO sense (pr.objsense)
run_benchmarks(i::Int) = run_benchmarks(valid_problems[i])
function run_benchmarks(f::String; threads::Bool=false, iterations::Int=100)
	pr = get_problem_file(f)
	Q,b,c,M,d = get_matrices(pr)
    p = form_problem(Q,b,M,d; threads=threads)
    Optimize(p, SC_Iterations(iterations))
    RAM.Resolve(p)
    @show p.statistics
    return p
end

function write_csv(test::String, solver::String, notes::String, overhead_times::Vector{T}, runtimes::Vector{T}, optimum; file::String = "result.csv") where {T<:Number}


    cols = ["Time", "Test", "Solver", "notes", "file load", "getting matrices", "JuMP build", "build (ram only)", "optimisation", "total time", "optimum value"]
    if solver == "RAM"
        row = [Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") test solver notes overhead_times... runtimes[1] runtimes[2] sum([overhead_times..., runtimes...]) optimum]
    else                                                                    
        row = [Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") test solver notes overhead_times... "" runtimes[1] sum([overhead_times..., runtimes...]) optimum]
    end
    file = "/home/ed/Documents/julia_benchmarks/RAMBenchmarks/results/"*file
    new = !isfile(file)

    open(file, new ? "w" : "a") do f
        CSV.write(f, Tables.table(row), header=cols, append=!new)
    end

end

run_osqp(f::String; notes::String="",file::String="results.csv", write=false) = run_jump(f, Model(OSQP.Optimizer), "OSQP", notes, file, write)
run_cosmo(f::String; notes::String="",file::String="results.csv", write=false) = run_jump(f, Model(COSMO.Optimizer), "COSMO", notes, file, write)
run_ipopt(f::String; notes::String="",file::String="results.csv", write=false) = run_jump(f, Model(Ipopt.Optimizer), "IPOPT", notes, file, write)

function run_ram(f::String; threads::Bool=false, notes::String="", file::String="results.csv", iterations::Int=100, write=false)
    p = Model(()->RAM.Optimizer("Hildreth"))
    set_optimizer_attribute(p, "iterations", iterations)
    set_optimizer_attribute(p, "threading", threads)

    combined_notes = threads ? "multi-threaded" : "single-threaded" 
    combined_notes *= "; $(iterations) iterations"
    if notes != ""
        combined_notes *= "; "
        combined_notes *= notes
    end

    run_jump(f, p, "RAM", combined_notes, file, write)
end

function run_jump(f::String, p, solver, notes, file, write)
    times = Vector{Float64}(); temp_t = time()

    pr = get_problem_file(f)
    push!(times, time() - temp_t); temp_t = time()

	Q,b,c,M,d = get_matrices(pr)
    try
        set_silent(p)
    catch 
        println("Can't set silent")
    end
    push!(times, time() - temp_t); temp_t = time()
    @variable(p, x[1:pr.nvar])
    @objective(p, Min, 0.5sum(x[i]Q[i,j]x[j] for i=1:pr.nvar,j=1:pr.nvar) + sum(x[i]b[i] for i=1:pr.nvar))
    for c in 1:pr.ncon
        @constraint(p, sum(x[i]M[c,i] for i=1:pr.nvar) <= d[c])
    end
    push!(times, time() - temp_t); temp_t = time()

    optimize!(p)

    if solver == "RAM"
        stats = p.moi_backend.optimizer.model.inner_model.statistics
        op_time = [stats.BuildTime, stats.OptimizeTime]
        @show stats
    else
        op_time = [time() - temp_t]
    end
    

    if write
        write_csv(f, solver, notes, times, op_time, objective_value(p); file=file)
    end

    return p, x
end


function get_problem_file(name::String; path::String="/home/ed/Documents/julia_benchmarks/CUTEst/marosmeszaros",ext::String=".SIF")
    problem_path = joinpath(path, name) * ext
    return readqps(problem_path)
end

function get_matrices(pr::QPSData)
	Q = sparse(pr.qcols, pr.qrows, pr.qvals)
	b = pr.c
	c = pr.c0
	
    le = Vector{Int}()
    ge = Vector{Int}()
	con = sparse(pr.arows, pr.acols, pr.avals, pr.ncon, pr.nvar)
	for (i, (l,u)) in enumerate(zip(pr.lcon, pr.ucon))
		if u != Inf
			push!(le, i)
		end
		if l != -Inf
			push!(ge, i)		
		end
	end
    
	M = spzeros(length(le)+length(ge), pr.nvar)		
    d = spzeros(length(le)+length(ge))
	
	for i in le
		M[i,:] = con[i,:]
        d[i] = pr.ucon[i]
	end

    for i in ge
        M[i+length(le),:] = -con[i,:]
        d[i+length(le)] = -pr.lcon[i]
    end
	
    return Q, b, c, M, d	
end

function form_problem(Q, b, M, d; method::String="Hildreth", threads=false)
    p = GetModel("Hildreth")
    SetThreads(p, threads=threads)
    Setup(p, Matrix(Q), b)
    for (i,l) in enumerate(d)
        AddConstraint(p, M[i,:], l)
    end
    return p
end
end


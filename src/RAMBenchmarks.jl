module RAMBenchmarks

export Benchmark, run_ram_range, run_OSQP, run_COSMO, run_Ipopt

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

const valid_problems = ["AUG2D", "AUG2DC", "AUG2DCQP", "AUG2DQP", "AUG3D", "AUG3DC",
			 "AUG3DCQP", "AUG3DQP", "BOYD1", "CONT-050", "CONT-100", "CONT-101",
			 "CONT-200", "CONT-201", "CONT-300", "HS118", "HS21", "HUES-MOD",
			 "HUESTIS", "KSIP", "LISWET1", "LISWET10", "LISWET11", "LISWET12",
			 "LISWET2", "3LISWET3", "LISWET4", "LISWET5", "LISWET6", "LISWET7",
			 "LISWET8", "LISWET9", "POWELL20", "QPCBLEND", "QPCBOEI1", "QPCBOEI2",
             "QPCSTAIR", "YAO"]

#known good: HS21, HS118, , LISWET1, LISWET2, AUG2DC, 
#poor performance: CONT-100
#running problems: AUG2D,

mutable struct Benchmark
    file::String

    base_problem::Model
    active_problem::Model

    jump_time::Float64


    #Stores the details of the last run
    build_time::Float64
    optimisation_time::Float64
    iterations::Int

    solver::String
    optimum::Float64

    threading::Bool
    threads::Int

    write_file::String
    notes::String

    function Benchmark(test::String, threading::Bool; out_file="results.csv", notes="")
        b = new()

        b.file = test
        b.threading = threading
        b.threads = Base.Threads.nthreads()
        b.write_file = out_file
        b.notes = notes

        t1 = time()
        b.base_problem = build_jump(test)
        b.jump_time = time() - t1

        return b
    end
end

function run_ram_range(bm::Benchmark, iterations::Vector{Int}, write::Bool)
    for i in iterations
        bm.iterations = i
        bm.solver = "RAM"
        bm.active_problem = copy(bm.base_problem)
        
        set_optimizer(bm.active_problem, ()->RAM.Optimizer("Hildreth"))
        optimize!(bm.active_problem)

        bm.build_time = bm.active_problem.moi_backend.optimizer.model.inner_model.statistics.BuildTime
        bm.optimisation_time = bm.active_problem.moi_backend.optimizer.model.inner_model.statistics.OptimizeTime
        bm.optimum = objective_value(bm.active_problem)

        if write
            write_csv(bm)
        end
    end
end

function run_OSQP(bm::Benchmark, write::Bool)
    bm.solver = "OSQP"
    bm.active_problem = copy(bm.base_problem)
    set_optimizer(bm.active_problem, OSQP.Optimizer)
    run_non_ram(bm, write)
end

function run_COSMO(bm::Benchmark, write::Bool)
    bm.solver = "COSMO"
    bm.active_problem = copy(bm.base_problem)
    set_optimizer(bm.active_problem, COSMO.Optimizer)
    run_non_ram(bm, write)
end

function run_Ipopt(bm::Benchmark, write::Bool)
    bm.solver = "Ipopt"
    bm.active_problem = copy(bm.base_problem)
    set_optimizer(bm.active_problem, Ipopt.Optimizer)
    run_non_ram(bm, write)
end

function run_non_ram(bm::Benchmark, write::Bool)
    bm.iterations = 0  

    t1 = time()
    optimize!(bm.active_problem)

    bm.optimisation_time = time() - t1
    bm.build_time = 0.0
    bm.optimum = objective_value(bm.active_problem)

    if write
        write_csv(bm)
    end
end

function write_csv(bm::Benchmark)
    cols = ["Time", "Test", "Solver", "Threading", "Threads", "Iterations", "Notes", 
            "JuMP build", "Build (RAM only)", "Optimisation", "Total", "Optimum Value"]

    if bm.solver == "RAM"
        row = [Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") bm.file bm.solver bm.threading bm.threads #=
               =# bm.iterations bm.notes bm.jump_time bm.build_time bm.optimisation_time #=
               =# sum([bm.jump_time, bm.build_time, bm.optimisation_time]) bm.optimum]
    else
        row = [Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") bm.file bm.solver "" "" "" #=
               =# bm.notes bm.jump_time "" bm.optimisation_time #=
               =# sum([bm.jump_time, bm.build_time, bm.optimisation_time]) bm.optimum]
    end


    file = "/home/eps116/FYP_Benchmarking/RAMBenchmarks/results/"*bm.write_file
    new = !isfile(file)

    open(file, new ? "w" : "a") do f
        CSV.write(f, Tables.table(row), header=cols, append=!new)
    end

end

function build_jump(f::String)
    pr = get_problem_file(f)
	Q,b,c,M,d = get_matrices(pr)
    p = Model()
    @variable(p, x[1:pr.nvar])
    expr = 0.5sum(x[i]Q[i,j]x[j] for i=1:pr.nvar,j=1:pr.nvar) + sum(x[i]b[i] for i=1:pr.nvar)
    @objective(p, Min, expr)
    for c in 1:pr.ncon
        @constraint(p, sum(x[i]M[c,i] for i=1:pr.nvar) <= d[c])
    end
    return p
end

function get_problem_file(name::String; 
                          path::String="/home/eps116/FYP_Benchmarking/problems/marosmeszaros",
                          ext::String=".SIF")
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

end

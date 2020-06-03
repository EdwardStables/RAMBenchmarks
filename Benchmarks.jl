using LinearAlgebra
using BenchmarkTools
using MathOptInterface
using JuMP

function single_benchmark_ram(n, m, i)
    @benchmark JuMP.optimize!(model) setup=(model=JuMP.Model(JuMP.with_optimizer(RAM.Optimizer, "Hildreth", iterations=$i));
                                            jump_generate_test(model, $n, $m))
end

#using Ipopt
#function single_benchmark_ipopt(n, m, i)
#    #Build problem build benchmark
#    @benchmark JuMP.optimize!(model) setup=(model=JuMP.Model(JuMP.with_optimizer(Ipopt.Optimizer));
#                                            jump_generate_test(model, $n, $m))
#end

"""
    function generate_random_mat(n::Int, ::Type{T}; entry_gen=rand, non_zero_gen=0.1)::Array{T,2} where T
    
Returns nxn matrix with non-zero entry selection defined by `entry_gen`, 
and probabilty of an entry being non-zero given by `non_zero_gen`. The first
entry in entry_gen should denote the return type.

"""
function generate_random_posdef_matrix(n::Int, ::Type{T}; entry_gen=rand, non_zero_gen::Float64=0.1)::Array{T,2} where T
    #FIXME: This function is terrible and needs to change
    M = zeros(T,n,n)

    while true
        for i=1:n,j=1:n
            if rand() < non_zero_gen
                M[i,j] = entry_gen(T)
            end
        end
        M = M'*M
        if isposdef(M)
            break
        end
        M = zeros(T, n, n)
    end

    return M
end

function generate_random_vector(n::Int, ::Type{T}; 
                                entry_gen=rand, 
                                non_zero_gen::Float64=0.1,
                                zero_vectors::Bool=true
                               )::Vector{T} where T
    V = zeros(T,n)

    while true
        for i=1:n
            if rand() < non_zero_gen
                V[i] = entry_gen(T)
            end
        end
        if zero_vectors || V != zeros(T,n)
            break
        end
    end

    return V
end


function moi_generate_test(model::Optimizer, n::Int, m::Int
                          )::Tuple{Vector{MOI.VariableIndex},Vector{MOI.ConstraintIndex}}
    #TODO make more generic to other solvers
    #TODO give specific value ranges (currently mostl limited to [0,1])
    
    x = MOI.add_variables(model, n)
    obj_quad_term = MOI.ScalarQuadraticTerm{Float64}[]
    obj_aff_term = MOI.ScalarAffineTerm{Float64}[]

    E = generate_random_posdef_matrix(n, Float64)
    F = generate_random_vector(n, Float64)

    for i=1:n, j=1:n
        push!(obj_quad_term, MOI.ScalarQuadraticTerm(E[i,j], x[i], x[j]))
    end

    for i=1:n
        push!(obj_aff_term, MOI.ScalarAffineTerm(F[i], x[i]))
    end

    obj = MOI.ScalarQuadraticFunction(obj_aff_term, obj_quad_term, 0.0)

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), obj) 

    #Constraint generation
    con = MOI.ConstraintIndex[]
    for i = 1:m
        F = 10 .* generate_random_vector(n, Float64, non_zero_gen=0.5, zero_vectors=false)
        func = MOI.ScalarAffineFunction(
               map((c, i) -> MOI.ScalarAffineTerm(c, i), F, x), 0.0)
        lim = MOI.LessThan(rand())
        push!(con, MOI.add_constraint(model, func, lim))
    end

    return x, con
end

"""
    jump_generate_test(model::JuMP.Model, n::Int, m::Int)

Generates a random QP minimisation problem in `model`. Sets constraints within the 
model with the name `con` and constraints with the name `x`. Access these variables 
using the syntax:
- con = model[:con]
- x = model[:x]
respectively.
"""
function jump_generate_test(model::JuMP.Model, n::Int, m::Int)
    JuMP.@variable(model, x[1:n])

    E = generate_random_posdef_matrix(n, Float64)
    F = generate_random_vector(n, Float64)

    JuMP.@objective(model, Min, 0.5sum(x[i]E[i,j]x[j] for i=1:n,j=1:n) + sum(x[i]F[i] for i=1:n))


    γ = generate_random_vector(m, Float64, non_zero_gen=1.0)
    M = zeros(m,n)
    for i=1:m
        M[i,:] = generate_random_vector(n, Float64, zero_vectors=false) 
    end
    
    JuMP.@constraint(model, con, M*x .<= γ)
    return x
end

function jump_generate_multiple(count::Int, n::Int, m::Int; identical::Bool=false)::Vector{JuMP.Model}
    models = JuMP.Model[]
    vars = []
    
    if identical == false
        for i in 1:count
            model = JuMP.Model()
            push!(vars, jump_generate_test(model, n, m))
            push!(models, model)
        end
    else
        model = JuMP.Model()
        jump_generate_test(model, n, m)
        for i in 1:count
            push!(models, copy(model))
        end
    end

    return models, vars
end

function jump_generate_pair(n, m; solver1=()->RAM.Optimizer("Hildreth"),solver2=()->RAM.Optimizer("Hildreth"))
    E = generate_random_posdef_matrix(n, Float64)
    F = generate_random_vector(n, Float64)

    γ = generate_random_vector(m, Float64, non_zero_gen=1.0)
    M = zeros(m,n)
    for i=1:m
        M[i,:] = generate_random_vector(n, Float64, zero_vectors=false) 
    end

    p1 = JuMP.Model()
    @variable(p1, x[1:n])  
    p2 = JuMP.Model()
    @variable(p2, y[1:n])  
    
    @objective(p1, Min, 0.5sum(x[i]E[i,j]x[j] for i=1:n,j=1:n) + sum(x[i]F[i] for i=1:n))
    @objective(p2, Min, 0.5sum(y[i]E[i,j]y[j] for i=1:n,j=1:n) + sum(y[i]F[i] for i=1:n))

    @constraint(p1, M*x .<= γ)
    @constraint(p2, M*y .<= γ)

    set_optimizer(p1, solver1)
    set_optimizer(p2, solver2)

    return p1, x, p2, y
end

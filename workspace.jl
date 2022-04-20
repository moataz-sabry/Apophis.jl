using Apophis, BenchmarkTools
using Symbolics: build_function, jacobian, derivative, Num, @variables
using LinearAlgebra
using SparseArrays: sparse, SparseMatrixCSC
using DifferentialEquations: ODEProblem, solve
using Sundials: CVODE_BDF

P₀ = only(@variables P₀)
T = only(@variables T)
Y = collect(only(@variables Y[1:53]));
T₀ = only(@variables T₀)
Y₀ = collect(only(@variables Y₀[1:53]));

#init(:H2, 1000.0, 1.59e6, H2 = 0.29, N2 = 0.56, O2 = 0.15)
#init(:GRI3, 1000.0, CH4 = 0.5, N2 = 0.75, O2 = 0.20)
#init(:H2, T₀, P₀, Y₀);
init(:GRI3, T₀, P₀, Y₀)

## Symbolically
#step!(gas, Y₀, T₀)
step!(gas, Y, T)

## Numerically
#step!(gas, gas.initial.mass_fractions, gas.initial.temperature)
#step!(gas, gas.current.mass_fractions, gas.current.temperature)

dṪ = only(gas.intermediate.temperature_change_rate)
dẎ = gas.intermediate.mass_change_rate
q = vcat(dẎ, dṪ)
u = vcat(Y, T)

A = jacobian(q, u);
_, A_expr = build_function(A, Y, T, Y₀, T₀, P₀)
A_function = eval(A_expr);

A = zeros(54, 54);
init(:GRI3, 1000.0, CH4=0.05, N2=0.75, O2=0.20)
step!(gas, gas.initial.mass_fractions, gas.initial.temperature)

tstart = time()
sol = backward(gas)
tend = time()
println(tend - tstart, " - ", sol.retcode)

open("./jac/GRI3.jl", "a") do io
   println(io, A_expr)
   #replace(io, "(ˍ₋out, ˍ₋arg1, T, ˍ₋arg3, T₀, P₀)" => "jac(ˍ₋out, ˍ₋arg1, T, ˍ₋arg3, T₀, P₀)")
end

#mohamed osman 
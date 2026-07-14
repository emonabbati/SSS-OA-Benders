using JuMP, CPLEX, Random, Combinatorics, Dates, Plots
using SparseArrays

num_Outsource, num_Supplier, num_AB_level, num_AD_level, num_RS_level, num_Periods, num_Scenarios = (3,4,2,2,2,15,5)
ξ = [1 2 3 4 5] # This permutation tells the model which part of generated data should be used for scenarios. For example if we use only scenario num. 3, we simply set num_Scenarios=1 and ξ[1] = 3
α1, α2 = 25,30
λ1, λ2 = 8e-3, 9e-3
include("SSS-OA-e-constraint-functions.jl")

include("generate-instance.jl")

global K, I, AB, AD, RS, T, Ω,
             fc_m, fc_bc, I_dis, I_safe,
             α, β, σ, γ, δ,
             tcm, tcb, ablc, adlc, rslc,
             cap, rcapc, sc, d, c, pc, inc, 
             n, m, M, t0, t_ab, t_ad, t_rs,
             T_1_t0, T_t0_t_ab, T_ab_ad, T_ad_rs,
             λ, μ_normal, μ_ad, μ_rs,
             FS_α, FS_η, FS_κ,
             π_vec, ϕ, ρ, itω_inds, CPT_SS_val = generate_data([λ1, λ2], [α1, α2])






w1 = 0.0
w2 = 1
ϵ = 2
SSSOA_RMP = solve_problem(w1, w2, ϵ)
ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
ϵ_min = value.(ret.SSSOA_SP[:MTR])

w1 = 0
w2 = -1
ϵ = 2
SSSOA_RMP = solve_problem(w1, w2, ϵ)
ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
ϵ_max = value.(ret.SSSOA_SP[:MTR])

ϵ_rng = range(ϵ_min,ϵ_max,50)

w1 = 1
w2 = 0

generate_figure5()
generate_figure6()
num_Scenarios0 = num_Scenarios
generate_figure7_8()
generate_figure9_10()
generate_figure11_12()
reset = () -> ( spzeros(num_Supplier), 
                spzeros(num_Supplier), 
                spzeros(num_Supplier, num_AB_level), 
                [spzeros(num_Supplier, num_AD_level) for ω in Ω], 
                [spzeros(num_Supplier, num_RS_level) for ω in Ω] )
function create_MP()

    SSSOA_RMP = Model(CPLEX.Optimizer)
    @variable(SSSOA_RMP, z[I], Bin)
    @variable(SSSOA_RMP, zp[I], Bin)
    @variable(SSSOA_RMP, ab[I,AB], Bin)
    @variable(SSSOA_RMP, ad[I,AD,Ω], Bin)
    @variable(SSSOA_RMP, rs[I,RS,Ω], Bin)

    @expression(SSSOA_RMP, C1, sum(fc_m[i]*z[i] for i ∈ I) + sum(fc_bc[i]*zp[i] for i ∈ I))
    @expression(SSSOA_RMP, C2, sum(ablc[i,a]*ab[i,a] for i ∈ I, a ∈ AB))
        
    @expression(SSSOA_RMP, C3[ω in Ω], (
        sum(adlc[i,l,ξ[ω]]*ad[i,l,ω] for i ∈ I, l ∈ AD) 
        + sum(rslc[i,r,ξ[ω]]*rs[i,r,ω] for i ∈ I, r ∈ RS) 
        )  )
        
    @expression(SSSOA_RMP, TC, C1 + C2 + sum(π_vec[ξ[ω]]*C3[ω] for ω in Ω))
    @expression(SSSOA_RMP, TC_sc[ω in Ω], C1 + C2 + C3[ω])
        
    @constraint(SSSOA_RMP,CON5[i ∈ I], sum(ab[i,AB]) == z[i] + zp[i]) # 5
    @constraint(SSSOA_RMP, CON9, sum(z) == n) # 9
    @constraint(SSSOA_RMP, CON91, sum(zp) == m) # 9-1
    @constraint(SSSOA_RMP,CON10[i ∈ I, ω ∈ Ω; i ∈ I_dis[ξ[ω]]], sum(ad[i,AD,ω]) == z[i] + zp[i]) # 10
    @constraint(SSSOA_RMP,CON101[i ∈ I, ω ∈ Ω; i ∈ I_safe[ξ[ω]]], sum(ad[i,AD,ω]) == 0) # 10-1

    @constraint(SSSOA_RMP,CON14[i ∈ I, ω ∈ Ω], sum(rs[i,RS,ω]) == sum(ad[i,AB,ω])) # 14
    @constraint(SSSOA_RMP,CON11[i ∈ I], z[i] + zp[i] ≤ 1) # 11
    @variable(SSSOA_RMP, θ ≥ -10_000_000)

 


    return SSSOA_RMP
end

function create_subproblem(z, zp, ab, ad, rs, w1, w2, ϵ)

    SSSOA_SP = Model(CPLEX.Optimizer)



    @variable(SSSOA_SP, pr[K,I,t ∈ T_1_t0] ≥ 0)
    @variable(SSSOA_SP, in_var[K,I,t ∈ t0:t_ab] ≥ 0)
    @variable(SSSOA_SP, x[K,i in I,T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs,ω in Ω; i in I_safe[ξ[ω]]] ≥ 0)
    @variable(SSSOA_SP, x′[K,i in I,t ∈ T_ad_rs,ω in Ω; i in I_dis[ξ[ω]]] ≥ 0)
    @variable(SSSOA_SP, y[K,i in I,j in I,t ∈ T_ab_ad, ω in Ω; j in I_dis[ξ[ω]]] ≥ 0)
    @variable(SSSOA_SP, h[K,i in I,t ∈ T_t0_t_ab,ω in Ω; i in I_dis[ξ[ω]]] ≥ 0)
    @variable(SSSOA_SP, ls[K,T,Ω] ≥ 0)
    @variable(SSSOA_SP, SL[K,T,Ω] ≥ 0)

    # Computing Ave

        Ave = Dict(k => 0.0 for k ∈ itω_inds)

        for i ∈ I, ω ∈ Ω
            Ave[(i,0,ω)] = 1.0
        end

        Ave_t0 = Dict(k => 0.0 for k ∈ itω_inds)

        for i ∈ I, ω ∈ Ω
            Ave_t0[(i,0,ω)] = 1.0
        end

        CPT_FS_val = [CPT_FS(i,ω,t,ab,ad,rs) for i in I, ω in Ω, t ∈ T]
        CPT_SSN = CPT_SS_val
        for i in I, ω in Ω
            for t in 1:t0
                CPT_SSN = CPT_SS_val[i,ω,t0]
                CPT_FSN = CPT_FS_val[i,ω,t0]
                CPT_SS = CPT_SSN 
                CPT_FS = CPT_FSN 
                Ave_t0[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave_t0[(i,t-1,ω)] + CPT_FS

            end    
            # println( [Ave_t0[(i,tp,ω)] for tp in 1:t0 ])
        end


        for ω in Ω
            for i in I_safe[ξ[ω]]
                for t in 1:t_rs
                    CPT_SS = CPT_SS_val[i,ω,t0]
                    CPT_FS = CPT_FS_val[i,ω,t0]
                    Ave[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave[(i,t-1,ω)] + CPT_FS
                end
            end

            for i in I_dis[ξ[ω]]
                for t in 1:t0
                    CPT_SS = CPT_SS_val[i,ω,t]
                    CPT_FS = CPT_FS_val[i,ω,t]
                    Ave[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave[(i,t-1,ω)] + CPT_FS
                end

                for t in t0+1:t_ab
                    CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,1], CPT_SS_val[i,ω,t] 
                    CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,1], 0
                    CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
                    CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
                    Ave[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave[(i,t-1,ω)] + CPT_FS
                end

                for t in t_ab+1:t_ad
                    CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,1], CPT_SS_val[i,ω,t]
                    CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,1], CPT_FS_val[i,ω,t]
                    CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
                    CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
                    Ave[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave[(i,t-1,ω)] + CPT_FS

                    
                end

                for t in t_ad+1:t_rs
                    CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,1], CPT_SS_val[i,ω,t]
                    CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,1], CPT_FS_val[i,ω,t]
                    CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
                    CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
                    Ave[(i,t,ω)] = (CPT_SS - CPT_FS)*Ave[(i,t-1,ω)] + CPT_FS

                end
            end
        end

    ∂Ave_val = ∂Ave(Ave,z,zp,ab,ad,rs)


    
    
    @expression(SSSOA_SP, C2, 
    sum(pc[k,i]*pr[k,i,t] for t ∈ T_1_t0, k ∈ K, i ∈ I) + 
    sum(inc[k,i]*in_var[k,i,t] for k ∈ K, i ∈ I, t ∈ T_t0_t_ab))
    
    @expression(SSSOA_SP, C3[ω in Ω], (
                                    sum(tcm[i,ξ[ω]]*h[k,i,t,ω] for k ∈ K, i ∈ I_dis[ξ[ω]] , t ∈ T if t0 ≤ t && t < t_ab) 
                                    + sum(tcm[i,ξ[ω]]*x[k,i,t,ω] for k ∈ K, i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs if i in I_safe[ξ[ω]]) 
                                    + sum(tcb[i,j,ξ[ω]]*y[k,i,j,t,ω] for k ∈ K, i ∈ I, j ∈ I, t ∈ T_ab_ad if j in I_dis[ξ[ω]] ) 
                                    + sum(rcapc[k,i,t,ξ[ω]]*x′[k,i,t,ω] for k ∈ K, i ∈ I_dis[ξ[ω]], t ∈ T_ad_rs)
                                ) )
    
    @expression(SSSOA_SP, TC, C2 + sum(π_vec[ξ[ω]]*C3[ω] for ω ∈ Ω))
    @expression(SSSOA_SP, TC_sc[ω ∈ Ω], C2 + C3[ω])


    @expression(SSSOA_SP, Res[k in K, ω in Ω], 1/2*(
        sum(
            1/(cap[i]*(t_rs-t0)*Ave_t0[(i,t0,ω)]*length(I_safe[ξ[ω]]))*sum(x[k,i,t,ω] for t in T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs)
                        for i in I_safe[ξ[ω]]
            ) 
        + sum(
            1/(cap[i]*(t_rs-t0)*Ave_t0[(i,t0,ω)]*length(I_dis[ξ[ω]]))*(sum(h[k,i,t,ω] for t in T_t0_t_ab)
                                                    + sum(y[k,j,i,t,ω] for j in I for t in T_ab_ad)
                                                    + sum(x′[k,i,t,ω] for t in T_ad_rs))
                                                    for i in I_dis[ξ[ω]]
                )
            )
        )
    @expression(SSSOA_SP, R[ω in Ω], sum(ϕ[k]*Res[k,ω] for k in K))
    @expression(SSSOA_SP, MTR, 1 - sum(π_vec[ξ[ω]]*R[ω]  for ω in Ω))
    # println(MTR)
    # return

    # CON4_RHS = [M*z[i] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω]  # 4
    # CON6_RHS = [sc[i]*z[i] for i ∈ I, t ∈ T_1_t0] # 6
    # CON7_RHS = [0 for k ∈ K, i ∈ I] # 7
    # CON8_RHS = [0 for t ∈ T_t0_t_ab, k ∈ K, i ∈ I, ω ∈ Ω] # 8
    # CON12_RHS = [M*zp[i]*sum(ad[j,AD,ω]) for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω]  # 12
    # CON13_RHS = [Ave[(i,t,ω)]*cap[i] for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω]  # 13
    # CON15_RHS = [ M*sum(rs[i,RS,ω]) for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω] # 15
    # CON16_RHS = [Ave[(i,t,ω)]*cap[i] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω] # 16
    # CON17_RHS = [-d[k,ω] for k ∈ K, ω ∈ Ω]  # 17

    @constraint(SSSOA_SP,CON4[i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω; i in I_safe[ξ[ω]]], sum(x[K,i,t,ω]) ≤ M*z[i])  # 4
    @constraint(SSSOA_SP,CON6[i ∈ I, t ∈ T_1_t0], sum(pr[k,i,t] for k ∈ K) ≤ sc[i]*z[i]) # 6
    @constraint(SSSOA_SP,CON7[k ∈ K, i ∈ I], in_var[k,i,t0] == sum(pr[k,i,t] for t ∈ T_1_t0))  # 7
    @constraint(SSSOA_SP,CON8[t ∈ T_t0_t_ab, k ∈ K, i ∈ I, ω ∈ Ω; i in I_dis[ξ[ω]]], in_var[k,i,t+1] <= in_var[k,i,t]-h[k,i,t,ω])  # 8
    @constraint(SSSOA_SP,CON12[i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω; j in I_dis[ξ[ω]]], sum(y[K,i,j,t,ω]) ≤ M*zp[i]*sum(ad[j,AD,ω]))  # 12
    @constraint(SSSOA_SP,CON13[i ∈ I, t ∈ T_ab_ad, ω ∈ Ω], sum(y[k,i,j,t,ω] for j ∈ I_dis[ξ[ω]], k ∈ K) ≤ Ave[(i,t,ω)]*cap[i]) # 13
    @constraint(SSSOA_SP,CON15[i ∈ I, t ∈ T_ad_rs, ω ∈ Ω; i in I_dis[ξ[ω]]], sum(x′[K,i,t,ω]) ≤ M*sum(rs[i,RS,ω])) # 15
    @constraint(SSSOA_SP, MTR ≤ ϵ)
    @constraint(SSSOA_SP,CON16[i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω; i in I_safe[ξ[ω]]], sum(x[k,i,t,ω] for k ∈ K)	≤ Ave[(i,t,ω)]*cap[i]) # 16

    
    @constraint(SSSOA_SP,CON161[i in I, t ∈ T_t0_t_ab, ω ∈ Ω; i in I_dis[ξ[ω]]], sum(h[k,i,t,ω] for k ∈ K)	≤ Ave[(i,t,ω)]*cap[i]) # 161
    @constraint(SSSOA_SP,CON162[i in I, t ∈ T_ad_rs, ω ∈ Ω; i in I_dis[ξ[ω]]], sum(x′[k,i,t,ω] for k ∈ K)	≤ Ave[(i,t,ω)]*cap[i]) # 162
    

    @constraint(SSSOA_SP,CON17[k ∈ K, ω ∈ Ω], sum(x[k,i,t,ω] for i in I_safe[ξ[ω]] for t in T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs) 
                                            + sum(
                                                sum(h[k,i,t,ω] for t ∈ T_t0_t_ab) 
                                                + sum(y[k,i,j,t,ω] for j in I_dis[ξ[ω]] for t in T_ab_ad) 
                                                + sum(x′[k,i,t,ω] for t in T_ad_rs)   for i in I_dis[ξ[ω]]) >= d[k,ξ[ω]])  # 17


    @objective(SSSOA_SP, Min, w1*TC+w2*MTR)


    ∂CON4 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)
    ∂CON6 = Array{Any}(undef, num_Supplier, num_Periods)
    ∂CON12 = Array{Any}(undef, num_Supplier, num_Supplier, num_Periods, num_Scenarios)
    ∂CON13 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)
    ∂CON15 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)
    ∂CON16 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)
    ∂CON161 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)
    ∂CON162 = Array{Any}(undef, num_Supplier, num_Periods, num_Scenarios)


    for i in I, t in T
        for ω in Ω
            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            z_val[i]=M
            ∂CON4[i,t,ω] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )

            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            ks = keys(∂Ave_val[(i,t,ω)])
            vl = .*( cap[i], values(∂Ave_val[(i,t,ω)]))
            ∂CON13[i,t,ω] = (; zip(ks, vl)...) 
            
            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            for r in RS
                rs_val[ω][i,r] = M
            end
            ∂CON15[i,t,ω] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )


            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            ∂CON16[i,t,ω] = (; zip(keys(∂Ave_val[(i,t,ω)]), .*( cap[i], values(∂Ave_val[(i,t,ω)])))...) 
            ∂CON161[i,t,ω] = (; zip(keys(∂Ave_val[(i,t,ω)]), .*( cap[i], values(∂Ave_val[(i,t,ω)])))...) 
            ∂CON162[i,t,ω] = (; zip(keys(∂Ave_val[(i,t,ω)]), .*( cap[i], values(∂Ave_val[(i,t,ω)])))...) 
            

            for j in I
                z_val, zp_val, ab_val, ad_val, rs_val = reset()
                zp_val[i] = -M*sum(ad[j,l,ω] for l ∈ AD)
                for l in AD
                    ad_val[ω][j,l] = M*zp[i]
                end
                ∂CON12[i,j,t,ω] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )            
            end
            


        end
        z_val, zp_val, ab_val, ad_val, rs_val = reset()
        z_val[i]=sc[i]
        ∂CON6[i,t] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )

    end
    set_optimizer_attribute(SSSOA_SP, "CPXPARAM_Preprocessing_Presolve", 0) 
    set_silent(SSSOA_SP)
    optimize!(SSSOA_SP)
    π_CON4 = dual.(CON4)
    π_CON6 = dual.(CON6)
    π_CON12 = dual.(CON12)
    π_CON13 = dual.(CON13)
    π_CON15 = dual.(CON15)
    π_CON16 = dual.(CON16)
    π_CON161 = dual.(CON161)
    π_CON162 = dual.(CON162)

    ∂f = (
        z = fc_m,
        zp = fc_bc,
        ab = ablc,
        ad = [π_vec[ξ[ω]]*adlc[i,l,ξ[ω]] for i in I, l in AD, ω in Ω],
        rs = [π_vec[ξ[ω]]*rslc[i,r,ξ[ω]] for i in I, r in RS, ω in Ω]
    )    


    return ( π_CON = (π_CON4 = π_CON4, π_CON6 = π_CON6, π_CON12 = π_CON12, π_CON13 = π_CON13, π_CON15 = π_CON15, π_CON16 = π_CON16, π_CON161 = π_CON161, π_CON162 = π_CON162 ), 
             ∂CON =  ( ∂CON4 =  ∂CON4,  ∂CON6 =  ∂CON6,  ∂CON12 =  ∂CON12,  ∂CON13 =  ∂CON13,  ∂CON15 =  ∂CON15,  ∂CON16 =  ∂CON16,  ∂CON161 =  ∂CON161,  ∂CON162 =  ∂CON162 ),
             ∂f = ∂f,
             SSSOA_SP = SSSOA_SP,
             Ave = Ave
    )

end


function print_var(io, model, sym)
    # print(keys(model[sym].data))
    for k in keys(model[sym].data)
        println(io,"$(String(sym))[$k] = $(value(model[sym][k])), ")
        # println(k)
    end


end

get_MP_solution = (SSSOA_RMP) -> (
    z = value.(SSSOA_RMP[:z]), 
    zp = value.(SSSOA_RMP[:zp]),
    ab = value.(SSSOA_RMP[:ab]),
    ad = value.(SSSOA_RMP[:ad]),
    rs = value.(SSSOA_RMP[:rs])
)

function CPT_FS(i,ω,t,ab,ad,rs)
    
    if (0 ≤ t ≤ t0)
        FS = sum(ab[i,a]*(1-exp(-μ_normal[i,a])) for a in AB)
    elseif ( t0+1 ≤ t ≤ t_ab)
        FS = 0.0        
    elseif ( t_ab+1 ≤ t ≤ t_ad)
        FS =  sum(ad[i,l,ω]*(1-exp(-μ_ad[i,l,ξ[ω]])) for l in AD)
    elseif ( t_ad+1 ≤ t ≤ t_rs)
        FS =  sum(rs[i,r,ω]*(1-exp(-μ_rs[i,r,ξ[ω]])) for r in RS) 
    end
    return FS
end

function ∂Ave(Ave,z,zp,ab,ad,rs)
    z_val, zp_val, ab_val, ad_val, rs_val = reset()
    res = Dict{Tuple{Int64, Int64, Int64}, Any}(itω => (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val ) for itω in itω_inds)
    CPT_FS_val = [CPT_FS(i,ω,t,ab,ad,rs) for i in I, ω in Ω, t ∈ T]

    for i in I, ω in Ω
        z_val, zp_val, ab_val, ad_val, rs_val = reset()
        res[(i,0,ω)] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )

        for t in 1:t0
            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,t], 0.0
            CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,t], 0.0
            CPT_SS = CPT_SSN + CPT_SSD
            CPT_FS = CPT_FSN + CPT_FSD        
            for a in AB
                ab_val[i,a] = -FS_α[i,a]*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].ab[i,a] + FS_α[i,a]
            end    
            res[(i,t,ω)] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )
        end

        for t in t0+1:t_ab
            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,t0], CPT_SS_val[i,ω,t]
            CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,t0], 0
            CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
            CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
            for a in AB
                ab_val[i,a] = -zp[i]*FS_α[i,a]*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].ab[i,a] + zp[i]*FS_α[i,a]
            end    
            zp_val[i] = (CPT_SSN - CPT_FSN)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].zp[i] + CPT_FSN
            z_val[i] = (CPT_SSD - CPT_FSD)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].z[i] + CPT_FSD
            res[(i,t,ω)] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )
        end


        for t in t_ab+1:t_ad
            z_val, zp_val, ab_val, ad_val, rs_val = reset()
            CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,t0], CPT_SS_val[i,ω,t]
            CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,t0], CPT_FS_val[i,ω,t]
            CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
            CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
            for a in AB
                ab_val[i,a] = -zp[i]*FS_α[i,a]*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].ab[i,a] + zp[i]*FS_α[i,a]
            end    
            for l in AD
                ad_val[ω][i,l] = -z[i]*FS_κ[i,l,ξ[ω]]*Ave[(i,t-1,ω)] + (CPT_SS - CPT_FS)*res[(i,t-1,ω)].ad[ω][i,l] + z[i]*FS_κ[i,l,ξ[ω]]
            end
            zp_val[i] = (CPT_SSN - CPT_FSN)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].zp[i] + CPT_FSN
            z_val[i] = (CPT_SSD - CPT_FSD)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].z[i] + CPT_FSD
            res[(i,t,ω)] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )
        end

        for t in t_ad+1:t_rs
            z_val, zp_val, ab_val, ad_val, rs_val = reset()

            CPT_SSN, CPT_SSD = CPT_SS_val[i,ω,t0], CPT_SS_val[i,ω,t]
            CPT_FSN, CPT_FSD = CPT_FS_val[i,ω,t0], CPT_FS_val[i,ω,t] 
            CPT_SS = zp[i]*CPT_SSN + z[i]*CPT_SSD
            CPT_FS = zp[i]*CPT_FSN + z[i]*CPT_FSD
            for a in AB
                ab_val[i,a] = -zp[i]*FS_α[i,a]*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].ab[i,a] + zp[i]*FS_α[i,a]
            end
            for l in AD
                ad_val[ω][i,l] = (CPT_SS - CPT_FS)*res[(i,t-1,ω)].ad[ω][i,l]
            end
            for r in RS
                rs_val[ω][i,r] = -z[i]*FS_η[i,r,ξ[ω]]*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].rs[ω][i,r] + z[i]*FS_η[i,r,ξ[ω]]
            end
            zp_val[i] = (CPT_SSN - CPT_FSN)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].zp[i] + CPT_FSN
            z_val[i] = (CPT_SSD - CPT_FSD)*Ave[(i,t-1,ω)] + (CPT_SS-CPT_FS)*res[(i,t-1,ω)].z[i] + CPT_FSD
            res[(i,t,ω)] = (z = z_val, zp = zp_val, ab = ab_val, ad = ad_val, rs = rs_val )
        end

    end
    return res

end




function solve_subproblem(fixed_variables, w1, w2, ϵ)
    z = fixed_variables.z
    zp = fixed_variables.zp
    ab = fixed_variables.ab
    ad = fixed_variables.ad
    rs = fixed_variables.rs

    ( π_CON , ∂CON, ∂f, SSSOA_SP, Ave) =  create_subproblem(z, zp, ab, ad, rs, w1, w2, ϵ)

    

    if is_solved_and_feasible(SSSOA_SP; dual = true)
        # println("the subproblem is solved and feasible. 👌")
        return (
            is_feasible = true,
            obj = objective_value(SSSOA_SP),
            π_CON = π_CON, 
            ∂CON = ∂CON, 
            ∂f = ∂f,
            SSSOA_SP = SSSOA_SP,
            Ave = Ave
        )
    end
    # println("the subproblem is infeasible. 🤷")
    # println("termination status = $(termination_status(SSSOA_SP))")
    # println("raw status = $(raw_status(SSSOA_SP))")
    # println("dual termination status = $(dual_status(SSSOA_SP))")
    return (
        is_feasible = false,
        obj = dual_objective_value(SSSOA_SP),
        π_CON = π_CON, 
        ∂CON = ∂CON, 
        ∂f = ∂f,
        SSSOA_SP = SSSOA_SP,
        Ave = Ave
    )    

end






function my_callback_param(cb_data, w1, w2, ϵ, SSSOA_RMP)
    status = callback_node_status(cb_data, SSSOA_RMP)
    if status != MOI.CALLBACK_NODE_STATUS_INTEGER
        ## Only add the constraint if `x` is an integer feasible solution
        return
    end
    z_k = callback_value.(cb_data, SSSOA_RMP[:z])
    zp_k = callback_value.(cb_data, SSSOA_RMP[:zp])
    ab_k = callback_value.(cb_data, SSSOA_RMP[:ab])
    ad_k = callback_value.(cb_data, SSSOA_RMP[:ad])
    rs_k = callback_value.(cb_data, SSSOA_RMP[:rs])
    fixed_variables = (z = z_k, zp = zp_k, ab = ab_k, ad = ad_k, rs = rs_k)
    θ_k = callback_value(cb_data, SSSOA_RMP[:θ])

    ret = solve_subproblem(fixed_variables, w1, w2, ϵ)


            ∂f = ret.∂f
            π_CON = ret.π_CON
            ∂CON = ret.∂CON
            ∂z =  sum(π_CON.π_CON4[i,t,ω]*∂CON.∂CON4[i,t,ω].z for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON6[i,t]*∂CON.∂CON6[i,t].z for i ∈ I, t ∈ T_1_t0; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON12[i,j,t,ω]*∂CON.∂CON12[i,j,t,ω].z for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω if j in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON13[i,t,ω]*∂CON.∂CON13[i,t,ω].z for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON15[i,t,ω]*∂CON.∂CON15[i,t,ω].z for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON16[i,t,ω]*∂CON.∂CON16[i,t,ω].z for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON161[i,t,ω]*∂CON.∂CON161[i,t,ω].z for i ∈ I, t ∈ T_t0_t_ab , ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON162[i,t,ω]*∂CON.∂CON162[i,t,ω].z for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  

            ∂zp =  sum(π_CON.π_CON4[i,t,ω]*∂CON.∂CON4[i,t,ω].zp for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON6[i,t]*∂CON.∂CON6[i,t].zp for i ∈ I, t ∈ T_1_t0; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON12[i,j,t,ω]*∂CON.∂CON12[i,j,t,ω].zp for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω if j in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON13[i,t,ω]*∂CON.∂CON13[i,t,ω].zp for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω)  
                      .+ sum(π_CON.π_CON15[i,t,ω]*∂CON.∂CON15[i,t,ω].zp for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON16[i,t,ω]*∂CON.∂CON16[i,t,ω].zp for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier) )  
                      .+ sum(π_CON.π_CON161[i,t,ω]*∂CON.∂CON161[i,t,ω].zp for i ∈ I, t ∈ T_t0_t_ab , ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  
                      .+ sum(π_CON.π_CON162[i,t,ω]*∂CON.∂CON162[i,t,ω].zp for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier))  

            ∂ab =  sum(π_CON.π_CON4[i,t,ω]*∂CON.∂CON4[i,t,ω].ab for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON6[i,t]*∂CON.∂CON6[i,t].ab for i ∈ I, t ∈ T_1_t0; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON12[i,j,t,ω]*∂CON.∂CON12[i,j,t,ω].ab for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω if j in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON13[i,t,ω]*∂CON.∂CON13[i,t,ω].ab for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON15[i,t,ω]*∂CON.∂CON15[i,t,ω].ab for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON16[i,t,ω]*∂CON.∂CON16[i,t,ω].ab for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON161[i,t,ω]*∂CON.∂CON161[i,t,ω].ab for i ∈ I, t ∈ T_t0_t_ab , ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  
                      .+ sum(π_CON.π_CON162[i,t,ω]*∂CON.∂CON162[i,t,ω].ab for i ∈ I, t ∈  T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AB_level))  

            ∂ad = zeros(num_Supplier, num_AB_level, num_Scenarios)
            for ω0 in Ω
                ∂ad[:,:,ω0] =  sum(π_CON.π_CON4[i,t,ω]*∂CON.∂CON4[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON6[i,t]*∂CON.∂CON6[i,t].ad[ω0] for i ∈ I, t ∈ T_1_t0; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON12[i,j,t,ω]*∂CON.∂CON12[i,j,t,ω].ad[ω0] for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω if j in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON13[i,t,ω]*∂CON.∂CON13[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON15[i,t,ω]*∂CON.∂CON15[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON16[i,t,ω]*∂CON.∂CON16[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON161[i,t,ω]*∂CON.∂CON161[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_t0_t_ab , ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
                        .+ sum(π_CON.π_CON162[i,t,ω]*∂CON.∂CON162[i,t,ω].ad[ω0] for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_AD_level))  
            end
            ∂rs = zeros(num_Supplier, num_AD_level, num_Scenarios)
            for ω0 in Ω
                ∂rs[:,:,ω0] =  sum(π_CON.π_CON4[i,t,ω]*∂CON.∂CON4[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON6[i,t]*∂CON.∂CON6[i,t].rs[ω0] for i ∈ I, t ∈ T_1_t0; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON12[i,j,t,ω]*∂CON.∂CON12[i,j,t,ω].rs[ω0] for i ∈ I, j ∈ I, t ∈ T_ab_ad, ω ∈ Ω if j in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON13[i,t,ω]*∂CON.∂CON13[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_ab_ad, ω ∈ Ω; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON15[i,t,ω]*∂CON.∂CON15[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω  if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON16[i,t,ω]*∂CON.∂CON16[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_t0_t_ab ∪ T_ab_ad ∪ T_ad_rs, ω ∈ Ω if i in I_safe[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON161[i,t,ω]*∂CON.∂CON161[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_t0_t_ab , ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
                        .+ sum(π_CON.π_CON162[i,t,ω]*∂CON.∂CON162[i,t,ω].rs[ω0] for i ∈ I, t ∈ T_ad_rs, ω ∈ Ω if i in I_dis[ξ[ω]]; init=zeros(num_Supplier,num_RS_level))  
            end
                
    if (ret.is_feasible)
        if θ_k < (ret.obj - 1e-6)

            ∂z = ∂z .+ ∂f.z
            ∂zp = ∂zp .+ ∂f.zp
            ∂ab = ∂ab .+ ∂f.ab
            ∂ad = ∂ad .+ ∂f.ad
            ∂rs = ∂rs .+ ∂f.rs

            
            cut = @build_constraint( 
                SSSOA_RMP[:θ] ≥ ret.obj + sum(∂z[i]*(SSSOA_RMP[:z][i]-z_k[i]) for i in I)
                            + sum(∂zp[i]*(SSSOA_RMP[:zp][i]-zp_k[i]) for i in I)
                            + sum(∂ab[i,a]*(SSSOA_RMP[:ab][i,a]-ab_k[i,a]) for i in I for a in AB)
                            + sum(∂ad[i,l,ω]*(SSSOA_RMP[:ad][i,l,ω]-ad_k[i,l,ω]) for i in I for l in AD for ω in Ω)
                            + sum(∂rs[i,r,ω]*(SSSOA_RMP[:rs][i,r,ω]-rs_k[i,r,ω]) for i in I for r in RS for ω in Ω)
            )

            MOI.submit(SSSOA_RMP, MOI.LazyConstraint(cb_data), cut)
        end
    else


            println("in infeasibility cut ret.obj = $(ret.obj)")
            cut = @build_constraint( 
                0 ≥ ret.obj + sum(∂z[i]*(SSSOA_RMP[:z][i]-z_k[i]) for i in I)
                            + sum(∂zp[i]*(SSSOA_RMP[:zp][i]-zp_k[i]) for i in I)
                            + sum(∂ab[i,a]*(SSSOA_RMP[:ab][i,a]-ab_k[i,a]) for i in I for a in AB)
                            + sum(∂ad[i,l,ω]*(SSSOA_RMP[:ad][i,l,ω]-ad_k[i,l,ω]) for i in I for l in AD for ω in Ω)
                            + sum(∂rs[i,r,ω]*(SSSOA_RMP[:rs][i,r,ω]-rs_k[i,r,ω]) for i in I for r in RS for ω in Ω)
            )

            MOI.submit(SSSOA_RMP, MOI.LazyConstraint(cb_data), cut)
    end
    return
end


function solve_problem(w1, w2, ϵ)




    SSSOA_RMP = create_MP()


    @objective(SSSOA_RMP, Min, w1*SSSOA_RMP[:TC] + SSSOA_RMP[:θ])

    my_callback = (cb_data) -> my_callback_param(cb_data, w1, w2, ϵ, SSSOA_RMP)
    set_attribute(SSSOA_RMP, MOI.LazyConstraintCallback(), my_callback)

    ϵ = 3
    set_silent(SSSOA_RMP)
    optimize!(SSSOA_RMP)
    return (
        SSSOA_RMP = SSSOA_RMP
    )
end



using XLSX


function generate_figure5()
    # # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Figure 1


    data_name = "Figure 5"
    if (~isdir(data_name))
        mkdir(data_name)
    end

    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Pareto Frtont")
        sheet[1,:] = ["ϵ","1-TRes","TCost"]
        k = 2
        for ϵ0 in ϵ_rng
            ϵ = ϵ0
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            MP_solution = get_MP_solution(SSSOA_RMP)
            ret = solve_subproblem(MP_solution, w1, w2, ϵ)

            sheet[k,1] = ϵ
            sheet[k,2] = value.(ret.SSSOA_SP[:MTR])
            sheet[k,3] = value.(ret.SSSOA_SP[:TC])+value.(SSSOA_RMP[:TC])
            k += 1
        end
    end
end

# # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Figure 6
function generate_figure6()
    data_name = "Figure 6"

    if (~isdir(data_name))
        mkdir(data_name)
    end
    ϵ = ϵ_rng[2]
    println("Generate Figure 6 for ϵ = $ϵ")

    SSSOA_RMP = solve_problem(w1, w2, ϵ)
    ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
        i = I[1]
        ω = Ω[1]
        println("create sheet for i=$i, ω=$ω")

        sheet = xf[1]
        XLSX.rename!(sheet, "i=$i, ω=$ω")
        sheet[1,:] = ["t", "Ave"]
        sheet["A2",dim=1] = collect(0:t_rs)
        sheet["B2", dim=1] = [ret.Ave[(i,t,ω)] for t in 0:t_rs] 
        for i in I
            for ω in Ω
                println("create sheet for i=$i, ω=$ω")
                (i==I[1] && ω==Ω[1]) && continue
                new_sheet = XLSX.addsheet!(xf, "i=$i, ω=$ω")
                new_sheet[1,:] = ["t", "Ave"]
                new_sheet["A2",dim=1] = collect(0:t_rs)
                new_sheet["B2", dim=1] = [ret.Ave[(i,t,ω)] for t in 0:t_rs]
            end
        end

        

    end
end


# # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Figure 7 & 8
function generate_figure7_8()
    global num_Scenarios, Ω, ξ
    data_name = "Figure 7, 8"
    if (~isdir(data_name))
        mkdir(data_name)
    end
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "Robust")
        sheet[1,:] = ["ϵ", ["Scenario $ω - TC" for ω in Ω]... , ["Scenario $ω - TR" for ω in Ω]...]
        sheet["A2", dim=1] = collect(ϵ_rng)


        
        for k in eachindex(ϵ_rng)
            ϵ0 = ϵ_rng[k]    
            ϵ = ϵ0
            println("Generate Figure 7 (Robust) for ϵ = $ϵ")
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
            for ω in Ω
                sheet[k+1,ω+1] = value(ret.SSSOA_SP[:TC_sc][ω]) + value.(SSSOA_RMP[:TC_sc][ω])
                sheet[k+1,num_Scenarios+ω+1] = value(ret.SSSOA_SP[:R][ω])
            end
        end


        for ω0 in Ω
            ξ = [ω0]
            num_Scenarios = 1
            Ω = 1:num_Scenarios
            new_sheet = XLSX.addsheet!(xf, "Deterministic - Scenario ω = $ω0")

            
            w1 = 0.0
            w2 = 1
            ϵ = 2
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
            ϵ_min0 = value.(ret.SSSOA_SP[:MTR])

            w1 = 0
            w2 = -1
            ϵ = 2
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
            ϵ_max0 = value.(ret.SSSOA_SP[:MTR])



            ϵ_rng0 = range(ϵ_min0,ϵ_max0,25)
            # println("ϵ_min0 = $ϵ_min0, ϵ_max0 = $ϵ_max0")


            for k in eachindex( ϵ_rng0)
                ϵ0 = ϵ_rng0[k]
                ϵ = ϵ0
                println("Generate Figure 7 data Deterministic for ω=$(ω0), ϵ = $ϵ")

                SSSOA_RMP = solve_problem(w1, w2, ϵ)
                ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)

                new_sheet[1,:] = ["ϵ", "Total Cost", "Total Resilience"]
                new_sheet[k+1,1] = ϵ
                new_sheet[k+1,2] = value(ret.SSSOA_SP[:TC]) + value.(SSSOA_RMP[:TC])
                new_sheet[k+1,3] = 1-value(ret.SSSOA_SP[:MTR])
            end
        end
    end
end


# # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Table 2    
function generate_table2_ver2()
    global num_Scenarios, Ω, ξ
    data_name = "Table 2"

    Table2 = Dict{Float64, Matrix{Vector{Float64}}}()
    for ϵ0 in ϵ_rng
        
        ϵ = ϵ0
        println("Generate Table 2 for ϵ = $ϵ")

        SSSOA_RMP = solve_problem(w1, w2, ϵ)
        ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
        # res_tb2 = Matrix{Float64}(undef,num_Scenarios,2)
        x_val = value.(ret.SSSOA_SP[:x])
        res_tb2 = Matrix{Vector{Float64}}(undef,num_Scenarios,2)

        for ω in Ω
            res_tb2[ω,1] = [ value(SSSOA_RMP[:C3][ω])+value(ret.SSSOA_SP[:C3][ω]), value(ret.SSSOA_SP[:MTR])]
        end
        Table2[ϵ] = res_tb2
        # println(Table2[ϵ])
    end

    for ω0 ∈ Ω
        ξ = [ω0]
        num_Scenarios = 1
        Ω = 1:num_Scenarios
        # ϵ = 2

        for ϵ0 in ϵ_rng
            ϵ = ϵ0
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)  
            x_val = value.(ret.SSSOA_SP[:x])

            Table2[ϵ][ω0,2] =  [ value(SSSOA_RMP[:C3][1])+value(ret.SSSOA_SP[:C3][1]), value(ret.SSSOA_SP[:MTR])]
        end
        # println("Table2[ϵ] = $(Table2[ϵ])")
    end
   num_Scenarios = num_Scenarios0
    Ω = 1:num_Scenarios
    ξ = Ω    
   if (~isdir(data_name))
        mkdir(data_name)
    end
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
        
        sheet = xf[1]
        XLSX.rename!(sheet, "ϵ")
        sheet[1,1] = ["ϵ"]
        sheet["A2", dim=1] = collect(ϵ_rng)
        for k in eachindex(ϵ_rng)
            new_sheet = XLSX.addsheet!(xf, "ϵ $k")
            new_sheet[1,:] = ["","Robust","","Deterministic",""]
            new_sheet["A2", dim=1] = ["Scenario $ω" for ω ∈ Ω]
            # [println(Table2[ϵ_rng[k]][ω,1][s]) for ω in Ω, s in 1:2]
            # [println(Table2[ϵ_rng[k]][ω][2]) for ω in Ω]
            # println(Table2[ϵ_rng[k]])

            new_sheet["B2"] = [Table2[ϵ_rng[k]][ω,1][s] for ω in Ω, s in 1:2]
            new_sheet["D2"] = [Table2[ϵ_rng[k]][ω,2][s] for ω in Ω, s in 1:2]
        end

    end    
end


# # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Figure 2--6
function generate_figure2_6()
    data_name = "Figure 2--6"

    Figures_2_6 = Dict{Int64, Dict{Float64, Dict{Tuple{Int64,Int64},Tuple{Any,Any}}}}()
    # for ω0 in Ω
    for ω0 in Ω
        Figure = Dict{Float64, Dict{Tuple{Int64,Int64},Tuple{Any,Any}}}()
        for ϵ0 in ϵ_rng
            
            ϵ = ϵ0
            println("Generate Figure $(ω0+1) for ϵ = $ϵ")
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            MP_solution = get_MP_solution(SSSOA_RMP)
            ret = solve_subproblem(MP_solution, w1, w2, ϵ)        

            res2 = Dict{Tuple{Int64,Int64},Tuple{Any,Any}}()
            for i in I
                # if (MP_solution.z[i]==1)
                    ω = ω0
                    res2[(i,ω)] = (plot(0:t_rs, [ret.Ave[(i,tp,ω)] for tp in 0:t_rs]),
                                    [ret.Ave[(i,tp,ω)] for tp in 0:t_rs]
                                    )
                # end
            end
            Figure[ϵ] = res2
        end
        Figures_2_6[ω0+1] = Figure
    end

   if (~isdir(data_name))
        mkdir(data_name)
    end
    # for ω in Ω
    for ω in Ω
        cdata_name = "Figure $(ω+1)" 
        XLSX.openxlsx("$(data_name)/Result-$(cdata_name).xlsx", mode = "w") do xf
            
            sheet = xf[1]
            XLSX.rename!(sheet, "ϵ")
            sheet[1,1] = ["ϵ"]
            sheet["A2", dim=1] = collect(ϵ_rng)

            for k in eachindex(ϵ_rng)
                new_sheet = XLSX.addsheet!(xf, "ϵ $k - Ave")
                cdata = Figures_2_6[ω+1][ϵ_rng[k]]
                # k = 1
                ind = 1
                for s in eachindex(cdata)
                    cAve = cdata[s][2]
                    plt = cdata[s][1]
                    savefig(plt,"$(data_name)/$(cdata_name)-$s-ϵ$k-ω-$ω.pdf")
                    println("$s")
                    new_sheet[1,ind] = "$s"
                    new_sheet[2:1+length(cAve),ind] = cAve
                    ind += 1
                end
            end

        end    
    end
end
# # # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Table 3
function generate_table3()

    data_name = "Table 3"

    Table3 = Dict{Float64,Matrix{Float64}}()
    for ϵ0 in ϵ_rng
        
        ϵ = ϵ0
        println("Generate Table 3 for ϵ = $ϵ")
        SSSOA_RMP = solve_problem(w1, w2, ϵ)
        MP_solution = get_MP_solution(SSSOA_RMP)
        ret = solve_subproblem(MP_solution, w1, w2, ϵ)    
        res_tb3 = Matrix{Float64}(undef,num_Scenarios,num_Supplier)
        Res_val = value.(ret.SSSOA_SP[:Res])
        for i in I
            if (MP_solution.z[i]==1)
                for ω in Ω

                    # Res_arr = [Res_val[k,ω] for k in K]
                    # res_tb3[ω,i] = [ minimum(Res_arr), maximum(Res_arr)]
                    res_tb3[ω,i] =  Res_val[1,ω]

                end
            else
                 for ω in Ω
                    res_tb3[ω,i] = 0
                end               
            end
        end
        Table3[ϵ] = res_tb3
    end
   if (~isdir(data_name))
        mkdir(data_name)
    end
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
        
        sheet = xf[1]
        XLSX.rename!(sheet, "ϵ")
        sheet[1,1] = ["ϵ"]
        sheet["A2", dim=1] = collect(ϵ_rng)
        for k in eachindex(ϵ_rng)
            println(Table3[ϵ_rng[k]])
            new_sheet = XLSX.addsheet!(xf, "ϵ $k")
            # headers = fill(" ", 1, 2*num_Supplier)
            headers = ["i$i" for i in I]
            new_sheet[1,:] = vec([" " headers...])
            new_sheet["A2", dim=1] = ["Scenario $ω" for ω ∈ Ω]
            new_sheet["B2"] = Table3[ϵ_rng[k]]
            # new_sheet["B2"] =  [Table3[ϵ_rng[k]][ω,1] for ω ∈ Ω]
            # new_sheet["D2"] = [Table3[ϵ_rng[k]][ω,2] for ω ∈ Ω]
            # new_sheet["F2"] = [Table3[ϵ_rng[k]][ω,3] for ω ∈ Ω]
            # new_sheet["H2"] = [Table3[ϵ_rng[k]][ω,4] for ω ∈ Ω]

        end

    end 
end

# # # # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Table 4
function generate_table4()
        global ϵ0
        data_name = "Table 4"
        α_bounds = [5 25;
                    25 30;
                    30 40
                    ]
        Table4 = Vector{NamedTuple{(:table,:Ave)}}(undef, num_Scenarios)
        
        
        for ω in Ω
            res_tb4 = Matrix{Float64}(undef,2,3)
            Ave = Matrix{Float64}(undef,num_Supplier,num_Periods)
            Table4[ω] = (table=res_tb4, Ave = Ave)
        end
        for i=1:3
            α_bnd = α_bounds[i,:]
            α = rand(α_bnd[1]:α_bnd[2],num_Supplier)

            ϵ = ϵ0
            println("Generate Table 4 for α in $α_bnd")
            SSSOA_RMP = solve_problem(w1, w2, ϵ)
            MP_solution = get_MP_solution(SSSOA_RMP)
            ret = solve_subproblem(MP_solution, w1, w2, ϵ)    

            for ω ∈ Ω            
                Table4[ω].table[1,i] = value(SSSOA_RMP[:C3][ω])+value(ret.SSSOA_SP[:C3][ω])
                Table4[ω].table[2,i] = value(ret.SSSOA_SP[:MTR])
                Table4[ω].Ave .= [ret.Ave[(ip,tp,ω)] for ip in I, tp in T]
            end
        end

        
        if (~isdir(data_name))
            mkdir(data_name)
        end
        XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
            
            sheet = xf[1]
            ω = 1
            XLSX.rename!(sheet, "ω $ω")
            sheet[1,:] = vec( [" " ["$(s[k,:])" for k in 1:3]...])
            sheet["A2", dim=1] = ["Z1", "Z2"]
            sheet["B2"] = Table4[ω].table

            # for case=2:3
            #     new_sheet = XLSX.addsheet!(xf, "Case $case")
            #     new_sheet[1,:] = ["1-TRes","TCost"]
            #     new_sheet["A2", dim=1] = figure7_data[case][:,1]
            #     new_sheet["B2", dim=1] = figure7_data[case][:,2]            
            #     plt = plot(figure7_data[case][:,1],figure7_data[case][:,2] )
            #     savefig(plt, "$(data_name)/$data_name - Case $(case).pdf")
            # end
        end


    # # # 📝📝📝📝📝📝📝📝📝📝📝📝📝📝📝                        Generate Figure 7
    #     data_name = "Figure 7"
        
    #     π_mat = [0.15 0.1 0.2 0.25 0.3
    #              0.05 0.05 0.1 0.3 0.5
    #             0.5 0.3 0.1 0.05 0.05]
    #     figure7_data = Vector{Matrix{Float64}}(undef,3)
    #     for case in 1:3
    #         π_vec = π_mat[case,:]
    #         MTR_vals = Vector{Float64}(undef,length(ϵ_rng))
    #         TC_vals = Vector{Float64}(undef,length(ϵ_rng))
    #         iter = 1
    #         for ϵ0 in ϵ_rng
                
    #             ϵ = ϵ0
    #             println("Generate Figure 7 in  case $case for ϵ = $ϵ")

    #             SSSOA_RMP = solve_problem(w1, w2, ϵ)
    #             ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
    #             MTR_vals[iter] = value(ret.SSSOA_SP[:MTR])
    #             TC_vals[iter] = value(ret.SSSOA_SP[:TC])
    #             iter += 1
    #         end
    #         res = Matrix{Float64}(undef,length(ϵ_rng),2)
    #         res[:,1] = MTR_vals
    #         res[:,2] = TC_vals

    #         figure7_data[case] = res

    #     end

    #     if (~isdir(data_name))
    #         mkdir(data_name)
    #     end
    #     XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf
            
    #         sheet = xf[1]
    #         case = 1
    #         XLSX.rename!(sheet, "Case $case")
    #         sheet[1,:] = ["1-TRes","TCost"]
    #         sheet["A2", dim=1] = figure7_data[case][:,1]
    #         sheet["B2", dim=1] = figure7_data[case][:,2]
    #         plt = plot(figure7_data[case][:,1],figure7_data[case][:,2] )
    #         savefig(plt, "$(data_name)/$data_name - Case $(case).pdf")

    #         for case=2:3
    #             new_sheet = XLSX.addsheet!(xf, "Case $case")
    #             new_sheet[1,:] = ["1-TRes","TCost"]
    #             new_sheet["A2", dim=1] = figure7_data[case][:,1]
    #             new_sheet["B2", dim=1] = figure7_data[case][:,2]            
    #             plt = plot(figure7_data[case][:,1],figure7_data[case][:,2] )
    #             savefig(plt, "$(data_name)/$data_name - Case $(case).pdf")
    #         end
    #     end
end


function generate_figure9_10()

    λ1, λ2 = 8e-3, 9e-3

    α_bounds = [5 25;
                25 30;
                30 40
                ]
    Ave_mat = Matrix{Float64}(undef,3,num_Periods)
    MTR_vals = Vector{Any}(undef,3)
    TC_vals = Vector{Any}(undef,3)
    ϵ_vals = Vector{Any}(undef,3)

    data_name = "Figure 9, 10"
    global K, I, AB, AD, RS, T, Ω,
                fc_m, fc_bc, I_dis, I_safe,
                α, β, σ, γ, δ,
                tcm, tcb, ablc, adlc, rslc,
                cap, rcapc, sc, d, c, pc, inc, 
                n, m, M, t0, t_ab, t_ad, t_rs,
                T_1_t0, T_t0_t_ab, T_ab_ad, T_ad_rs,
                λ, μ_normal, μ_ad, μ_rs,
                FS_α, FS_η, FS_κ,
                π_vec, ϕ, ρ, itω_inds, CPT_SS_val
    if (~isdir(data_name))
        mkdir(data_name)
    end
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf

        for case=1:3
            global ϵ, ϵ_min, ϵ_max, w1, w2


            α1, α2 = α_bounds[case,:]
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

            MTR_vals[case] = []
            TC_vals[case] = []
            ϵ_vals[case] = []
            for k in eachindex(ϵ_rng)
                ϵ0 = ϵ_rng[k]
                ϵ = ϵ0
                println("Generate figure 9, 10 in case $case for ϵ = $ϵ, α in $(α_bounds[case,:])")
                SSSOA_RMP = solve_problem(w1, w2, ϵ)
                MP_solution = get_MP_solution(SSSOA_RMP)

                ret = solve_subproblem(MP_solution, w1, w2, ϵ)
                i = 1
                ω = 5
                for t in T
                    Ave_mat[case,t] = ret.Ave[(i,t,ω)]
                end
                
                push!(ϵ_vals[case],ϵ)
                push!(MTR_vals[case],value.(ret.SSSOA_SP[:MTR]))
                push!(TC_vals[case], value.(ret.SSSOA_SP[:TC])+value.(SSSOA_RMP[:TC]))
            end

            if (case==1)
                sheet = xf[1]
                XLSX.rename!(sheet, "Pareto Frtont (Case $case)")
                sheet[1,:] = ["1-TRes","TCost", "ϵ"]
                

                sheet["A2", dim=1] = MTR_vals[case]
                sheet["B2", dim=1] = TC_vals[case]
                sheet["C2", dim=1] = ϵ_vals[case]

                new_sheet = XLSX.addsheet!(xf, "Case $case - Ave_(1,t,5)")
                new_sheet["A1"] = [Ave_mat[case,:]';  T']        
            else
                new_sheet = XLSX.addsheet!(xf, "Pareto Frtont (Case $case)")
                new_sheet[1,:] = ["1-TRes","TCost"]
                

                new_sheet["A2", dim=1] = MTR_vals[case]
                new_sheet["B2", dim=1] = TC_vals[case]
                new_sheet["C2", dim=1] = ϵ_vals[case]

                new_sheet = XLSX.addsheet!(xf, "Case $case - Ave_(1,t,5)")
                new_sheet["A1"] = [Ave_mat[case,:]';  T']            
            end
        end




    end

end


function generate_figure11_12()

    α1, α2 = 25, 35

    λ_bounds = [8e-3 9e-3
                1e-2 2e-2
                8e-2 9e-2
                ]
    global K, I, AB, AD, RS, T, Ω,
                fc_m, fc_bc, I_dis, I_safe,
                α, β, σ, γ, δ,
                tcm, tcb, ablc, adlc, rslc,
                cap, rcapc, sc, d, c, pc, inc, 
                n, m, M, t0, t_ab, t_ad, t_rs,
                T_1_t0, T_t0_t_ab, T_ab_ad, T_ad_rs,
                λ, μ_normal, μ_ad, μ_rs,
                FS_α, FS_η, FS_κ,
                π_vec, ϕ, ρ, itω_inds, CPT_SS_val
    Ave_mat = Matrix{Float64}(undef,3,num_Periods)
    MTR_vals = Vector{Any}(undef,3)
    TC_vals = Vector{Any}(undef,3)
    data_name = "Figure 11, 12"
    if (~isdir(data_name))
        mkdir(data_name)
    end
    XLSX.openxlsx("$(data_name)/Result-$(data_name).xlsx", mode = "w") do xf

        for case=1:3
            global λ1, λ2, ϵ, ϵ_min, ϵ_max, w1, w2


            λ1, λ2 = λ_bounds[case,:]

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

            MTR_vals[case] = []
            TC_vals[case] = []
            for k in eachindex(ϵ_rng)
                ϵ0 = ϵ_rng[k]
                ϵ = ϵ0
                println("Generate figure 11, 12 in case $case for ϵ = $ϵ, λ in $(λ_bounds[case,:])")
                SSSOA_RMP = solve_problem(w1, w2, ϵ)
                ret = solve_subproblem(get_MP_solution(SSSOA_RMP), w1, w2, ϵ)
                i = 1
                ω = 5
                for t in T
                    Ave_mat[case,t] = ret.Ave[(i,t,ω)]
                end
                
                push!(MTR_vals[case],value.(ret.SSSOA_SP[:MTR]))
                push!(TC_vals[case], value.(ret.SSSOA_SP[:TC])+value.(SSSOA_RMP[:TC]))
            end

            if (case==1)
                sheet = xf[1]
                XLSX.rename!(sheet, "Pareto Frtont (Case $case)")
                sheet[1,:] = ["1-TRes","TCost"]
                

                sheet["A2", dim=1] = MTR_vals[case]
                sheet["B2", dim=1] = TC_vals[case]

                new_sheet = XLSX.addsheet!(xf, "Case $case - Ave_(1,t,5)")
                new_sheet["A1"] = [Ave_mat[case,:]';  T']        
            else
                new_sheet = XLSX.addsheet!(xf, "Pareto Frtont (Case $case)")
                new_sheet[1,:] = ["1-TRes","TCost"]
                

                new_sheet["A2", dim=1] = MTR_vals[case]
                new_sheet["B2", dim=1] = TC_vals[case]

                new_sheet = XLSX.addsheet!(xf, "Case $case - Ave_(1,t,5)")
                new_sheet["A1"] = [Ave_mat[case,:]';  T']            
            end
        end




    end



end
using  Random

function generate_data(λbound, αbound)

    Random.seed!(1)



    K = 1:num_Outsource
    I = 1:num_Supplier

    AB = 1:num_AB_level
    AD = 1:num_AD_level
    RS = 1:num_RS_level
    T = 1:num_Periods
    Ω = 1:num_Scenarios

    fc_m = rand(400:1000,1,num_Supplier)
    fc_bc = rand(700:1200,1,num_Supplier)

    I_dis = [
        [],
        [1],
        [2,3],
        [1,3,4],
        [1,2,3,4]
    ]
    I_safe = [
        [1,2,3,4],
        [2,3,4],
        [1,4],
        [2],
        []
    ]

    randx = (a,b,args...) -> a .+ rand(Float32,args...)*(b-a)
    α = rand(αbound[1]:αbound[2],num_Supplier)
    β = randx(0.25,0.3,num_Supplier,num_Scenarios)
    σ = randx(0.25,0.3,num_Supplier)
    γ = randx(3.5,5,num_Supplier)
    δ = randx(2.5,3,num_Outsource)

    tcm = randx(5,20,num_Supplier, num_Scenarios)
    for ω in Ω
        for i in I_dis[ω]
            tcm[i,ω] *= γ[i]
        end
    end

    tcb = randx(10,25,num_Supplier, num_Supplier, num_Scenarios)

    for ω in Ω
        for j in I_dis[ω]
            @. tcb[:,j,ω] = tcb[:,j,ω]*γ[j]
        end
    end





    ablc = zeros(num_Supplier, num_AB_level)
    ablc[:, 1] = rand(100:300, num_Supplier,1)
    ablc[:, 2] = ablc[:,1] + rand(300:500, num_Supplier,1)

    adlc = zeros(num_Supplier, num_AD_level, num_Scenarios)
    adlc[:, 1, :] = rand(100:300, num_Supplier,1,num_Scenarios)
    adlc[:, 2, :] = adlc[:,1, :].*(rand(num_Supplier,num_Scenarios).+2)

    rslc = zeros(num_Supplier, num_RS_level, num_Scenarios)
    rslc[:, 1, :] = rand(300:900, num_Supplier, num_Scenarios)
    rslc[:, 2, :] = rslc[:,1,:].*(rand(num_Supplier, num_Scenarios).+2)

    cap = rand(400:1000,1,num_Supplier)
    # cap = rand(200:400,1,num_Supplier)

    # cap = rand(0:1,1,num_Supplier)

    rcapc = convert.(Float64,rand(10:25, num_Outsource, num_Supplier, num_Periods, num_Scenarios))
    for k in K, i in I, t in T, ω in Ω
        rcapc[k,i,t,ω] = Float64(rcapc[k,i,t,ω])*γ[i]
    end
    sc = cap
    pc = rand(1:2, num_Outsource, num_Supplier)

    inc = rand(1:2, num_Outsource, num_Supplier)
    d = convert.(Float64,rand(100:400, num_Outsource, num_Periods))
    for k in K, ω in Ω
        d[k,ω] = d[k,ω]*δ[k]
    end
    # d = rand(1:4, num_Outsource, num_Periods).*δ


    c = rand(1:2, num_Outsource, num_Supplier)
    n = floor(num_Supplier/2)+1
    m = num_Supplier - n
    M = 1_000_000
    t0 = 3
    t_ab = 5
    t_ad = 7
    t_rs = num_Periods


    T_1_t0 = 1:(t0-1)
    T_t0_t_ab = t0:(t_ab-1)
    T_ab_ad = t_ab:(t_ad-1)
    T_ad_rs = t_ad:t_rs



    λ = randx(λbound[1], λbound[2],  num_Supplier, num_Scenarios)  
    for ω in Ω
        for i in I_dis[ω]
            λ[i,ω] *= α[i]
        end
    end

    μ_normal = zeros(num_Supplier, num_AB_level)
    μ_normal[:,1] = randx(0.06, 0.08,  num_Supplier, 1)
    μ_normal[:,2] = randx(0.09, 0.12, num_Supplier, 1)


    μ_ad = zeros(num_Supplier, num_AD_level, num_Scenarios) 
    μ_ad[:,1,:] = randx(0.125, 0.3,  num_Supplier,  num_Scenarios)   
    μ_ad[:,2,:] = randx(0.35, 0.6,  num_Supplier, num_Scenarios )  

    μ_rs = zeros(num_Supplier, num_RS_level, num_Scenarios)
    μ_rs[:,1,:] = randx(0.65, 0.8, num_Supplier, num_Scenarios) 
    μ_rs[:,2,:] = randx(0.85, 1,  num_Supplier, num_Scenarios) 


    FS_α = 1 .- exp.(μ_normal)
    FS_κ = 1 .-exp.(μ_ad)
    FS_η = 1 .- exp.(μ_rs)

    π_vec = [0.15 0.1 0.2 0.25 0.3]

    ϕ = rand(1,num_Supplier)
    ϕ /= sum(ϕ)


    ρ = rand(num_Supplier, num_Scenarios)
    for ω in Ω
        ρ[I_safe[ω], ω] .= 0
        ρ[I_dis[ω], ω] .= ρ[I_dis[ω], ω]./sum(ρ[I_dis[ω], ω])

    end

    itω_inds = [(i,t,ω) for i ∈ I for t ∈  [0 T'] for ω ∈ Ω]

    function CPT_SS(i,ω,t)
        if (0 ≤ t ≤ t0)
            SS = exp(-λ[i,ξ[ω]])
        elseif ( t0+1 ≤ t ≤ t_ab)
            SS = exp(-β[i,ξ[ω]])            
        elseif ( t_ab+1 ≤ t ≤ t_ad)
            SS= 1.0
        elseif ( t_ad+1 ≤ t ≤ t_rs)
            SS = 1.0           
        end
        return SS
    end

    CPT_SS_val = [CPT_SS(i,ω,t) for i in I, ω in Ω, t in T]

    return K, I, AB, AD, RS, T, Ω,
             fc_m, fc_bc, I_dis, I_safe,
             α, β, σ, γ, δ,
             tcm, tcb, ablc, adlc, rslc,
             cap, rcapc, sc, d, c, pc, inc, 
             n, m, M, t0, t_ab, t_ad, t_rs,
             T_1_t0, T_t0_t_ab, T_ab_ad, T_ad_rs,
             λ, μ_normal, μ_ad, μ_rs,
             FS_α, FS_η, FS_κ,
             π_vec, ϕ, ρ, itω_inds, CPT_SS_val

end

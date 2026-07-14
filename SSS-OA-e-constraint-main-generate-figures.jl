using Plots, LaTeXStrings, XLSX, Printf, Measures  

num_Outsource, num_Supplier, num_AB_level, num_AD_level, num_RS_level, num_Periods, num_Scenarios = (3,4,2,2,2,15,5)
ξ = [1 2 3 4 5] # This permutation tells the model which part of generated data should be used for scenarios. For example if we use only scenario num. 3, we simply set num_Scenarios=1 and ξ[1] = 3
α1, α2 = 25,30
λ1, λ2 = 8e-3, 9e-3
include("generate-data-V8.jl")


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




function pareto_set_min_min(Z1,Z2)
    pZ1 = []
    pZ2 = []
    for i in eachindex(Z1)
        # println("new point = $((Z1[i], Z2[i]))")
        local is_dominated = false
        for k in eachindex(pZ1)
            # println("Current Pareto = $((pZ1[k], pZ2[k]))")

            if ((pZ1[k]<= Z1[i] && pZ2[k]<= Z2[i]) && ( pZ1[k]<Z1[i] || pZ2[k]<Z2[i]))
                is_dominated = true
                break
            end

        end
        if (~is_dominated)
            push!(pZ2, Z2[i])
            push!(pZ1, Z1[i])
        end
    end
    return pZ1, pZ2
end


function grouped_range_bar_advanced_avg(categories, data1, data2;
    labels=("Attribute A", "Attribute B"),
    colors=(:blue, :red),
    average_color=:white,
    title="Grouped Range Bar Chart (Min, Avg, Max)",
    xlabel="Categories",
    ylabel="Values",
    bar_width=0.3,
    box_alpha=0.7,
    average_size=6,
    show_average_legend=true,
    )

    # Unpack the data
    min1, avg1, max1 = data1
    min2, avg2, max2 = data2
    n_categories = length(categories)

    # Validate data
    @assert length(min1) == length(avg1) == length(max1) == length(categories) "Data1 arrays must have same length as categories"
    @assert length(min2) == length(avg2) == length(max2) == length(categories) "Data2 arrays must have same length as categories"
    @assert all(min1 .<= avg1 .<= max1) "Data1 must satisfy min ≤ avg ≤ max"
    @assert all(min2 .<= avg2 .<= max2) "Data2 must satisfy min ≤ avg ≤ max"

    # Create plot
    p = plot(size=plot_size,
            legend=legend_position,
            left_margin=left_margin,
            right_margin=right_margin,
            top_margin=top_margin,
            bottom_margin=bottom_margin)

    # Plot first attribute (boxes + average markers)
    for i in 1:n_categories
        x_left = i - bar_width
        x_right = i - bar_width/3  # Slightly narrower for visual separation
        
        # Draw box from min to max
        rect_x = [x_left, x_right, x_right, x_left, x_left]
        rect_y = [min1[i], min1[i], max1[i], max1[i], min1[i]]
        
        plot!(p, rect_x, rect_y, 
            seriestype=:shape, 
            fillcolor=colors[1], 
            linecolor=darken_color(colors[1]),
            linewidth=2,
            label=(i == 1 ? labels[1] : ""),
            alpha=box_alpha)
        
        # Draw average line inside the box
        avg_line_x = [x_left, x_right]
        avg_line_y = [avg1[i], avg1[i]]
        
        plot!(p, avg_line_x, avg_line_y,
            seriestype=:path,
            linecolor=average_color,
            linewidth=3,
            label=(i == 1 && show_average_legend ? "Average" : ""))
    end

    # Plot second attribute (boxes + average markers)
    for i in 1:n_categories
        x_left = i + bar_width/3  # Start after first bar
        x_right = i + bar_width
        
        # Draw box from min to max
        rect_x = [x_left, x_right, x_right, x_left, x_left]
        rect_y = [min2[i], min2[i], max2[i], max2[i], min2[i]]
        
        plot!(p, rect_x, rect_y, 
            seriestype=:shape, 
            fillcolor=colors[2], 
            linecolor=darken_color(colors[2]),
            linewidth=2,
            label=(i == 1 ? labels[2] : ""),
            alpha=box_alpha)
        
        # Draw average line inside the box
        avg_line_x = [x_left, x_right]
        avg_line_y = [avg2[i], avg2[i]]
        
        plot!(p, avg_line_x, avg_line_y,
            seriestype=:path,
            linecolor=average_color,
            linewidth=3,
            label=false)  # Don't repeat average in legend
    end

    # Add average points for better visibility (optional)
    scatter!(p, 1:n_categories .- bar_width/2, avg1,
            markerstrokecolor=average_color,
            markercolor=average_color,
            markersize=average_size,
            label=false)

    scatter!(p, 1:n_categories .+ bar_width/2, avg2,
            markerstrokecolor=average_color,
            markercolor=average_color,
            markersize=average_size,
            label=false)

    # Customize plot
    xticks!(1:n_categories, categories)
    xlabel!(xlabel)
    ylabel!(ylabel)
    title!(title)
    xlims!(0.5, n_categories + 0.5)

    return p

end

function grouped_range_bar_advanced(categories, data1, data2;
    labels=("Attribute A", "Attribute B"),
    colors=(:blue, :red),
    title="Grouped Range Bar Chart",
    xlabel="Categories",
    ylabel="Value Range",
    bar_width=0.35,
    alpha=0.7,
    legend_position=:topleft)
    
    # Data validation
    min1, max1 = data1
    min2, max2 = data2
    
    @assert length(categories) == length(min1) == length(max1) == length(min2) == length(max2) "All input arrays must have the same length"
    @assert all(min1 .<= max1) "min1 values must be <= max1 values"
    @assert all(min2 .<= max2) "min2 values must be <= max2 values"
    
    n_categories = length(categories)
    
    # Create plot with error handling
    try
        p = plot(size=(800, 500), legend=legend_position)
        
        # Helper function to plot one set of bars
        function plot_bars!(min_vals, max_vals, x_offset, color, label)
            for i in 1:n_categories
                x_left = i + x_offset - bar_width/2
                x_right = i + x_offset + bar_width/2
                
                rect_x = [x_left, x_right, x_right, x_left, x_left]
                rect_y = [min_vals[i], min_vals[i], max_vals[i], max_vals[i], min_vals[i]]
                
                plot!(p, rect_x, rect_y, 
                seriestype=:shape, 
                fillcolor=color, 
                linecolor=darken_color(color),
                linewidth=2,
                label=(i == 1 ? label : ""),
                alpha=alpha)
            end
        end
        
        # Plot both sets of bars
        plot_bars!(min1, max1, -bar_width/2, colors[1], labels[1])
        plot_bars!(min2, max2, bar_width/2, colors[2], labels[2])
        
        # Customize plot
        xticks!(collect( 1:n_categories), categories)
        xlabel!(xlabel)
        ylabel!(ylabel)
        title!(title)
        xlims!(0.5, n_categories + 0.5)
        
        return p
        
    catch e
        println("Error creating plot: ", e)
        rethrow(e)
    end
end


function darken_color(color)
    if color == :blue
        return :darkblue
    elseif color == :red
        return :darkred
    elseif color == :green
        return :darkgreen
    elseif color == :orange
        return :darkorange
    else
        return color # Fallback
    end
end

gr()


    # plt = scatter(MTR_vals, TC_vals/10.0^4, 
    #         xlabel = "\${Z_{resilience}} = \$ 1 - Total Resilience",
    #         ylabel = "\$ Z_{cost} \\;(\\!\\times 10^4 ) \$",
    #         yformatter = y -> y >= 10000 ? string(round(y/10000, digits=2)) : string(y),
    #         legend = false
    #         #  label="\$(Z_{resilience},Z_{cost})\$"
    #         )


# Figure 5: Pareto front

    data_path = "Figure 5/Result-Figure 5.xlsx"
    global Z1_values, Z2_values
    XLSX.openxlsx(data_path, mode = "r") do xf
        global Z1_values, Z2_values
        sheet = xf[1]
        Z1_values =  Matrix{Float64}(sheet[2:end,2])
        Z2_values =  Matrix{Float64}(sheet[2:end,3] )
        
    end
    Z1_values, Z2_values = pareto_set_min_min(Z1_values, Z2_values)
    plt = plot(1 .- Z1_values, Z2_values, legend=:none,
    xlabel = "Total Resilience",
    ylabel = "Total Cost",
    xlims = (0.0, 0.3),
    ylims = (1e2, 9e5),
    xguidefontsize = 14,
    lw = 2,
    yguidefontsize = 14 )

    xticks_props = collect(0:0.1:0.3)
    yticks_props = [3e4, collect((2:2:10).*1e5)...]
    xticks!(xticks_props,latexstring.(xticks_props),xtickfont = font(12))

    math_string = replace.([@sprintf("%.2e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
    yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))


    if (~isdir("Figures"))
        mkdir("Figures")
    end

    savefig(plt, "Figures/Figure 5.pdf")


# Figure 6:  the availability curves of safe (center 4) and disrupted suppliers (center 1) under scenario 2.


    data_path = "Figure 6/Result-Figure 6.xlsx"
    scen = 2

    global Ave1_values, Ave4_values
    XLSX.openxlsx(data_path, mode = "r") do xf
        global Ave1_values, Ave4_values
        sheet = xf["i=1, ω=$scen"]
        Ave1_values =  Matrix{Float64}(sheet[2:end,2])
        sheet = xf["i=4, ω=$scen"]
        Ave4_values =  Matrix{Float64}(sheet[2:end,2] )
        
    end

    plt = plot(1:length(Ave4_values), [Ave1_values,Ave4_values], legend=:bottom,
    xlabel = "Time",
    ylabel = "Availability",
    # xlims = (-0.1, 1),
    ylims = (0.1, 1.1),
    xguidefontsize = 14,
    yguidefontsize = 14,
    lw = 2.5,
    label = ["Supplier 1" "Supplier 4"] )


    xticks_props = collect(1:16)
    yticks_props = collect(0.1:0.1:1)
    xticks!(xticks_props,latexstring.(collect(0:15)),xtickfont = font(12))

    yticks!(yticks_props,latexstring.(yticks_props),ytickfont = font(12))
    # display(plt)


    # math_string = replace.([@sprintf("%.2e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
    # yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))


    if (~isdir("Figures"))
        mkdir("Figures")
    end

    savefig(plt, "Figures/Figure 6.pdf")

    data_path = "Figure 7, 8/Result-Figure 7, 8.xlsx"

# Figure 7, 8
    global TC_det_mins, TC_det_maxs
    global TC_rob_maxs, TC_rob_mins
    global scn_text

    global TR_det_mins, TR_det_maxs
    global TR_rob_maxs, TR_rob_mins


    XLSX.openxlsx(data_path, mode = "r") do xf
        global TC_det_mins, TC_det_maxs
        global TC_rob_maxs, TC_rob_mins
        global scn_text
        global TR_det_mins, TR_det_maxs
        global TR_rob_maxs, TR_rob_mins

        scn_text = []
        TC_det_mins = []
        TC_det_maxs = []
        TC_rob_mins = []
        TC_rob_maxs = []

        TR_det_mins = []
        TR_det_maxs = []
        TR_rob_mins = []
        TR_rob_maxs = []

        Rb = xf["Robust"]
        for ω in Ω
            push!(scn_text, "Scenario $ω")
            TC_mx = maximum( Rb[2:end,ω+1])
            TC_mn = minimum( Rb[2:end,ω+1])
            push!(TC_rob_maxs, TC_mx)
            push!(TC_rob_mins, TC_mn)


            TR_mx = maximum( 1 .- Rb[2:end,num_Scenarios+ω+1])
            TR_mn = minimum( 1 .- Rb[2:end,num_Scenarios+ω+1])
            push!(TR_rob_maxs, TR_mx)
            push!(TR_rob_mins, TR_mn)



            Dt = xf["Deterministic - Scenario ω = $ω"]
            TC_Dt_mx = maximum( Dt[2:end,2])
            TC_Dt_mn = minimum( Dt[2:end,2])
            push!(TC_det_maxs, TC_Dt_mx)
            push!(TC_det_mins, TC_Dt_mn)


            TR_Dt_mx = maximum( Dt[2:end,3])
            TR_Dt_mn = minimum( Dt[2:end,3])
            push!(TR_det_maxs, TR_Dt_mx)
            push!(TR_det_mins, TR_Dt_mn)

        end
    end

    p = grouped_range_bar_advanced(string.(scn_text), (TC_det_mins,TC_det_maxs), (TC_rob_mins,TC_rob_maxs),
    labels=("Deterministic", "Robust"),
    colors=(:purple, :gold),
    title=" ",
    xlabel="Scenarios",
    ylabel="Total Cost Range",
    bar_width=0.3)

    yticks_props = [3e4, collect((2:2:10).*1e5)...]
    math_string = replace.([@sprintf("%.2e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
    yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))


    savefig(p, "Figures/Figure 7.pdf")



    p = grouped_range_bar_advanced(string.(scn_text), (TR_det_mins,TR_det_maxs), (TR_rob_mins,TR_rob_maxs),
    labels=("Deterministic", "Robust"),
    legend_position=:outertopright,
    colors=(:purple, :gold),
    title=" ",
    xlabel="Scenarios",
    ylabel="Total Resilience Range",
    bar_width=0.3)

    yticks_props = collect(0:0.1:1)
    # math_string = replace.([@sprintf("%.2e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
    yticks!(yticks_props,latexstring.(yticks_props),ytickfont = font(12))
    plot!(p, left_margin=2mm)

    savefig(p, "Figures/Figure 8.pdf")


# # Figure 9, 10
    global plt
    global plt_avl
    data_path = "Figure 9, 10/Result-Figure 9, 10.xlsx"
    XLSX.openxlsx(data_path, mode = "r") do xf
        global plt
        global plt_avl
        case = 1
        sheet = xf["Pareto Frtont (Case $case)"]
        Z1_values =  Matrix{Float64}(sheet[2:end,1])
        Z2_values =  Matrix{Float64}(sheet[2:end,2])
        Z1_values, Z2_values = pareto_set_min_min(Z1_values, Z2_values)
        # println("Z1_values = $Z1_values, Z2_values = $Z2_values")
        # plt = plot(1 .- Z1_values, Z2_values,
        # xlabel = "Total Resilience",
        # ylabel = "Total Cost",
        # xlims = (0.0, 0.3),
        # ylims = (1e2, 10e5),
        # xguidefontsize = 14,
        # yguidefontsize = 14,
        # label="Severity level $case",
        # legend_position=:bottomright,
        # linewidth=2 )

        # xticks_props = collect(0:0.1:0.3)
        # yticks_props = [3e4, collect((1:2:9).*1e5)...]
        # xticks!(xticks_props,latexstring.(xticks_props),xtickfont = font(12))

        # math_string = replace.([@sprintf("%.1e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
        # yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))
        plt = plot(1 .- Z1_values, Z2_values, legend=:none,
        xlabel = "Total Resilience",
        ylabel = "Total Cost",
        xlims = (0.0, 0.3),
        ylims = (1e2, 9e5),
        xguidefontsize = 14,
        lw = 2,
        yguidefontsize = 14 )

        xticks_props = collect(0:0.1:0.3)
        yticks_props = [3e4, collect((2:2:10).*1e5)...]
        xticks!(xticks_props,latexstring.(xticks_props),xtickfont = font(12))

        math_string = replace.([@sprintf("%.2e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
        yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))

        sheet = xf["Case $case - Ave_(1,t,5)"]
        avl_array = [1 Matrix{Float64}(sheet[1,:])]
        plt_avl = plot(1:length(avl_array), avl_array', legend=:bottomright,
        xlabel = "Time",
        ylabel = "Availability",
        ylims = (0.1, 1.1),
        xguidefontsize = 14,
        yguidefontsize = 14,
        lw = 2.5,
        label="Severity level $case")



        xticks_props = collect(1:16)
        yticks_props = collect(0.1:0.1:1)
        xticks!(xticks_props,latexstring.(collect(0:15)),xtickfont = font(12))

        yticks!(yticks_props,latexstring.(yticks_props),ytickfont = font(12))    


        for case = 2:3
            sheet = xf["Pareto Frtont (Case $case)"]
            Z1_values =  Matrix{Float64}(sheet[2:end,1])
            Z2_values =  Matrix{Float64}(sheet[2:end,2])
            Z1_values, Z2_values = pareto_set_min_min(Z1_values, Z2_values)

            plot!(plt, 1 .- Z1_values, Z2_values, label="Severity level $case",linewidth=2 )


            sheet = xf["Case $case - Ave_(1,t,5)"]
            avl_array = [1 Matrix{Float64}(sheet[1,:])]
            plot!(plt_avl, 1:length(avl_array), avl_array',
            lw = 2.5,
            label="Severity level $case")

        end

    end
    savefig(plt, "Figures/Figure 9.pdf")
    savefig(plt_avl, "Figures/Figure 10.pdf")


# # Figure 11 , 12
    global plt
    global plt_avl
    data_path = "Figure 11, 12/Result-Figure 11, 12.xlsx"
    labels = ["Low failure rate", "Medium Failure rate", "High Failure rate"]
    XLSX.openxlsx(data_path, mode = "r") do xf
        global plt
        global plt_avl
        case = 1
        sheet = xf["Pareto Frtont (Case $case)"]
        Z1_values =  Matrix{Float64}(sheet[2:end,1])
        Z2_values =  Matrix{Float64}(sheet[2:end,2])
        Z1_values, Z2_values = pareto_set_min_min(Z1_values, Z2_values)
        # println("Z1_values = $Z1_values, Z2_values = $Z2_values")
        plt = plot(1 .- Z1_values, Z2_values,
        xlabel = "Total Resilience",
        ylabel = "Total Cost",
        xlims = (0.0, 1),
        ylims = (1e2, 10e5),
        xguidefontsize = 14,
        yguidefontsize = 14,
        label=labels[case],
        legend_position=:topright,
        linewidth=2 )

        xticks_props = collect(0:0.2:1)
        yticks_props = [3e4, collect((1:2:9).*1e5)...]
        xticks!(xticks_props,latexstring.(xticks_props),xtickfont = font(12))

        math_string = replace.([@sprintf("%.1e", val) for val in yticks_props], r"([0-9.]+)[Ee]\+?0*([1-9][0-9]*|0)" => s"\1 \\times 10^\2")
        yticks!(yticks_props,latexstring.(math_string),ytickfont = font(12))

        sheet = xf["Case $case - Ave_(1,t,5)"]
        avl_array = [1 Matrix{Float64}(sheet[1,:])]
        plt_avl = plot(1:length(avl_array), avl_array', legend=:bottomright,
        xlabel = "Time",
        ylabel = "Availability",
        ylims = (0, 1.1),
        xguidefontsize = 14,
        yguidefontsize = 14,
        lw = 2.5,
        label=labels[case]
        )


        xticks_props = collect(1:16)
        yticks_props = collect(0.1:0.1:1)
        xticks!(xticks_props,latexstring.(collect(0:15)),xtickfont = font(12))

        yticks!(yticks_props,latexstring.(yticks_props),ytickfont = font(12))    


        for case = 2:3
            sheet = xf["Pareto Frtont (Case $case)"]
            Z1_values =  Matrix{Float64}(sheet[2:end,1])
            Z2_values =  Matrix{Float64}(sheet[2:end,2])
            Z1_values, Z2_values = pareto_set_min_min(Z1_values, Z2_values)

            plot!(plt, 1 .- Z1_values, Z2_values, label=labels[case] ,linewidth=2 )


            sheet = xf["Case $case - Ave_(1,t,5)"]
            avl_array = [1 Matrix{Float64}(sheet[1,:])]
            plot!(plt_avl, 1:length(avl_array), avl_array',
            lw = 2.5,
            label=labels[case])

        end

    end
    savefig(plt, "Figures/Figure 11.pdf")
    savefig(plt_avl, "Figures/Figure 12.pdf")

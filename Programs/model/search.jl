# Some simple dynamic simulations with random search

using Plots, Distributions, LaTeXStrings; gr(border = :box, grid = true, minorgrid = true, gridalpha = 0.2,
xguidefontsize = 13, yguidefontsize = 13, xtickfontsize = 13, ytickfontsize=13,
linewidth = 2, gridstyle = :dash, gridlinewidth = 1.2, margin = 10* Plots.px,legendfontsize =12)

include("static_model.jl")

# Standard calibration (unemployment rate of 0.4)
J       = 100
χ       = 0.9
w       = 1
m       = Mkt(J = J, w = w, χ = χ)
nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = true)

# Function that simulates job movements
function simple_random_search(m, pgrid, q, T, r; d = 5)
    # Define some starting values
    # Employment at each firm with full staffing
    p_all_staffed = pgrid[1,2:end]
    J = size(pgrid)[1]-1

    # Employment where ratio r of firms have a negative shock so emp = Lbar
    nfirms = Int(floor(J*r))
    emp_nfirms = 1 .-pgrid[nfirms, 1] # Employment level at nfirms ratio

    p_start = [repeat([m.Lbar], nfirms); 
               repeat([(emp_nfirms - (nfirms * m.Lbar))/nfirms], J - nfirms)]
    p_start = round.(p_start, digits = d)

    # Employment at equilibrium with nfirms understaffed
    p_end = round.(pgrid[nfirms + 1, 2:end], digits = d)

    # Some grids to track flows
    p_flows = [p_all_staffed p_start]
    p_current = p_start
    incr = 1.0 * (10.0^-d)
    rounded_lbar = round(m.Lbar, digits = d)
    t = 1
    switch_draws = Binomial(1, q) # Probability of switching jobs
    while minimum(p_current) >= minimum(p_end) && 
    maximum(p_current) <= maximum(p_end) && t <= T
        # Simulate some market movements
        # Which understaffed firms are still not equilibriated
        under_non_eq_ind = findall(rounded_lbar .>= p_current .> p_end)
        # Which overstaffed firms are still not equilibriated
        over_non_eq_ind = findall(rounded_lbar .<= p_current .< p_end)
        # Which understaffed firms lose people?
        nlosses = rand(switch_draws, length(under_non_eq_ind))
        losses = [i for i in nlosses .* under_non_eq_ind if i != 0]
        # Which overstaffed firms gain people?
        gains = sample(over_non_eq_ind, length(losses), replace = false)

        # Update the market
        p_current[losses] .-= incr
        p_current[gains] .+= incr
        p_flows = [p_flows p_current]
        t += 1
    end
    return p_flows, t
end

# Path to equilibrium with a negative labor shock to a ratio r of firms
nshock_path_05, nshock_t_05 = simple_random_search(m, pgrid, 0.2, 1000, 0.5)
nshock_plot_05 = plot(legend=:bottomright)
xlabel!("Firm index")
ylabel!("Employment")
plot!(nshock_plot_05, nshock_path_05[:,2], label = "Path: t = 1")
plot!(nshock_plot_05, nshock_path_05[:,Int(floor(nshock_t_05/2))], 
        label = "Path: t = " * string(Int(floor(nshock_t_05/2))))
plot!(nshock_plot_05, nshock_path_05[:,end], label = "Path: t = " * string(nshock_t_05))
plot!(nshock_plot_05, pgrid[50, 2:end], label = "Equilibrium")




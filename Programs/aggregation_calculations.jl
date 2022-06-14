using Plots; gr(border = :box, grid = true, minorgrid = true, gridalpha = 0.2,
xguidefontsize = 15, yguidefontsize = 15, xtickfontsize = 13, ytickfontsize=13,
linewidth = 2, gridstyle = :dash, gridlinewidth = 1.2, margin = 10* Plots.px,legendfontsize =12)

include("static_model.jl")

# Some potential utility functions
u1(w)           = w 
u2(w; σ = 0.5)  = (w^(1-σ) - 1)/(1-σ)
u3(w)           = log(w)
u4(w)           = -exp(-w)


## Preliminary results
J       = 100
χ       = 0.9
w       = 1
m       = Mkt(J = J, w = w, χ = χ)
p, s    = choiceProbs(m, false)

nonemp, pgrid, sgrid, shares = checkProbs(m)

# Plot employment for different J
# Nested logit for under-staffed/over-staffed (independence of irrelvant alternatives)?
# Temporary fix: add normalization for share of firms under-staffed = 0
p1 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for j = J:-10:20
    local m = Mkt(J = j, w = w, χ = χ)
    local nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = false)
    plot!(p1, shares, 1 .-nonemp, label = string(j)*" Firms")
end
p1
savefig("plots/vary_J.pdf")

# Plot employment for different J
# With normalization of unemployment rate at 0 understaffed firms = 0.05.
# Note: normalization gets rid of J differences in slope
p1_normalized = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for j = J:-10:20
    local m =Mkt(J = j, w = w, χ = χ)
    local nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = true)
    plot!(p1_normalized, shares, 1 .-nonemp, label = string(j)*" Firms")
end
p1_normalized
savefig("plots/vary_J_normalized.pdf")

# Plot employment for different A(=w)
# Note: Normalization gets rid of wage differences in slope.
p2 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for w = 0.5:0.1:1
    local m = Mkt(J = J, w = w, χ = χ)
    local nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = false)
    plot!(p2, shares, 1 .-nonemp, label ="w = "*string(w))
end
p2
savefig("plots/vary_w.pdf")

# Plot employment for different Lbar
# Normalization just changes value of unemployment -- affects slope.
# Shows how result is sensitive to value of unempoyment.
p3 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for l = 1.0:-0.1:0.5
    local lbar  = l/(J+1)
    local m     = Mkt(J = J, w = w, χ = χ, Lbar = lbar)
    local nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = false)
    plot!(p3, shares, 1 .-nonemp, 
        label = "Lbar = "*string(round(lbar, digits = 3)))
end
p3
savefig("plots/vary_lbar.pdf")

# Plot employment for different χ
# Here, χ differences will still matter.
p4 = plot(legend=:topright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for chi in 0.5:0.1:1
    local m = Mkt(J = J, w = w, χ = chi)
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = true)
    plot!(p4, shares, 1 .-nonemp, label = "χ="*string(chi))
end
p4
savefig("plots/vary_chi.pdf")

# Plot employment for different utility functions
# No difference in employment for different u with normalization
p5 = plot(legend=:topright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for u in [u1, u2, u3, u4]
    local m = Mkt(J = J, w = w, χ = χ, u = u)
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = false)
    plot!(p5, shares, 1 .-nonemp, 
        label = "Utility: "*string(u))
end
p5
savefig("plots/vary_u.pdf")

# Plot employment for different non-emp rates
p6 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for unrate in -0.030:0.01:0.03
    local m = Mkt(J = J, w = w, χ = χ, unrate = (0.4 + unrate))
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = true)
    plot!(p6, shares, 1 .-nonemp, 
        label = "UB: 60% + "*string(-unrate))
end
p6
savefig("plots/vary_ub.pdf")

# Plot employment for different non-emp rates
p7 = plot(legend=:bottomleft)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for unrate in [0.05 0.4]
    local m = Mkt(J = J, w = w, χ = χ, unrate = unrate)
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = true)
    plot!(p7, shares, 1 .-nonemp, 
        label = "Unemployment:  "*string(unrate))
end
p7
savefig("plots/ub_industry_firm.pdf")

# Elasticity calculations
e_len = 1000 # Market size (large)
elasticityMkt = Mkt(J = e_len, w = 1, χ = 0.9) # Define market
e_nonemp, e_pgrid, e_sgrid, e_shares = checkProbs(elasticityMkt)

# For each entry where the understaffed share leads to defined employment,
# Calculate the elasticity by subtracting the employment of the
# Last understaffed firm from the first fully staffed firm
first_defined = findmax(findall(isnan, e_shares))[1] + 1
elasticity = zeros(e_len - first_defined + 1)
for f = first_defined:e_len
    understaffed_emp = e_pgrid[f, f]
    staffed_emp = e_pgrid[f+1, f+1]
    elasticity[f - first_defined + 1] = 
        (staffed_emp - understaffed_emp)/understaffed_emp * 100
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("Elasticity (% ↑ Employment)")
plot!(e_plot, e_shares[first_defined:e_len], elasticity, 
    label = string(e_len)*" Firms")
savefig("plots/elasticity.pdf")



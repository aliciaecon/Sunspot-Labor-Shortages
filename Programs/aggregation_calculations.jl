using Plots, LaTeXStrings; gr(border = :box, grid = true, minorgrid = true, gridalpha = 0.2,
xguidefontsize = 13, yguidefontsize = 13, xtickfontsize = 13, ytickfontsize=13,
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

nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = true)

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
@inbounds for util in [u1, u2, u3, u4]
    local m = Mkt(J = J, w = w, χ = χ, u = util)
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = false)
    plot!(p5, shares, 1 .-nonemp, 
        label = "Utility: "*string(u))
end
p5
savefig("plots/vary_u.pdf")

# Plot employment for different non-emp rates (EPOP = 0.6)
p6 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for u = -0.030:0.01:0.03
    local m = Mkt(J = J, w = w, χ = χ, unrate = (0.4 + u))
    local nonemp, pgrid, sgrid, shares   = checkProbs(m, normalization = true)
    plot!(p6, shares, 1 .-nonemp, label = "Nonemp = "*string(round(0.4 + u, digits = 3)))
end
p6
savefig("plots/vary_ub.pdf")

#= Plot employment for different non-emp rates. 
These correspond to industry (vs aggregate) calculations.
So we make the value of the outside option very high. 
Adjust Lbar to to reflect this. =#
m = Mkt(J = J, w = w, χ = χ, unrate = 0.96, Lbar = 0.7*(1 - 0.96)/J)
nonemp, pgrid, sgrid, shares = checkProbs(m, normalization = true)
plot(shares, 1 .-nonemp, legend =:false)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
savefig("plots/ub_industry_firm.pdf")

# Now, we compute some different "elasticity" measures.
# Define the Market.
e_len   = 1000 # Market size (large)
e_Mkt   = Mkt(J = e_len, w = 1, χ = 0.9) # Aggregate
#e_Mkt   = Mkt(J = e_len, w = w, χ = χ, unrate = 0.96, Lbar = 0.7*(1 - 0.96)/e_len) # Industry
e_nonemp, e_pgrid, e_sgrid, e_shares = checkProbs(e_Mkt)

#= For each entry where the understaffed share leads to a defined eq,
Calculate the % change in the employment of an understaffed firm as 
the share of understaffed firms increases in the economy. =#

first_defined = isempty(findall(isnan, e_shares)) ? 2 : findmax(findall(isnan, e_shares))[1] + 1
pchange       = zeros(e_len - first_defined + 1)
for f = first_defined:e_len
    emp1      = e_pgrid[f, f]
    emp2      = e_pgrid[f+1, f+1]
    pchange[f - first_defined + 1] = (emp2 - emp1)/emp1 * 100
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("% ↑ in Understaffed Employment)")
plot!(e_plot, e_shares[first_defined:e_len], pchange, 
    label = string(e_len)*" Firms")
savefig("plots/pchange_under.pdf")

#=
Now compute an elasticity with respect to this share.
Use forward differences for now.
=#
pchange       = zeros(e_len - first_defined + 1)
for f = first_defined:e_len
    emp1    = e_pgrid[f, f]
    emp2    = e_pgrid[f+1, f+1]
    pchange[f - first_defined + 1] = 
    ((emp2 - emp1)/(e_shares[f+1] - e_shares[f]))*e_shares[f]/emp1
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("Elasticity: Overstaffed Employment")
plot!(e_plot, e_shares[first_defined:e_len], pchange, 
    label = string(e_len)*" Firms")
savefig("plots/elasticity_under.pdf")

#= For each entry where the understaffed share leads to a defined eq,
Calculate the % change in the employment of an overstaffed firm as 
the share of understaffed firms increases in the economy. =#
pchange       = zeros(e_len - first_defined)
for f = first_defined:e_len-1
    emp1  = e_pgrid[f, f+1]
    emp2  = e_pgrid[f+1, f+2]
    pchange[f - first_defined + 1] = (emp2 - emp1)/emp1 * 100
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("% ↑ in Overstaffed Employment")
plot!(e_plot, e_shares[first_defined:e_len-1], pchange, 
    label = string(e_len)*" Firms")
savefig("plots/pchange_over.pdf")

# Now compute an elasticity with respect to this share.
elasticity    = zeros(e_len - first_defined)
for f = first_defined:e_len-1
    emp1  = e_pgrid[f, f+1]
    emp2       = e_pgrid[f+1, f+2]
    elasticity[f - first_defined + 1] = 
        ((emp2 - emp1)/(e_shares[f+1] - e_shares[f]))*e_shares[f]/emp1
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("Elasticity: Overstaffed Employment")
plot!(e_plot, e_shares[first_defined:e_len-1], elasticity, 
    label = string(e_len)*" Firms")
savefig("plots/elasticity_over.pdf")

#= Say you are in a market with a given level of understaffing. 
What is the difference in employment between being understaffed/overstaffed? =#
change        = zeros(e_len - first_defined+1)
for f = first_defined:e_len
    understaffed_emp  = e_pgrid[f, f]
    staffed_emp       = e_pgrid[f, f+1]
    change[f - first_defined + 1] =  (staffed_emp - understaffed_emp)
end

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("Pr((Over) - Pr(Under)")
plot!(e_plot, e_shares[first_defined:e_len], change, 
    label = string(e_len)*" Firms")
savefig("plots/change_over_under.pdf")

# What is the elasticity?
shares     = e_shares[first_defined:e_len]
elasticity = (change[2:end] - change[1:end-1])./(shares[2:end]-shares[1:end-1])
elasticity = elasticity.*shares[1:end-1]./change[1:end-1]

e_plot = plot(legend=:bottomright)
xlabel!("Share of Firms Understaffed")
ylabel!("Elasticity: Pr(Over) - Pr(Under)")
plot!(e_plot, shares[1:end-1], elasticity, 
    label = string(e_len)*" Firms")
savefig("plots/elasticity_over_under.pdf")
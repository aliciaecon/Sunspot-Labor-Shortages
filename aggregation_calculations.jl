"""
Simplest case: w = A, Lbar, χ identical across firms
"""

using Parameters, LinearAlgebra

using Plots; gr(border = :box, grid = true, minorgrid = true, gridalpha = 0.2,
xguidefontsize = 15, yguidefontsize = 15, xtickfontsize = 13, ytickfontsize=13,
linewidth = 2, gridstyle = :dash, gridlinewidth = 1.2, margin = 10* Plots.px,legendfontsize =12)

## Structure that holds all details about the labor market 
struct Mkt
    J       ::Int64
    w       ::Float64
    Lbar    ::Float64
    χ       ::Float64
    u       ::Function
end

# Some potential utility functions
u1(w)           = w 
u2(w; σ = 0.5)  = (w^(1-σ) - 1)/(1-σ)
u3(w)           = log(w)
u4(w)           = -exp(-w)

# Define a constructor for default Lbar, u settings
Mkt(J, w, χ)                 = Mkt(J, w, 0.7/(J+1), χ, u1)
Mkt(J, w, χ, Lbar::Float64)  = Mkt(J, w, Lbar, χ, u1)
Mkt(J, w, χ, u::Function)    = Mkt(J, w, 0.7/(J+1), χ, u)

"""
Get all under/over-staffed combinations of a market of size J,
with heterogenous firms (so which firms in the market
are under-staffed would matter for total employment).
1 = under-staffed, 0 = over-staffed 
"""
function staffCombosDiff(J)
    return reverse.(digits.(0:2^J - 1, base = 2, pad = J))
end

"""
Get under/over-staffed combinations of a market of size J,
for markets with identical firms (since markets
with the same number of under-staffed firms will have
identical employment shares).
1 = under-staffed, 0 = over-staffed.
"""
function staffCombosIdent(J)
    # match staffCombosDiff's output structure
    c = LowerTriangular(ones(Bool, J, J))
    x = [zeros(Bool, 1, J); c]
    s = [x[j,:] for j = 1:J+1]
    return s
end

"""
Calculate the value of the unemployment benefit b such that
nonemployment is normalized to 5% when there are 
no understaffed firms.
(Matches the US natural rate of unemployment of 4.4%)
"""
function findUB(u, w, J; unemp = 0.05)
    denom_sum = exp(u(w))*J
    ub        = (unemp * denom_sum)/(1 - unemp)
    return ub
end

"""
Calculate the choice probability of working at firm j,
given other firms' employment status and wages. 
Note: j = 0 denotes non-employment (w = 0, s = 0),
and b denotes unemployment benefit.
"""
function probWork(u, j, w, χ, sgrid, ub)
    if j == 0 # non-employment
        num      = ub
    else
        num      = exp(u(w) - (χ * sgrid[j]))
    end
    denom        = sum(exp.(u(w) .- (χ .* sgrid))) + ub
    return num/denom
end

"""
For every combination of under/over-staffed firms in a mkt,
compute the choice probabilities p(i -> j) for firm j in each mkt,
where j == 0 denotes non-employment. Each mkt is a row.
"""
function choiceProbs(m::Mkt, normalization; ident = true, b = 0.4)
    @unpack J, w, χ, Lbar, u = m

    if ident == true
        sgrid = staffCombosIdent(J)
        probs = zeros(J+1, J+1)
    else
        sgrid = staffCombosDiff(J)
        probs = zeros(2^J, J+1)
    end
    
    if normalization == true
        ub = findUB(u, w, J) # Find the normalized value of b
    else
        ub = exp(u(b))
    end

    @inbounds for combo = 1:length(sgrid) # combo refers to a version of our market

        probs[combo, 1] = probWork(u, 0, w, χ, sgrid[combo], ub) # non-employment
        
        @inbounds for firm = 1:J
            probs[combo, firm + 1] = probWork(u, firm, w, χ, sgrid[combo], ub) # firm j
        end
    end

    @assert maximum(abs.(sum(probs,dims = 2) .- 1)) < 10^-10

    return probs, sgrid
end

"""
For the simplest case, where all firms are identical,
compute non-employment shares and then check
whether the over/under-staffed assignments agree with Lbar.
"""
function checkProbs(m; normalization = true)
    @unpack J, w, χ, Lbar, u = m

    pgrid, sgrid = choiceProbs(m, normalization)
    shares       = [sum(sgrid[i]) for i = 1:length(sgrid)]./J

    # Check that our under/over-staffed assignments makes sense
    nonemp = pgrid[:, 1]
    pu     = Vector{Any}(undef, J+1)
    po     = Vector{Any}(undef, J+1)

    #= Check P(i -> j | understaffed) < Lbar & P(i -> j | overstaffed) >= Lbar.
    If this does not hold, then we don't compute nonemp shares for that combination.=#
    for j  = 1:J+1
        
        sj = sgrid[j]
        pj = pgrid[j, 2:end]
        
        pu[j] = isempty(pj[sj.==1]) ? missing : first(pj[sj.==1]) 
        po[j] = isempty(pj[sj.==0]) ? missing : first(pj[sj.==0])

        # ismissing if all are understaffed or all are overstaffed
        if (!ismissing(po[j]) && po[j] < m.Lbar) || (!ismissing(pu[j]) && pu[j] >=  m.Lbar)
           nonemp[j] = NaN 
           shares[j] = NaN
        end

    end

    return nonemp, pgrid, sgrid, shares
end

## Preliminary results
J       = 100
χ       = 0.9
w       = 1
m       = Mkt(J, w, χ)
p, s    = choiceProbs(m, false)

nonemp, pgrid, sgrid, shares = checkProbs(m)

# Plot employment for different J
# Nested probit for under-staffed/over-staffed (independence of irrelvant alternatives)?
# Temporary fix: add normalization for share of firms under-staffed = 0
p1 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for j = J:-10:20
    local m = Mkt(j, w, χ)
    local nonemp, pnew, sgrid, shares = checkProbs(m, normalization = false)
    plot!(p1, shares, 1 .-nonemp, label = string(j)*" Firms")
end
p1
savefig("plots/vary_J.pdf")

# Plot employment for different J
# With normalization of unemployment rate at 0 understaffed firms = 0.05
p1_normalized = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for j = J:-10:20
    local m = Mkt(j, w, χ)
    local nonemp, pnew, sgrid, shares = checkProbs(m, normalization = true)
    plot!(p1_normalized, shares, 1 .-nonemp, label = string(j)*" Firms")
end
p1_normalized
savefig("plots/vary_J_normalized.pdf")

# Plot employment for different A(=w)
p2 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for w = 0.5:0.1:1
    local m = Mkt(J, w, χ)
    local nonemp, pnew, sgrid, shares = checkProbs(m)
    plot!(p2, shares, 1 .-nonemp, label ="w = "*string(w))
end
p2
savefig("plots/vary_w.pdf")

# Plot employment for different Lbar
p3 = plot(legend=:outertopright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for l = 1.0:-0.1:0.5
    local lbar  = l/(J+1)
    local m     = Mkt(J, w, χ, lbar)
    local nonemp, pnew, sgrid, shares = checkProbs(m)
    plot!(p3, shares, 1 .-nonemp, 
        label = "Lbar="*string(round(lbar, digits = 3)))
end
p3
savefig("plots/vary_lbar.pdf")

# Plot employment for different χ
p4 = plot(legend=:topright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for chi in 0.5:0.1:1
    local m = Mkt(J, w, chi)
    local nonemp, pnew, sgrid, shares   = checkProbs(m)
    plot!(p4, shares, 1 .-nonemp, label = "χ="*string(chi))
end
p4
savefig("plots/vary_chi.pdf")

# Plot employment for different u
p5 = plot(legend=:topright)
xlabel!("Share of Firms Understaffed")
ylabel!("Employment")
@inbounds for u in [u1, u2, u3, u4]
    local m = Mkt(J, w, χ, u)
    local nonemp, pnew, sgrid, shares   = checkProbs(m)
    plot!(p5, shares, 1 .-nonemp, 
        label = "Utility: "*string(u))
end
p5
savefig("plots/vary_u.pdf")

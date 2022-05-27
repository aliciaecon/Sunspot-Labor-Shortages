"""
Simplest case: w = A, Lbar, χ identical across firms
"""

using Parameters, LinearAlgebra

using Plots; gr(border = :box, grid = true, minorgrid=true, gridalpha=0.2,
xguidefontsize =15, yguidefontsize=15, xtickfontsize=13, ytickfontsize=13,
linewidth = 2, gridstyle = :dash, gridlinewidth = 1.2, margin = 10* Plots.px,legendfontsize =12)

## Structure that holds all details about market 
struct Mkt
    J       ::Int64
    w       ::Float64
    Lbar    ::Float64
    χ       ::Float64
end

# Constructor for default Lbar setting
Mkt(J, w, χ) = Mkt(J, w, 0.7/(J+1), χ)

# Utility for type i working at firm j with wage w
u(w) = w 

"""
Calculate the choice probability of working at firm j,
given other firms' employment status and wages. 
Note: j = 0 denotes non-employment (w = 0, s = 0)
"""
function probWork(u, j, w, χ, sgrid)
    if j == 0 # non-employment
        num      = 1
    else
        num      = exp(u(w) - (χ * sgrid[j]))
    end
    denom        = sum(exp.(u(w) .- (χ .* sgrid))) + 1
    return(num/denom)
end

"""
Get all under/over-staffed combinations of a market of size J,
with heterogenous firms (so which firms in the market
are under-staffed would matter for non-employment).
1 = under-staffed, 0 = over-staffed 
"""
function staffCombosDiff(J)
    return(reverse.(digits.(0:2^J - 1, base = 2, pad = J)))
end

"""
Get under/over-staffed combinations of a market of size J,
for markets with identifcal firms (since markets
with the same number of under-staffed firms will have
identical non-employment shares).
1 = under-staffed, 0 = over-staffed.
"""
function staffCombosIdent(J)
    c = LowerTriangular(ones(Bool, J,J))
    x = [zeros(Bool, 1,J); c]
    s = [x[j,:] for j = 1:J+1]
    return s
end

"""
For every combination of under/over-staffed firms in a mkt,
compute the choice probabilities p(i -> j) for firm j in each mkt,
where j == 0 denotes non-employment. Each mkt is a row.
"""
function choiceProbs(m::Mkt; ident = true)

    @unpack J, w, χ, Lbar = m

    if ident == true
        sgrid = staffCombosIdent(J)
        probs = zeros(J+1, J+1)
    else
        sgrid = staffCombosDiff(J)
        probs = zeros(2^J, J+1)
    end

    @inbounds for combo = 1:length(sgrid) # combo refers to a version of our market

        probs[combo, 1] = probWork(u, 0, w, χ, sgrid[combo]) # non-employment
        
        @inbounds for firm = 1:J
            probs[combo, firm + 1] = probWork(u, firm, w, χ, sgrid[combo]) # firm j
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
function checkProbs(m)

    @unpack J, w, χ, Lbar = m
    pgrid, sgrid          = choiceProbs(m)

    s      = [sum(sgrid[i]) for i = 1:length(sgrid)]./J

    # check that our under/over-staffed assignments makes sense
    nonemp = pgrid[:, 1]
    pu     = Vector{Any}(undef, J+1)
    po     = Vector{Any}(undef, J+1)

    for j  = 1:J+1
        
        sj = sgrid[j]
        pj = pgrid[j, 2:end]
        
        pu[j] = isempty(pj[sj.==1]) ? missing : first( pj[sj.==1]) 
        po[j] = isempty(pj[sj.==0]) ? missing : first( pj[sj.==0])

        if (!ismissing(po[j]) && po[j] < m.Lbar) || (!ismissing(pu[j]) && pu[j] >=  m.Lbar)
           nonemp[j]= NaN 
        end

    end

    return nonemp, pgrid, sgrid, s

end

## Preliminary results
J       = 20
χ       = 0.9
w       = 1
m       = Mkt(J, w, χ)
p, s    = choiceProbs(m)

nonemp, pgrid, sgrid, shares = checkProbs(m)

# Plot non-employment for different J
p1 = plot(legend=:topright, nrows=2)
xlabel!("Share of Firms Understaffed")
ylabel!("Share of Non-employment")
@inbounds for j = 2:4:J
    local m = Mkt(j, w, χ)
    local nonemp, pnew, sgrid, shares     = checkProbs(m)
    plot!(p1, shares, nonemp, label = string(j)*" Firms")
end
p1

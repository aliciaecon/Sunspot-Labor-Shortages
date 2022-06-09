#=
Compute choice probabilities and non-employment shares
for the simplest static model: w = A, Lbar, χ identical across firms.
Eventually can set up as package.
=#

using Parameters, LinearAlgebra

## Structure that holds all details about the labor market 
@with_kw struct Mkt
    J             ::Int64
    w             ::Float64
    Lbar          ::Float64 = 0.7/(J+1)
    χ             ::Float64
    u             ::Function = u(w) = w
    unrate        ::Float64 = 0.4
end

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
Calculate the value of the unemployment benefit b 
based on the unemployment rate with no understaffed firms.
(Default constructor normalizes unemployment to 40% when there are 
no understaffed firms, which matches the 
US Employment-Population ratio of ~60%)
"""
function findUB(u, w, J, unemp)
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
    @unpack J, w, χ, Lbar, u, unrate = m

    if ident == true
        sgrid = staffCombosIdent(J)
        probs = zeros(J+1, J+1)
    else
        sgrid = staffCombosDiff(J)
        probs = zeros(2^J, J+1)
    end
    
    if normalization == true
        ub = findUB(u, w, J, unrate) # Find the normalized value of b
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
    @unpack J, w, χ, Lbar, u, unrate = m

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


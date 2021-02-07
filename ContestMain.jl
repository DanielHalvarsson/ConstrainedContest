# Solving the budget constrained lottery contest with asymmetric valuation, with homogenous budget constraints
# Find "m" given "k"
# Daniel Halvarsson

# Comments:
# if k_star >= 2 then S > 2b and as M increases S such that k_star is binding


using Plots

function generate_budget(n)
    budget = Vector{Any}(undef, n);
    #srand(433421)
    for i = 1:n
        budget[i] = 2/(2-1)
        #*rand()
    end
    # Unsorted valuations
   return budget
end

function generate_utility(n)
    utility = Vector{Any}(undef, n);
    #srand(484777777765391659721973)
    for i = 1:n
        utility[i] = 1/(1-rand()^(1/2))
        # 1-log(rand())/1
        #1/(1-rand()^(1/10))
        #1-log(rand())/0.04
        #(1+rand())
    end
    # Unsorted valuations
    return utility
end

# Number of potential players

n = 7

# Generate budget and utility vector

M1 = [generate_utility(n) generate_budget(n)]

M2 = sortslices(M1, dims=1, rev = true)
display(M2)

# Check for real roots, i.e. the possibility of there being constrained players
RealRoot = zeros(n)

for i = 1:n
    RealRoot[i]= (-M2[i,1])^2-4*M2[i,1]*M2[i,2]
end

M3 = [RealRoot M2] # New vector with real roots >0

# M3 and M4 the same for fixed budget
M4 = sortslices(M3, dims = 1, rev = true) # Sort players collecting constrained players first, however unordered


 # Find the number of potentially constrained players - first k rows in the M matrix
# Vector{Any}(undef, n)
# Matrix{Float64}(undef, 2, 3)
NrPotConst = Vector{Int8}(undef, 1)

NrPotConst = length(filter(x -> x >= 0, RealRoot))

# Calculate the roots to the polynomial for all potentially constrained players
Roots = zeros(NrPotConst,2)

for i = 1:NrPotConst
    Roots[i,1] = (M4[i,2]-sqrt((-M4[i,2]).^2-4*M4[i,2].*M4[i,3]))/2
    Roots[i,2] = (M4[i,2]+sqrt((-M4[i,2]).^2-4*M4[i,2].*M4[i,3]))/2
end

SortRoots = sortslices([Roots M4[1:NrPotConst,:]], dims = 1) ## Sort according to lambda1>lambda2...

# Check whether the nesting property is satisfied

nonested = zeros(Int8,1)
println(" ")
println("Number of potential constrained (k): ", NrPotConst, ", ")

for i = 1:NrPotConst-1
    if SortRoots[i+1,2]>= SortRoots[i,2]
        nonested = 1
        warn("Not nested")
        #println("not nested",", " )
        break
    else
        continue
    end
end
    if NrPotConst > 0
        M4_aug = [zeros(n-NrPotConst,2) M4[NrPotConst+1:n,:]]
        M5 = [vcat(cumsum(ones(NrPotConst)),zeros(n-NrPotConst)) vcat(SortRoots,M4_aug)] # Matrix with ordered potentially constrained
    else

    M4_aug = M4
    M5     = [zeros(n,3) M4_aug]
end

# end module

###################################################################################################################

# M5:
# c1: Nr. potentially constrained
# c2-c3:  Left and right lambda root
# c4: Discrimand
# c5: Valuation
# c6: Budget


LR0 = zeros(n)
Kmat0  = zeros(Int8, n)
harm0  = zeros(n)
value0 = zeros(n)
X0 = zeros(n)
mPair0 = zeros(Int8,n, 1)
Pair0  = zeros(Int8,n, 2)
NrPotConstM = Vector{Any}(undef, 1);
InvVal = zeros(n)
# Module 1. "Find m*"

# Reciprocal valuation
InvVal = 1 ./ M2[:,1]

# Solving for eq.
# j = 0;
S = zeros(n,n)
k_star = 0
m_star = 0 # Vector{Any}(undef, 1);
for j = 0:n-1
    global i = max(j,2)
    while M2[i,1] > S[i-1,j+1] && M2[i,1] > j*M2[1,2]
        # Aggegate expenditures with K = 0
        if i == j
            S[i,j+1] = j*M2[1,2]
        else
            S[i,j+1] = ((i - 1 - j) + sqrt((i - 1 - j)^2 + 4*j*M2[i,2]*sum(InvVal[j+1:i,1])))/(2*sum(InvVal[j+1:i,1]))
        end
        global i += 1
    end
    # Eq. nr. of constrained & active players
    global k_star = j
    global m_star = i - 1
    if M2[j+1,1]-S[i-1,j+1]^2/(S[i-1,j+1]-M2[1,2])<=0 || S[i-1,j+1] < M2[1,2]
        break
    end
end

display(S[1:m_star,1:k_star+1])
display(NrPotConst)
display(k_star)
display(m_star)

# Control eq. conditions
A = M2[m_star,1] > S[m_star,k_star+1]

B = M2[m_star+1,1] <= S[m_star,k_star+1]

if k_star > 0
    C = M2[k_star,1] > S[m_star,k_star+1]^2/(S[m_star,k_star+1] - M2[1,2])

    D = M2[k_star+1,1] <= S[m_star,k_star+1]^2/(S[m_star,k_star+1] - M2[1,2])
else
    C = 1>0;
    D = 1>0;
end
display(NrPotConst>=k_star && A && B && C && D)

# Solving the budget constrained Tullock contest with asymmetric valuation
# Find "m" given "k"
using Plots

if k_star > 0
    equation(x) = x^2/(x-M2[1,2])
    S_k = filter(x -> x > 0, maximum(S, dims = 1))'
    yS_k = zeros(1,k_star+1)
    for i=1:k_star+1
        yS_k[i] = equation(S_k[i])
    end
    trim(equation;val=equation(yS_k[1]))= x -> abs(equation(x))>val ? NaN : equation(x)
    plot(trim(equation), 0, 1+S_k[1])
    scatter!(S_k, yS_k)
end





# IF there are k>1 then it must be the case that S²/(S-b) is decreasing over k,max(m) above 2b

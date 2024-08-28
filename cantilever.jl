
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, ∫vᵢgᵢds, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_cantilever.jl")

ndiv = 4
elements, nodes = import_fem("./msh/cantilever_tri3_"*string(ndiv)*".msh")

nₚ = length(nodes)

L = 48.0
D = 12.0
P = 1000
E = 3e6
# ν = 0.3
ν = 0.4999999
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I

u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν̄)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν̄*y^2*(L-x) + (4+5*ν̄)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν̄)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν̄*y^2 + (4+5*ν̄)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν̄

σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)

# Debug
# n = 3
# u(x,y) = (1+2*x+3*y)^n
# v(x,y) = (4+5*x+6*y)^n
# ∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
# ∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
# ∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
# ∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
# ∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
# ∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
# ∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

# ε₁₁(x,y) = ∂u∂x(x,y)
# ε₂₂(x,y) = ∂v∂y(x,y)
# ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
# ∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
# ∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
# ∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
# ∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
# ∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
# ∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))
# σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
# σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
# σ₃₃(x,y) = E*ν/(1+ν̄)/(1-2*ν)*(ε₁₁(x,y) + ε₂₂(x,y))
# σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
# 𝑝(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3
# ∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
# ∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
# ∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
# ∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
# ∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
# ∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
# b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
# b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
# End debug

prescribe!(elements["Ω"],:E=>(x,y,z)->Ē,index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν̄,index=:𝑔)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E,index=:𝑔)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν,index=:𝑔)
# prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->Ē,index=:𝑔)
# prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν̄,index=:𝑔)
# prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
# prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e9*E,index=:𝑔)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

ops = [
    ∫∫εᵢⱼσᵢⱼdxdy=>elements["Ω"],
    ∫∫vᵢbᵢdxdy=>elements["Ω"],
    ∫vᵢtᵢds=>elements["Γᵗ"]∪elements["Γʳ"],
    ∫vᵢgᵢds=>elements["Γᵍ"],
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](k)
# ops[2](f)
ops[3](f)
ops[4](k,f)

d = k\f

push!(nodes,:d₁=>d[1:2:2*nₚ],:d₂=>d[2:2:2*nₚ])

Hₑ, L₂_𝒖 = Hₑ_PlaneStress(elements["Ωᵍ"])

Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍ"])
# log10(Hₑ)

println(log10(L₂_𝒖))
println(log10(Hₑ))
println(log10(Hₑ_dev))

using LinearAlgebra
using ApproxOperator
using ApproxOperator.Elasticity: ∫εᵢⱼσᵢⱼdΩ, ∫vᵢbᵢdΩ, ∫vᵢtᵢdΓ, ∫vᵢgᵢdΓ, L₂, Hₑ

include("import_block.jl")

ndiv = 4
elements, nodes = import_fem("./msh/block_"*string(ndiv)*".msh")
# elements, nodes = import_fem("./msh/block_hex8_"*string(ndiv)*".msh")

nₚ = length(nodes)

E = 1.0
ν = 0.3

n = 2
# u(x,y,z) = x^n
# v(x,y,z) = 0.0
# w(x,y,z) = 0.0
# ∂u∂x(x,y,z) = n*x^abs(n-1)
# ∂u∂y(x,y,z) = 0.0
# ∂u∂z(x,y,z) = 0.0
# ∂v∂x(x,y,z) = 0.0
# ∂v∂y(x,y,z) = 0.0
# ∂v∂z(x,y,z) = 0.0
# ∂w∂x(x,y,z) = 0.0
# ∂w∂y(x,y,z) = 0.0
# ∂w∂z(x,y,z) = 0.0
# ∂²u∂x²(x,y,z)  = n*(n-1)*x^abs(n-2)
# ∂²u∂y²(x,y,z)  = 0.0
# ∂²u∂z²(x,y,z)  = 0.0
# ∂²u∂x∂y(x,y,z) = 0.0
# ∂²u∂x∂z(x,y,z) = 0.0
# ∂²u∂y∂z(x,y,z) = 0.0
# ∂²v∂x²(x,y,z)  = 0.0
# ∂²v∂y²(x,y,z)  = 0.0
# ∂²v∂z²(x,y,z)  = 0.0
# ∂²v∂x∂y(x,y,z) = 0.0
# ∂²v∂x∂z(x,y,z) = 0.0
# ∂²v∂y∂z(x,y,z) = 0.0
# ∂²w∂x²(x,y,z)  = 0.0
# ∂²w∂y²(x,y,z)  = 0.0
# ∂²w∂z²(x,y,z)  = 0.0
# ∂²w∂x∂y(x,y,z) = 0.0
# ∂²w∂x∂z(x,y,z) = 0.0
# ∂²w∂y∂z(x,y,z) = 0.0

u(x,y,z) = (x+y+z)^n
v(x,y,z) = (x+y+z)^n
w(x,y,z) = (x+y+z)^n
∂u∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂u∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂u∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂v∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂x(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂y(x,y,z) = n*(x+y+z)^abs(n-1)
∂w∂z(x,y,z) = n*(x+y+z)^abs(n-1)
∂²u∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²u∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²v∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
∂²w∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)

ε₁₁(x,y,z) = ∂u∂x(x,y,z)
ε₂₂(x,y,z) = ∂v∂y(x,y,z)
ε₃₃(x,y,z) = ∂w∂z(x,y,z)
ε₁₂(x,y,z) = 0.5*(∂u∂y(x,y,z) + ∂v∂x(x,y,z))
ε₁₃(x,y,z) = 0.5*(∂u∂z(x,y,z) + ∂w∂x(x,y,z))
ε₂₃(x,y,z) = 0.5*(∂v∂z(x,y,z) + ∂w∂y(x,y,z))
∂ε₁₁∂x(x,y,z) = ∂²u∂x²(x,y,z)
∂ε₁₁∂y(x,y,z) = ∂²u∂x∂y(x,y,z)
∂ε₁₁∂z(x,y,z) = ∂²u∂x∂z(x,y,z)
∂ε₂₂∂x(x,y,z) = ∂²v∂x∂y(x,y,z)
∂ε₂₂∂y(x,y,z) = ∂²v∂y²(x,y,z)
∂ε₂₂∂z(x,y,z) = ∂²v∂y∂z(x,y,z)
∂ε₃₃∂x(x,y,z) = ∂²w∂x∂z(x,y,z)
∂ε₃₃∂y(x,y,z) = ∂²w∂y∂z(x,y,z)
∂ε₃₃∂z(x,y,z) = ∂²w∂z²(x,y,z)
∂ε₁₂∂x(x,y,z) = 0.5*(∂²u∂x∂y(x,y,z) + ∂²v∂x²(x,y,z))
∂ε₁₂∂y(x,y,z) = 0.5*(∂²u∂y²(x,y,z) + ∂²v∂x∂y(x,y,z))
∂ε₁₂∂z(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²v∂x∂z(x,y,z))
∂ε₁₃∂x(x,y,z) = 0.5*(∂²u∂x∂z(x,y,z) + ∂²w∂x²(x,y,z))
∂ε₁₃∂y(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
∂ε₁₃∂z(x,y,z) = 0.5*(∂²u∂z²(x,y,z) + ∂²w∂x∂z(x,y,z))
∂ε₂₃∂x(x,y,z) = 0.5*(∂²v∂x∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
∂ε₂₃∂y(x,y,z) = 0.5*(∂²v∂y∂z(x,y,z) + ∂²w∂y²(x,y,z))
∂ε₂₃∂z(x,y,z) = 0.5*(∂²v∂z²(x,y,z) + ∂²w∂y∂z(x,y,z))
σ₁₁(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
σ₂₂(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + (1-ν)*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
σ₃₃(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + (1-ν)*ε₃₃(x,y,z))
σ₁₂(x,y,z) = E/(1+ν)*ε₁₂(x,y,z)
σ₁₃(x,y,z) = E/(1+ν)*ε₁₃(x,y,z)
σ₂₃(x,y,z) = E/(1+ν)*ε₂₃(x,y,z)
𝑝(x,y,z) = (σ₁₁(x,y,z)+σ₂₂(x,y,z)+σ₃₃(x,y,z))/3
∂σ₁₁∂x(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
∂σ₁₁∂y(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
∂σ₁₁∂z(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
∂σ₂₂∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + (1-ν)*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
∂σ₂₂∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + (1-ν)*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
∂σ₂₂∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + (1-ν)*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
∂σ₃₃∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + (1-ν)*∂ε₃₃∂x(x,y,z))
∂σ₃₃∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + (1-ν)*∂ε₃₃∂y(x,y,z))
∂σ₃₃∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + (1-ν)*∂ε₃₃∂z(x,y,z))
∂σ₁₂∂x(x,y,z) = E/(1+ν)*∂ε₁₂∂x(x,y,z)
∂σ₁₂∂y(x,y,z) = E/(1+ν)*∂ε₁₂∂y(x,y,z)
∂σ₁₂∂z(x,y,z) = E/(1+ν)*∂ε₁₂∂z(x,y,z)
∂σ₁₃∂x(x,y,z) = E/(1+ν)*∂ε₁₃∂x(x,y,z)
∂σ₁₃∂y(x,y,z) = E/(1+ν)*∂ε₁₃∂y(x,y,z)
∂σ₁₃∂z(x,y,z) = E/(1+ν)*∂ε₁₃∂z(x,y,z)
∂σ₂₃∂x(x,y,z) = E/(1+ν)*∂ε₂₃∂x(x,y,z)
∂σ₂₃∂y(x,y,z) = E/(1+ν)*∂ε₂₃∂y(x,y,z)
∂σ₂₃∂z(x,y,z) = E/(1+ν)*∂ε₂₃∂z(x,y,z)
b₁(x,y,z) = - ∂σ₁₁∂x(x,y,z) - ∂σ₁₂∂y(x,y,z) - ∂σ₁₃∂z(x,y,z)
b₂(x,y,z) = - ∂σ₁₂∂x(x,y,z) - ∂σ₂₂∂y(x,y,z) - ∂σ₂₃∂z(x,y,z)
b₃(x,y,z) = - ∂σ₁₃∂x(x,y,z) - ∂σ₂₃∂y(x,y,z) - ∂σ₃₃∂z(x,y,z)

prescribe!(elements["Ω"],:E=>(x,y,z)->E)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵍ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ω"],:b₁=>b₁)
prescribe!(elements["Ω"],:b₂=>b₂)
prescribe!(elements["Ω"],:b₃=>b₃)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
prescribe!(elements["Γᵗ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)
prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
prescribe!(elements["Γʳ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)
prescribe!(elements["Γᵍ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
prescribe!(elements["Γᵍ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
prescribe!(elements["Γᵍ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)
prescribe!(elements["Γʳ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γʳ"],:g₁=>u)
prescribe!(elements["Γʳ"],:g₂=>v)
prescribe!(elements["Γʳ"],:g₃=>w)
prescribe!(elements["Γʳ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γʳ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γʳ"],:n₃₃=>(x,y,z)->1.0)
prescribe!(elements["Γʳ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γʳ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γʳ"],:n₂₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γᵗ"],:g₁=>u)
prescribe!(elements["Γᵗ"],:g₂=>v)
prescribe!(elements["Γᵗ"],:g₃=>w)
prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵗ"],:n₃₃=>(x,y,z)->1.0)
prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:n₂₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γᵍ"],:g₁=>u)
prescribe!(elements["Γᵍ"],:g₂=>v)
prescribe!(elements["Γᵍ"],:g₃=>w)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₃₃=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₃=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍ"],:u₁=>u)
prescribe!(elements["Ωᵍ"],:u₂=>v)
prescribe!(elements["Ωᵍ"],:u₃=>w)
prescribe!(elements["Ωᵍ"],:∂u₁∂x=>∂u∂x)
prescribe!(elements["Ωᵍ"],:∂u₁∂y=>∂u∂y)
prescribe!(elements["Ωᵍ"],:∂u₁∂z=>∂u∂z)
prescribe!(elements["Ωᵍ"],:∂u₂∂x=>∂v∂x)
prescribe!(elements["Ωᵍ"],:∂u₂∂y=>∂v∂y)
prescribe!(elements["Ωᵍ"],:∂u₂∂z=>∂v∂z)
prescribe!(elements["Ωᵍ"],:∂u₃∂x=>∂w∂x)
prescribe!(elements["Ωᵍ"],:∂u₃∂y=>∂w∂y)
prescribe!(elements["Ωᵍ"],:∂u₃∂z=>∂w∂z)

𝑎 = ∫εᵢⱼσᵢⱼdΩ=>elements["Ω"]
𝑓 = [
    ∫vᵢbᵢdΩ=>elements["Ω"],
    ∫vᵢtᵢdΓ=>elements["Γᵗ"]∪elements["Γʳ"],
    # ∫vᵢtᵢdΓ=>elements["Γᵗ"]∪elements["Γʳ"]∪elements["Γᵍ"],
]
𝑎ᵅ = ∫vᵢgᵢdΓ=>elements["Γᵍ"]
# 𝑎ᵅ = ∫vᵢgᵢdΓ=>elements["Γᵍ"]∪elements["Γᵗ"]∪elements["Γʳ"]

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(k,f)

d = k\f

push!(nodes,:d₁=>d[1:3:end],:d₂=>d[2:3:end],:d₃=>d[3:3:end])

# L₂_𝒖 = L₂(elements["Ωᵍ"])
Hₑ_𝒖, L₂_𝒖 = Hₑ(elements["Ωᵍ"])

println(log10(L₂_𝒖))

# d = zeros(3*nₚ)
# for i in 1:nₚ
#     x = nodes[i].x
#     y = nodes[i].y
#     z = nodes[i].z
#     d[3*i-2] = u(x,y,z)
#     d[3*i-1] = v(x,y,z)
#     d[3*i]   = w(x,y,z)
# end

# err = k*d .- f
# norm(err)
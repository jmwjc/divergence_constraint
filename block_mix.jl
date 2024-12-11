
using TimerOutputs 
using Pardiso
using SparseArrays, LinearAlgebra
using SharedArrays, Distributed
using WriteVTK
using ApproxOperator
using ApproxOperator.Elasticity: ∫qpdΩ, ∫εᵈᵢⱼσᵈᵢⱼdΩ, ∫p∇udΩ, ∫vᵢbᵢdΩ, ∫vᵢtᵢdΓ, ∫vᵢgᵢdΓ, Hₑ

# addprocs(3)
# println(nprocs())
println(Threads.nthreads())

include("import_block.jl")

const to = TimerOutput()
ps = MKLPardisoSolver()

ndiv = 4
ndiv_p = 4
# poly = "tet4"
poly = "hex8"
@timeit to "import data" begin
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/block_"*string(ndiv)*".msh","./msh/block_"*string(ndiv_p)*".msh",ndiv_p)
elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/block_hex8_"*string(ndiv)*".msh","./msh/block_"*string(ndiv_p)*".msh",ndiv_p)
end

nᵤ = length(nodes)
nₚ = length(nodes_p)

E = 240.56839
ν = 0.5-1e-8
P = 80.0

n₁₁(n₁,n₂,n₃) = n₃ ≈ 1.0 || n₁ ≈ -1.0 ? 1.0 : 0.0
n₂₂(n₁,n₂,n₃) = n₃ ≈ 1.0 || n₂ ≈ -1.0 ? 1.0 : 0.0
n₃₃(n₁,n₂,n₃) = n₃ ≈ -1.0 ? 1.0 : 0.0
prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:t₃=>(x,y,z)->-P)
prescribe!(elements["Γᵗ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γᵗ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:g₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z,n₁,n₂,n₃)->1.0)
prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z,n₁,n₂,n₃)->1.0)
prescribe!(elements["Γᵗ"],:n₃₃=>(x,y,z,n₁,n₂,n₃)->0.0)
prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵗ"],:n₂₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e12*E)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z,n₁,n₂,n₃)->n₁₁(n₁,n₂,n₃))
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z,n₁,n₂,n₃)->n₂₂(n₁,n₂,n₃))
prescribe!(elements["Γᵍ"],:n₃₃=>(x,y,z,n₁,n₂,n₃)->n₃₃(n₁,n₂,n₃))
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₃=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₃=>(x,y,z)->0.0)

# n = 2
# u(x,y,z) = (x+y+z)^n
# v(x,y,z) = (x+y+z)^n
# w(x,y,z) = (x+y+z)^n
# ∂u∂x(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂u∂y(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂u∂z(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂v∂x(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂v∂y(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂v∂z(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂w∂x(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂w∂y(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂w∂z(x,y,z) = n*(x+y+z)^abs(n-1)
# ∂²u∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²u∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²u∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²u∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²u∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²u∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²v∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂x²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂y²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂z²(x,y,z)  = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂x∂y(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂x∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)
# ∂²w∂y∂z(x,y,z) = n*(n-1)*(x+y+z)^abs(n-2)

# ε₁₁(x,y,z) = ∂u∂x(x,y,z)
# ε₂₂(x,y,z) = ∂v∂y(x,y,z)
# ε₃₃(x,y,z) = ∂w∂z(x,y,z)
# ε₁₂(x,y,z) = 0.5*(∂u∂y(x,y,z) + ∂v∂x(x,y,z))
# ε₁₃(x,y,z) = 0.5*(∂u∂z(x,y,z) + ∂w∂x(x,y,z))
# ε₂₃(x,y,z) = 0.5*(∂v∂z(x,y,z) + ∂w∂y(x,y,z))
# ∂ε₁₁∂x(x,y,z) = ∂²u∂x²(x,y,z)
# ∂ε₁₁∂y(x,y,z) = ∂²u∂x∂y(x,y,z)
# ∂ε₁₁∂z(x,y,z) = ∂²u∂x∂z(x,y,z)
# ∂ε₂₂∂x(x,y,z) = ∂²v∂x∂y(x,y,z)
# ∂ε₂₂∂y(x,y,z) = ∂²v∂y²(x,y,z)
# ∂ε₂₂∂z(x,y,z) = ∂²v∂y∂z(x,y,z)
# ∂ε₃₃∂x(x,y,z) = ∂²w∂x∂z(x,y,z)
# ∂ε₃₃∂y(x,y,z) = ∂²w∂y∂z(x,y,z)
# ∂ε₃₃∂z(x,y,z) = ∂²w∂z²(x,y,z)
# ∂ε₁₂∂x(x,y,z) = 0.5*(∂²u∂x∂y(x,y,z) + ∂²v∂x²(x,y,z))
# ∂ε₁₂∂y(x,y,z) = 0.5*(∂²u∂y²(x,y,z) + ∂²v∂x∂y(x,y,z))
# ∂ε₁₂∂z(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²v∂x∂z(x,y,z))
# ∂ε₁₃∂x(x,y,z) = 0.5*(∂²u∂x∂z(x,y,z) + ∂²w∂x²(x,y,z))
# ∂ε₁₃∂y(x,y,z) = 0.5*(∂²u∂y∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
# ∂ε₁₃∂z(x,y,z) = 0.5*(∂²u∂z²(x,y,z) + ∂²w∂x∂z(x,y,z))
# ∂ε₂₃∂x(x,y,z) = 0.5*(∂²v∂x∂z(x,y,z) + ∂²w∂x∂y(x,y,z))
# ∂ε₂₃∂y(x,y,z) = 0.5*(∂²v∂y∂z(x,y,z) + ∂²w∂y²(x,y,z))
# ∂ε₂₃∂z(x,y,z) = 0.5*(∂²v∂z²(x,y,z) + ∂²w∂y∂z(x,y,z))
# σ₁₁(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
# σ₂₂(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + (1-ν)*ε₂₂(x,y,z) + ν*ε₃₃(x,y,z))
# σ₃₃(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y,z) + ν*ε₂₂(x,y,z) + (1-ν)*ε₃₃(x,y,z))
# σ₁₂(x,y,z) = E/(1+ν)*ε₁₂(x,y,z)
# σ₁₃(x,y,z) = E/(1+ν)*ε₁₃(x,y,z)
# σ₂₃(x,y,z) = E/(1+ν)*ε₂₃(x,y,z)
# 𝑝(x,y,z) = (σ₁₁(x,y,z)+σ₂₂(x,y,z)+σ₃₃(x,y,z))/3
# ∂σ₁₁∂x(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
# ∂σ₁₁∂y(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
# ∂σ₁₁∂z(x,y,z) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
# ∂σ₂₂∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + (1-ν)*∂ε₂₂∂x(x,y,z) + ν*∂ε₃₃∂x(x,y,z))
# ∂σ₂₂∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + (1-ν)*∂ε₂₂∂y(x,y,z) + ν*∂ε₃₃∂y(x,y,z))
# ∂σ₂₂∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + (1-ν)*∂ε₂₂∂z(x,y,z) + ν*∂ε₃₃∂z(x,y,z))
# ∂σ₃₃∂x(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y,z) + ν*∂ε₂₂∂x(x,y,z) + (1-ν)*∂ε₃₃∂x(x,y,z))
# ∂σ₃₃∂y(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y,z) + ν*∂ε₂₂∂y(x,y,z) + (1-ν)*∂ε₃₃∂y(x,y,z))
# ∂σ₃₃∂z(x,y,z) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂z(x,y,z) + ν*∂ε₂₂∂z(x,y,z) + (1-ν)*∂ε₃₃∂z(x,y,z))
# ∂σ₁₂∂x(x,y,z) = E/(1+ν)*∂ε₁₂∂x(x,y,z)
# ∂σ₁₂∂y(x,y,z) = E/(1+ν)*∂ε₁₂∂y(x,y,z)
# ∂σ₁₂∂z(x,y,z) = E/(1+ν)*∂ε₁₂∂z(x,y,z)
# ∂σ₁₃∂x(x,y,z) = E/(1+ν)*∂ε₁₃∂x(x,y,z)
# ∂σ₁₃∂y(x,y,z) = E/(1+ν)*∂ε₁₃∂y(x,y,z)
# ∂σ₁₃∂z(x,y,z) = E/(1+ν)*∂ε₁₃∂z(x,y,z)
# ∂σ₂₃∂x(x,y,z) = E/(1+ν)*∂ε₂₃∂x(x,y,z)
# ∂σ₂₃∂y(x,y,z) = E/(1+ν)*∂ε₂₃∂y(x,y,z)
# ∂σ₂₃∂z(x,y,z) = E/(1+ν)*∂ε₂₃∂z(x,y,z)
# b₁(x,y,z) = - ∂σ₁₁∂x(x,y,z) - ∂σ₁₂∂y(x,y,z) - ∂σ₁₃∂z(x,y,z)
# b₂(x,y,z) = - ∂σ₁₂∂x(x,y,z) - ∂σ₂₂∂y(x,y,z) - ∂σ₂₃∂z(x,y,z)
# b₃(x,y,z) = - ∂σ₁₃∂x(x,y,z) - ∂σ₂₃∂y(x,y,z) - ∂σ₃₃∂z(x,y,z)

# prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
# prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
# prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
# prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
# prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E)
# prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν)
# prescribe!(elements["Ωᵘ"],:b₁=>b₁)
# prescribe!(elements["Ωᵘ"],:b₂=>b₂)
# prescribe!(elements["Ωᵘ"],:b₃=>b₃)
# prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
# prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
# prescribe!(elements["Γᵗ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)
# prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂,n₃)->σ₁₁(x,y,z)*n₁+σ₁₂(x,y,z)*n₂+σ₁₃(x,y,z)*n₃)
# prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂,n₃)->σ₁₂(x,y,z)*n₁+σ₂₂(x,y,z)*n₂+σ₂₃(x,y,z)*n₃)
# prescribe!(elements["Γʳ"],:t₃=>(x,y,z,n₁,n₂,n₃)->σ₁₃(x,y,z)*n₁+σ₂₃(x,y,z)*n₂+σ₃₃(x,y,z)*n₃)
# prescribe!(elements["Γʳ"],:α=>(x,y,z)->1e12*E)
# prescribe!(elements["Γʳ"],:g₁=>u)
# prescribe!(elements["Γʳ"],:g₂=>v)
# prescribe!(elements["Γʳ"],:g₃=>w)
# prescribe!(elements["Γʳ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γʳ"],:n₂₂=>(x,y,z)->1.0)
# prescribe!(elements["Γʳ"],:n₃₃=>(x,y,z)->1.0)
# prescribe!(elements["Γʳ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Γʳ"],:n₁₃=>(x,y,z)->0.0)
# prescribe!(elements["Γʳ"],:n₂₃=>(x,y,z)->0.0)
# prescribe!(elements["Γᵗ"],:α=>(x,y,z)->1e9*E)
# prescribe!(elements["Γᵗ"],:g₁=>u)
# prescribe!(elements["Γᵗ"],:g₂=>v)
# prescribe!(elements["Γᵗ"],:g₃=>w)
# prescribe!(elements["Γᵗ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γᵗ"],:n₂₂=>(x,y,z)->1.0)
# prescribe!(elements["Γᵗ"],:n₃₃=>(x,y,z)->1.0)
# prescribe!(elements["Γᵗ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Γᵗ"],:n₁₃=>(x,y,z)->0.0)
# prescribe!(elements["Γᵗ"],:n₂₃=>(x,y,z)->0.0)
# prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e12*E)
# prescribe!(elements["Γᵍ"],:g₁=>u)
# prescribe!(elements["Γᵍ"],:g₂=>v)
# prescribe!(elements["Γᵍ"],:g₃=>w)
# prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍ"],:n₃₃=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Γᵍ"],:n₁₃=>(x,y,z)->0.0)
# prescribe!(elements["Γᵍ"],:n₂₃=>(x,y,z)->0.0)
# prescribe!(elements["Ωᵍᵘ"],:u₁=>u)
# prescribe!(elements["Ωᵍᵘ"],:u₂=>v)
# prescribe!(elements["Ωᵍᵘ"],:u₃=>w)
# prescribe!(elements["Ωᵍᵘ"],:∂u₁∂x=>∂u∂x)
# prescribe!(elements["Ωᵍᵘ"],:∂u₁∂y=>∂u∂y)
# prescribe!(elements["Ωᵍᵘ"],:∂u₁∂z=>∂u∂z)
# prescribe!(elements["Ωᵍᵘ"],:∂u₂∂x=>∂v∂x)
# prescribe!(elements["Ωᵍᵘ"],:∂u₂∂y=>∂v∂y)
# prescribe!(elements["Ωᵍᵘ"],:∂u₂∂z=>∂v∂z)
# prescribe!(elements["Ωᵍᵘ"],:∂u₃∂x=>∂w∂x)
# prescribe!(elements["Ωᵍᵘ"],:∂u₃∂y=>∂w∂y)
# prescribe!(elements["Ωᵍᵘ"],:∂u₃∂z=>∂w∂z)

𝑎ᵘ = ∫εᵈᵢⱼσᵈᵢⱼdΩ=>elements["Ωᵘ"]
𝑎ᵖ = ∫qpdΩ=>elements["Ωᵖ"]
𝑏ᵖ = ∫p∇udΩ=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑓 = [
    # ∫vᵢbᵢdΩ=>elements["Ωᵘ"],
    ∫vᵢtᵢdΓ=>elements["Γᵗ"],
    # ∫vᵢtᵢdΓ=>elements["Γᵗ"]∪elements["Γʳ"],
]
𝑎ᵅ = ∫vᵢgᵢdΓ=>elements["Γᵍ"]∪elements["Γᵗ"]
# 𝑎ᵅ = ∫vᵢgᵢdΓ=>elements["Γᵍ"]∪elements["Γᵗ"]∪elements["Γʳ"]

kᵘᵘ = zeros(3*nᵤ,3*nᵤ)
kᵖᵖ = zeros(nₚ,nₚ)
kᵖᵘ = zeros(nₚ,3*nᵤ)
fᵖ = zeros(nₚ)
fᵘ = zeros(3*nᵤ)

# kᵘᵘ = SharedMatrix{Float64}(3*nᵤ,3*nᵤ)
# kᵖᵖ = SharedMatrix{Float64}(nₚ,nₚ)
# kᵖᵘ = SharedMatrix{Float64}(nₚ,3*nᵤ)
# fᵖ  = SharedVector{Float64}(nₚ)
# fᵘ  = SharedVector{Float64}(3*nᵤ)

@timeit to "assembly" begin
𝑎ᵘ(kᵘᵘ)
𝑎ᵖ(kᵖᵖ)
𝑏ᵖ(kᵖᵘ)
𝑎ᵅ(kᵘᵘ,fᵘ)
𝑓(fᵘ)
end

k =sparse([-kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ])
# k = [-kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [-fᵘ;fᵖ]
d = zeros(3*nᵤ+nₚ)

set_matrixtype!(ps, -2)
k = get_matrix(ps,k,:N)
@timeit to "solve" pardiso(ps,d,k,f)
# d = k\f

𝑢₁ = d[1:3:3*nᵤ]
𝑢₂ = d[2:3:3*nᵤ]
𝑢₃ = d[3:3:3*nᵤ]
𝑝 = d[3*nᵤ+1:3*nᵤ+nₚ]
push!(nodes,:u₁=>𝑢₁,:u₂=>𝑢₂,:u₃=>𝑢₃)
push!(nodes_p,:p=>𝑝)
# Hₑ_𝒖, L₂_𝒖 = Hₑ(elements["Ωᵍᵘ"])


# println(log10(L₂_𝒖))
# println(log10(Hₑ_𝒖))

colors = zeros(nᵤ)
𝗠 = zeros(10)
for (i,node) in enumerate(nodes)
    x = node.x
    y = node.y
    z = node.z
    indices = sp(x,y,z)
    ni = length(indices)
    𝓒 = [nodes_p[i] for i in indices]
    data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[z]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
    ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
    𝓖 = [ξ]
    a = type(𝓒,𝓖)
    set𝝭!(a)
    p = 0.0
    N = ξ[:𝝭]
    for (k,xₖ) in enumerate(𝓒)
        p += N[k]*xₖ.p
    end
    colors[i] = p
end
α = 1.0
points = [[node.x+α*node.u₁ for node in nodes]';[node.y+α*node.u₂ for node in nodes]';[node.z+α*node.u₃ for node in nodes]']
# cells = [MeshCell(VTKCellTypes.VTK_TETRA,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON,[xᵢ.𝐼 for xᵢ in elm.𝓒]) for elm in elements["Ωᵘ"]]
vtk_grid("./vtk/block_"*poly*"_"*string(ndiv)*"_"*string(nₚ),points,cells) do vtk
    vtk["u"] = (𝑢₁,𝑢₂,𝑢₃)
    vtk["𝑝"] = colors
end

println(nodes[5])

show(to)
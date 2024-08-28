
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵢⱼσᵢⱼdxdy, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, ∫vᵢgᵢds, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_patchtest.jl")

ndiv = 2
elements, nodes = import_patchtest_fem("./msh/patchtest_"*string(ndiv)*".msh")

nₚ = length(nodes)

n = 1
u(x,y) = (1+2*x+3*y)^n
v(x,y) = (4+5*x+6*y)^n
∂u∂x(x,y) = 2*n*(1+2*x+3*y)^abs(n-1)
∂u∂y(x,y) = 3*n*(1+2*x+3*y)^abs(n-1)
∂v∂x(x,y) = 5*n*(4+5*x+6*y)^abs(n-1)
∂v∂y(x,y) = 6*n*(4+5*x+6*y)^abs(n-1)
∂²u∂x²(x,y)  = 4*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂x∂y(x,y) = 6*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²u∂y²(x,y)  = 9*n*(n-1)*(1+2*x+3*y)^abs(n-2)
∂²v∂x²(x,y)  = 25*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂x∂y(x,y) = 30*n*(n-1)*(4+5*x+6*y)^abs(n-2)
∂²v∂y²(x,y)  = 36*n*(n-1)*(4+5*x+6*y)^abs(n-2)

Ē = 1.0
# ν̄ = 0.4999999
ν̄ = 0.3
E = Ē/(1.0-ν̄^2)
ν = ν̄/(1.0-ν̄)

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))
σ₁₁(x,y) = Ē/(1+ν̄)/(1-2*ν̄)*((1-ν̄)*ε₁₁(x,y) + ν̄*ε₂₂(x,y))
σ₂₂(x,y) = Ē/(1+ν̄)/(1-2*ν̄)*(ν̄*ε₁₁(x,y) + (1-ν̄)*ε₂₂(x,y))
σ₃₃(x,y) = Ē*ν̄/(1+ν̄)/(1-2*ν̄)*(ε₁₁(x,y) + ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
𝑝(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3
∂σ₁₁∂x(x,y) = E/(1-ν^2)*(∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1-ν^2)*(∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂x(x,y) + ∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1-ν^2)*(ν*∂ε₁₁∂y(x,y) + ∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)

prescribe!(elements["Ω"],:E=>(x,y,z)->E,index=:𝑔)
prescribe!(elements["Ω"],:ν=>(x,y,z)->ν,index=:𝑔)
prescribe!(elements["Ω"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ω"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ¹"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γ¹"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γ²"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γ²"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γ³"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γ³"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
prescribe!(elements["Γ⁴"],:α=>(x,y,z)->1e9,index=:𝑔)
prescribe!(elements["Γ⁴"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ⁴"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

ops = [
    ∫∫εᵢⱼσᵢⱼdxdy=>elements["Ω"],
    ∫∫vᵢbᵢdxdy=>elements["Ω"],
    ∫vᵢtᵢds=>elements["Γ¹"],
    ∫vᵢtᵢds=>elements["Γ²"],
    ∫vᵢtᵢds=>elements["Γ³"],
    ∫vᵢgᵢds=>elements["Γ⁴"],
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](k)
ops[2:5](f)
ops[6](k,f)

d = k\f

push!(nodes,:d₁=>d[1:2:2*nₚ],:d₂=>d[2:2:2*nₚ])

Hₑ, L₂_𝒖 = Hₑ_PlaneStress(elements["Ω"])

Hₑ_dev = Hₑ_PlaneStrain_Deviatoric(elements["Ω"])

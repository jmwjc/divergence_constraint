
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫∫qpdxdy, ∫∫δsᵢⱼsᵢⱼdxdy, ∫∫p∇udxdy, ∫∫sᵢⱼεᵢⱼdxdy, ∫pnᵢgᵢds, ∫sᵢⱼnⱼgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_patchtest.jl")

ndiv = 8
nₚ = 49
elements, nodes, nodes_p = import_patchtest_elasticity_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_c_"*string(nₚ)*".msh")
# elements, nodes, nodes_p = import_patchtest_elasticity_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_"*string(ndiv)*".msh")

nₑ = length(elements["Ωᵘ"])
nₛ = 3
nᵤ = length(nodes)

E = 1.0
ν = 0.3

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

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = 0.5*(∂u∂y(x,y) + ∂v∂x(x,y))
σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

prescribe!(elements["Ωᵇ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵇ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Ωᵇ"],:b₁=>(x,y,z)->b₁(x,y))
prescribe!(elements["Ωᵇ"],:b₂=>(x,y,z)->b₂(x,y))
prescribe!(elements["Γ¹ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ¹ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ²ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ³ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ⁴ᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴ᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γ¹ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ¹ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ¹ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ²ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ²ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ³ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ³ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ³ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γ⁴ᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->p(x,y))

𝑎ᵇ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>elements["Ωᵇ"]
𝑎ˢ = ∫∫δsᵢⱼsᵢⱼdxdy=>elements["Ωˢ"]
𝑎ᵖ = ∫∫qpdxdy=>elements["Ωᵖ"]
𝑏ˢ = ∫∫sᵢⱼεᵢⱼdxdy=>(elements["Ωˢ"],elements["Ωᵘ"])
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑏ˢᵇ = ∫∫sᵢⱼεᵢⱼdxdy=>(elements["Ωˢ"],elements["Ωᵇ"])
𝑏ᵖᵇ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵇ"])
𝑏ᵇᵘ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>(elements["Ωᵇ"],elements["Ωᵘ"])
𝑏ˢᵅ = ∫sᵢⱼnⱼgᵢds=>(elements["Γˢ"],elements["Γᵘ"])
𝑏ᵖᵅ = ∫pnᵢgᵢds=>(elements["Γᵖ"],elements["Γᵘ"])
𝑓 = ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
𝑓ᵇ = ∫∫vᵢbᵢdxdy=>elements["Ωᵇ"]

kˢˢ = zeros(4*nₛ*nₑ,4*nₛ*nₑ)
kᵖᵖ = zeros(nₚ,nₚ)
kᵇᵇ = zeros(2*nₑ,2*nₑ)
kˢᵘ = zeros(4*nₛ*nₑ,2*nᵤ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
kˢᵇ = zeros(4*nₛ*nₑ,2*nₑ)
kᵖᵇ = zeros(nₚ,2*nₑ)
kᵇᵘ = zeros(2*nₑ,2*nᵤ)
fˢ = zeros(4*nₛ*nₑ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)
fᵇ = zeros(2*nₑ)

𝑎ˢ(kˢˢ)
𝑎ᵖ(kᵖᵖ)
𝑎ᵇ(kᵇᵇ)
𝑏ˢ(kˢᵘ)
𝑏ᵖ(kᵖᵘ)
# 𝑏ˢ(kˢᵇ)
𝑏ᵖᵇ(kᵖᵇ)
# 𝑏ᵇ(kᵇᵘ)
𝑏ˢᵅ(kˢᵘ,fˢ)
𝑏ᵖᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)
𝑓ᵇ(fᵇ)

d = [
    zeros(2*nᵤ,2*nᵤ)            kᵖᵘ'             kˢᵘ' kᵇᵘ';
                kᵖᵘ             kᵖᵖ  zeros(nₚ,4*nₛ*nₑ) kᵖᵇ;
                kˢᵘ zeros(4*nₛ*nₑ,nₚ)             kˢˢ  kˢᵇ;
                kᵇᵘ             kᵖᵇ'             kˢᵇ' kᵇᵇ
    ]\[fᵘ;fᵖ;fˢ;fᵇ]
# d = [zeros(2*nᵤ,2*nᵤ) kᵖᵘ' kˢᵘ';kᵖᵘ kᵖᵖ zeros(nₚ,4*nₛ*nₑ);kˢᵘ zeros(4*nₛ*nₑ,nₚ) kˢˢ]\[fᵘ;fᵖ;fˢ]

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

L₂_𝑢 = L₂(elements["Ωᵍᵘ"])
L₂_𝒑 = L₂𝑝(elements["Ωᵍᵖ"])

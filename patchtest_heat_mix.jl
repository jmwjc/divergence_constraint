
using BenchmarkTools
using ApproxOperator
using ApproxOperator.Heat: ∫∫qᵢpᵢdxdy, ∫pᵢnᵢuds, ∫∫∇𝒑udxdy, ∫pᵢnᵢgⱼds, ∫vbdΩ, ∫vgdΓ, L₂, L₂𝒑, H₁

include("import_patchtest.jl")

ndiv = 32
nᵤ = 1885
elements, nodes, nodes_u = import_patchtest_mix("./msh/patchtest_"*string(ndiv)*".msh","./msh/patchtest_c_"*string(nᵤ)*".msh")
# elements, nodes, nodes_u = import_patchtest_mix("./msh/patchtest_quad_"*string(ndiv)*".msh","./msh/patchtest_c_"*string(nᵤ)*".msh")

nₚ = length(nodes)

n = 5
u(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

prescribe!(elements["Ωᵘ"],:b=>b)
prescribe!(elements["Γ¹ᵘ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²ᵘ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³ᵘ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴ᵘ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂z=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍᵖ"],:𝑝₁=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵖ"],:𝑝₂=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵖ"],:𝑝₃=>(x,y,z)->0.0)

𝑎 = ∫∫qᵢpᵢdxdy=>elements["Ωᵖ"]
𝑏 = [
    ∫pᵢnᵢuds=>(elements["∂Ωᵖ"],elements["∂Ωᵘ"]),
    ∫∫∇𝒑udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"]),
]
𝑏ᵅ = ∫pᵢnᵢgⱼds=>(elements["Γᵖ"],elements["Γᵘ"])
𝑓 = ∫vbdΩ=>elements["Ωᵘ"]

kᵖᵖ = zeros(2*nₚ,2*nₚ)
fᵖ = zeros(2*nₚ)
kᵖᵘ = zeros(2*nₚ,nᵤ)
fᵘ = zeros(nᵤ)

𝑎(kᵖᵖ)
𝑏(kᵖᵘ)
𝑏ᵅ(kᵖᵘ,fᵖ)
𝑓(fᵘ)

d = [kᵖᵖ kᵖᵘ;kᵖᵘ' zeros(nᵤ,nᵤ)]\[fᵖ;-fᵘ]

𝑝₁ = d[1:2:2*nₚ]
𝑝₂ = d[2:2:2*nₚ]
𝑢 = d[2*nₚ+1:end]
push!(nodes,:p₁=>𝑝₁)
push!(nodes,:p₂=>𝑝₂)
push!(nodes,:p₃=>zeros(nₚ))
push!(nodes_u,:d=>𝑢)

L₂_𝒑 = L₂𝒑(elements["Ωᵍᵖ"])
H₁_𝑢, L₂_𝑢 = H₁(elements["Ωᵍᵘ"])
println(log10(L₂_𝑢))
println(log10(H₁_𝑢))
println(log10(L₂_𝒑))

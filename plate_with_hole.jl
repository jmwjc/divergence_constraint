
using BenchmarkTools
using ApproxOperator
using ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫vtdΓ, ∫vbdΩ, ∫vgdΓ, H₁

include("import_plate_with_hole.jl")

ndiv = 2
poly = "tri3"
elements, nodes = import_fem("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh")

nₚ = length(nodes)

n = 1
u(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

prescribe!(elements["Ω"],:k=>(x,y,z)->1.0)
prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e9)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵗ"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ω"],:∂u∂z=>(x,y,z)->0.0)

𝑎 = ∫∫∇v∇udxdy=>elements["Ω"]
𝑓 = [
    ∫vbdΩ=>elements["Ω"],
    ∫vtdΓ=>elements["Γᵗ"],
]
𝑎ᵅ = ∫vgdΓ=>elements["Γᵍ"]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(k,f)

d = k\f

push!(nodes,:d=>d)

Hₑ, L₂ = H₁(elements["Ω"][50:50])

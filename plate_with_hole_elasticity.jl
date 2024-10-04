
using BenchmarkTools
using ApproxOperator
using ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫vtdΓ, ∫vbdΩ, ∫vgdΓ, H₁

include("import_plate_with_hole.jl")

ndiv = 4
poly = "tri3"
elements, nodes = import_fem("./msh/plate_with_hole_"*poly*"_"*string(ndiv)*".msh")

nₚ = length(nodes)

# n = 1
# u(x,y) = (x+y)^n
# ∂u∂x(x,y) = n*(x+y)^abs(n-1)
# ∂u∂y(x,y) = n*(x+y)^abs(n-1)
# ∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
# ∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
# b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

r(x,y) = (x^2+y^2)^0.5
θ(x,y) = atan(y/x)
u(x,y) = (r(x,y) + 1/r(x,y))cos(θ(x,y))
∂u∂x(x,y) = 1 - 1/r(x,y)^2 + 2*1/r(x,y)^2 * sin(θ(x,y))^2
∂u∂y(x,y) = - 2/r(x,y)^2 * sin(θ(x,y))*cos(θ(x,y))

prescribe!(elements["Ω"],:k=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:α=>(x,y,z)->1e9)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵗ"],:t=>(x,y,z,n₁,n₂)->∂u∂x(x,y)*n₁ + ∂u∂y(x,y)*n₂)
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ω"],:∂u∂z=>(x,y,z)->0.0)

𝑎 = ∫∫∇v∇udxdy=>elements["Ω"]
𝑓 = ∫vtdΓ=>elements["Γᵗ"]
𝑎ᵅ = ∫vgdΓ=>elements["Γᵍ"]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

𝑎(k)
𝑓(f)
𝑎ᵅ(k,f)

d = k\f

push!(nodes,:d=>d)

H₁_𝑢, L₂_𝑢 = H₁(elements["Ω"])
println(log10(L₂_𝑢))
println(log10(H₁_𝑢))

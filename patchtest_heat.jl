
using BenchmarkTools
using ApproxOperator
using ApproxOperator.Heat: ∫∫∇v∇udxdy, ∫vbdΩ, ∫vgds, H₁

include("import_patchtest.jl")

ndiv = 2
elements, nodes = import_patchtest_fem("./msh/patchtest_"*string(ndiv)*".msh")

nₚ = length(nodes)

n = 2
u(x,y) = (x+y)^n
v(x,y) = (x+y)^n
∂u∂x(x,y) = n*(x+y)^abs(n-1)
∂u∂y(x,y) = n*(x+y)^abs(n-1)
∂²u∂x²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
∂²u∂x∂y(x,y) = n*(n-1)*(x+y)^abs(n-2)
∂²u∂y²(x,y)  = n*(n-1)*(x+y)^abs(n-2)
b(x,y,z) = -∂²u∂x²(x,y)-∂²u∂y²(x,y)

prescribe!(elements["Ω"],:k=>(x,y,z)->1.0,index=:𝑔)
prescribe!(elements["Ω"],:b=>b)
prescribe!(elements["Γ¹"],:α=>(x,y,z)->1e9,index=:𝑔)
prescribe!(elements["Γ²"],:α=>(x,y,z)->1e9,index=:𝑔)
prescribe!(elements["Γ³"],:α=>(x,y,z)->1e9,index=:𝑔)
prescribe!(elements["Γ⁴"],:α=>(x,y,z)->1e9,index=:𝑔)
prescribe!(elements["Γ¹"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ²"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ³"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Γ⁴"],:g=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))

ops = [
    ∫∫∇v∇udxdy=>elements["Ω"],
    ∫vgds=>elements["Γ"],
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)

ops[1](k)
ops[2](k,f)

d = k\f

push!(nodes,:d=>d)

# H₁, L₂ = H₁2D(elements["Ω"])

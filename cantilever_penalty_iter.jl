using Revise
using TimerOutputs 
using XLSX
using LinearSolve
using ApproxOperator
using ApproxOperator.Elasticity: ∫∫εᵈᵢⱼσᵈᵢⱼdxdy, ∫qpdΩ, ∫∫p∇udxdy, ∫vᵢgᵢds, ∫∫vᵢbᵢdxdy, ∫vᵢtᵢds, L₂, L₂𝑝, Hₑ_PlaneStress, Hₑ_PlaneStrain_Deviatoric

include("import_cantilever.jl")

ndiv = 4

indices = 2:3
nₜ = length(indices)
L₂_𝒖   = zeros(nₜ)
Hₑ_𝒖   = zeros(nₜ)
Hₑ_dev = zeros(nₜ)
L₂_𝑝   = zeros(nₜ)

for (i,n) in enumerate(indices)

# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_tri3_"*string(ndiv)*".msh","./msh/cantilever_"*string(n)*".msh",4*n,n)
# elements, nodes, nodes_p, sp, type = import_quadratic_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_"*string(n)*".msh",4*n,n)
elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_"*string(n)*".msh",4*n,n)
# elements, nodes, nodes_p, sp, type = import_quadratic_mix("./msh/cantilever_quad8_"*string(ndiv)*".msh","./msh/cantilever_"*string(n)*".msh",4*n,n)
# nx = n
# ny = 31
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_tri3_"*string(ndiv)*".msh","./msh/cantilever_"*string(ny)*"_"*string(nx)*".msh",nx,ny)
# elements, nodes, nodes_p, sp, type = import_linear_mix("./msh/cantilever_quad_"*string(ndiv)*".msh","./msh/cantilever_"*string(ny)*"_"*string(nx)*".msh",nx,ny)
# elements, nodes, nodes_p, sp, type = import_quadratic_mix("./msh/cantilever_tri6_"*string(ndiv)*".msh","./msh/cantilever_"*string(ny)*"_"*string(nx)*".msh",nx,ny)
# elements, nodes, nodes_p, sp, type = import_quadratic_mix("./msh/cantilever_quad8_"*string(ndiv)*".msh","./msh/cantilever_"*string(ny)*"_"*string(nx)*".msh",nx,ny)
nₚ = length(nodes_p)

nₑ = length(elements["Ωᵘ"])
nᵤ = length(nodes)

L = 48.0
D = 12.0
P = 1000
E = 3e6
# E = 1.0
ν = 0.5-1e-8
# ν = 0.3
Ē = E/(1.0-ν^2)
ν̄ = ν/(1.0-ν)
I = D^3/12
EI = Ē*I
Cᵢᵢᵢᵢ = E/(1+ν)/(1-2*ν)*(1-ν)
Cᵢᵢⱼⱼ = E/(1+ν)/(1-2*ν)*ν
Cᵢⱼᵢⱼ = E/(1+ν)/2

u(x,y) = -P*y/6/EI*((6*L-3*x)*x + (2+ν̄)*(y^2-D^2/4))
v(x,y) = P/6/EI*(3*ν̄*y^2*(L-x) + (4+5*ν̄)*D^2*x/4 + (3*L-x)*x^2)
∂u∂x(x,y) = -P/EI*(L-x)*y
∂u∂y(x,y) = -P/6/EI*((6*L-3*x)*x + (2+ν̄)*(3*y^2-D^2/4))
∂v∂x(x,y) = P/6/EI*((6*L-3*x)*x - 3*ν̄*y^2 + (4+5*ν̄)*D^2/4)
∂v∂y(x,y) = P/EI*(L-x)*y*ν̄

ε₁₁(x,y) = ∂u∂x(x,y)
ε₂₂(x,y) = ∂v∂y(x,y)
ε₁₂(x,y) = ∂u∂y(x,y) + ∂v∂x(x,y)
σ₁₁(x,y) = -P*(L-x)*y/I
σ₂₂(x,y) = 0.0
σ₃₃(x,y) = Cᵢᵢⱼⱼ*ε₁₁(x,y) + Cᵢᵢⱼⱼ*ε₂₂(x,y)
σ₁₂(x,y) = P/2/I*(D^2/4-y^2)
prescribe!(elements["Ωᵘ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵘ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E)
prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν)
prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂) 
prescribe!(elements["Γᵍᵘ"],:α=>(x,y,z)->1e12)
prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->(σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3)

## Debug
# n = 1
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
# σ₁₁(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*ε₁₁(x,y) + ν*ε₂₂(x,y))
# σ₂₂(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + (1-ν)*ε₂₂(x,y))
# σ₃₃(x,y) = E/(1+ν)/(1-2*ν)*(ν*ε₁₁(x,y) + ν*ε₂₂(x,y))
# σ₁₂(x,y) = E/(1+ν)*ε₁₂(x,y)
# ∂ε₁₁∂x(x,y) = ∂²u∂x²(x,y)
# ∂ε₁₁∂y(x,y) = ∂²u∂x∂y(x,y)
# ∂ε₂₂∂x(x,y) = ∂²v∂x∂y(x,y)
# ∂ε₂₂∂y(x,y) = ∂²v∂y²(x,y)
# ∂ε₁₂∂x(x,y) = 0.5*(∂²u∂x∂y(x,y) + ∂²v∂x²(x,y))
# ∂ε₁₂∂y(x,y) = 0.5*(∂²u∂y²(x,y) + ∂²v∂x∂y(x,y))

# ∂σ₁₁∂x(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂x(x,y) + ν*∂ε₂₂∂x(x,y))
# ∂σ₁₁∂y(x,y) = E/(1+ν)/(1-2*ν)*((1-ν)*∂ε₁₁∂y(x,y) + ν*∂ε₂₂∂y(x,y))
# ∂σ₂₂∂x(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂x(x,y) + (1-ν)*∂ε₂₂∂x(x,y))
# ∂σ₂₂∂y(x,y) = E/(1+ν)/(1-2*ν)*(ν*∂ε₁₁∂y(x,y) + (1-ν)*∂ε₂₂∂y(x,y))
# ∂σ₁₂∂x(x,y) = E/(1+ν)*∂ε₁₂∂x(x,y)
# ∂σ₁₂∂y(x,y) = E/(1+ν)*∂ε₁₂∂y(x,y)
# b₁(x,y) = -∂σ₁₁∂x(x,y) - ∂σ₁₂∂y(x,y)
# b₂(x,y) = -∂σ₁₂∂x(x,y) - ∂σ₂₂∂y(x,y)
# p(x,y) = (σ₁₁(x,y)+σ₂₂(x,y)+σ₃₃(x,y))/3

# prescribe!(elements["Ωˢ"],:E=>(x,y,z)->E, index=:𝑔)
# prescribe!(elements["Ωˢ"],:ν=>(x,y,z)->ν, index=:𝑔)
# prescribe!(elements["Ωᵖ"],:E=>(x,y,z)->E, index=:𝑔)
# prescribe!(elements["Ωᵖ"],:ν=>(x,y,z)->ν, index=:𝑔)
# prescribe!(elements["Ωᵍᵘ"],:E=>(x,y,z)->E, index=:𝑔)
# prescribe!(elements["Ωᵍᵘ"],:ν=>(x,y,z)->ν, index=:𝑔)
# prescribe!(elements["Ωᵘ"],:b₁=>(x,y,z)->b₁(x,y))
# prescribe!(elements["Ωᵘ"],:b₂=>(x,y,z)->b₂(x,y))
# prescribe!(elements["Γᵍᵘ"],:g₁=>(x,y,z)->u(x,y))
# prescribe!(elements["Γᵍᵘ"],:g₂=>(x,y,z)->v(x,y))
# prescribe!(elements["Γᵍᵘ"],:n₁₁=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍᵘ"],:n₂₂=>(x,y,z)->1.0)
# prescribe!(elements["Γᵍᵘ"],:n₁₂=>(x,y,z)->0.0)
# prescribe!(elements["Γᵗ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γᵗ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
# prescribe!(elements["Γʳ"],:t₁=>(x,y,z,n₁,n₂)->σ₁₁(x,y)*n₁+σ₁₂(x,y)*n₂)
# prescribe!(elements["Γʳ"],:t₂=>(x,y,z,n₁,n₂)->σ₁₂(x,y)*n₁+σ₂₂(x,y)*n₂)
# prescribe!(elements["Ωᵍᵘ"],:u=>(x,y,z)->u(x,y))
# prescribe!(elements["Ωᵍᵘ"],:v=>(x,y,z)->v(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
# prescribe!(elements["Ωᵍᵘ"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))
# prescribe!(elements["Ωᵍᵖ"],:p=>(x,y,z)->p(x,y))
## End debug

𝑎ᵘ = ∫∫εᵈᵢⱼσᵈᵢⱼdxdy=>elements["Ωᵘ"]
𝑎ᵖ = ∫qpdΩ=>elements["Ωᵖ"]
𝑏ᵖ = ∫∫p∇udxdy=>(elements["Ωᵖ"],elements["Ωᵘ"])
𝑎ᵘᵅ = ∫vᵢgᵢds=>elements["Γᵍᵘ"]
𝑓 = ∫vᵢtᵢds=>elements["Γᵗ"]
# 𝑓 = [
#     ∫vᵢtᵢds=>elements["Γᵗ"]∪elements["Γʳ"],
#     ∫∫vᵢbᵢdxdy=>elements["Ωᵘ"]
# ]

kᵘᵘ = zeros(2*nᵤ,2*nᵤ)
kᵖᵖ = zeros(nₚ,nₚ)
kᵖᵘ = zeros(nₚ,2*nᵤ)
fᵖ = zeros(nₚ)
fᵘ = zeros(2*nᵤ)

# kˢˢ = spzeros(4*nₛ*nₑ,4*nₛ*nₑ)
# kᵖᵖ = spzeros(nₚ,nₚ)
# kˢᵘ = spzeros(4*nₛ*nₑ,2*nᵤ)
# kᵖᵘ = spzeros(nₚ,2*nᵤ)
# fˢ = spzeros(4*nₛ*nₑ)
# fᵖ = spzeros(nₚ)
# fᵘ = spzeros(2*nᵤ)

𝑎ᵘ(kᵘᵘ)
𝑎ᵖ(kᵖᵖ)
𝑏ᵖ(kᵖᵘ)
𝑎ᵘᵅ(kᵘᵘ,fᵘ)
𝑓(fᵘ)

k =[-kᵘᵘ kᵖᵘ';kᵖᵘ kᵖᵖ]
f = [-fᵘ;fᵖ]
d = zeros(2*nᵤ+nₚ)
# d = k\f
prob = LinearProblem(k,f)
sol = solve(prob)
d = sol.u

𝑢₁ = d[1:2:2*nᵤ]
𝑢₂ = d[2:2:2*nᵤ]
𝑝 = d[2*nᵤ+1:2*nᵤ+nₚ]
push!(nodes,:d₁=>𝑢₁)
push!(nodes,:d₂=>𝑢₂)
push!(nodes_p,:p=>𝑝)

Hₑ_𝒖_, L₂_𝒖_ = Hₑ_PlaneStress(elements["Ωᵍᵘ"])
Hₑ_dev_ = Hₑ_PlaneStrain_Deviatoric(elements["Ωᵍᵘ"])
L₂_𝑝_ = L₂𝑝(elements["Ωᵍᵖ"])

L₂_𝒖[i] = log10(L₂_𝒖_)
Hₑ_𝒖[i] = log10(Hₑ_𝒖_)
Hₑ_dev[i] = log10(Hₑ_dev_)
L₂_𝑝[i] = log10(L₂_𝑝_)

println("n = $n, L₂_𝒖 = $L₂_𝒖_, Hₑ_𝒖 = $Hₑ_𝒖_, Hₑ_dev = $Hₑ_dev_, L₂_𝑝 = $L₂_𝑝_")

end

XLSX.openxlsx("./xlsx/plate_with_hole_linear_mix.xlsx", mode = "rw") do xf
    sheet = xf[1]
    row = "A"
    row_L₂_𝒖 = "B"
    row_Hₑ_𝒖 = "C"
    row_Hₑ_dev = "D"
    row_L₂_𝑝 = "E"
    for (n,L₂_𝒖_,Hₑ_𝒖_,Hₑ_dev_,L₂_𝑝_) in zip(indices,L₂_𝒖,Hₑ_𝒖,Hₑ_dev,L₂_𝑝)
        sheet[row*string(n)] = n
        sheet[row_L₂_𝒖*string(n)] = L₂_𝒖_
        sheet[row_Hₑ_𝒖*string(n)] = Hₑ_𝒖_
        sheet[row_Hₑ_dev*string(n)] = Hₑ_dev_
        sheet[row_L₂_𝑝*string(n)] = L₂_𝑝_
    end
end
